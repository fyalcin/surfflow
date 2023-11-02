import itertools
import pickle
from datetime import datetime
from typing import Union, Optional

import numpy as np
from pymatgen.analysis.wulff import WulffShape
from pymatgen.core import Structure
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from htflow_utils.caching import cache_results
from htflow_utils.misc_tools import generate_input_dict, make_calculation_hash
from htflow_utils.misc_tools import parse_miller
from htflow_utils.shaper import Shaper
from surfflow.phys.surfen import calculate_surface_energy_gen
from surfflow.phys.wulff import get_wulff_shape
from surfflow.utils.db_tools import VaspDB
from surfflow.utils.structure_manipulation import get_conv_bulk_from_db


def get_surfen_inputs_from_mpid(
    mpid: str,
    bulk_coll: str,
    sg_params: dict,
    sg_filter: dict,
    comp_params: dict,
    db_file: str,
    high_level: Union[str, bool],
    custom_id: Optional[str] = None,
) -> list:
    """
    Generates an input dictionary that contains all the necessary information about
    the slabs that are consistent with the given parameters for which surface energy
    calculations will be performed for.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param bulk_coll: Collection to query bulk structures from in the database.
    :type bulk_coll: str

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict

    :param sg_filter: Parameters for filtering generated slabs depending on their BVS values.
    :type sg_filter: dict

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param custom_id: Unique ID to use for the surface energy workflow. This will replace
        the unique ID generated from the hash of computational parameters.
        Use this for debugging or one-off surface energy calculations to easily
        find the results in the database.
    :type custom_id: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :raises ValueError: If the bulk structure is not found in the database.

    :return: List of input dictionaries.
    :rtype: list
    """
    tol = sg_params.get("tol")
    candidates = generate_candidate_slabs_from_mpid(
        mpid, sg_params, sg_filter, bulk_coll, db_file, high_level
    )
    inputs_list = [
        get_surfen_inputs_from_slab(c[0], c[1], tol, custom_id) for c in candidates
    ]
    inputs_list = update_inputs_list(inputs_list, comp_params)

    return inputs_list


def put_surfen_inputs_into_db(
    slab_dict: dict,
    sg_params: dict,
    comp_params: dict,
    fltr: dict,
    coll: str,
    metadata: dict,
    db_file: str = "auto",
    high_level: Union[str, bool] = True,
) -> None:
    """
    Writes all the inputs and parameters for surface energy calculations described in
    the inputs_list into the database to be later updated with the surface energies.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be performed.
    :type slab_dict: dict

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param metadata: Metadata to be added to the database.
    :type metadata: dict

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :return: None
    """
    nav = VaspDB(db_file=db_file, high_level=high_level)

    inputs = slab_dict.get("inputs")
    slab_params = slab_dict.get("slab_params")
    slab = slab_dict.get("struct")
    sg = slab_dict.get("slab_generator")

    for calc in inputs:
        struct = calc["struct"]
        uid_input = calc["uid"]
        tag = calc["tag"]
        loc = calc["loc"]
        loc_input = ".".join(loc)
        nav.update_data(
            collection=coll,
            fltr=fltr,
            new_values={
                "$set": {
                    loc_input + ".structure": struct.as_dict(),
                    loc_input + ".uid": uid_input,
                    loc_input + ".task_label": tag,
                }
            },
            upsert=True,
        )

    tol = sg_params.get("tol")
    layers = Shaper.get_layers(slab, tol)
    top_layer = [str(slab[site].species) for site in layers[max(layers)]]
    bot_layer = [str(slab[site].species) for site in layers[min(layers)]]
    terminations = {"top": top_layer, "bottom": bot_layer}

    sg_params.update(
        {
            "primitive": slab_dict["primitive"],
            "lll_reduce": slab_dict["lll_reduce"],
            "max_normal_search": slab_dict["max_normal_search"],
            "param_modified": slab_dict["sg_modified"],
        }
    )

    nav.update_data(
        collection=coll,
        fltr=fltr,
        new_values={
            "$set": {
                "structure": slab.as_dict(),
                "slab_generator": pickle.dumps(sg),
                "slab_params": slab_params,
                "sg_params": sg_params,
                "comp_params": comp_params,
                "terminations": terminations,
                "created_on": datetime.now(),
            }
        },
        upsert=True,
    )

    nav.update_data(
        collection=coll, fltr=fltr, new_values={"$set": {**metadata}}, upsert=True
    )


def generate_candidate_slabs_from_mpid(
    mpid: str,
    sg_params: dict,
    sg_filter: dict,
    coll: str = "PBE.bulk_data",
    db_file: str = "auto",
    high_level: Union[str, bool] = True,
) -> list:
    """
    Generates slabs within certain constraints from a given bulk structure.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict

    :param sg_filter: Parameters for filtering generated slabs depending on their BVS values.
    :type sg_filter: dict

    :param coll: Name of the bulk collection to load the bulk structure from. Defaults to "PBE.bulk_data".
    :type coll: str, optional

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :return: List of dictionaries containing the slabs and the inputs for the calculations needed.
    :rtype: list
    """
    bulk_conv = get_conv_bulk_from_db(mpid, coll, db_file, high_level)
    if sg_params.get("minimize_structures", False):
        params_dict = {}
        # it's difficult to predict which combinations of primitive, lll_reduce, and max_normal_search
        # will result in the smallest structures, so we generate slabs for all possible combinations
        slabs_dict = None
        for prim, lll, mns in itertools.product(
            *[[True, False], [True, False], ["max", None]]
        ):
            sg_params.update(
                {"primitive": prim, "lll_reduce": lll, "max_normal_search": mns}
            )
            slabs_dict, sg_dict = Shaper.generate_slabs(bulk_conv, sg_params)
            params_dict[(prim, lll, mns)] = {}
            for hkl, slabs in slabs_dict.items():
                # we assign costs to each combination of parameters using number of sites in structures
                weight = sum(
                    [
                        1.5 * slab.num_sites + slab.oriented_unit_cell.num_sites
                        for slab in slabs
                    ]
                )
                params_dict[(prim, lll, mns)][hkl] = {
                    "weight": round(weight),
                    "slabs": slabs,
                    "slabgen": sg_dict[hkl],
                }

        # ws = {hkl: {k: v[hkl]['weight'] for k, v in params_dict.items()} for hkl in sg_dict}
        opt_slabs_dict = {}
        opt_sg_dict = {}
        for hkl in slabs_dict:
            ws = {k: v[hkl]["weight"] for k, v in params_dict.items()}
            opt_param = min(ws, key=ws.get)
            # we choose the parameters that minimize the structures for each miller index
            # In the end, we may have different parameters for different hkl
            opt_slabs_dict[hkl] = params_dict[opt_param][hkl]["slabs"]
            opt_sg_dict[hkl] = params_dict[opt_param][hkl]["slabgen"]
        slabs_dict = opt_slabs_dict
        sg_dict = opt_sg_dict
    else:
        slabs_dict, sg_dict = Shaper.generate_slabs(bulk_conv, sg_params)

    if slabs_dict is None:
        return []

    max_nsites = sg_params.get("max_nsites", None)
    if max_nsites:
        slabs_dict = {
            hkl: [slab for slab in slabs if slab.num_sites <= max_nsites]
            for hkl, slabs in slabs_dict.items()
        }

    for hkl, slabs in slabs_dict.items():
        for slab in slabs:
            # polarity checks need oxidation state. We only add them to slabs
            # that are non-polar to not mess up hash compatibility from previous
            # calculations that did not have these oxidation states.
            tmp = slab.copy()
            tmp.add_oxidation_state_by_guess()
            if tmp.is_polar():
                slab.add_oxidation_state_by_guess()

        if sg_params["filter_polar"]:
            slabs_dict[hkl] = [slab for slab in slabs if not slab.is_polar()]

    # this is where further filtering with respect to bond valence sum happens
    method = sg_filter.get("method")
    if method:
        slabs_list = [
            slab for slabs_by_hkl in slabs_dict.values() for slab in slabs_by_hkl
        ]
        bvs_list = [slab.energy["bvs_per_area"] for slab in slabs_list]
        if method == "bvs_threshold":
            # here, we take all slabs with BVS values within a threshold of the minimum one encountered
            bvs_tol = sg_filter.get("bvs_param")
            min_bvs = min(bvs_list)
            filtered_slabs = [
                slab
                for (index, slab) in enumerate(slabs_list)
                if bvs_list[index] / min_bvs - 1 < bvs_tol
            ]
        elif method == "bvs_min_N":
            # here, we take N slabs with the lowest bond valence sums
            n = sg_filter.get("bvs_param")
            if n < len(slabs_list):
                sorted_ind = np.argsort(bvs_list)
                filtered_slabs = [slabs_list[i] for i in sorted_ind[:n]]
            else:
                filtered_slabs = slabs_list
        elif method == "bvs_min_N_hkl":
            # if we instead want at least some slabs for every hkl, we can use this to get the N
            # slabs with the lowest BVS for each hkl
            n = sg_filter.get("bvs_param")
            filtered_slabs = []
            for hkl, slabs in slabs_dict.items():
                if n < len(slabs):
                    bvs_list_hkl = [slab.energy["bvs_per_area"] for slab in slabs]
                    sorted_ind = np.argsort(bvs_list_hkl)
                    filtered_slabs += [slabs[i] for i in sorted_ind[:n]]
                else:
                    filtered_slabs += slabs
        else:
            filtered_slabs = slabs_list
        slabs_list = filtered_slabs
    else:
        slabs_list = [slab for slab_list in slabs_dict.values() for slab in slab_list]
    return [(slab, sg_dict[slab.miller_index]) for slab in slabs_list]


def get_surfen_inputs_from_slab(
    slab: Slab, sg: SlabGenerator = None, tol: float = 0.1, custom_id: str = None
) -> dict:
    """
    Generates a dictionary containing all the sub structures needed to compute the
    surface energy of the given slab.

    :param slab: Pymatgen Slab object.
    :type slab: pymatgen.core.surface.Slab

    :param sg: Pymatgen SlabGenerator object. Necessary to generate substructures for
        non-stoichiometric slabs with a specific termination. Defaults to None.
    :type sg: pymatgen.core.surface.SlabGenerator, optional

    :param tol: Tolerance parameters used in the identification of layers, given in units of Angstroms. Defaults to 0.1.
    :type tol: float, optional

    :param custom_id: Unique id to identify the database entry linked to the surface energy
        calculation of the slab. Defaults to None.
    :type custom_id: str, optional

    :return: Dictionary containing a summary of the parameters of the input slab,
        the slab itself, and the substructures needed in surface energy calculation.
    :rtype: dict

    """
    # symmetry and stoichiometry of the input slab is identified in order to add
    # the correct calculations to the inputs_dict
    id_slab = Shaper.identify_slab(slab)
    sym = id_slab["symmetric"]
    sto = id_slab["stoichiometric"]

    # oriented unit cell is used for the reference bulk energies
    ouc = slab.oriented_unit_cell
    slab_layers = len(Shaper.get_layers(slab, tol))
    slab_thickness = Shaper.get_proj_height(slab, "slab")
    vac_thickness = round(Shaper.get_proj_height(slab, "vacuum"), 3)
    ouc_input = generate_input_dict(ouc, "static", "ouc")
    millerstr = "".join([str(i) for i in slab.miller_index])
    inputs_dict = {
        "struct": slab,
        "slab_generator": sg,
        "primitive": sg.primitive,
        "lll_reduce": sg.lll_reduce,
        "max_normal_search": sg.max_normal_search,
        "sg_modified": getattr(slab, "param_modified", False),
        "custom_id": custom_id,
        "inputs": [ouc_input],
        "slab_params": {
            "sym": sym,
            "sto": sto,
            "layer_tol": tol,
            "thickness_layers": slab_layers,
            "thickness_A": np.round(slab_thickness, 6),
            "vac_thickness_A": np.round(vac_thickness, 6),
            "hkl": millerstr,
            # "bvs": slab.energy,
            "area": np.round(slab.surface_area, 6),
        },
    }

    try:
        ouc_layers = len(Shaper.get_layers(ouc, tol))
    except ValueError:
        ouc_layers = 1

    try:
        pmg_layer_size = slab.pmg_layer_size
    except AttributeError:
        print(
            'Your slab does not have a "pmg_layer_size" attribute which is needed to\n'
            "determine if a slab has complementary terminations. The workflow will proceed\n"
            "assuming that your slab has complementary terminations, which could lead to\n"
            "incorrect surface energies."
        )
        comp = True
        pmg_layer_size = ouc_layers
    else:
        comp = True if slab_layers % pmg_layer_size == 0 else False
        inputs_dict["slab_params"].update({"comp": comp})

    slab_relax_input = generate_input_dict(slab, "relax", "slab_relax")
    inputs_dict["inputs"] += [slab_relax_input]
    if sym:
        if not sto:
            # for non-stoichiometric slabs, we need the cleavage energy for which we need
            # a stoichiometric version of the slab with terminations that are complementary
            # to each other, this is done by SlabGenerator.get_slab() which always creates
            # a stoichiometric slab.
            sto_slab = sg.get_slab(slab.shift, tol)
            sto_slab_layers = len(Shaper.get_layers(sto_slab, tol))
            # Since SlabGenerator creates a larger slab than we want, we remove layers
            # and the number of layers removed is an integer multiple of the number of layers
            # in the d_hkl portion of the oriented unit cell to preserve stoichiometry
            layers_to_remove = int(
                pmg_layer_size * np.floor((sto_slab_layers - slab_layers) / ouc_layers)
            )
            target_layers = sto_slab_layers - layers_to_remove
            sto_slab = Shaper.resize(sto_slab, target_layers, vac_thickness, tol)
            # sto_slab = Shaper.remove_layers(sto_slab, layers_to_remove, tol, method='layers')
            sto_slab_input = generate_input_dict(sto_slab, "static", "sto_slab")
            slab_static_input = generate_input_dict(slab, "static", "slab_static")
            inputs_dict["inputs"] += [slab_static_input, sto_slab_input]
    else:
        # Asymmetric slabs have different surface energies on the top and the bottom,
        # which means we need to relax those regions separately.
        slab_tf = Shaper.fix_regions(slab, tol, fix_type="top_half")
        slab_bf = Shaper.fix_regions(slab, tol, fix_type="bottom_half")

        slab_tf_input = generate_input_dict(slab_tf, "relax", "slab_top_fixed_relax")
        slab_bf_input = generate_input_dict(slab_bf, "relax", "slab_bot_fixed_relax")
        slab_static_input = generate_input_dict(slab, "static", "slab_static")

        inputs_dict["inputs"] += [slab_tf_input, slab_bf_input, slab_static_input]

        if not comp:
            # For non-stoichiometric slabs, we need the stoichiometric versions in order
            # to calculate the cleavage energy which is a component in the surface energy
            # calculation.
            all_shifts = sg._calculate_possible_shifts(tol)
            top_shift_index = all_shifts.index(slab.shift)
            bot_shift_index = (top_shift_index - slab_layers) % len(all_shifts)
            sto_slab_top = sg.get_slab(slab.shift, tol)
            sto_slab_bot = sg.get_slab(all_shifts[bot_shift_index], tol)
            sto_slab_layers = len(Shaper.get_layers(sto_slab_top, tol))
            layers_to_remove = int(
                pmg_layer_size * np.floor((sto_slab_layers - slab_layers) / ouc_layers)
            )
            target_layers = sto_slab_layers - layers_to_remove
            sto_slab_top = Shaper.resize(
                sto_slab_top, target_layers, vac_thickness, tol
            )
            sto_slab_bot = Shaper.resize(
                sto_slab_bot, target_layers, vac_thickness, tol
            )
            sto_slab_top_input = generate_input_dict(
                sto_slab_top, "static", "sto_slab_top"
            )
            sto_slab_bot_input = generate_input_dict(
                sto_slab_bot, "static", "sto_slab_bot"
            )
            inputs_dict["inputs"] += [sto_slab_top_input, sto_slab_bot_input]

    return inputs_dict


def generate_wulff_data_from_mpid(
    mpid: str, coll: str, db_file="auto", high_level: Union[str, bool] = True
) -> None:
    """
    Goes through all the calculated surface energies in the database and finds the minimum
    of each orientation, creating a new entry in the same document with an easy-to-navigate
    list of surface energies.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :return: None
    """
    nav = VaspDB(db_file=db_file, high_level=high_level)
    data = nav.find_data(coll, {"mpid": mpid})
    surfen_dict = {}
    for hkl, miller_data in data["miller_list"].items():
        tmp_hkl = {}
        rel_str = {}
        for uid, slab_data in miller_data.items():
            try:
                surfen = slab_data["surface_energy"]
            except KeyError:
                continue
            else:
                tmp_hkl[uid] = [surfen["top"], surfen["bottom"]]
                rel_str[uid] = {
                    k: v["output"]["structure"]
                    for k, v in slab_data["calcs"].items()
                    if k.endswith("relax")
                }
        if not tmp_hkl:
            continue
        min_uid = min(tmp_hkl, key=tmp_hkl.get)
        min_surfen = min(tmp_hkl.values())
        params = ["structure", "comp_params", "slab_params"]
        surfen_dict[hkl] = {"surface_energy": min(min_surfen), "uid": min_uid}
        params_dict = {k: v for k, v in miller_data[min_uid].items() if k in params}
        relaxed_structure = rel_str[min_uid]
        surfen_dict[hkl].update({"relaxed_structure": relaxed_structure}, **params_dict)

    nav.update_data(
        collection=coll,
        fltr={"mpid": mpid},
        new_values={"$set": {"surfen_list": surfen_dict}},
        upsert=True,
    )


def update_inputs_list(inputs_list: list, comp_params: dict) -> list:
    """
    Generates Workflows for the calculations needed for surface energy and
    updates the input dictionary with them.

    :param inputs_list: List containing information about the inputs necessary for surface
        energy calculations for a given material represented by its Materials Project ID.
        Generated by the function get_surfen_inputs_from_slab by inputting a Slab and
        a SlabGenerator corresponding to this Slab object.
    :type inputs_list: list

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :return: inputs_list with each element updated with the Workflows, unique IDs, locations, and tags.
    :rtype: list
    """
    for slab_data in inputs_list:
        slab = slab_data.get("struct")
        slab_params = slab_data.get("slab_params")
        # we grab the fractional coords and round them to get the unique ID
        # we also handle the PBC by settings all frac coords that are 1.0 to 0.0
        # because depending on the platform the code is running on, the sites at
        # the boundaries sometimes have a fractional coord of 1.0 and sometimes
        # 0.0.
        # if there is a custom_id, we use it, otherwise we use the unique ID
        # uid_slab is only used to define the location of the entry in the database,
        # it is not used in the calculations.
        uid_slab = slab_data.get("custom_id")
        if not uid_slab:
            uid_slab = make_calculation_hash(
                structure=slab, comp_params=comp_params, slab_params=slab_params
            )
        # slab_params.update({"uid": uid_slab})
        slab_data.update({"uid": uid_slab})
        millerstr = slab_params.get("hkl")
        inputs = slab_data.get("inputs")
        for calc in inputs:
            struct = calc.get("struct")
            calc_tag = calc.get("calc_tag")
            calc_type = calc.get("calc_type")
            loc = ["calcs", calc_tag]
            # we do the same thing for the substructures, so this way we have uids
            # for the substructures as well, in case any calculation needs to use
            # the same reference structure, i.e. the same miller index but different
            # termination.
            struct_frac_coords = np.round(struct.frac_coords, 6)
            struct_frac_coords[struct_frac_coords == 1.0] = 0.0
            # we assign a different uid to each calculation, separate from slab_uid which
            # is only used to determine the location in the DB, so that even if
            # we pass a custom_id to the calculation for any reason, the workflow
            # will still find any previous calculation if the same inputs are used.
            uid_input = make_calculation_hash(structure=struct, comp_params=comp_params)
            formula = struct.formula
            tag = f"{formula}_{millerstr}_{calc_tag}_{uid_input}"
            tag_dict = {
                "tag": tag,
                "loc": loc,
                "calc_type": calc_type,
                "uid": uid_input,
            }
            calc.update(tag_dict)

    return inputs_list


def write_surface_energies_to_db(
    slab_dict: dict,
    fltr: dict,
    coll: str,
    db_file: str = "auto",
    high_level: Union[str, bool] = True,
) -> None:
    """
    Calculates the surface energies of the slabs given in the inputs_list and writes
    the surface energies to the database under the relevant entries.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be performed.
    :type slab_dict: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :return: None
    """
    nav = VaspDB(db_file, high_level)

    surf_en, ens = calculate_surface_energy_gen(
        slab_dict, fltr, coll, db_file, high_level
    )

    nav.update_data(
        collection=coll,
        fltr=fltr,
        new_values={"$set": {"surface_energy": surf_en, "surfen_components": ens}},
        upsert=True,
    )


def generate_wulff_shape(
    mpid: str,
    surfen_coll: str,
    bulk_coll: str,
    wulff_coll: str,
    db_file: str = "auto",
    high_level: Union[str, bool] = True,
) -> Union[None, WulffShape]:
    """
    Calculates the Wulff shape of a given material.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param surfen_coll: Collection to query the results from in the database.
    :type surfen_coll: str

    :param bulk_coll: Collection to query bulk structures from in the database.
    :type bulk_coll: str

    :param wulff_coll: Collection to write the Wulff shape to in the database.
    :type wulff_coll: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :return: WulffShape object containing the Wulff shape of the material.
    :rtype: WulffShape
    """

    nav = VaspDB(db_file=db_file, high_level=high_level)
    bulk = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid})
    try:
        bulk_struct = Structure.from_dict(bulk["structure_fromMP_conv"])
    except KeyError:
        try:
            bulk_struct = Structure.from_dict(bulk["structure_fromMP"])
        except KeyError:
            print("No bulk structure found for this material.")
            return None
        else:
            bulk_struct = SpacegroupAnalyzer(
                bulk_struct
            ).get_conventional_standard_structure(keep_site_properties=True)

    lattice = bulk_struct.lattice

    data = list(
        nav.find_many_data(
            collection=surfen_coll,
            fltr={"mpid": mpid},
            projection={"mpid": 1, "hkl": 1, "surface_energy": 1, "uid": 1},
        )
    )
    if not data:
        print(f"No entry found for mpid {mpid}")
        return None

    miller_surfen_dict = {
        f'{d["hkl"]}-{d["uid"][:4]}': min(d["surface_energy"].values()) for d in data
    }

    nav.update_data(
        collection=wulff_coll,
        fltr={"mpid": mpid},
        new_values={
            "$set": {
                "miller_surfen_dict": miller_surfen_dict,
            }
        },
        upsert=True,
    )

    miller_list = [v["hkl"] for v in data]
    surfen_list = [min(v["surface_energy"].values()) for v in data]
    miller_parsed = [parse_miller(m) for m in miller_list]
    try:
        ws = get_wulff_shape(lattice, miller_parsed, surfen_list)
    except ValueError:
        print(f"Wulff shape calculation failed for mpid {mpid}")
        return None
    else:
        nav.update_data(
            collection=wulff_coll,
            fltr={"mpid": mpid},
            new_values={
                "$set": {
                    "wulff_shape": pickle.dumps(ws),
                    "miller_energy_dict": {
                        "".join([str(m) for m in k]): v
                        for k, v in ws.miller_energy_dict.items()
                    },
                    "area_fraction_dict": {
                        "".join([str(m) for m in k]): v
                        for k, v in ws.area_fraction_dict.items()
                    },
                }
            },
            upsert=True,
        )

    return ws
