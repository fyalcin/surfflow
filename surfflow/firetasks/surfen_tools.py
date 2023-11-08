import os
from typing import Optional, List, Union

from fireworks import explicit_serialize, FiretaskBase, FWAction, Workflow, Firework

from htflow_utils.misc_tools import make_calculation_hash
from htflow_utils.optimade import optimadeify_document
from htflow_utils.vasp_tools import get_vis
from htflow_utils.workflows import get_calc_wf, use_fake_vasp
from surfflow.firetasks.core import SurfenFT
from surfflow.firetasks.database import MoveResults
from surfflow.utils.db_tools import get_entry_by_loc, VaspDB
from surfflow.utils.misc_tools import check_input
from surfflow.utils.structure_manipulation import find_slabs
from surfflow.utils.surfen_tools import (
    put_surfen_inputs_into_db,
    write_surface_energies_to_db,
    generate_wulff_shape,
)

module_dir = os.path.dirname(os.path.abspath(__file__))
fake_vasp_dir = os.path.join(module_dir, "..", "utils", "fake_vasp")


@explicit_serialize
class WriteSurfenInputsToDB(SurfenFT):
    """
    Firetask that writes all the inputs and parameters for surface energy calculations described in
    the inputs_list into the database to be later updated with the surface energies.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be
    performed.
    :type slab_dict: dict

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
    task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.
    """

    _fw_name = "Put the structures and computational parameters for the surface energy calculation into the DB"
    required_params = ["slab_dict", "sg_params", "comp_params", "fltr", "coll"]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        put_surfen_inputs_into_db(**self._inp)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class RelaxSurfaceEnergyInputs(FiretaskBase):
    """Perform VASP calculations on the input structures in order to calculate
    surface energy.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be
    performed.
    :type slab_dict: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
    task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    :param fake_vasp: Whether to use the fake vasp or not.
    :type fake_vasp: bool, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.

    """

    _fw_name = "Perform various VASP calculations on the structures in order to calculate surface energy."
    required_params = ["slab_dict", "fltr", "coll", "comp_params"]
    optional_params = ["db_file", "high_level", "fake_vasp"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}

        wf_list = get_surface_energy_wfs(**inp)
        return FWAction(detours=wf_list, update_spec=fw_spec)


def calculate_and_move(
        input_dict: dict,
        fltr: dict,
        coll: str,
        comp_params: dict,
        hkl: str,
        db_file: str = "auto",
        high_level: Union[str, bool] = True,
        fake_vasp: bool = False,
) -> Optional[Workflow]:
    """
    Perform a VASP calculation on the input structure and move the results to the database.

    :param input_dict: Dictionary containing all the necessary information about the structures and calculations to
    be performed.
    :type input_dict: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
    task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param hkl: Miller index of the main slab. Only used in the workflow and firework names.
    :type hkl: str

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    Defaults to True, in which case the
        value in the db.json file is used.
    :type high_level: bool, str, optional

    :param fake_vasp: Whether to use the fake vasp or not. Defaults to False.
    :type fake_vasp: bool, optional

    :return: Workflow containing the Fireworks that will perform the VASP calculations.
    :rtype: Workflow
    """
    nav_high = VaspDB(db_file=db_file, high_level=high_level)

    struct = input_dict.get("struct")
    calc_type = input_dict.get("calc_type")
    calc_tag = input_dict.get("calc_tag")
    tag = input_dict.get("tag")
    loc = input_dict.get("loc")

    # If the calculation is already in the high_level database, we do not return
    # a new workflow for this structure.
    result_high = get_entry_by_loc(nav_high, fltr, coll, loc).get("output")
    if result_high:
        return None

    comp = str(struct.composition)
    wf_move = Workflow(
        [
            Firework(
                MoveResults(
                    tag=tag,
                    fltr=fltr,
                    coll=coll,
                    to_loc=loc,
                    from_loc=[],
                    db_file=db_file,
                    high_level=high_level,
                    fake_calc=fake_vasp,
                ),
                name=f"Move {comp}-{hkl}-{calc_tag} results FW.",
            )
        ],
        name=f"Move {comp}-{hkl}-{calc_tag} results WF.",
    )

    vis = get_vis(struct, comp_params, calc_type)
    wf_calc = get_calc_wf(
        struct=struct, vis=vis, tag=tag, add_static=True, db_file=db_file
    )
    if fake_vasp and wf_calc:
        wf_calc = use_fake_vasp(
            wf_calc,
            ref_dir=fake_vasp_dir,
            clear_inputs=False,
            check_incar=False,
            check_poscar=False,
            check_potcar=False,
            check_kpoints=False,
        )

    if wf_calc:
        wf_calc.append_wf(wf_move, wf_calc.leaf_fw_ids)
        return wf_calc
    else:
        return wf_move


def get_surface_energy_wfs(
        slab_dict, fltr, coll, comp_params, db_file="auto", high_level=True, fake_vasp=False
) -> List[Workflow]:
    """Perform VASP calculations on the input structures in order to calculate
    surface energy.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be
    performed.
    :type slab_dict: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
    task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
        is checked from the environment variable FW_CONFIG_FILE.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :param fake_vasp: Whether to use the fake vasp or not. Defaults to False.
    :type fake_vasp: bool, optional

    :return: List of Workflows containing the Fireworks that will perform the VASP calculations.
    :rtype: List[Workflow]
    """
    wf_list = []
    hkl = slab_dict["slab_params"]["hkl"]
    for entry in slab_dict.get("inputs"):
        new_wf = calculate_and_move(
            entry, fltr, coll, comp_params, hkl, db_file, high_level, fake_vasp
        )
        if new_wf:
            wf_list.append(new_wf)
    return wf_list


@explicit_serialize
class GenerateOptimadeEntry(FiretaskBase):
    """Firetask that generates a structure from the given slab dictionary and
    adds it to the OPTIMADE database.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be
    performed.
    :type slab_dict: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
    task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional
    """

    _fw_name = "Generate an OPTIMADE entry from the slab dictionary"
    required_params = ["fltr", "surfen_coll"]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        fltr = inp["fltr"]
        surfen_coll = inp["surfen_coll"]
        db_file = inp["db_file"]
        high_level = inp["high_level"]

        nav = VaspDB(db_file=db_file, high_level=high_level)
        slab_entry = nav.find_data(collection=surfen_coll, fltr=fltr)
        optimade_entry = optimadeify_document(slab_entry, "surface")
        nav.update_data(
            collection="optimade",
            fltr=fltr,
            new_values={"$set": {**optimade_entry}},
            upsert=True,
        )

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class WriteSurfaceEnergiesToDB(FiretaskBase):
    """Calculates the surface energies and updates the LEO collection entries
    with them.

    :param slab_dict: Dictionary containing all the necessary info on the structures and the calculations to be
    performed.
    :type slab_dict: dict

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
    task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.
    """

    _fw_name = "Calculate the surface energies and write to DB"
    required_params = ["slab_dict", "fltr", "coll"]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        write_surface_energies_to_db(**inp)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class GenerateWulffShape(FiretaskBase):
    """Firetask that calculates the Wulff shape of a given structure.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param surfen_coll: Collection to query the surface energy results from in the database.
    :type surfen_coll: str

    :param bulk_coll: Name of the bulk structure collection in the database.
    :type bulk_coll: str

    :param wulff_coll: Name of the Wulff shape collection in the database.
    :type wulff_coll: str

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.
    """

    _fw_name = "Calculate the Wulff shape of a given structure"
    required_params = ["mpid", "surfen_coll", "bulk_coll", "wulff_coll"]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        generate_wulff_shape(**inp)
        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FullyRelaxSlab(FiretaskBase):
    """Perform a full relaxation on a slab efficiently if not already present.

    :param input_dict: Dictionary containing the information about the slab to be relaxed.
    :type input_dict: dict

     :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or
     task_label.
    :type fltr: dict

    :param coll: Collection to query the results from in the database.
    :type coll: str

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    :param fake_vasp: Whether to use fake_vasp or not.
    :type fake_vasp: bool, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.
    """

    _fw_name = "Perform various VASP calculations on the structures in order to calculate surface energy."
    required_params = ["input_dict", "fltr", "coll", "comp_params"]
    optional_params = ["db_file", "high_level", "fake_vasp"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}

        fltr = inp["fltr"]
        collection = inp["coll"]
        db_file = inp["db_file"]
        high_level = inp["high_level"]
        comp_params = inp["comp_params"]
        target_loc = []

        is_relaxed, slab = find_slabs(
            fltr=fltr,
            collection=collection,
            db_file=db_file,
            high_level=high_level,
        )

        if is_relaxed != "fully_relaxed":
            vis = get_vis(struct=slab, comp_params=comp_params, calc_type="relax")
            tag = make_calculation_hash(structure=slab, comp_params=comp_params)
            wf = get_calc_wf(
                struct=slab, vis=vis, tag=tag, add_static=False, db_file=db_file
            )
            move_fw_wf = Workflow(
                [
                    Firework(
                        MoveResults(
                            tag=tag,
                            fltr=fltr,
                            coll=collection,
                            to_loc=target_loc,
                            from_loc=["output"],
                            project_fields=["structure"],
                            rename=["fully_relaxed_slab"],
                        )
                    )
                ]
            )

            wf.append_wf(move_fw_wf, wf.leaf_fw_ids)
            return FWAction(detours=wf, update_spec=fw_spec)
        else:
            nav = VaspDB(db_file=db_file, high_level=high_level)
            nav.update_data(
                collection=collection,
                fltr=fltr,
                new_values={"$set": {"fully_relaxed_slab": slab.as_dict()}},
            )
            return
