from typing import Union, Any

from surfflow.utils.db_tools import get_entry_by_loc, VaspDB


def calculate_surface_energy_gen(
        slab_dict: dict,
        fltr: dict,
        coll: str,
        db_file: str = "auto",
        high_level: Union[str, bool] = True,
) -> tuple[dict[str, float | Any], dict[str, float | Any]]:
    """
    Calculates the surface energies of top and bottom terminations of the slab given as a dictionary.

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

    :return: Surface energies of the top and bottom terminations of the Slab described by the slab_dict.
    :rtype: dict
    """
    nav = VaspDB(db_file, high_level)
    slab_params = slab_dict.get("slab_params")

    sym = slab_params.get("sym")
    sto = slab_params.get("sto")
    comp = slab_params.get("comp")
    area = slab_params.get("area")
    en_dict = {}

    inputs = slab_dict.get("inputs")
    for entry in inputs:
        loc = entry["loc"]
        calc_tag = entry["calc_tag"]
        result = get_entry_by_loc(nav, fltr, coll, loc)
        output = result["output"]
        en_dict[calc_tag] = {
            "nsites": len(output["structure"]["sites"]),
            "energy": output["energy"],
            "energy_per_atom": output["energy_per_atom"],
        }

    ens = {}
    bulk_en = en_dict["ouc"]["energy_per_atom"]
    if sym:
        slab_relax_en = en_dict["slab_relax"]["energy"]
        if sto:
            nsites = en_dict["slab_relax"]["nsites"]
            surf_en_top = (slab_relax_en - nsites * bulk_en) / 2
            surf_en_bot = surf_en_top
            ens["surf_en_top"] = surf_en_top
            ens["surf_en_bot"] = surf_en_bot

        else:
            slab_relax_en = en_dict["slab_relax"]["energy"]
            slab_static_en = en_dict["slab_static"]["energy"]
            sto_slab_en = en_dict["sto_slab"]["energy"]

            sto_slab_nsites = en_dict["sto_slab"]["nsites"]

            e_cle = (sto_slab_en - sto_slab_nsites * bulk_en) / 2
            e_rel = (slab_relax_en - slab_static_en) / 2

            surf_en_top = e_cle + e_rel
            surf_en_bot = surf_en_top
            ens["E_cle"] = e_cle
            ens["E_rel"] = e_rel

    else:
        slab_tf_relax_en = en_dict["slab_top_fixed_relax"]["energy"]
        slab_bf_relax_en = en_dict["slab_bot_fixed_relax"]["energy"]
        slab_static_en = en_dict["slab_static"]["energy"]

        e_rel_top = slab_bf_relax_en - slab_static_en
        e_rel_bot = slab_tf_relax_en - slab_static_en
        ens["E_rel_top"] = e_rel_top
        ens["E_rel_bot"] = e_rel_bot

        if comp:
            nsites = en_dict["slab_bot_fixed_relax"]["nsites"]
            e_cle = (slab_static_en - nsites * bulk_en) / 2
            surf_en_top = e_cle + e_rel_top
            surf_en_bot = e_cle + e_rel_bot
            ens["E_cle"] = e_cle
        else:
            sto_slab_top_en = en_dict["sto_slab_top"]["energy"]
            sto_slab_bot_en = en_dict["sto_slab_bot"]["energy"]

            sto_slab_top_nsites = en_dict["sto_slab_top"]["nsites"]
            sto_slab_bot_nsites = en_dict["sto_slab_bot"]["nsites"]

            e_cle_top = (sto_slab_top_en - sto_slab_top_nsites * bulk_en) / 2
            e_cle_bot = (sto_slab_bot_en - sto_slab_bot_nsites * bulk_en) / 2

            surf_en_top = e_cle_top + e_rel_top
            surf_en_bot = e_cle_bot + e_rel_bot
            ens["E_cle_top_term"] = e_cle_top
            ens["E_cle_bot_term"] = e_cle_bot

    surf_en_top *= 16.02176565 / area
    surf_en_bot *= 16.02176565 / area

    return {"top": surf_en_top, "bottom": surf_en_bot}, ens
