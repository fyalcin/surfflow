from typing import List, Union, Optional

from fireworks import Workflow, Firework
from pymatgen.core.surface import Slab
from pymatgen.core.surface import SlabGenerator

from htflow_utils.misc_tools import get_pmg_sg_params
from htflow_utils.misc_tools import make_calculation_hash
from htflow_utils.shaper import Shaper
from htflow_utils.vasp_tools import get_vis
from htflow_utils.workflows import get_calc_wf
from surfflow.firetasks.database import MoveResults
from surfflow.utils.structure_manipulation import get_conv_bulk_from_db


def run_conv_calc_and_move_results(
    comp_params: dict,
    conv_coll: str,
    fltr: dict,
    slab: Slab,
    db_file: str,
    high_level: Union[bool, str],
) -> List[Workflow]:
    """
    Run a single relaxation on a slab and move the results to the database.

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param conv_coll: Name of the collection to store the convergence results in.
    :type conv_coll: str

    :param fltr: Filter to use when looking up results in the database. We use the MPID for convergence.
    :type fltr: dict

    :param slab: Slab to calculate the surface energy for.
    :type slab: Slab

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: Union[bool, str], optional

    :return: List of Workflows containing the Fireworks that will perform the VASP calculations.
    :rtype: List[Workflow]
    """
    vis = get_vis(struct=slab, comp_params=comp_params, calc_type="relax")
    formula = slab.formula
    hklstr = "".join([str(i) for i in slab.miller_index])

    slab_thick = len(Shaper.get_layers(slab, 0.05))

    uid = make_calculation_hash(structure=slab, comp_params=comp_params)
    tag = f"{formula}-{hklstr}-{slab_thick}-{uid}"
    wf_list = []
    wf_calc = get_calc_wf(struct=slab, vis=vis, tag=tag, db_file=db_file)
    wf_move = Workflow(
        [
            Firework(
                MoveResults(
                    tag=tag,
                    fltr=fltr,
                    coll=conv_coll,
                    from_loc=[],
                    to_loc=["conv_data", str(slab_thick)],
                    project_fields=["output"],
                    db_file=db_file,
                    high_level=high_level,
                )
            )
        ]
    )
    if wf_calc:
        wf_calc.append_wf(wf_move, wf_calc.leaf_fw_ids)
        wf_list.append(wf_calc)
    else:
        wf_list.append(wf_move)
    return wf_list


def generate_slab_for_conv(
    mpid: str,
    miller: tuple,
    shift: float,
    conv_params: dict,
    db_file: str,
    en_dict: dict,
    high_level: Union[bool, str],
    sg_params: dict,
    bulk_coll: str,
) -> Optional[tuple]:
    """
    Generate a single slab for convergence calculations.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param miller: Miller indices of the slab.
    :type miller: tuple

    :param shift: Shift of the slab.
    :type shift: float

    :param conv_params: Dictionary containing the convergence parameters.
    :type conv_params: dict

    :param db_file: Full path of the db.json file.
    :type db_file: str

    :param en_dict: Dictionary containing the energies of the slabs.
    :type en_dict: dict

    :param high_level: Name of the high level database to use. If False, the
        low-level database will be used. If True, the high-level database will be
        queried from the db.json file.
    :type high_level: Union[bool, str]

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict

    :param bulk_coll: Name of the bulk structure collection in the database.
    :type bulk_coll: str

    :return: Slab and slab thickness.
    :rtype: tuple
    """
    miller = tuple(miller)
    sg_params["miller"] = miller
    bulk_conv = get_conv_bulk_from_db(
        mpid=mpid, coll=bulk_coll, db_file=db_file, high_level=high_level
    )
    pmg_sg_params = get_pmg_sg_params(
        bulk_conv=bulk_conv, miller=miller, sg_params=sg_params
    )
    sg = SlabGenerator(**pmg_sg_params)
    _, pmg_layer_size = Shaper.get_pmg_layer_size(
        bulk_conv=bulk_conv, miller=miller, sg=sg, tol=sg_params["tol"]
    )

    slab_thick_arr = [int(a) for a in en_dict]
    if not slab_thick_arr:
        slab_thick = conv_params["thick_start"]
    else:
        slab_thick = slab_thick_arr[-1] + pmg_layer_size
        if slab_thick > conv_params["thick_end"]:
            return None
    sg_params["slab_thick"] = slab_thick
    slabs_dict, _ = Shaper.generate_slabs(bulk_conv=bulk_conv, sg_params=sg_params)
    slab = [slab for slab in slabs_dict[miller] if round(slab.shift, 6) == shift][0]

    actual_slab_thick = len(Shaper.get_layers(struct=slab, tol=sg_params["tol"]))

    return slab, actual_slab_thick


def is_list_converged(
    list_to_check: list, nelm: int = 4, threshold: float = 0.01
) -> bool:
    """
    Check if the values in a list are converged.

    :param list_to_check: List of values to check.
    :type list_to_check: list

    :param nelm: Number of elements to check.
    :type nelm: int

    :param threshold: Threshold for convergence.
    :type threshold: float

    :return: True if the list is converged, False otherwise.
    :rtype: bool
    """
    if len(list_to_check) < nelm:
        return False

    for i in range(nelm):
        if (
            abs(list_to_check[-i] - list_to_check[-i - 1])
            > threshold * list_to_check[-1]
        ):
            return False
    return True
