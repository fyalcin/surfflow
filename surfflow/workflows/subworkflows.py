import logging

from fireworks import Firework, Workflow

from surfflow.firetasks.surfen_tools import (
    RelaxSurfaceEnergyInputs,
    WriteSurfaceEnergiesToDB,
    FullyRelaxSlab,
    GenerateOptimadeEntry
)
from surfflow.utils.db_tools import VaspDB
from surfflow.utils.surfen_tools import (
    get_surfen_inputs_from_mpid,
    put_surfen_inputs_into_db,
)

logger = logging.getLogger(__name__)


def surface_energy_swf(
        mpid: str,
        comp_params: dict,
        sg_params: dict,
        sg_filter: dict,
        database_params: dict,
) -> tuple[list[Workflow], list]:
    """
    Create a workflow for calculating surface energies.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict

    :param sg_filter: Parameters for filtering generated slabs depending on their BVS values.
    :type sg_filter: dict

    :param db_file: Full path of the db.json file.
    :type db_file: str

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
        Defaults to True, in which case the value in the db.json file is used.
    :type high_level: bool, str, optional

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict

    :param custom_id: Custom ID for the workflow.
    :type custom_id: str

    :param bulk_coll: Name of the collection where the bulk structures are stored.
    :type bulk_coll: str

    :param surfen_coll: Name of the collection where the surface energies are stored.
    :type surfen_coll: str

    :param add_full_relax: If true, we relax the slab and copy it to the high_level
        database, even if the symmetry of the slab does not require a full
        relaxation for the computation of the surface energy. Defaults to True.
    :type add_full_relax: bool

    :return: List of workflows.
    :rtype: list[Workflow]
    """
    db_file = database_params["db_file"]
    high_level = database_params["high_level"]
    bulk_coll = database_params["bulk_coll"]
    surfen_coll = database_params["surfen_coll"]
    custom_id = database_params.get("custom_id", None)
    add_full_relax = database_params.get("add_full_relax", True)
    inputs_list = get_surfen_inputs_from_mpid(
        mpid=mpid,
        bulk_coll=bulk_coll,
        sg_params=sg_params,
        sg_filter=sg_filter,
        comp_params=comp_params,
        db_file=db_file,
        high_level=high_level,
        custom_id=custom_id,
    )

    wf_list = []
    for slab_dict in inputs_list:
        wf = surface_energy_swf_from_slab_dict(
            mpid=mpid,
            slab_dict=slab_dict,
            surfen_coll=surfen_coll,
            db_file=db_file,
            high_level=high_level,
            sg_params=sg_params,
            comp_params=comp_params,
            add_full_relax=add_full_relax,
        )
        if wf:
            wf_list.append(wf)

    return wf_list, inputs_list


def surface_energy_swf_from_slab_dict(
        mpid,
        slab_dict,
        surfen_coll,
        db_file,
        high_level,
        sg_params,
        comp_params,
        add_full_relax,
):
    functional = comp_params["functional"]
    metadata = {"mpid": mpid, "functional": functional}
    nav = VaspDB(db_file=db_file, high_level=high_level)
    uid = slab_dict["uid"]
    slab = slab_dict["struct"]
    result = nav.find_data(collection=surfen_coll, fltr={"uid": uid})
    if result:
        surface_energy = result.get("surface_energy")
        if surface_energy:
            print(f"Surface energy for {mpid}-{uid} already in DB.")
            logger.info(f"Surface energy for {mpid}-{uid} already in DB.")
            return

    hkl = slab_dict["slab_params"]["hkl"]
    uid_short = uid[:4]
    fltr = {"uid": uid}
    metadata.update({"hkl": hkl, "formula": slab.formula})
    put_surfen_inputs_into_db(
        slab_dict=slab_dict,
        sg_params=sg_params,
        comp_params=comp_params,
        fltr=fltr,
        coll=surfen_coll,
        db_file=db_file,
        high_level=high_level,
        metadata=metadata,
    )

    fw1 = Firework(
        RelaxSurfaceEnergyInputs(
            slab_dict=slab_dict,
            fltr=fltr,
            coll=surfen_coll,
            comp_params=comp_params,
            db_file=db_file,
            high_level=high_level,
        ),
        name=f"Generate and relax surface energy inputs for {mpid}-{hkl}-{uid_short} with {functional}.",
    )

    fw2 = Firework(
        WriteSurfaceEnergiesToDB(
            slab_dict=slab_dict,
            fltr=fltr,
            coll=surfen_coll,
            db_file=db_file,
            high_level=high_level,
        ),
        parents=fw1,
        name=f"Calculate the surface energies for {mpid}-{hkl}-{uid_short} with {functional} and put into DB.",
    )

    fw3 = Firework(GenerateOptimadeEntry(fltr=fltr,
                                         surfen_coll=surfen_coll,
                                         db_file=db_file,
                                         high_level=high_level),
                   parents=fw2,
                   name=f"Generate OPTIMADE entry for {mpid}-{hkl}-{uid_short} with {functional} and put into DB.",
                   )

    if add_full_relax:
        fw4 = Firework(
            FullyRelaxSlab(
                input_dict=slab_dict,
                fltr=fltr,
                coll=surfen_coll,
                comp_params=comp_params,
                db_file=db_file,
                high_level=high_level,
            ),
            parents=fw2,
            name=f"Fully relax the slab for {mpid}-{hkl}-{uid_short} with {functional}",
        )

        wf = Workflow(
            fireworks=[fw1, fw2, fw3, fw4],
            name=f"Surface energy SWF for {mpid}-{hkl}-{uid_short} with {functional}.",
        )
    else:
        wf = Workflow(
            fireworks=[fw1, fw2, fw3],
            name=f"Surface energy SWF for {mpid}-{hkl}-{uid_short} with {functional}.",
        )

    return wf
