import copy
from fireworks import Workflow, Firework
from types import SimpleNamespace

from surfflow.firetasks.start_swfs import StartSurfaceEnergy
from surfflow.firetasks.surfen_tools import GenerateWulffShape
from surfflow.utils.db_tools import VaspDB
from surfflow.utils.misc_tools import check_input
from surfflow.utils.structure_manipulation import get_conv_bulk_from_db


def surface_energy_wf(mpid: str, **kwargs) -> Workflow:
    """Create a workflow for calculating surface energies.

    Args:
        mpid (str): Materials Project ID.
        **kwargs (dict): Keyword arguments for the workflow.

    Returns:
        Workflow: A workflow for calculating surface energies.
    """
    # Parse the input
    input_params_dict = check_input(
        kwargs,
        [
            "comp_params",
            "sg_params",
            "sg_filter",
            "database_params",
        ],
    )

    # Get the computational parameters
    database_params = input_params_dict["database_params"]
    db_file = database_params["db_file"]
    high_level = database_params["high_level"]
    bulk_coll = database_params["bulk_coll"]

    nav = VaspDB(db_file=db_file, high_level=high_level)
    conv_bulk = get_conv_bulk_from_db(
        mpid=mpid, coll=bulk_coll, db_file=db_file, high_level=high_level
    )

    if kwargs.get("comp_params_loc") == "bulk_db":
        print(
            f'comp_params_loc supplied by the user as "bulk_db". Querying the DB for computational params...'
        )
        bulk_data = nav.find_data(collection=bulk_coll, fltr={"mpid": mpid})
        comp_params_db = bulk_data.get("comp_parameters")
        try:
            input_params_dict["comp_params"].update(comp_params_db)
        except TypeError:
            print(
                "No computational parameters in the database! Continuing with the default comp_params."
            )

    sg_params = input_params_dict["sg_params"]
    miller = sg_params.get("miller")
    max_index = sg_params.get("max_index")
    name_suffix = ""
    if miller and max_index:
        raise ValueError(
            "Both miller and max_index are provided! Please provide only one of them."
        )
    elif miller and not max_index:
        millerstr = "-".join(
            [a for a in ["".join([str(i) for i in j]) for j in miller]]
        )
        name_suffix = f"({millerstr})"
    elif not miller and max_index:
        name_suffix = f"MMI {max_index}"
    elif not miller and not max_index:
        raise ValueError(
            "Neither miller nor max_index is provided! Please provide one of them."
        )

    custom_id = kwargs.get("custom_id", None)
    input_params_dict["custom_id"] = custom_id

    req_params, opt_params = (
        StartSurfaceEnergy.required_params,
        StartSurfaceEnergy.optional_params,
    )
    input_params_dict = {
        k: v for k, v in input_params_dict.items() if k in req_params + opt_params
    }
    formula = conv_bulk.composition.reduced_formula
    wf_list = []
    try:
        input_params_dict["sg_params"]["slab_thick"][0]
    except TypeError:
        surface_energy = Firework(
            StartSurfaceEnergy(mpid=mpid, **input_params_dict),
            name=f"Start surface energy SWF for {formula} ({mpid}) - {name_suffix}",
        )
        wf_list.append(surface_energy)
    else:
        slab_thick_list = sorted(input_params_dict["sg_params"]["slab_thick"])
        for slab_thick in slab_thick_list:
            tmp = copy.deepcopy(input_params_dict)
            tmp["sg_params"]["slab_thick"] = slab_thick
            tmp["custom_id"] = f"{custom_id}-{slab_thick}" if custom_id else None
            surface_energy = Firework(
                StartSurfaceEnergy(mpid=mpid, **tmp),
                name=f"Start surface energy SWF for {formula} ({mpid}) - {name_suffix} "
                f"- {slab_thick} layers",
            )
            wf_list.append(surface_energy)

    wf = Workflow(
        wf_list,
        name=f"Surface energy workflow for {formula} ({mpid}) - {name_suffix}",
    )
    return wf


def wulff_shape_wf(mpid, **kwargs) -> Workflow:
    input_params_dict = check_input(
        kwargs,
        [
            "comp_params",
            "high_level",
            "db_file",
            "surfen_coll",
            "bulk_coll",
            "wulff_coll",
        ],
    )

    # Create a namespace for the input parameters
    inp = SimpleNamespace(**input_params_dict)
    functional = inp.comp_params.get("functional")
    wulff_shape = Firework(
        GenerateWulffShape(
            mpid=mpid,
            surfen_coll=inp.surfen_coll,
            bulk_coll=inp.bulk_coll,
            wulff_coll=inp.wulff_coll,
            db_file=inp.db_file,
            high_level=inp.high_level,
        ),
        name=f"Start Wulff shape SWF for {mpid}",
        spec={"_allow_fizzled_parents": True},
    )
    wf = Workflow(
        [wulff_shape],
        name=f"Wulff shape workflow for {mpid}",
    )
    return wf
