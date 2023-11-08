from fireworks import explicit_serialize, FiretaskBase, FWAction, Workflow, Firework

from surfflow.utils.misc_tools import check_input
from surfflow.firetasks.convergence import SurfenConvergence
from surfflow.utils.db_tools import VaspDB, get_entry_by_loc
from surfflow.workflows.subworkflows import surface_energy_swf


@explicit_serialize
class StartSurfaceEnergy(FiretaskBase):
    """Starts a detour workflow for calculating surface energies.

    :param mpid: Unique MaterialsProject ID of the material.
    :type mpid: str

    :param comp_params: Computational parameters for the VASP calculations.
    :type comp_params: dict, optional

    :param sg_params: Parameters used in slab generation.
    :type sg_params: dict, optional

    :param sg_filter: Parameters for filtering generated slabs depending on their BVS values.
    :type sg_filter: dict, optional

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    :param custom_id: Custom ID for the workflow.
    :type custom_id: str, optional

    :param bulk_coll: Name of the collection where the bulk structures are stored.
    :type bulk_coll: str, optional

    :param surfen_coll: Name of the collection where the surface energies are stored.
    :type surfen_coll: str, optional

    :param add_full_relax: If true, the slab is fully relaxed and copied to the high_level
        database, even if the symmetry of the slab does not require a full
        relaxation for the computation of the surface energy.
    :type add_full_relax: bool, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.
    """

    _fw_name = "Starts a sub-workflow that calculates surface energies as detour"
    required_params = ["mpid"]
    optional_params = [
        "comp_params",
        "sg_params",
        "sg_filter",
        "database_params",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        wfs, inputs_list = surface_energy_swf(**inp)
        uids = [slab_dict["uid"] for slab_dict in inputs_list]
        fw_spec.update({"StartSurfaceEnergy_info": {'mpid': self['mpid'], 'uids': uids}})

        return FWAction(detours=wfs, update_spec=fw_spec)


@explicit_serialize
class StartSurfenConvergence(FiretaskBase):
    """Firetask that starts the convergence workflow for the surface energies.
    =====================================================================

    :param mpid: Materials Project ID of the structure to calculate the surface energies for.
    :type mpid: str

    :param sg_params: Parameters for generating slabs.
    :type sg_params: dict

    :param sg_filter: Parameters for filtering slabs.
    :type sg_filter: dict

    :param db_file: Full path to the db.json file which holds the location and access credentials to the database. The default is 'auto'.
    :type db_file: str, optional

    :param high_level: Whether to query the results from the high level database or not. The default is True.
    :type high_level: str, optional

    :param comp_params: Parameters for the VASP calculations. The default is None.
    :type comp_params: dict, optional

    :param custom_id: Custom ID for the calculation. The default is None.
    :type custom_id: str, optional

    :param surfen_coll: Name of the collection to store the surface energies in. The default is '.surfen'.
    :type surfen_coll: str, optional

    :param bulk_coll: Name of the collection to store the bulk structures in. The default is '.bulk'.
    :type bulk_coll: str, optional

    :param add_full_relax: If true, the slab is fully relaxed and copied to the high_level database, even if the symmetry of the slab does not require a full relaxation for the computation of the surface energy.
    :type add_full_relax: bool

    :returns: FWAction that detours to a Workflow containing the Fireworks that will perform the VASP calculations.
    :rtype: FWAction
    """

    _fw_name = "Start the convergence workflow for the surface energies"
    required_params = ["mpid", "miller", "shift"]
    optional_params = [
        "comp_params",
        "sg_params",
        "conv_params",
        "db_file",
        "high_level",
        "conv_coll",
        "bulk_coll",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        converged = self.is_surfen_converged(inp)
        if converged:
            return FWAction(update_spec=fw_spec)
        else:
            wf = Workflow([Firework(SurfenConvergence(**inp))])
            return FWAction(detours=wf, update_spec=fw_spec)

    @staticmethod
    def is_surfen_converged(inp):
        mpid = inp["mpid"]
        miller = inp["sg_params"]["miller"]
        functional = inp["comp_params"]["functional"]
        nav = VaspDB(db_file=inp["db_file"], high_level=inp["high_level"])
        fltr = {"mpid": mpid, "miller": miller, "functional": functional, "shift": inp["shift"]}
        conv_loc = ["converged"]
        converged = get_entry_by_loc(
            nav=nav,
            fltr=fltr,
            coll=inp["conv_coll"],
            loc=conv_loc,
        )
        return True if converged else False
