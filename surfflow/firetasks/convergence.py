import warnings
from types import SimpleNamespace

import numpy as np
from fireworks import explicit_serialize, FiretaskBase, FWAction, Firework, Workflow
from scipy.stats import linregress

from surfflow.utils.convergence import is_list_converged
from surfflow.utils.convergence import (
    run_conv_calc_and_move_results,
    generate_slab_for_conv,
)
from surfflow.utils.db_tools import VaspDB, get_entry_by_loc
from surfflow.utils.misc_tools import check_input


@explicit_serialize
class RunConvCalcMoveResults(FiretaskBase):
    """Firetask that runs VASP and moves the results to the database."""

    _fw_name = "Run VASP and move the results."
    required_params = [
        "conv_coll",
        "comp_params",
        "fltr",
        "slab",
        "slab_thick",
    ]
    optional_params = ["db_file", "high_level"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        wf = run_conv_calc_and_move_results(**inp)
        return FWAction(detours=wf)


@explicit_serialize
class SurfenConvergence(FiretaskBase):
    """Workflow that performs the convergence calculations for the surface energies.

    Parameters
    ----------
    comp_params : dict
        Dictionary containing the computational parameters to use for the calculations.
    thick_start : float
        Thickness of the slab to start the convergence calculations with.
    thick_end : float
        Thickness of the slab to end the convergence calculations with.
    mpid : str
        Materials Project ID of the structure to calculate the surface energies for.
    sg_params : dict
        Dictionary containing the surface generation parameters.
    tol : float
        Tolerance to use when comparing the surface energies.
    Returns
    -------
        Workflow containing the Fireworks that will perform the VASP calculations.
    """

    _fw_name = "Perform a series of calculations for surface energy convergence."
    required_params = ["mpid", "miller", "shift"]
    optional_params = [
        "comp_params",
        "sg_params",
        "conv_params",
        "db_file",
        "high_level",
        "functional",
        "bulk_coll",
        "conv_coll",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        inp = SimpleNamespace(**inp)
        mpid = inp.mpid
        db_file = inp.db_file
        high_level = inp.high_level
        sg_params = inp.sg_params
        comp_params = inp.comp_params
        conv_params = inp.conv_params
        functional = comp_params["functional"]
        bulk_coll = inp.bulk_coll
        conv_coll = inp.conv_coll

        miller = inp.miller

        nav = VaspDB(db_file=db_file, high_level=high_level)
        millerstr = "".join([str(i) for i in miller])
        fltr = {
            "mpid": mpid,
            "functional": functional,
            "hkl": millerstr,
            "shift": inp.shift,
        }

        en_dict_loc = ["energies_by_layer"]
        en_dict = get_entry_by_loc(nav=nav, fltr=fltr, coll=conv_coll, loc=en_dict_loc)
        en_dict = en_dict if en_dict else {}
        surf_en_arr = [v["surf_en"] for v in en_dict.values()]
        surf_en_arr = [v for v in surf_en_arr if not np.isnan(v)]
        bulk_en_arr = [v["bulk_en"] for v in en_dict.values()]
        bulk_en_arr = [v for v in bulk_en_arr if not np.isnan(v)]

        surf_en_conv = is_list_converged(
            list_to_check=surf_en_arr,
            nelm=conv_params["nelm"],
            threshold=conv_params["conv_threshold"],
        )

        if surf_en_conv:
            thick_arr = list(en_dict.keys())
            conv_thick = thick_arr[-1]
            nav.update_data(
                collection=conv_coll,
                fltr=fltr,
                new_values={
                    "$set": {
                        "bulk_en": bulk_en_arr[-1],
                        "surf_en": surf_en_arr[-1],
                        "conv_thick": conv_thick,
                        "converged": True,
                    }
                },
                upsert=True,
            )

            return FWAction(update_spec=fw_spec)

        try:
            slab, slab_thick = generate_slab_for_conv(
                mpid=mpid,
                miller=miller,
                shift=inp.shift,
                conv_params=conv_params,
                db_file=db_file,
                en_dict=en_dict,
                high_level=high_level,
                sg_params=sg_params,
                bulk_coll=bulk_coll,
            )
        except TypeError:
            warnings.warn("Convergence can not be achieved with the given parameters.")
            nav.update_data(
                collection=conv_coll,
                fltr=fltr,
                new_values={"$set": {"converged": False}},
                upsert=True,
            )
            return FWAction(update_spec=fw_spec)

        fw1 = Firework(
            RunConvCalcMoveResults(
                comp_params=comp_params,
                fltr=fltr,
                conv_coll=conv_coll,
                slab=slab,
                db_file=db_file,
                high_level=high_level,
            )
        )
        fw2 = Firework(
            UpdateSurfenConvData(
                mpid=mpid,
                miller=miller,
                slab=slab,
                conv_coll=conv_coll,
                db_file=db_file,
                high_level=high_level,
            )
        )
        fw3 = Firework(
            SurfenConvergence(
                mpid=mpid,
                comp_params=comp_params,
                conv_params=conv_params,
                db_file=db_file,
                high_level=high_level,
                functional=functional,
                bulk_coll=inp.bulk_coll,
                conv_coll=inp.conv_coll,
                sg_params=sg_params,
            )
        )

        wf = Workflow(
            [fw1, fw2, fw3], links_dict={fw1.fw_id: [fw2.fw_id], fw2.fw_id: [fw3.fw_id]}
        )
        return FWAction(detours=wf, update_spec=fw_spec)


@explicit_serialize
class UpdateSurfenConvData(FiretaskBase):
    _fw_name = "Update the slab energies in the database"
    required_params = ["mpid", "slab"]
    optional_params = ["conv_coll", "db_file", "high_level"]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        self.update_surfen_conv_data(**inp)
        return FWAction(update_spec=fw_spec)

    @staticmethod
    def update_surfen_conv_data(mpid, slab, conv_coll, db_file, high_level):
        nav = VaspDB(db_file=db_file, high_level=high_level)
        fltr = {"mpid": mpid}
        conv_data_loc = ["conv_data"]
        conv_data = get_entry_by_loc(
            nav=nav, fltr=fltr, coll=conv_coll, loc=conv_data_loc
        )
        conv_data = conv_data if conv_data else {}
        slab_en_dict = {
            k: {
                "energy": v["output"]["energy"],
                "sites": len(v["output"]["structure"]["sites"]),
            }
            for k, v in conv_data.items()
        }
        area = slab.surface_area
        if slab_en_dict:
            thick_arr = list(slab_en_dict.keys())
            slab_en_arr = [v["energy"] for v in slab_en_dict.values()]
            sites_arr = [v["sites"] for v in slab_en_dict.values()]

            # if len(slab_en_arr) > 2:
            #     slab_en_arr = slab_en_arr[-3:]
            #     sites_arr = sites_arr[-3:]

            linreg = linregress(sites_arr, slab_en_arr)
            slope = linreg.slope

            en_dict_loc = ["energies_by_layer"]
            en_dict = get_entry_by_loc(
                nav=nav, fltr=fltr, coll=conv_coll, loc=en_dict_loc
            )
            en_dict = en_dict if en_dict else {}
            if en_dict:
                for index, (k, v) in enumerate(en_dict.items()):
                    try:
                        v["surf_en"] = (
                            16.02176565
                            * (slab_en_arr[index] - slope * sites_arr[index])
                            / (2 * area)
                        )
                    except TypeError:
                        print("TypeError IN UpdateSurfenConvData")
                        continue

            en_dict[thick_arr[-1]] = {
                "bulk_en": slope,
                "surf_en": 16.02176565
                * (slab_en_arr[-1] - slope * sites_arr[-1])
                / (2 * area),
            }
            nav.update_data(
                collection=conv_coll,
                fltr=fltr,
                new_values={"$set": {"energies_by_layer": en_dict}},
                upsert=True,
            )
