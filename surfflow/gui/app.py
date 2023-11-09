import os
from collections import defaultdict
from datetime import datetime
from time import perf_counter

import numpy as np
import pandas as pd
import streamlit as st
from atomate.vasp.powerups import add_modify_incar
from fireworks import LaunchPad
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from htflow_utils.db_tools import VaspDB
from htflow_utils.vasp_tools import (
    get_custom_vasp_relax_settings,
    get_custom_vasp_static_settings,
)
from htflow_utils.workflows import get_calc_wf
from surfflow.workflows.main import surface_energy_wf, wulff_shape_wf

st.set_page_config(
    layout="wide",
    page_title="Surfgen Control Center",
    page_icon="curly_loop",
    initial_sidebar_state="expanded",
)
today = datetime.today().strftime("%Y-%m-%d")

start = perf_counter()

bucket_name = "streamlit-general"
comment_file_path = "comments.txt"


def update_fw(fw_id, update_dict):
    """
    Update a Firework in the database.
    Parameters
    ----------
    id_list : [int]
        List of fw_ids
    update : dict
        Dictionary of VASP INCAR flags and corresponding values.
    rerun : bool
        Selects if the fireworks only get updated or also set to rerun.

    Returns
    -------
    None.

    """
    fw = lpad.get_fw_by_id(fw_id)
    try:
        metagga = fw.spec["_tasks"][1]["vasp_input_set_params"]["user_incar_settings"][
            "METAGGA"
        ]
    except KeyError:
        metagga = False
    if metagga:
        for k, v in update_dict.items():
            lpad.update_spec(
                [fw_id],
                {"_tasks.1.vasp_input_set_params.user_incar_settings.{}".format(k): v},
            )
    else:
        for k, v in update_dict.items():
            lpad.update_spec(
                [fw_id], {"_tasks.0.vasp_input_set.user_incar_settings.{}".format(k): v}
            )


def parse_miller(miller: str) -> tuple:
    """
    Parse a miller index from a string to a tuple.
    Parameters
    ----------
    miller : str
        Miller index in string format.

    Returns
    -------
    miller : tuple
    """
    emp = []
    flag = False
    for i, m in enumerate(miller):
        try:
            tmp = (1 - 2 * int(flag)) * int(m)
            emp.append(tmp)
        except ValueError:
            flag = True
        else:
            flag = False
    return tuple(emp)


default_inp = {
    "comp_params": {
        "use_vdw": False,
        "functional": "PBE",
        "use_spin": True,
        "encut": 400,
        "k_dens": 5.0,
        "epsilon": 1.0,
    },
    "sg_params": {
        "symmetrize": False,
        "slab_thick": 8,
        "vac_thick": 30.0,
        "primitive": True,
        "lll_reduce": True,
        "tol": 0.05,
        "max_normal_search": None,
        "resize": True,
        "preserve_terminations": True,
        "min_thick_A": 10.0,
        "minimize_structures": False,
        "match_ouc_lattice": True,
        "calculate_bonds": True,
        "max_nsites": 100,
        "nn_method": "all",
        "weight": "BVS",
        "filter_polar": True,
    },
    "sg_filter": {"method": "all"},
    "database_params": {
        "db_file": "auto",
        "high_level": True,
        "bulk_coll": "bulk_data",
        "surfen_coll": "surfen_data",
        "conv_coll": "conv_data",
        "wulff_coll": "wulff_data",
    }}

options_dict = {
    "comp_params": {
        "functional": ["PBE", "LDA"],
        "use_vdw": [True, False],
        "use_spin": [True, False],
        "encut": 400,
        "k_dens": 5.0,
        "epsilon": 1.0,
    },
    "sg_params": {
        "symmetrize": [True, False],
        "slab_thick": 8,
        "vac_thick": 30,
        "primitive": [True, False],
        "lll_reduce": [True, False],
        "tol": 0.05,
        "resize": [True, False],
        "preserve_terminations": [True, False],
        "min_thick_A": 10.0,
        "minimize_structures": [True, False],
        "match_ouc_lattice": [True, False],
        "calculate_bonds": [True, False],
        "filter_polar": [True, False],
        "max_nsites": 100,
    },
}


def st_set_fireworks_db():
    fireworks_db = st.text_input("Fireworks database", "fireworks")
    if fireworks_db:
        st.session_state.fireworks_db = fireworks_db
    else:
        st.session_state.fireworks_db = "fireworks"


st_set_fireworks_db()

lpad = LaunchPad.auto_load()


def streamlit_refresh_wfs(num_wfs=None):
    if num_wfs:
        st.session_state.wfs_to_show = num_wfs
    else:
        num_wfs = st.session_state.get("wfs_to_show", 10)

    if num_wfs == -1:
        wfs = list(
            lpad.workflows.find({}, {"name": 1, "state": 1, "nodes": 1, "fw_states": 1})
        )[::-1]
    else:
        wfs = list(
            lpad.workflows.find({}, {"name": 1, "state": 1, "nodes": 1, "fw_states": 1})
        )[::-1][:num_wfs]

    wfs_dict = {wf["nodes"][0]: wf for wf in wfs}
    main_fws = {
        fw_id: lpad.fireworks.find_one({"fw_id": fw_id}, {"spec._priority": 1})
        for fw_id in wfs_dict
    }
    st.session_state.wfs = wfs_dict
    st.session_state.main_fws = main_fws


streamlit_refresh_wfs()

lostruns, reservations, _, _ = st.columns(4)
with lostruns:
    if st.button("Check for and rerun lost runs", key="lostruns"):
        lpad.detect_lostruns(rerun=True)
        st.success("Done!")

with reservations:
    if st.button("Check for and rerun unreserved runs"):
        lpad.detect_unreserved(rerun=True)
        st.success("Done!")

inp = {}

sb = st.sidebar
sb.title("Surfgen Control Center")


def grab_input_from_text_area(label, placeholder, key):
    txt = st.text_area(label=label, value="", placeholder=placeholder, key=key)
    custom_input_dict = {}
    if txt:
        txt = txt.splitlines()
        for line in txt:
            line = line.split("=")
            line = [el.strip() for el in line]
            custom_input_dict[line[0]] = line[-1]
    return custom_input_dict


with sb.expander("**Submit a new surface energy workflow**"):
    st.markdown("## MPID")
    mpids = st.text_input(
        label="mpid", value="mp-66", key="mpid", label_visibility="collapsed"
    ).split()

    st.markdown("## Miller Index Selection")

    surfaces_option = st.radio(
        label="Surface selection",
        index=1,
        options=["max_index", "miller"],
        key="surfaces_option",
        label_visibility="collapsed",
    )

    if surfaces_option == "max_index":
        max_index = int(
            st.selectbox(
                label="Maximum Miller index", options=np.arange(1, 4), key="max_index"
            )
        )
    elif surfaces_option == "miller":
        miller = st.text_input(
            label="Miller indices (separated by spaces)",
            key="miller",
            placeholder="11-1 110 100",
        )

    st.markdown("## Surface filtering (Broken bonds)")
    sg_filter = st.radio(
        label="Surface filtering (broken bond)",
        options=["all (default)", "N minimum", "N minimum per hkl"],
        key="sg_filter",
        label_visibility="collapsed",
    )

    bvs_options_map = {
        "all (default)": "all",
        "N minimum": "bvs_min_N",
        "N minimum per hkl": "bvs_min_N_hkl",
    }

    inp["sg_filter"] = {"method": bvs_options_map[sg_filter]}
    if sg_filter != "all (default)":
        sg_filter_n = st.number_input(label=sg_filter, value=1, key="sg_filter_n")
        inp["sg_filter"]["bvs_param"] = sg_filter_n

    for key, value in options_dict.items():
        try:
            value.items()
        except AttributeError:
            continue
        inp[key] = {}
        st.markdown(f"## {key}")
        for k, v in value.items():
            if isinstance(v, list):
                user_inp = st.selectbox(
                    label=f"{k} (default: {default_inp[key][k]})",
                    options=v,
                    index=v.index(default_inp[key][k]),
                    key=f"surface_wf_{k}",
                )
            else:
                inp_type = type(v)
                user_inp = inp_type(
                    st.number_input(
                        label=f"{k} (default: {default_inp[key][k]})",
                        value=default_inp[key][k],
                        key=f"surface_wf_{k}",
                    )
                )
            inp[key][k] = user_inp

    if surfaces_option == "max_index":
        inp["sg_params"]["max_index"] = max_index
    elif surfaces_option == "miller":
        inp["sg_params"]["miller"] = [parse_miller(m) for m in miller.split()]

    fake_vasp = st.checkbox(label="Fake VASP", key="fake_vasp")

    inp["fake_vasp"] = fake_vasp

    st.markdown("## Submit workflow")
    for key, value in default_inp.items():
        if isinstance(value, dict):
            for key2, value2 in value.items():
                if not inp.get(key):
                    inp[key] = {}

                if inp[key].get(key2) is None:
                    inp[key][key2] = default_inp[key][key2]
        else:
            if not inp.get(key):
                inp[key] = default_inp[key]

    if st.checkbox(label="Show input", key="show_input"):
        st.write(inp)

    calculate_wulff = st.checkbox(label="Calculate Wulff shape", key="wulff")
    st.checkbox(label="Set ISMEAR to 0", key="set_ismear_0", value=False)
    if st.button(label="Submit", key="submit_surfgen"):
        if surfaces_option == "miller" and not miller:
            st.error("Please enter miller indices!")
        else:
            for mpid in mpids:
                inp["mpid"] = mpid
                wf_surfen = surface_energy_wf(**inp)
                if st.session_state.get("set_ismear_0"):
                    wf_surfen = add_modify_incar(
                        wf_surfen, modify_incar_params={"incar_update": {"ISMEAR": 0}}
                    )
                if calculate_wulff:
                    wf_wulff = wulff_shape_wf(**inp)
                    wf_surfen.append_wf(wf_wulff, wf_surfen.leaf_fw_ids)
                try:
                    lpad.add_wf(wf_surfen)
                except ValueError:
                    st.error(
                        f"Workflow could not be submitted! Check your launchpad config!"
                    )
                st.success(f"Workflow for {mpid} submitted!")

with sb.expander("**Submit a custom workflow**"):
    structure = st.file_uploader(label="**Structure file (POSCAR):**", key="structure")
    if structure:
        try:
            struct = Structure.from_str(
                structure.getvalue().decode("utf-8"), fmt="poscar"
            )
        except (IndexError, ValueError):
            try:
                struct = Structure.from_str(
                    structure.getvalue().decode("utf-8"), fmt="cif"
                )
            except (IndexError, ValueError):
                struct = None
                st.error("Invalid structure file!")
            else:
                st.success("Structure loaded from cif file successfully!")
        else:
            st.success("Structure loaded from POSCAR file successfully!")
        if struct:
            # write in markdown "structure loaded successfully" in green
            calc_type = st.selectbox(
                label="**Calculation type**",
                options=["relaxation", "static"],
                key="calc_type",
            )

            input_set = MPRelaxSet if calc_type == "relaxation" else MPStaticSet

            allowed_relax_types = [
                "bulk_full_relax",
                "bulk_vol_relax",
                "bulk_pos_relax",
                "bulk_shape_relax",
                "bulk_pos_shape_relax",
                "slab_shape_relax",
                "slab_pos_relax",
                "interface_shape_relax",
                "interface_pos_relax",
                "interface_z_relax",
            ]

            allowed_static_types = [
                "bulk_from_scratch",
                "bulk_follow_up",
                "bulk_nscf",
                "slab_from_scratch",
                "slab_follow_up",
                "slab_nscf",
                "bulk_epsilon_from_scratch",
                "bulk_epsilon_follow_up",
            ]

            calc_subtype = st.selectbox(
                label="**Calculation subtype**",
                options=allowed_relax_types
                if calc_type == "relaxation"
                else allowed_static_types,
                key="calc_subtype",
            )

            vis = (
                get_custom_vasp_relax_settings
                if calc_type == "relaxation"
                else get_custom_vasp_static_settings
            )

            st.markdown("## Computational Parameters")
            user_comp_params = {}
            for k, v in options_dict["comp_params"].items():
                if isinstance(v, list):
                    user_inp = st.selectbox(
                        label=f'{k} (default: {default_inp["comp_params"][k]})',
                        options=v,
                        index=v.index(default_inp["comp_params"][k]),
                        key=f"custom_wf_{k}",
                    )
                else:
                    user_inp = float(
                        st.text_input(
                            label=f'{k} (default: {default_inp["comp_params"][k]})',
                            value=default_inp["comp_params"][k],
                            key=f"custom_wf_{k}",
                        )
                    )
                user_comp_params[k] = user_inp

            custom_uis = grab_input_from_text_area(
                label="**Custom VASP input parameters:**",
                placeholder="LREAL = .FALSE.\nPREC = Accurate",
                key="vasp_input_custom",
            )

            if st.button(
                    label="Reset custom VASP input parameters", key="reset_custom_uis"
            ):
                custom_uis = {}

            vis = vis(struct, user_comp_params, calc_subtype, custom_uis=custom_uis)
            show_incar = st.checkbox(label="Show INCAR", key="show_incar")
            if show_incar:
                st.write(vis.incar)

            add_static = st.checkbox(
                label="Add static calculation after relaxation",
                value=True,
                key="add_static",
            )

            task_label = st.text_input(
                label="**Custom tag:**",
                value="",
                key="custom_tag",
                placeholder="Task label for calculation (leave blank for default)",
            )

            if not task_label:
                task_label = (
                    f"{struct.composition.reduced_formula}_{calc_type}_{calc_subtype}"
                )

            wf = get_calc_wf(struct, vis, task_label, add_static)
            wf.name = f"{task_label}_custom_wf"

            st.markdown("## Submit workflow")

            if st.button(label="Submit", key="submit_custom"):
                try:
                    lpad.add_wf(wf)
                except ValueError:
                    st.error(
                        f"Workflow could not be submitted! Check your launchpad config!"
                    )
                else:
                    st.success(f"Workflow submitted!")

lpad_tab, control_tab, results_tab = st.tabs(
    ["Launchpad", "Fireworks Control", "Results"]
)

state_color_template = {
    "FIZZLED": "red",
    "RUNNING": "green",
    "COMPLETED": "green",
    "DEFUSED": "orange",
    "READY": "green",
    "WAITING": "orange",
    "RESERVED": "orange",
    "ARCHIVED": "violet",
    "PAUSED": "blue",
}

state_color = {k: f":{v}[{k}]" for k, v in state_color_template.items()}


def add_to_selected_wf(fw_id):
    if "selected_wfs" not in st.session_state:
        st.session_state.selected_wfs = []
    if fw_id not in st.session_state.selected_wfs:
        st.session_state.selected_wfs.append(fw_id)


def delete_selected_wfs():
    selected_wfs = st.session_state.get("selected_wfs", [])
    for fw_id in selected_wfs:
        lpad.delete_wf(fw_id=fw_id)
        st.success(f"Workflow {fw_id} deleted!")
    st.session_state.selected_wfs = []


def st_render_launchpad():
    wfs = st.session_state.wfs
    main_fws = st.session_state.main_fws
    wf_id_map = {}
    fw_states_grouped_all = {}

    col_widths = [2.5, 0.2, 0.5, 0.6, 0.2, 0.5, 1.5, 0.5]
    a1, _, b1, c1, d1, _, e1, f1 = st.columns(col_widths)
    with a1:
        st.markdown("### Workflow (most recent first)")
    with b1:
        st.markdown("### ID")
    with c1:
        st.markdown("### State")
    with d1:
        st.markdown("### Priority")
    with e1:
        st.markdown("### Control Panel")
    with f1:
        st.button(
            label="Delete", key="delete_selected_wf", on_click=delete_selected_wfs
        )

    for main_fw_id, wf in wfs.items():
        name = wf["name"]
        state = wf["state"]
        fw_states_grouped = defaultdict(list)
        for fw_id, state in wf["fw_states"].items():
            fw_states_grouped[state].append(int(fw_id))
            wf_id_map[int(fw_id)] = main_fw_id
        fw_states_grouped_all[main_fw_id] = fw_states_grouped
        a, _, b, c, d, _, e, f = st.columns(col_widths)

        with a:
            with st.expander(f"**{name}**"):
                fw_states_grouped_tabs = st.tabs(list(fw_states_grouped.keys()))
                for state, tab in zip(
                        list(fw_states_grouped.keys()), fw_states_grouped_tabs
                ):
                    with tab:
                        st.write(fw_states_grouped[state])
        with b:
            st.markdown(f"**{main_fw_id}**")
        with c:
            state_md = state_color[state]
            if state != "RUNNING" and "RUNNING" in fw_states_grouped:
                state_md += " - :green[RUN]"
            if state != "FIZZLED" and "FIZZLED" in fw_states_grouped:
                state_md += " - :red[FIZ]"
            st.markdown(f"**{state_md}**")

        with d:
            priority = main_fws[main_fw_id]["spec"].get("_priority")
            if priority is not None:
                st.markdown(f"**:green[{priority}]**")

        with e:
            with st.expander("Show control panel"):
                if state == "PAUSED":
                    if st.button(label="Resume", key=f"resume_wf_{main_fw_id}"):
                        lpad.resume_wf(fw_id=main_fw_id)
                        st.success("Workflow resumed!")
                elif state == "RUNNING":
                    if st.button(label="Pause", key=f"pause_wf_{main_fw_id}"):
                        lpad.pause_wf(fw_id=main_fw_id)
                        st.success("Workflow paused!")

                    if st.button(label="Defuse", key=f"defuse_wf_{main_fw_id}"):
                        lpad.defuse_wf(fw_id=main_fw_id)
                        st.success("Workflow defused!")

                elif state == "DEFUSED":
                    if st.button(label="Reignite", key=f"reignite_wf_{main_fw_id}"):
                        lpad.reignite_wf(fw_id=main_fw_id)
                        st.success("Workflow reignited!")

                fizzled_fw_ids = fw_states_grouped.get("FIZZLED")
                if fizzled_fw_ids:
                    custom_uis_wf = grab_input_from_text_area(
                        label="Custom Input Params",
                        key=f"custom_uis_{main_fw_id}",
                        placeholder="LREAL = .FALSE.\nPREC = Accurate",
                    )
                    if st.button(
                            label="Update and Rerun fizzled", key=f"rerun_wf_{main_fw_id}"
                    ):
                        for i in fizzled_fw_ids:
                            lpad.rerun_fw(fw_id=i)
                            if custom_uis_wf:
                                update_fw(i, custom_uis_wf)
                            st.success(f"Rerun Firework with id {i}")
                        st.success("Workflow rerun!")
                    st.markdown("---")
                if st.button(label="Delete", key=f"delete_wf_{main_fw_id}"):
                    lpad.delete_wf(fw_id=main_fw_id)
                    st.success("Workflow deleted!")

                all_fw_ids = wf["nodes"]
                if state != "COMPLETED":
                    priority = st.text_input(
                        label="**Priority (higher executes first)**",
                        value="0",
                        key=f"set_priority_{main_fw_id}",
                    )
                    try:
                        priority = int(priority)
                    except ValueError:
                        st.error("Priority must be an integer!")
                    else:
                        if st.button(
                                label="Set priority", key=f"set_priority_wf_{main_fw_id}"
                        ):
                            for i in all_fw_ids:
                                lpad.set_priority(fw_id=i, priority=priority)
                            st.success("Priority set!")
        st.text("")
        with f:
            if st.checkbox(
                    label="Delete WF",
                    label_visibility="collapsed",
                    value=st.session_state.get(f"selected_wf_{main_fw_id}", False),
                    key=f"selected_wf_{main_fw_id}",
            ):
                add_to_selected_wf(main_fw_id)
    st.session_state["wf_id_map"] = wf_id_map
    st.session_state["fw_states_grouped_all"] = fw_states_grouped_all


def st_load_more_lpad():
    wfs = st.session_state.get("wfs", [])
    if st.button(label="Load 10 more", key="load_more"):
        streamlit_refresh_wfs(num_wfs=len(wfs) + 50)
        st.experimental_rerun()
    if st.button(label="Load all", key="load_all"):
        streamlit_refresh_wfs(num_wfs=-1)
        st.experimental_rerun()


with lpad_tab:
    st_render_launchpad()
    st.divider()
    st_load_more_lpad()


def text_area_to_rerun_single_fw(fw_id, widget_key_append):
    key_append = f"{fw_id}_{widget_key_append}"
    custom_uis = grab_input_from_text_area(
        label="**Update INCAR and rerun**",
        placeholder="LREAL = .FALSE.\nPREC = Accurate",
        key=f"update_incar_{fw_id}_{widget_key_append}",
    )
    if st.button(label="Rerun", key=f"rerun_{key_append}"):
        lpad.rerun_fw(fw_id=fw_id)
        st.success("Firework rerun!")
        if custom_uis:
            update_fw(fw_id=fw_id, update_dict=custom_uis)
            st.success("Firework updated!")


def print_wf(fw_id):
    wf_id_map = st.session_state.get("wf_id_map", {})
    wfs = st.session_state.get("wfs", {})
    fw_states_grouped_all = st.session_state.get("fw_states_grouped_all", {})
    main_id = wf_id_map.get(fw_id)
    wf_query = wfs.get(main_id)
    if not wf_query:
        st.error("Workflow not found!")
        return
    name = wf_query["name"]
    state = wf_query["state"]
    fws = {
        result["fw_id"]: result
        for result in lpad.fireworks.find(
            {"fw_id": {"$in": wf_query["nodes"]}}, {"fw_id": 1, "state": 1, "name": 1}
        )
    }
    st.markdown(f"### **{name} - {state_color[state]}**")
    a, b = st.columns(2)
    with a:
        st.markdown(f"**Workflow details**")
        for k, v in wf_query.items():
            with st.expander(label=k):
                st.write(v)
    with b:
        st.markdown(f"**Fireworks overview by state**")
        fw_states_grouped = fw_states_grouped_all.get(main_id, None)
        if fw_states_grouped:
            for state, fw_ids in fw_states_grouped.items():
                with st.expander(label=f"{state}"):
                    for fw_id in fw_ids:
                        fw = fws.get(fw_id, None)
                        st.write(
                            f"**:{state_color_template[state]}[{fw_id}] -"
                            f" {fw['name']}**"
                        )


def print_single_fw_controls(fw_id, widget_key_append):
    fw = lpad.fireworks.find_one(filter={"fw_id": fw_id}, projection={"state": 1})
    if not fw:
        return
    fwstate = fw["state"]
    st.markdown(f"**Firework controls**")
    key_append = f"{fw_id}_{widget_key_append}"
    if fwstate == "PAUSED":
        if st.button(label="Resume", key=f"resume_{key_append}"):
            lpad.resume_fw(fw_id=fw_id)
            st.success("Firework resumed!")
    elif fwstate == "RUNNING":
        if st.button(label="Pause", key=f"pause_{key_append}"):
            lpad.pause_fw(fw_id=fw_id)
            st.success("Firework paused!")

        if st.button(label="Defuse", key=f"defuse_{key_append}"):
            lpad.defuse_fw(fw_id=fw_id)
            st.success("Firework defused!")

    elif fwstate == "DEFUSED":
        if st.button(label="Reignite", key=f"reignite_{key_append}"):
            lpad.reignite_fw(fw_id=fw_id)
            st.success("Firework reignited!")

    elif fwstate == "FIZZLED" or fwstate == "COMPLETED":
        text_area_to_rerun_single_fw(fw_id=fw_id, widget_key_append=widget_key_append)


def print_fw(fw_id):
    fw = lpad.fireworks.find_one(filter={"fw_id": fw_id})
    if not fw:
        st.error("Firework not found!")
        return
    name = fw["name"]
    fwstate = fw["state"]
    st.markdown(
        f"**:{state_color_template[fwstate]}[{fw_id}]: "
        f"{name} - {state_color[fwstate]}**"
    )
    for k, v in fw.items():
        label = k if k != "spec" else "spec (might take a while to load)"
        with st.expander(label=label):
            st.write(v)
    print_single_fw_controls(fw_id=fw_id, widget_key_append="print_fw")


def print_launch(fw_id):
    fw = lpad.fireworks.find_one(
        filter={"fw_id": fw_id},
        projection={"launches": 1, "archived_launches": 1, "name": 1, "state": 1},
    )
    if not fw:
        st.error("Firework not found!")
        return

    fw_name = fw["name"]
    fw_state = fw["state"]
    fw_md = f"**:{state_color_template[fw_state]}[{fw_id}] - {fw_name} - {state_color[fw_state]}**"
    launch_ids = fw["launches"] + fw["archived_launches"]
    launches_found = list(
        lpad.launches.find(
            filter={"launch_id": {"$in": launch_ids}},
            projection={
                "launch_dir": 1,
                "launch_id": 1,
                "state": 1,
                "action": 1,
                "fworker": 1,
                "state_history": 1,
            },
        )
    )
    if not launches_found:
        st.error("Launches not found!")
        return
    st.markdown(fw_md)
    for launch in launches_found:
        state = launch["state"]
        md_append = (
            f" - :blue[ARCHIVED]"
            if launch["launch_id"] in fw["archived_launches"]
            else ""
        )
        st.markdown(f"**{launch['launch_id']} - {state_color[state]}{md_append}**")
        for k, v in launch.items():
            with st.expander(label=k):
                st.write(v)


def print_fw_by_state():
    st.markdown(f"## **Fireworks by state**")
    states = ["FIZZLED", "RUNNING", "DEFUSED", "READY", "RESERVED", "PAUSED"]

    tabs = st.tabs(states)
    for index, state in enumerate(states):
        tab = tabs[index]
        with tab:
            fws_query = list(
                lpad.fireworks.find(
                    filter={"state": state},
                    projection={
                        "fw_id": 1,
                        "name": 1,
                        "launches": 1,
                        "archived_launches": 1,
                    },
                )
            )

            if not fws_query:
                st.error(f"No Fireworks found in state {state}")
                continue
            a1, b1 = st.columns([2, 1])
            with a1:
                st.markdown("### **ID  -  Name**")
            with b1:
                st.markdown("### **Controls**")
            for fw in fws_query:
                fw_id = fw["fw_id"]
                fw_name = fw["name"]
                launch_ids = fw["launches"] + fw["archived_launches"]
                launches_found = list(
                    lpad.launches.find(
                        filter={"launch_id": {"$in": launch_ids}},
                        projection={
                            "launch_id": 1,
                            "state": 1,
                            "action": 1,
                            "fworker": 1,
                            "state_history": 1,
                        },
                    )
                )
                a, b = st.columns([2, 1])
                with a:
                    st.markdown(
                        f"**:{state_color_template[state]}[{fw_id}]  -  {fw_name}**"
                    )
                    st.markdown(f"**Launches**")
                    for launch in launches_found:
                        md_append = (
                            f" - ARCHIVED"
                            if launch["launch_id"] in fw["archived_launches"]
                            else ""
                        )
                        with st.expander(
                                label=f"**{launch['launch_id']} - {launch['state']}{md_append}**"
                        ):
                            st.markdown(
                                f"**:{state_color_template[state]}[{launch['launch_id']}]**"
                            )
                            launch_dict = {
                                k: launch.get(k)
                                for k in ["launch_dir", "action", "fworker"]
                            }
                            st.write(launch_dict)
                with b:
                    print_single_fw_controls(
                        fw_id=fw_id, widget_key_append="print_fw_by_state"
                    )
                if len(fws_query) > 1 and fw_id != fws_query[-1]:
                    st.markdown("---")


def print_fw_by_state_short():
    st.markdown(f"## **Fireworks by state**")
    states = ["FIZZLED", "RUNNING", "DEFUSED", "READY", "RESERVED", "PAUSED"]

    for index, state in enumerate(states):
        if st.button(label=f"Show {state} Fireworks"):
            fws_query = list(
                lpad.fireworks.find(
                    filter={"state": state},
                    projection={
                        "fw_id": 1,
                        "name": 1,
                        "launches": 1,
                        "archived_launches": 1,
                    },
                )
            )

            if not fws_query:
                st.error(f"No Fireworks found in state {state}")
                continue
            a1, b1 = st.columns([2, 1])
            with a1:
                st.markdown("### **ID  -  Name**")
            with b1:
                st.markdown("### **Controls**")
            for fw in fws_query:
                fw_id = fw["fw_id"]
                fw_name = fw["name"]
                launch_ids = fw["launches"] + fw["archived_launches"]
                launches_found = list(
                    lpad.launches.find(
                        filter={"launch_id": {"$in": launch_ids}},
                        projection={
                            "launch_id": 1,
                            "host": 1,
                            "fworker": 1,
                            "launch_dir": 1,
                            "state": 1,
                            "action": 1,
                            "state_history": 1,
                        },
                    )
                )
                a, b = st.columns([2, 1])
                with a:
                    st.markdown(
                        f"**:{state_color_template[state]}[{fw_id}]  -  {fw_name}**"
                    )
                    st.markdown(f"**Launches**")
                    for launch in launches_found:
                        md_append = (
                            f" - ARCHIVED"
                            if launch["launch_id"] in fw["archived_launches"]
                            else ""
                        )
                        with st.expander(
                                label=f"**{launch['launch_id']} - {launch['state']}{md_append}**"
                        ):
                            st.markdown(
                                f"**:{state_color_template[state]}[{launch['launch_id']}]**"
                            )
                            launch_dict = {
                                k: launch.get(k)
                                for k in ["launch_dir", "action", "fworker", "host"]
                            }
                            st.write(launch_dict)
                with b:
                    print_single_fw_controls(
                        fw_id=fw_id, widget_key_append="print_fw_by_state"
                    )
                if len(fws_query) > 1 and fw_id != fws_query[-1]:
                    st.markdown("---")


def print_query_result():
    st.markdown(f"## **Query by Firework ID**")
    query_types = ["Workflow", "Firework", "Launch"]
    tabs = st.tabs(query_types)

    for index, query_type in enumerate(query_types):
        if query_type == "Workflow":
            with tabs[index]:
                query = st.text_input(label="Firework ID belonging to workflow")
                if query:
                    try:
                        query = int(query)
                    except ValueError:
                        st.error("Please enter a single valid integer!")
                    else:
                        print_wf(query)
        elif query_type == "Firework":
            with tabs[index]:
                query = st.text_input(label="Firework IDs (space separated)")
                if query:
                    st.write(query)
                    try:
                        query = [
                            int(d)
                            for d in [a.strip(",") for a in query.split()]
                            if d.isdigit()
                        ]
                    except ValueError:
                        st.error("Please enter a single integer or list of integers!")
                    else:
                        for q in query:
                            a, b = st.columns(2)
                            with a:
                                print_fw(q)
                            with b:
                                print_launch(q)
                            if len(query) > 1 and q != query[-1]:
                                st.markdown("---")
        elif query_type == "Launch":
            with tabs[index]:
                query = st.text_input(
                    label="Firework IDs containing launches (space separated)"
                )
                if query:
                    try:
                        query = [
                            int(d)
                            for d in [a.strip(",") for a in query.split()]
                            if d.isdigit()
                        ]
                    except ValueError:
                        st.error("Please enter a single integer or list of integers!")
                    else:
                        for q in query:
                            print_launch(q)
                            if len(query) > 1 and q != query[-1]:
                                st.markdown("---")


def print_batch_fw_rerunner():
    st.markdown(f"**BATCH FIREWORK RERUNNER**")
    batch_rerun_fws = st.text_input(
        label="Enter Firework IDs to rerun (separated by space)"
    )
    custom_uis = grab_input_from_text_area(
        label="**Update INCAR and rerun**",
        placeholder="LREAL = .FALSE.\nPREC = Accurate",
        key=f"custom_input_batch_rerun_{batch_rerun_fws}",
    )
    if batch_rerun_fws:
        try:
            batch_rerun_fws = [
                int(d)
                for d in [a.strip(",") for a in batch_rerun_fws.split()]
                if d.isdigit()
            ]
        except ValueError:
            st.error("Please enter valid integers!")
            return
        else:
            st.markdown(f"**Fireworks to rerun:**")
            st.write(batch_rerun_fws)

    if st.button(label="Rerun", key=f"batch_rerun_{batch_rerun_fws}"):
        for fw_id in batch_rerun_fws:
            lpad.rerun_fw(fw_id=fw_id)
            if custom_uis:
                update_fw(fw_id=fw_id, update_dict=custom_uis)
            st.success(f"Firework with ID {fw_id} has been rerun!")


def print_batch_wf_deleter():
    st.markdown(f"**BATCH WORKFLOW DELETER**")
    batch_delete_fws = st.text_input(
        label="Enter Firework IDs to delete (separated by space)"
    )
    if batch_delete_fws:
        try:
            batch_delete_fws = [
                int(d)
                for d in [a.strip(",") for a in batch_delete_fws.split()]
                if d.isdigit()
            ]
        except ValueError:
            st.error("Please enter valid integers!")
            return
        else:
            st.markdown(f"**Fireworks to delete:**")
            st.write(batch_delete_fws)

    if st.button(label="Delete", key=f"batch_delete_{batch_delete_fws}"):
        for fw_id in batch_delete_fws:
            lpad.delete_wf(fw_id=fw_id)
            st.success(f"Firework with ID {fw_id} has been deleted!")


def print_fizzled_fw_finder(traceback_only=False):
    st.markdown(f"**FIZZLED FIREWORK FINDER**")
    fizzled_fws = list(
        lpad.fireworks.find(
            filter={"state": "FIZZLED"},
            projection={"fw_id": 1, "name": 1, "launches": 1, "archived_launches": 1},
        )
    )
    if not fizzled_fws:
        st.error(f"No Fizzled Fireworks found!")
        return

    st.markdown(f"**Found {len(fizzled_fws)} Fizzled Fireworks**")
    st.write([fw["fw_id"] for fw in fizzled_fws])
    if st.button(label="Find Fizzled Fireworks", key="find_fizzled_fws"):
        for fw in fizzled_fws:
            fw_name = fw["name"]
            fw_id = fw["fw_id"]
            fw_md = f"**:{state_color_template['FIZZLED']}[{fw_id}] - {fw_name} - {state_color['FIZZLED']}**"
            launch_ids = fw["launches"]
            launches_found = list(
                lpad.launches.find(
                    filter={"launch_id": {"$in": launch_ids}},
                    projection={
                        "launch_dir": 1,
                        "launch_id": 1,
                        "action.stored_data._message": 1,
                        "action.stored_data._exception": 1,
                    },
                )
            )
            st.markdown(fw_md)
            for launch in launches_found:
                st.markdown(f"{launch['launch_id']}")
                if traceback_only:
                    st.write(
                        launch["action"]["stored_data"]["_exception"]["_stacktrace"][
                        -200:
                        ]
                    )
                else:
                    for k, v in launch.items():
                        with st.expander(label=f"**{k}**"):
                            st.write(v)


with control_tab:
    print_query_result()

    # st.markdown("---")

    print_fw_by_state_short()

    st.markdown("---")

    print_fizzled_fw_finder(traceback_only=True)

    st.markdown("---")

    print_batch_fw_rerunner()
#


with results_tab:
    st.markdown(f"## **Results**")
    query = st.text_input(
        label="MPID or formula (space separated)", placeholder="mp-1234 Fe2O3"
    )
    if query:
        try:
            query = [q for q in query.split()]
        except ValueError:
            st.error("Please enter values separated by space!")
        else:
            cwd = os.getcwd()
            file = os.path.join(cwd, "db.json")
            db = VaspDB(db_file=file, high_level="surfgen_db_new")
            all_results = []
            for q in query:
                results = list(
                    db.find_many_data(
                        collection="PBE.surfen_data",
                        fltr={
                            "surface_energy": {"$exists": True},
                            "$or": [{"mpid": q}, {"formula": q}],
                        },
                        projection={
                            "hkl": 1,
                            "surface_energy": 1,
                            "terminations": 1,
                            "mpid": 1,
                            "formula": 1,
                            "uid": 1,
                            "comp_params": 1,
                            "_id": 0,
                        },
                    )
                )

                if results:
                    key_order = [
                        "mpid",
                        "formula",
                        "hkl",
                        "surface_energy",
                        "terminations",
                        "comp_params",
                        "uid",
                    ]
                    flattened_results = []
                    for result in results:
                        reordered_result = {k: result.get(k) for k in key_order}
                        flattened_result = {}
                        for k, v in reordered_result.items():
                            print(k, v)
                            if k in [
                                "comp_params",
                                "sg_params",
                                "slab_params",
                                "surface_energy",
                                "terminations",
                            ]:
                                for k2, v2 in v.items():
                                    flattened_result[f"{k}.{k2}"] = v2
                            else:
                                flattened_result[k] = v
                        flattened_results.append(flattened_result)
                    all_results.extend(flattened_results)

            if results:
                df = pd.DataFrame(all_results)
                for q in query:
                    if q.startswith("mp-"):
                        filtered_df = df[df["mpid"] == q]
                    else:
                        filtered_df = df[df["formula"] == q]
                    st.markdown(f"### **{q}**")
                    # filter df by formula or mpid having a value of q
                    st.dataframe(filtered_df)

end = perf_counter()
print(f"Time taken: {end - start}")
