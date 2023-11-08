from fireworks import LaunchPad

from htflow_utils.misc_tools import load_defaults
from surfflow.workflows.main import surface_energy_wf

defaults = load_defaults("surfflow")

sg_filter = defaults["sg_filter"]
sg_params = defaults["sg_params"]
comp_params = defaults["comp_params"]
database_params = defaults["database_params"]

sg_params["miller"] = [(0, 0, 1)]

mpid = "mp-149"

wf = surface_energy_wf(mpid=mpid,
                       sg_params=sg_params,
                       sg_filter=sg_filter,
                       comp_params=comp_params,
                       database_params=database_params)

lpad = LaunchPad.auto_load()

lpad.add_wf(wf)
