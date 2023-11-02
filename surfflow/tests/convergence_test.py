from fireworks import Firework, Workflow
from fireworks import LaunchPad

from surfflow.firetasks.start_swfs import StartSurfenConvergence

lpad = LaunchPad.auto_load()
mpids = ["mp-149"]
shift = 0.125
miller = (1, 0, 0)
conv_params = {"nelm": 5, "thick_start": 6, "thick_end": 60}
sg_params = {
    "min_thick_A": 0,
    "match_ouc_lattice": False,
    "calculate_bonds": False,
}

for mpid in mpids:
    wf = Workflow.from_Firework(
        Firework(
            StartSurfenConvergence(
                mpid=mpid,
                miller=miller,
                shift=shift,
                sg_params=sg_params,
                conv_params=conv_params,
            )
        )
    )
    lpad.add_wf(wf)
