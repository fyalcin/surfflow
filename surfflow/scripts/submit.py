import json
from argparse import Namespace

import yaml
from fireworks import LaunchPad

from htflow_utils.misc_tools import parse_miller
from surfflow.workflows.main import surface_energy_wf


def submit(args: Namespace) -> None:
    """
    Submit a workflow to the FireWorks database.

    :param args: Command line arguments.
    :type args: Namespace

    :return: None
    """
    lpad = LaunchPad.auto_load()

    if args.filename is not None:
        try:
            with open(args.filename, "r") as f:
                inp = yaml.safe_load(f)
            print(f"Loaded YAML input from {args.filename}")
        except FileNotFoundError or yaml.YAMLError:
            try:
                with open(args.filename, "r") as f:
                    inp = json.load(f)
                print(f"Loaded JSON input from {args.filename}")
            except FileNotFoundError or json.JSONDecodeError:
                raise ValueError(f"Could not load file {args.filename}")
        else:
            print(f"Loaded input from {args.filename}")
            try:
                inp["sg_params"]["miller"] = [
                    parse_miller(str(m)) for m in inp["sg_params"]["miller"]
                ]
            except KeyError:
                pass
    else:
        inp = vars(args)
        inp["sg_params"] = {} if inp["sg_params"] is None else inp["sg_params"]
        inp["sg_params"]["miller"] = inp.get("miller")
        if inp["sg_params"]["miller"] is not None:
            try:
                inp["sg_params"]["miller"] = [
                    parse_miller(m) for m in inp["sg_params"]["miller"]
                ]
            except TypeError:
                inp["sg_params"]["miller"] = [parse_miller(inp["sg_params"]["miller"])]
    wf = surface_energy_wf(**inp)
    lpad.add_wf(wf)
