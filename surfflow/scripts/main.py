import json
import os
import subprocess
import sys
from argparse import ArgumentParser
from typing import Dict, Optional, Sequence, Union

import yaml
from fireworks.core.launchpad import LaunchPad

from surfflow.scripts.submit import submit

DEFAULT_CONF_PATH = os.path.expanduser("~") + "/config"
DEFAULT_LPAD_YAML = "my_launchpad.yaml"

if sys.version_info < (3, 8):
    import importlib_metadata as metadata
else:
    from importlib import metadata


def set_default_config_path() -> str:
    fw_config_file = os.environ.get("FW_CONFIG_FILE")
    if fw_config_file:
        print(
            f"WARN: Existing configuration directory found at {os.path.dirname(fw_config_file)}. "
            f"You have the option to change this directory."
        )
        default_path = os.path.dirname(fw_config_file)
    else:
        default_path = DEFAULT_CONF_PATH
    return default_path


def get_user_inputs(fields: Sequence[tuple]) -> Dict[str, Union[str, int, bool, None]]:
    """
    Get user inputs for the fields specified.
    :param fields: List of tuples of the form (key, default, help_text)
    :type fields: List[Tuple[str, Union[str, int, bool, None], str]]

    :return:  Dictionary of the form {key: value}
    :rtype: Dict[str, Union[str, int, bool, None]]
    """
    doc: Dict[str, Union[str, int, bool, None]] = {}
    for k, default, help_text in fields:
        val = input(f"Enter {k} parameter. (default: {default}). {help_text}: ")
        doc[k] = val or default
    if "port" in doc:
        doc["port"] = int(doc["port"])  # enforce the port as an int
    return doc


class Color:
    """
    Class to color the output of the CLI.
    """

    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


def safe_open_w(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, "w")


def cli(argv: Optional[Sequence[str]] = None) -> int:
    """
    Command line interface to SurfFlow.
    :param argv: Arguments passed to the CLI.
    :type argv: Optional[Sequence[str]]

    :return: 0 if successful.
    :rtype: int
    """
    m_description = "A command line interface to SurfFlow. For more help on a specific command, type 'surfflow " \
                    "<command> -h'."

    parser = ArgumentParser("surfflow", description=m_description)
    surfflow_version = metadata.version("surfflow")
    parser.add_argument("-v", "--version", action="version", version=surfflow_version)

    subparsers = parser.add_subparsers(help="command", dest="command")

    config_parser = subparsers.add_parser(
        "config", help="Generate config files required for running SurfFlow"
    )
    config_parser.add_argument(
        "init",
        choices=["init"],
        help=f"Start the configuration wizard to generate files in the chosen directory.",
    )
    config_parser.set_defaults(func=generate_config)

    submit_parser = subparsers.add_parser(
        "submit", help="Submit a workflow to the launchpad"
    )
    submit_parser.add_argument(
        "-m", "--mpid", dest="mpid", type=str, help="MPID of the material to submit."
    )
    submit_parser.add_argument(
        "--miller",
        dest="miller",
        type=list,
        nargs="+",
        help="Miller indices of the surface to submit.",
    )
    submit_parser.add_argument(
        "-f",
        "--file",
        dest="filename",
        type=str,
        help="Path to the file containing the input parameters to use.",
    )
    submit_parser.add_argument(
        "-sg",
        "--sg-params",
        dest="sg_params",
        type=json.loads,
        help="SlabGenerator parameters given in a dict format.",
    )
    submit_parser.set_defaults(func=submit)

    args = parser.parse_args(argv)
    if args.command is None:
        # if no command supplied, print help
        parser.print_help()
    elif args.command == "submit":
        print(args)
        args.func(args)
    elif args.command == "config" and args.init:
        generate_config()

    return 0


def generate_init_config() -> None:
    """
    Generate the initialization configuration files.
    """
    fw_config_file = os.environ.get("FW_CONFIG_FILE")
    if fw_config_file:
        prompt = input(
            f"Existing configuration directory found at {os.path.dirname(fw_config_file)}. "
            f"You can define another configuration folder and overwrite this variable. "
            f"Do you want to continue? (y/n) "
        )
    else:
        prompt = "y"

    if prompt == "y":
        path = input(
            f"Enter path to where you want to save the configuration files."
            f" If you already have a configuration folder, you can input that to set it "
            f"as a conda environment variable instead. (default: {DEFAULT_CONF_PATH}): "
        )
        path = path or DEFAULT_CONF_PATH
        with safe_open_w(os.path.join(path, "FW_config.yaml")) as f:
            f.write(f"CONFIG_FILE_DIR: {path}")
        print(f"Configuration directory initialized at {path}!")
        print(
            "Would you like to set the conda environment variable for the config file?"
        )
        while True:
            inp = input("[Y/n]: ")
            if inp in ["y", "Y", ""]:
                print("Setting conda environment variable...")
                subprocess.run(
                    f"conda env config vars set FW_CONFIG_FILE={os.path.join(path, 'FW_config.yaml')}",
                    shell=True,
                )
                break
            elif inp in ["n", "N"]:
                break
            else:
                print("Please enter either Y or n.")
    else:
        print("Aborting configuration wizard.")
    print("\n")


def generate_lpad_config() -> None:
    """
    Generate the launchpad configuration files.
    """
    fields = (
        (
            "host",
            "localhost",
            "Example: 'localhost' or 'mongodb+srv://CLUSTERNAME.mongodb.net'",
        ),
        ("port", 27017, ""),
        (
            "name",
            "fireworks",
            "Database under which to store the fireworks collections",
        ),
        ("username", None, "Username for MongoDB authentication"),
        ("password", None, "Password for MongoDB authentication"),
        (
            "ssl_ca_file",
            None,
            "Path to any client certificate to be used for Mongodb connection",
        ),
        (
            "authsource",
            None,
            "Database used for authentication, if not connection db. e.g., for MongoDB Atlas this is sometimes "
            "'admin'.",
        ),
    )

    print("Please supply the following configuration values")
    print("(press Enter if you want to accept the defaults)\n")
    default_path = set_default_config_path()
    path = input(
        f"Enter path to where you want to save the launchpad configuration file. (default: {default_path}): "
    )
    path = path or default_path
    filename = input(f"Enter filename to write to. (default: {DEFAULT_LPAD_YAML}): ")
    filename = filename or DEFAULT_LPAD_YAML

    doc = get_user_inputs(fields)
    lp = LaunchPad.from_dict(doc)
    lp.to_file(os.path.join(path, filename))
    print(f"\nLaunchpad configuration file written to {os.path.join(path, filename)}!\n")


def generate_db_config() -> None:
    """
    Generate the database configuration files.
    """
    fields = (
        (
            "host",
            "localhost",
            "Example: 'localhost' or 'mongodb+srv://CLUSTERNAME.mongodb.net'",
        ),
        ("port", 27017, ""),
        (
            "database",
            "fireworks",
            "Database under which to store the fireworks collections",
        ),
        ("high_level", "surfflow", "Database under which to store the high-level data"),
        ("admin_user", None, "Admin username for MongoDB authentication"),
        ("admin_password", None, "Admin password for MongoDB authentication"),
        ("readonly_user", None, "MongoDB username for readonly operations"),
        ("admin_password", None, "MongoDB password for readonly operations"),
        (
            "authsource",
            None,
            "Database used for authentication, if not connection db. e.g., for MongoDB Atlas this is sometimes "
            "'admin'.",
        ),
    )

    print("Please supply the following configuration values")
    print("(press Enter if you want to accept the defaults)\n")
    default_path = set_default_config_path()
    path = input(
        f"Enter path to where you want to save the database configuration file. (default: {default_path}): "
    )
    path = path or default_path

    doc = get_user_inputs(fields)
    with safe_open_w(os.path.join(path, "db.json")) as f:
        json.dump(doc, f, indent=4)

    print(f"\nConfiguration written to {os.path.join(path, 'db.json')}!\n")


def generate_fworker_config() -> None:
    """
    Generate the fworker configuration files.
    """
    fields = (
        ("name", "fworker", "Name of the fireworker, must be unique"),
        (
            "db_file",
            os.path.join(DEFAULT_CONF_PATH, "db.json"),
            "Path to the db.json file",
        ),
        ("vasp_cmd", "mpirun -n x vasp_std", "Command to run VASP"),
        (
            "scratch_dir",
            os.path.dirname(DEFAULT_CONF_PATH),
            "Directory for scratch data",
        ),
        ("vdw_kernel_dir", DEFAULT_CONF_PATH, "Directory for vdw kernel data"),
        ("KPAR", 2, "Number of k-points per node"),
        ("NCORE", 2, "Number of cores per node"),
    )

    print("Please supply the following configuration values")
    print("(press Enter if you want to accept the defaults)\n")
    default_path = set_default_config_path()
    path = input(
        f"Enter path to where you want to save the fworker configuration file. (default: {default_path}): "
    )
    path = path or default_path
    doc = get_user_inputs(fields)
    doc["category"] = ""
    doc["query"] = "{}"

    fworker_doc = {
        "name": doc["name"],
        "query": doc["query"],
        "category": doc["category"],
        "env": {
            "db_file": doc["db_file"],
            "scratch_dir": doc["scratch_dir"],
            "vdw_kernel_dir": doc["vdw_kernel_dir"],
            "vasp_cmd": doc["vasp_cmd"],
            "incar_update": {"KPAR": doc["KPAR"], "NCORE": doc["NCORE"]},
        },
    }

    with safe_open_w(os.path.join(path, "my_fworker.yaml")) as outfile:
        yaml.dump(fworker_doc, outfile, default_flow_style=False)

    print(f"\nConfiguration written to {os.path.join(path, 'my_fworker.yaml')}!\n")


def generate_config() -> int:
    """
    Generate the configuration files for SurfFlow.
    """
    # check if the init argument it True, otherwise print help
    print(f"{Color.BOLD}Welcome to the SurfFlow configuration generator!{Color.END}")
    print(
        f"{Color.BOLD}This script will help you generate the necessary configuration files for SurfFlow.{Color.END}"
    )
    print(f"{Color.BOLD}Please choose your action:{Color.END}\n")

    while True:
        choice = input(
            f"{Color.BOLD}(1/init)     initialize/set configuration directory{Color.END}\n"
            f"{Color.BOLD}(2/lpad)     generate my_launchpad.yaml file{Color.END}\n"
            f"{Color.BOLD}(3/db)       generate db.json file{Color.END}\n"
            f"{Color.BOLD}(4/fworker)  generate my_fworker.yaml file{Color.END}\n"
            f"{Color.BOLD}(a)          all{Color.END}\n"
            f"{Color.BOLD}(q)          quit{Color.END}\n"
        )
        if choice in ["1", "init"]:
            generate_init_config()
        elif choice in ["2", "lpad"]:
            generate_lpad_config()
        elif choice in ["3", "db"]:
            generate_db_config()
        elif choice in ["4", "fworker"]:
            generate_fworker_config()
        elif choice in ["a", "all"]:
            generate_lpad_config()
            generate_db_config()
            generate_fworker_config()
        elif choice in ["q", "quit"]:
            break
        else:
            print("Invalid choice. Please try again.\n")
            continue
        # ask the user if they want to generate other config files
        restart = restart_config_wizard()
        if restart:
            continue
        else:
            break
    return 0


# function that asks the user if they want to generate other config files
def restart_config_wizard():
    restart = False
    while True:
        inp = input(
            f"{Color.BOLD}Would you like to generate other configuration files? (y/n){Color.END}: "
        )
        if inp in ["y", "Y", ""]:
            restart = True
            break
        elif inp in ["n", "N"]:
            break
        else:
            print("Please enter either Y or n.")
    return restart
