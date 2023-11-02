import json
import os
from types import SimpleNamespace
from typing import Union, Any

from atomate.vasp.database import VaspCalcDb

from htflow_utils.misc_tools import get_by_path


def move_result(
        from_dest: dict,
        to_dest: dict,
        project_fields: list,
        custom_dict: dict = None,
        rename: list = None,
) -> None:
    """
    Move the result of a VASP calculation from the Fireworks database to the destination.
    Parameters

    :param from_dest: Dictionary containing the information to connect to the source database.
    :type from_dest: dict

    :param to_dest: Dictionary containing the information to connect to the destination database.
    :type to_dest: dict

    :param project_fields: List of fields to project from the source database.
    :type project_fields: list

    :param custom_dict: Custom dictionary to write into the destination. Defaults to {}.
    :type custom_dict: dict, optional

    :param rename: List of names corresponding to the project fields so that the fields are renamed
        on the destination.
    :type rename: list, optional

    :return: None.
    :rtype: None
    """
    from_dest = SimpleNamespace(**from_dest)
    nav_from = VaspDB(db_file=from_dest.db_file, high_level=from_dest.high_level)
    calc = nav_from.find_data(from_dest.coll, from_dest.fltr)
    print(f'calc: {calc}')
    # out is the whole output of the calculation including everything,
    # fltr_out = ['input', 'output'] if fake_calc else ['input', 'output', 'orig_inputs', 'custodian']
    out = get_by_path(calc, from_dest.path)
    out = {f: out[f] for f in project_fields}
    if custom_dict:
        out.update(custom_dict)

    to_dest = SimpleNamespace(**to_dest)
    nav_to = VaspDB(db_file=to_dest.db_file, high_level=to_dest.high_level)

    path = ".".join(to_dest.path)

    if rename and len(rename) == len(project_fields):
        nav_to.update_data(
            collection=to_dest.coll,
            fltr=to_dest.fltr,
            new_values={
                "$set": {path + f".{rename[i]}": v for i, v in enumerate(out.values())}
            },
            upsert=True,
        )
    else:
        nav_to.update_data(
            collection=to_dest.coll,
            fltr=to_dest.fltr,
            new_values={"$set": {path + f".{k}": v for k, v in out.items()}},
            upsert=True,
        )


class VaspDB:
    """
    Class to interact with the VASP database.
    """

    def __init__(self, db_file: str = "auto", high_level: Union[str, bool] = False):
        """
        Class to interact with the VASP database.

        :param db_file: Full path of the db.json file. Defaults to 'auto', in which case the path
            is checked from the environment variable FW_CONFIG_FILE.
        :type db_file: str, optional

        :param high_level: Name of the high level database to use. If set to False, the low level database is used.
            Defaults to True, in which case the value in the db.json file is used.
        :type high_level: bool, str, optional
        """
        if db_file == "auto":
            config_path = os.environ.get("FW_CONFIG_FILE")
            if not config_path:
                raise Exception(
                    "The environmental variable FW_CONFIG_FILE "
                    "is not set. Please set it to the path of "
                    "your db.json file."
                )

            db_file = config_path.replace("FW_config.yaml", "db.json")

        self.db_file = db_file
        self.high_level = high_level

        with open(db_file, "r") as f:
            db_dict = json.load(f)
        self.db_dict = db_dict

        self.vasp_calc_db = VaspCalcDb.from_db_file(db_file, admin=True)

        try:
            db = self.vasp_calc_db.connection[high_level]
        except TypeError:
            if high_level:
                db = self.vasp_calc_db.connection[db_dict.get("high_level")]
            else:
                db = self.vasp_calc_db.db
        self.db = db

    def find_data(self, collection: str, fltr: dict, projection: dict = None):
        return self.db[collection].find_one(filter=fltr, projection=projection)

    def find_many_data(self, collection: str, fltr: dict, projection: dict = None):
        return self.db[collection].find(filter=fltr, projection=projection)

    def update_data(
            self, collection: str, fltr: dict, new_values: dict, upsert: bool = True
    ):
        self.db[collection].update_one(fltr, new_values, upsert=upsert)

    def insert_data(self, collection: str, data: dict):
        self.db[collection].insert_one(data)

    def delete_data(self, collection: str, fltr: dict):
        self.db[collection].delete_one(fltr)


def get_entry_by_loc(
        nav: VaspDB, fltr: dict, coll: str, loc: list[str]
) -> Union[dict[Any, Any], object]:
    """
    Locates and returns the entry at the given location.

    :param nav: VaspDB object for the database being queried.
    :type nav: VaspDB

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or task_label.
    :type fltr: dict

    :param coll: MongoDB collection in which to search for.
    :type coll: str

    :param loc: Location in which to search the entry at.
    :type loc: list

    :return: The database entry queried or None.
    :rtype: Union[dict[Any, Any], object]
    """
    data = nav.find_data(coll, fltr)
    if data:
        try:
            entry = get_by_path(data, loc)
        except KeyError:
            print(
                f'No entry found at {".".join([l for l in loc])} in {coll} for {fltr}'
            )
            return None
        if entry:
            return entry
