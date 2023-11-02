from fireworks import explicit_serialize, FiretaskBase, FWAction

from surfflow.utils.db_tools import move_result
from surfflow.utils.misc_tools import check_input


@explicit_serialize
class MoveResults(FiretaskBase):
    """
    Firetask to move the result of a VASP calculation from the Fireworks database to the destination.

    :param tag: Task label assigned to the calculation.
    :type tag: str

    :param fltr: Filter to use when looking up results in the database. Generally involves either the mpid or task_label.
    :type fltr: dict

    :param coll: Collection to move the results into in the destination database.
    :type coll: str

    :param to_loc: Location in the collection the results should be written to.
    :type to_loc: List[str]

    :param from_loc: Location in the collection the results should be read from.
    :type from_loc: List[str]

    :param custom_dict: Custom dictionary to write into the destination.
    :type custom_dict: dict, optional

    :param db_file: Full path of the db.json file.
    :type db_file: str, optional

    :param high_level: Name of the high level database to use. If set to False, the low level database is used.
    :type high_level: bool, str, optional

    :param fake_calc: Whether the calculation is a fake calculation.
    :type fake_calc: bool, optional

    :param project_fields: Fields to project from the Fireworks database.
    :type project_fields: List[str], optional

    :param rename: Dictionary to rename the fields in the database.
    :type rename: dict, optional

    .. note::
        The default values of the optional parameters are populated from the defaults.json file if not provided.
    """

    _fw_name = "Move the results from the tag to the entry found by flag."
    required_params = ["tag", "fltr", "coll", "to_loc", "from_loc"]
    optional_params = [
        "custom_dict",
        "db_file",
        "high_level",
        "fake_calc",
        "project_fields",
        "rename",
    ]

    def run_task(self, fw_spec):
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        from_dest = {
            "db_file": inp["db_file"],
            "high_level": False,
            "coll": "tasks",
            "fltr": {"task_label": inp["tag"]},
            "path": inp["from_loc"],
        }
        to_dest = {
            "db_file": inp["db_file"],
            "high_level": inp["high_level"],
            "coll": inp["coll"],
            "fltr": inp["fltr"],
            "path": inp["to_loc"],
        }

        project_fields = inp.get("project_fields", None)
        if not project_fields:
            project_fields = (
                ["input", "output"]
                if inp["fake_calc"]
                else ["input", "output", "orig_inputs", "custodian", "run_stats"]
            )
        move_result(from_dest=from_dest, to_dest=to_dest, project_fields=project_fields)

        return FWAction(update_spec=fw_spec)
