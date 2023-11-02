import warnings

from htflow_utils.misc_tools import load_defaults

defaults = load_defaults("surfflow")


def check_input(input_dict: dict, required_keys: list[str]) -> dict:
    """
    Checks if the input dictionary contains all the required keys and replaces the ones
    missing with the default values.

    :param input_dict: Dictionary containing all the input parameters needed
        to start a subworkflow.
    :type input_dict: dict

    :param required_keys: List of required keys.
    :type required_keys: list[str]

    :return: Input dictionary with the required keys.
    :rtype: dict
    """
    input_params_dict = {}
    for key in required_keys:
        default = defaults.get(key, None)
        user_input = input_dict.get(key)
        try:
            input_params_dict[key] = {**default, **user_input}
        except TypeError:
            input_params_dict[key] = user_input or default
    return input_params_dict


def check_input_dict(input_dict):
    """
    Checks if the input dictionary contains all the required keys and replaces the ones
    missing with the default values.

    :param input_dict: Dictionary containing all the input parameters needed
        to start a subworkflow.
    :type input_dict: dict

    :return: Input dictionary with the required keys.
    :rtype: dict
    """
    final_dict = {}
    all_extra_keys = {}
    all_missing_keys = {}
    for k, v in input_dict.items():
        default_val = defaults.get(k, None)
        if default_val is None:
            raise KeyError(f"Key {k} not found in defaults.")
        # check if v is a dict
        if isinstance(v, dict):
            # check if all the keys in v are in default_dict
            extras = set(v.keys()) - set(default_val.keys())
            if extras:
                all_extra_keys[k] = list(extras)
            updated_v = {**default_val, **v}
            # check if "NO_DEFAULT" is in any of the values
            missing_required = [k for k, v in updated_v.items() if v == "NO_DEFAULT"]
            if missing_required:
                all_missing_keys[k] = missing_required
            final_dict[k] = updated_v
        else:
            final_dict[k] = v or default_val
            if final_dict[k] == "NO_DEFAULT":
                all_missing_keys[k] = [k]
    if all_extra_keys:
        warnings.warn(f"Extra keys found in the input dictionary: {all_extra_keys}")
    if all_missing_keys:
        warnings.warn(f"Required keys not found in the input dictionary: {all_missing_keys}"
        )

    return final_dict