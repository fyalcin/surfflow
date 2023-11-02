import abc

from fireworks import FiretaskBase

from surfflow.utils.misc_tools import check_input


class SurfenFT(FiretaskBase, metaclass=abc.ABCMeta):
    """
    Base class for all Firetasks in the surfflow package. Currently, not in use.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._inp = self.load_params()

    def load_params(self):
        """
        Loads the input parameters for the Firetask from the Firework spec.
        :return: Dictionary of input parameters.
        :rtype: dict
        """
        req_inp = {k: self.get(k) for k in self.required_params}
        opt_inp = {k: self.get(k) for k in self.optional_params}
        opt_inp = check_input(opt_inp, self.optional_params)
        inp = {**req_inp, **opt_inp}
        return inp
