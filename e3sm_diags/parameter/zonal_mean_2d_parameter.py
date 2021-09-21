import numpy

from e3sm_diags.driver.utils.general import monotonic

from .core_parameter import CoreParameter


class ZonalMean2dParameter(CoreParameter):
    def __init__(self):
        super(ZonalMean2dParameter, self).__init__()
        self.plevs = numpy.linspace(50, 1000, 20).tolist()
        # self.plevs = numpy.logspace(2.0, 3.0, num=17).tolist()
        # self.plevs = [
        #    30.0,
        #    50.0,
        #    75.0,
        #    100.0,
        #    150.0,
        #    200.0,
        #    250.0,
        #    300.0,
        #    350.0,
        #    400.0,
        #    450.0,
        #    500.0,
        #    550.0,
        #    600.0,
        #    650.0,
        #    700.0,
        #    750.0,
        #    800.0,
        #    850.0,
        #    875.0,
        #    900.0,
        #    925.0,
        #    950.0,
        #    975.0,
        #    1000.0,
        # ]
        self.plot_log_plevs = False
        self.plot_plevs = False
        # Granulating plevs causes duplicate plots in this case.
        # So keep all of the default values except plevs.
        self.granulate.remove("plevs")

    def check_values(self):
        plevs = self.plevs
        if not isinstance(plevs, list):
            msg = "plevs needs to be a list"
            raise RuntimeError(msg)

        if len(plevs) > 1:
            if monotonic(plevs):
                if plevs[0] > plevs[1]:
                    plevs = plevs[
                        ::-1
                    ]  # Force plevs to be monotonically increasing by reversing the list.
                    self.plevs = plevs
            else:
                msg = "plevs should be monotonically increasing or decreasing"
                raise RuntimeError(msg)
        else:
            msg = "At least 2 plevs needed"
            raise RuntimeError(msg)
