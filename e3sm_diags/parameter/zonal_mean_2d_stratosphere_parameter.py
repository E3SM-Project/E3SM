import numpy

from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter


class ZonalMean2dStratosphereParameter(ZonalMean2dParameter):
    def __init__(self):
        super(ZonalMean2dStratosphereParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.plevs = numpy.logspace(0, 2.0, num=10).tolist()
        self.plot_log_plevs = True
