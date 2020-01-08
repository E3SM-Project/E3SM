import numpy
from .core_parameter import CoreParameter
from acme_diags.driver.utils.general import monotonic

class ZonalMean2dParameter(CoreParameter):
    def __init__(self):
        super(ZonalMean2dParameter, self).__init__()
        self.plevs = numpy.logspace(2.0, 3.0, num=17).tolist()
        self.plot_log_plevs = False 
        self.plot_plevs = False 
        # Granulating plevs causes duplicate plots in this case.
        # So keep all of the default values except plevs.
        self.granulate.remove('plevs')

    def check_values(self):
        plevs = self.plevs
        if not isinstance(plevs,list):
           msg = 'plevs needs to be a list'
           raise RuntimeError(msg)


        if len(plevs) > 1:
            if monotonic(plevs):
                if plevs[0] > plevs[1]:
                    plevs = plevs[::-1]  # Force plevs to be monotonically increasing by reversing the list. 
                    self.plevs = plevs
            else:
                msg = 'plevs should be monotonically increasing or decreasing'
                raise RuntimeError(msg)
        else:
            msg = "At least 2 plevs needed"
            raise RuntimeError(msg)
             
                   

