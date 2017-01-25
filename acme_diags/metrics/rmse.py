import numpy
import genutil
import cdp.cdp_metric

class RMSE(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(RMSE, self).__init__(metric_path)

    def compute(self, model, obs):
        rmse = -numpy.infty
        try:
            rmse = float(genutil.statistics.rms(model, obs, axis='xy', weights='generate'))
        except Exception, err:
            print err
        return rmse
