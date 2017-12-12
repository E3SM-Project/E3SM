import genutil
import numpy
import cdp.cdp_metric


class STD(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(STD, self).__init__(metric_path)

    def compute(self, variable, axis='xy'):
        std = -numpy.infty
        try:
            std = float(genutil.statistics.std(variable, axis=axis, weights='generate'))
        except Exception as err:
            print(err)

        return std
