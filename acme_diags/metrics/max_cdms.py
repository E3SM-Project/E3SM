import cdp.cdp_metric

class Max(cdp.cdp_metric.CDPMetric):
    def __init__(self):
        metric_path = __file__
        super(Max, self).__init__(metric_path)

    def compute(self, variable):
        return float(variable.max())
