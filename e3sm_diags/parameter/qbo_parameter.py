from .time_series_parameter import TimeSeriesParameter


class QboParameter(TimeSeriesParameter):
    def __init__(self):
        super(QboParameter, self).__init__()
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True
        self.granulate.remove("seasons")
