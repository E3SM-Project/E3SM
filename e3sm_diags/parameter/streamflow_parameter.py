from .time_series_parameter import TimeSeriesParameter


class StreamflowParameter(TimeSeriesParameter):
    def __init__(self):
        super(StreamflowParameter, self).__init__()
        self.gauges_path = None
        self.main_title_seasonality_map = "Seasonality Map"
        self.main_title_annual_map = "Mean Annual Streamflow Map"
        self.main_title_annual_scatter = "Mean Annual Streamflow Scatter Plot"
        self.max_num_gauges = None
        self.output_file_seasonality_map = "seasonality_map"
        self.output_file_annual_map = "annual_map"
        self.output_file_annual_scatter = "annual_scatter"
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True
        self.granulate.remove("seasons")
