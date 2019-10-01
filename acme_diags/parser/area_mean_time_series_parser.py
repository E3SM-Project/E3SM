from .core_parser import CoreParser
from acme_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter


class AreaMeanTimeSeriesParser(CoreParser):
    def __init__(self, *args, **kwargs):
        if 'parameter_cls' in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(parameter_cls=AreaMeanTimeSeriesParameter, *args, **kwargs)


    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)

        self.add_argument(
            '--ref_names',
            type=str,
            nargs='+',
            dest='ref_names',
            help='List of reference names.',
            required=False)

        self.add_argument(
            '--ref_timeseries_input',
            dest='ref_timeseries_input',
            help='The input reference data are timeseries files.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--test_timeseries_input',
            dest='test_timeseries_input',
            help='The input test data are timeseries files.',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--start_yr',
            dest='start_yr',
            help="Start year for the timeseries files.",
            required=False)

        self.add_argument(
            '--end_yr',
            dest='end_yr',
            help="End year for the timeseries files.",
            required=False)




