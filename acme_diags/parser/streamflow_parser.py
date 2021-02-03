from .core_parser import CoreParser
from acme_diags.parameter.streamflow_parameter import StreamflowParameter


class StreamflowParser(CoreParser):
    def __init__(self, *args, **kwargs):
        if "parameter_cls" in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(
                parameter_cls=StreamflowParameter, *args, **kwargs
            )

    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)

        self.add_argument(
            "--gauges_path",
            dest="gauges_path",
            help="The file containing the gauge data.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--max_num_gauges",
            dest="max_num_gauges",
            help="The maximum number of gauges that should be processed.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--print_statements",
            dest="print_statements",
            help="Print information useful for debugging.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--ref_timeseries_input",
            dest="ref_timeseries_input",
            help="The input reference data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--test_timeseries_input",
            dest="test_timeseries_input",
            help="The input test data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--start_yr",
            dest="start_yr",
            help="Start year for the timeseries files.",
            required=False,
        )

        self.add_argument(
            "--end_yr",
            dest="end_yr",
            help="End year for the timeseries files.",
            required=False,
        )
