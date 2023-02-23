from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.parser.core_parser import CoreParser


class ZonalMean2dParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=ZonalMean2dParameter, *args, **kwargs)  # type: ignore

    def add_arguments(self):
        super().add_arguments()

        self.parser.add_argument(
            "--plevs",
            type=float,
            nargs="+",
            dest="plevs",
            help="Selected pressure level.[take list as input]",
            required=False,
        )

        self.parser.add_argument(
            "--plot_plevs",
            dest="plot_plevs",
            help="plot specified plevs",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--plot_log_plevs",
            dest="plot_log_plevs",
            help="plot plevs on log-scale",
            action="store_const",
            const=True,
            required=False,
        )
