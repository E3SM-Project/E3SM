from e3sm_diags.parser.core_parser import CoreParser


class MeridionalMean2dParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=MeridionalMean2dParser, *args, **kwargs)  # type: ignore

    def add_arguments(self):
        super().add_arguments()

        # The parameters unique to MeridionalMean2dParser are added here.
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
