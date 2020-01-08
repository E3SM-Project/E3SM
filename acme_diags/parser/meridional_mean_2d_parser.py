from .core_parser import CoreParser
from acme_diags.parameter.meridional_mean_2d_parameter import MeridionalMean2dParameter


class MeridionalMean2dParser(CoreParser):
    def __init__(self, *args, **kwargs):
        if 'parameter_cls' in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(parameter_cls=MeridionalMean2dParser, *args, **kwargs)


    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)

        # The parameters unique to MeridionalMean2dParser are added here.
        self.add_argument(
            '--plevs',
            type=float,
            nargs='+',
            dest='plevs',
            help='Selected pressure level.[take list as input]',
            required=False)

        self.add_argument(
            '--plot_plevs',
            dest='plot_plevs',
            help='plot specified plevs',
            action='store_const',
            const=True,
            required=False)

        self.add_argument(
            '--plot_log_plevs',
            dest='plot_log_plevs',
            help='plot plevs on log-scale',
            action='store_const',
            const=True,
            required=False)
