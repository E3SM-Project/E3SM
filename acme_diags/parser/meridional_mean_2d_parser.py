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
            help='Selected pressure level.',
            required=False)

