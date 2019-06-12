from .core_parser import CoreParser
from acme_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter


class ZonalMean2dParser(CoreParser):
    def __init__(self, *args, **kwargs):
        if 'parameter_cls' in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(parameter_cls=ZonalMean2dParameter, *args, **kwargs)


    def load_default_args(self, files=[]):
        # This has '-p' and '--parameter' reserved.
        super().load_default_args(files)

        # The parameters unique to ZonalMean2dParameter are added here.
        self.add_argument(
            '--plevs',
            type=float,
            nargs='+',
            dest='plevs',
            help='Selected pressure level.',
            required=False)

