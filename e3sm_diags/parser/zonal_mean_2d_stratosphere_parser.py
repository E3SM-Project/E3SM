from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    ZonalMean2dStratosphereParameter,
)
from e3sm_diags.parser.zonal_mean_2d_parser import ZonalMean2dParser


class ZonalMean2dStratosphereParser(ZonalMean2dParser):
    def __init__(self, *args, **kwargs):
        if "parameter_cls" in kwargs:
            super().__init__(*args, **kwargs)
        else:
            super().__init__(
                parameter_cls=ZonalMean2dStratosphereParameter, *args, **kwargs
            )
