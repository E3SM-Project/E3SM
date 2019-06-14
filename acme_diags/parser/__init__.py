from .core_parser import CoreParser
from .zonal_mean_2d_parser import ZonalMean2dParser
from .meridional_mean_2d_parser import MeridionalMean2dParser


SET_TO_PARSER = {
    'zonal_mean_xy': CoreParser,
    'zonal_mean_2d': ZonalMean2dParser,
    'meridional_mean_2d': MeridionalMean2dParser,
    'lat_lon': CoreParser,
    'polar': CoreParser,
    'cosp_histogram': CoreParser,
    'area_mean_time_series': CoreParser,
}

