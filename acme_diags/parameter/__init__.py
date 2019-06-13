from .core_parameter import CoreParameter
from .zonal_mean_2d_parameter import ZonalMean2dParameter

SET_TO_PARAMETERS = {
    'zonal_mean_xy': CoreParameter,
    'zonal_mean_2d': ZonalMean2dParameter,
    'meridional_mean_2d': CoreParameter,
    'lat_lon': CoreParameter,
    'polar': CoreParameter,
    'cosp_histogram': CoreParameter,
    'area_mean_time_series': CoreParameter,
}
