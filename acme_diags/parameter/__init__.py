from .core_parameter import CoreParameter
from .zonal_mean_2d_parameter import ZonalMean2dParameter
from .meridional_mean_2d_parameter import MeridionalMean2dParameter
from .area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
from .enso_diags_parameter import EnsoDiagsParameter

SET_TO_PARAMETERS = {
    'zonal_mean_xy': CoreParameter,
    'zonal_mean_2d': ZonalMean2dParameter,
    'meridional_mean_2d': MeridionalMean2dParameter,
    'lat_lon': CoreParameter,
    'polar': CoreParameter,
    'cosp_histogram': CoreParameter,
    'area_mean_time_series': AreaMeanTimeSeriesParameter,
    'enso_diags': EnsoDiagsParameter,
}
