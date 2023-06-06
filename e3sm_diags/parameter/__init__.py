from .annual_cycle_zonal_mean_parameter import ACzonalmeanParameter
from .area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
from .arm_diags_parameter import ARMDiagsParameter
from .core_parameter import CoreParameter
from .diurnal_cycle_parameter import DiurnalCycleParameter
from .enso_diags_parameter import EnsoDiagsParameter
from .lat_lon_land_parameter import LatLonLandParameter
from .lat_lon_river_parameter import LatLonRiverParameter
from .meridional_mean_2d_parameter import MeridionalMean2dParameter
from .mp_partition_parameter import MPpartitionParameter
from .qbo_parameter import QboParameter
from .streamflow_parameter import StreamflowParameter
from .tc_analysis_parameter import TCAnalysisParameter
from .zonal_mean_2d_parameter import ZonalMean2dParameter
from .zonal_mean_2d_stratosphere_parameter import ZonalMean2dStratosphereParameter

SET_TO_PARAMETERS = {
    "zonal_mean_xy": CoreParameter,
    "zonal_mean_2d": ZonalMean2dParameter,
    "zonal_mean_2d_stratosphere": ZonalMean2dStratosphereParameter,
    "meridional_mean_2d": MeridionalMean2dParameter,
    "lat_lon": CoreParameter,
    "polar": CoreParameter,
    "cosp_histogram": CoreParameter,
    "area_mean_time_series": AreaMeanTimeSeriesParameter,
    "enso_diags": EnsoDiagsParameter,
    "qbo": QboParameter,
    "streamflow": StreamflowParameter,
    "diurnal_cycle": DiurnalCycleParameter,
    "arm_diags": ARMDiagsParameter,
    "tc_analysis": TCAnalysisParameter,
    "annual_cycle_zonal_mean": ACzonalmeanParameter,
    "lat_lon_land": LatLonLandParameter,
    "lat_lon_river": LatLonRiverParameter,
    "aerosol_aeronet": CoreParameter,
    "aerosol_budget": CoreParameter,
    "mp_partition": MPpartitionParameter,
}
