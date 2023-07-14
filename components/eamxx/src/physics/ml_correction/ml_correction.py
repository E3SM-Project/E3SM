import numpy as np
import xarray as xr
import datetime
from vcm import cos_zenith_angle
from scream_run.steppers.machine_learning import (
    MultiModelAdapter,
    predict,
)


def get_ML_correction(model_path, T_mid, qv):
    pass


def update_fields(T_mid, qv, u, v, lat, lon, Ncol, Nlev, model_path, current_time):
    T_mid = np.reshape(T_mid, (-1, Nlev))
    cos_zenith = cos_zenith_angle(
        datetime.datetime.strptime(current_time, "%Y-%m-%d %H:%M:%S"),
        lon,
        lat,
    )
    print("[Python] update fields completed without any changes")
