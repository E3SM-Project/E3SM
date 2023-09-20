import numpy as np
import xarray as xr
import datetime
from vcm import cos_zenith_angle
from scream_run.steppers.machine_learning import (
    MachineLearningConfig,
    open_model,
    predict,
)

def get_ML_model(model_path):
    if model_path == "NONE":
        return None
    config = MachineLearningConfig(models=[model_path])
    model = open_model(config)
    return model    

def get_ML_correction(model, T_mid, qv, cos_zenith, dt):
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    return predict(model, ds, dt)


def update_fields(
    T_mid,
    qv,
    u,
    v,
    lat,
    lon,
    Ncol,
    Nlev,
    num_tracers,
    dt,
    model,
    current_time,
):
    T_mid = np.reshape(T_mid, (-1, Nlev))
    # qv is a 3D array of shape (Ncol, num_tracers, Nlev)
    # by default, qv is the frist tracer variable
    qv = np.reshape(qv, (-1, num_tracers, Nlev))
    current_datetime = datetime.datetime.strptime(current_time, "%Y-%m-%d %H:%M:%S")
    cos_zenith = cos_zenith_angle(
        current_datetime,
        lon,
        lat,
    )
    correction = get_ML_correction(model, T_mid, qv[:, 0, :], cos_zenith, dt)
    T_mid[:, :] += correction["dQ1"].values * dt
    qv[:, 0, :] += correction["dQ2"].values * dt
