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


def ensure_correction_ordering(correction):
    for key in correction:
        correction[key] = correction[key].transpose("ncol", "z")
    return correction


def get_ML_correction_dQ1_dQ2(model, T_mid, qv, cos_zenith, dt):
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    return ensure_correction_ordering(predict(model, ds, dt))

def get_ML_correction_dQu_dQv(model, T_mid, qv, cos_zenith, lat, phis, u, v, dt):
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            U=(["ncol", "z"], u),
            V=(["ncol", "z"], v),
            lat=(["ncol"], lat),
            surface_geopotential=(["ncol"], phis),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    output = ensure_correction_ordering(predict(model, ds, dt))
    # rename dQxwind and dQywind to dQu and dQv if needed
    if "dQxwind" in output.keys():
        output["dQu"] = output.pop("dQxwind")
    if "dQywind" in output.keys():
        output["dQv"] = output.pop("dQywind")
            
    return output

def update_fields(
    T_mid,
    qv,
    u,
    v,
    lat,
    lon,
    phis,
    Ncol,
    Nlev,
    num_tracers,
    dt,
    model_tq,
    model_uv,
    current_time,
):
    """
    T_mid: temperature
    qv: specific humidity
    u: x-component of wind
    v: y-component of wind
    lat: latitude
    lon: longitude
    phis: surface geopotential
    Ncol: number of columns
    Nlev: number of levels
    num_tracers: number of tracers
    dt: time step (s)
    model_tq: path to the ML model for temperature and specific humidity
    model_uv: path to the ML model for u and v
    current_time: current time in the format "YYYY-MM-DD HH:MM:SS"
    """
    T_mid = np.reshape(T_mid, (-1, Nlev))
    u = np.reshape(u, (-1, Nlev))
    v = np.reshape(v, (-1, Nlev))
    # qv is a 3D array of shape (Ncol, num_tracers, Nlev)
    # by default, qv is the frist tracer variable
    qv = np.reshape(qv, (-1, num_tracers, Nlev))
    current_datetime = datetime.datetime.strptime(current_time, "%Y-%m-%d %H:%M:%S")
    cos_zenith = cos_zenith_angle(
        current_datetime,
        lon,
        lat,
    )
    if model_tq is not None:
        correction_tq = get_ML_correction_dQ1_dQ2(model_tq, T_mid, qv[:, 0, :], cos_zenith, dt)
        T_mid[:, :] += correction_tq["dQ1"].values * dt
        qv[:, 0, :] += correction_tq["dQ2"].values * dt
    if model_uv is not None:
        correction_uv = get_ML_correction_dQu_dQv(model_uv, T_mid, qv[:, 0, :], cos_zenith, lat, phis, u, v, dt)
        u[:, :] += correction_uv["dQu"].values * dt
        v[:, :] += correction_uv["dQv"].values * dt        