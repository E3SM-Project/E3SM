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
    """Ensure that the ordering of the correction is always (ncol, z)"""
    for key in correction:
        if "z" in correction[key].dims:
            correction[key] = correction[key].transpose("ncol", "z")
    return correction


def get_ML_correction_dQ1_dQ2(model, T_mid, qv, cos_zenith, dt):
    """Get ML correction for air temperature (dQ1) and specific humidity (dQ2)

    Args:
        model: pre-trained ML model for dQ1 and dQ2
        T_mid: air temperature
        qv: specific humidity
        cos_zenith: cosine zenith angle
        dt: time step (s)
    """
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    return ensure_correction_ordering(predict(model, ds, dt))


def get_ML_correction_dQu_dQv(model, T_mid, qv, cos_zenith, lat, phis, u, v, dt):
    """Get ML correction for eastward wind (dQu or dQxwind) and northward wind (dQv or dQywind)

    Args:
        model: pre-trained ML model for dQu and dQv
        T_mid: air  temperature
        qv: specific humidity
        cos_zenith: cosine zenith angle
        lat: latitude
        phis: surface geopotential
        u: horizontal wind in x-direction
        v: horizontal wind in y-direction
        dt: time step (s)
    """
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


def get_ML_correction_sfc_fluxes(
    model,
    T_mid,
    qv,
    cos_zenith,
    lat,
    phis,
    sfc_alb_dif_vis,
    sw_flux_dn,
    dt,
):    
    """Get ML correction for overriding surface fluxes (net shortwave and downward longwave)
    ML model should have the following output variables:
        net_shortwave_sfc_flux_via_transmissivity
        override_for_time_adjusted_total_sky_downward_longwave_flux_at_surface

    Args:
        model: pre-trained ML model for radiative fluxes
        T_mid: air  temperature
        qv: specific humidity
        cos_zenith: cosine zenith angle
        lat: latitude
        phis: surface geopotential
        sfc_alb_dif_vis: surface albedo for diffuse shortwave radiation
        sw_flux_dn: downward shortwave flux
        dt: time step (s)
    """    
    SW_flux_dn_at_model_top = sw_flux_dn[:, 0]
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            lat=(["ncol"], lat),
            surface_geopotential=(["ncol"], phis),
            cos_zenith_angle=(["ncol"], cos_zenith),
            surface_diffused_shortwave_albedo=(["ncol"], sfc_alb_dif_vis),
            total_sky_downward_shortwave_flux_at_top_of_atmosphere=(
                ["ncol"],
                SW_flux_dn_at_model_top,
            ),
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
    phis,
    sw_flux_dn,
    sfc_alb_dif_vis,
    sfc_flux_sw_net,
    sfc_flux_lw_dn,
    Ncol,
    Nlev,
    num_tracers,
    dt,
    model_tq,
    model_uv,
    model_sfc_fluxes,
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
    SW_flux_dn_at_model_top: downwelling shortwave flux at the top of the model
    sfc_alb_dif_vis: surface diffuse shortwave albedo
    sfc_flux_sw_net
    sfc_flux_lw_dn
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
        correction_tq = get_ML_correction_dQ1_dQ2(
            model_tq, T_mid, qv[:, 0, :], cos_zenith, dt
        )
        T_mid[:, :] += correction_tq["dQ1"].values * dt
        qv[:, 0, :] += correction_tq["dQ2"].values * dt
    if model_uv is not None:
        correction_uv = get_ML_correction_dQu_dQv(
            model_uv, T_mid, qv[:, 0, :], cos_zenith, lat, phis, u, v, dt
        )
        u[:, :] += correction_uv["dQu"].values * dt
        v[:, :] += correction_uv["dQv"].values * dt
    if model_sfc_fluxes is not None:
        sw_flux_dn = np.reshape(sw_flux_dn, (-1, Nlev+1))
        correction_sfc_fluxes = get_ML_correction_sfc_fluxes(
            model_sfc_fluxes,
            T_mid,
            qv[:, 0, :],
            cos_zenith,
            lat,
            phis,
            sfc_alb_dif_vis,
            sw_flux_dn,
            dt,
        )
        sfc_flux_sw_net[:] = correction_sfc_fluxes["net_shortwave_sfc_flux_via_transmissivity"].values
        sfc_flux_lw_dn[:] = correction_sfc_fluxes["override_for_time_adjusted_total_sky_downward_longwave_flux_at_surface"].values
