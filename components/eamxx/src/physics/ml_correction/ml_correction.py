import numpy as np
import xarray as xr
import datetime
from vcm import cos_zenith_angle
from scream_run.steppers.machine_learning import (
    MachineLearningConfig,
    open_model,
    predict,
)


def get_ML_correction(model_path, T_mid, qv, cos_zenith, current_time):
    config = MachineLearningConfig(models=[model_path])
    model = open_model(config)
    ds = xr.Dataset(
        data_vars=dict(
            T_mid=(["ncol", "z"], T_mid),
            qv=(["ncol", "z"], qv),
            cos_zenith_angle=(["ncol"], cos_zenith),
        )
    )
    return predict(model, ds)


def update_fields(T_mid, qv, u, v, lat, lon, Ncol, Nlev, model_path, current_time):
    T_mid = np.reshape(T_mid, (-1, Nlev))
    qv = np.reshape(qv, (-1, Nlev))
    current_datetime = datetime.datetime.strptime(current_time, "%Y-%m-%d %H:%M:%S")
    cos_zenith = cos_zenith_angle(
        current_datetime,
        lon,
        lat,
    )
    correction = get_ML_correction(model_path, T_mid, qv, cos_zenith, current_datetime)
    print(f"[Python] prediction for dQ1 is of shape {correction['dQ1'].shape}")
    print(f"[Python] prediction for dQ2 is of shape {correction['dQ2'].shape}")
    print("[Python] update fields completed without any changes")
