import numpy as np
import h5py
import xarray as xr
import fv3fit
from scream_run.steppers.machine_learning import (
    MultiModelAdapter,
    predict,
)


def get_ML_model(model_path):
    if model_path == "NONE":
        return None
    config = MachineLearningConfig(models=[model_path])
    model = open_model(config)
    return model


def sample_ML_prediction(
    nz: int, input_data: np.ndarray, ML_model_tq: str, ML_model_uv: str
):
    """
    This function is used to generate a sample ML prediction for the given input data.
    We use a constant output predictor to generate the prediction.
    """
    output_variables = ["qv"]
    outputs = {
        "qv": np.full(nz, 1e-4),
    }
    predictor = fv3fit.testing.ConstantOutputPredictor(
        input_variables=["qv"],
        output_variables=output_variables,
    )
    predictor.set_outputs(**outputs)
    model = MultiModelAdapter([predictor])
    if len(input_data.shape) < 2:
        input_data = input_data[np.newaxis, :]
    input_data = xr.Dataset({"qv": xr.DataArray(data=input_data, dims=["ncol", "z"])})
    output = predict(model, input_data, dt=1.0)
    return output["qv"].values


def modify_view(data, Ncol, Nlev, model_tq, model_uv):
    data = np.reshape(data, (-1, Nlev))
    prediction = sample_ML_prediction(Nlev, data[1, :], model_tq, model_uv)
    data[1, :] = prediction
