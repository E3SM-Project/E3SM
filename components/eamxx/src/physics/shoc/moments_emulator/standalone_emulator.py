from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn

"""Standalone emulator exported from a trained SHOC experiment.

This wrapper loads its trained weights and normalization from the adjacent
`standalone_emulator.pt` file.
"""

MODEL_FILENAME = 'standalone_emulator.pt'
SUPPORTED_FORMAT = 1

_BLOB_CACHE = None
_MODEL_CACHE = None


class ExportedMlp(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_dims: list[int],
        activation: str,
        dropout: float,
    ) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        current_dim = input_dim
        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(current_dim, int(hidden_dim)))
            layers.append(_activation_module(activation))
            if dropout > 0.0:
                layers.append(nn.Dropout(dropout))
            current_dim = int(hidden_dim)
        layers.append(nn.Linear(current_dim, output_dim))
        self.network = nn.Sequential(*layers)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.network(inputs)


class ExportedColumnConv(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_channels: list[int],
        kernel_size: int,
        activation: str,
        dropout: float,
    ) -> None:
        super().__init__()
        padding = kernel_size // 2
        layers: list[nn.Module] = []
        in_channels = input_dim
        for hidden_channel in hidden_channels:
            layers.append(nn.Conv1d(in_channels, int(hidden_channel), kernel_size=kernel_size, padding=padding))
            layers.append(_activation_module(activation))
            if dropout > 0.0:
                layers.append(nn.Dropout(dropout))
            in_channels = int(hidden_channel)
        layers.append(nn.Conv1d(in_channels, output_dim, kernel_size=1))
        self.network = nn.Sequential(*layers)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.network(inputs.transpose(1, 2)).transpose(1, 2)


def _activation_module(name: str) -> nn.Module:
    registry = {
        "relu": nn.ReLU,
        "gelu": nn.GELU,
        "silu": nn.SiLU,
        "tanh": nn.Tanh,
    }
    if name not in registry:
        raise KeyError(f"Unsupported activation {name!r}.")
    return registry[name]()


def _load_blob(model_path: str | Path | None = None) -> dict[str, Any]:
    global _BLOB_CACHE
    resolved = Path(model_path) if model_path is not None else Path(__file__).with_name(MODEL_FILENAME)
    if _BLOB_CACHE is None or Path(_BLOB_CACHE["_model_path"]) != resolved:
        blob = torch.load(resolved, map_location="cpu")
        if int(blob["format_version"]) != SUPPORTED_FORMAT:
            raise RuntimeError(
                f"Unsupported standalone emulator format {blob['format_version']}; "
                f"expected {SUPPORTED_FORMAT}."
            )
        blob["_model_path"] = str(resolved)
        _BLOB_CACHE = blob
    return _BLOB_CACHE


def _build_model(blob: dict[str, Any]) -> nn.Module:
    model_cfg = blob["model"]
    input_dim = len(blob["feature_names"])
    output_dim = len(blob["target_names"])
    model_type = model_cfg["type"]
    if model_type == "mlp":
        model = ExportedMlp(
            input_dim=input_dim,
            output_dim=output_dim,
            hidden_dims=list(model_cfg["hidden_dims"]),
            activation=str(model_cfg["activation"]),
            dropout=float(model_cfg["dropout"]),
        )
    elif model_type == "column_conv":
        model = ExportedColumnConv(
            input_dim=input_dim,
            output_dim=output_dim,
            hidden_channels=list(model_cfg["hidden_channels"]),
            kernel_size=int(model_cfg["kernel_size"]),
            activation=str(model_cfg["activation"]),
            dropout=float(model_cfg["dropout"]),
        )
    else:
        raise KeyError(f"Unsupported model type {model_type!r}.")
    model.load_state_dict(blob["state_dict"])
    model.eval()
    model.requires_grad_(False)
    return model


def _get_model(model_path: str | Path | None = None) -> tuple[dict[str, Any], nn.Module]:
    global _MODEL_CACHE
    blob = _load_blob(model_path)
    resolved = blob["_model_path"]
    if _MODEL_CACHE is None or _MODEL_CACHE[0] != resolved:
        _MODEL_CACHE = (resolved, _build_model(blob))
    return blob, _MODEL_CACHE[1]


def feature_names(model_path: str | Path | None = None) -> list[str]:
    return list(_load_blob(model_path)["feature_names"])


def target_names(model_path: str | Path | None = None) -> list[str]:
    return list(_load_blob(model_path)["target_names"])


def _standardize_features(features: np.ndarray, blob: dict[str, Any]) -> np.ndarray:
    mean = np.asarray(blob["feature_scaler"]["mean"], dtype=np.float32)
    std = np.asarray(blob["feature_scaler"]["std"], dtype=np.float32)
    return ((features - mean) / std).astype(np.float32)


def _restore_targets(normalized: np.ndarray, blob: dict[str, Any]) -> np.ndarray:
    restored = normalized.astype(np.float32, copy=True)
    transforms = blob["target_transforms"]
    for index, transform in enumerate(transforms):
        restored[..., index] = restored[..., index] * float(transform["std"]) + float(transform["mean"])
        if transform["transform"] == "standardize":
            continue
        if transform["transform"] == "log1p_standardize":
            restored[..., index] = np.expm1(restored[..., index])
            continue
        raise KeyError(f"Unsupported target transform {transform['transform']!r}.")
    return restored.astype(np.float32)


def predict_emulator(features: np.ndarray, model_path: str | Path | None = None) -> np.ndarray:
    blob, model = _get_model(model_path)
    values = np.asarray(features, dtype=np.float32)
    model_type = blob["model"]["type"]

    if model_type == "mlp":
        single_sample = values.ndim == 1
        if single_sample:
            values = values[None, :]
        if values.shape[-1] != len(blob["feature_names"]):
            raise ValueError(
                f"Expected the last feature dimension to have length {len(blob['feature_names'])}, "
                f"got {values.shape[-1]}."
            )
    elif model_type == "column_conv":
        single_sample = values.ndim == 2
        if single_sample:
            values = values[None, :, :]
        if values.ndim != 3:
            raise ValueError(
                "Column-convolution export expects shape (z, feature) or (batch, z, feature)."
            )
        if values.shape[-1] != len(blob["feature_names"]):
            raise ValueError(
                f"Expected the last feature dimension to have length {len(blob['feature_names'])}, "
                f"got {values.shape[-1]}."
            )
    else:
        raise KeyError(f"Unsupported model type {model_type!r}.")

    standardized = _standardize_features(values, blob)
    with torch.no_grad():
        prediction = model(torch.from_numpy(standardized)).cpu().numpy()
    restored = _restore_targets(prediction, blob)
    if single_sample:
        return restored[0]
    return restored


def predict_from_named_features(
    feature_mapping: dict[str, float | np.ndarray],
    model_path: str | Path | None = None,
) -> dict[str, float]:
    blob = _load_blob(model_path)
    ordered = np.asarray([feature_mapping[name] for name in blob["feature_names"]], dtype=np.float32)
    prediction = predict_emulator(ordered, model_path=model_path)
    return {name: float(prediction[index]) for index, name in enumerate(blob["target_names"])}


def predict_column_from_named_features(
    feature_mapping: dict[str, np.ndarray],
    model_path: str | Path | None = None,
) -> dict[str, np.ndarray]:
    blob = _load_blob(model_path)
    ordered = np.stack(
        [np.asarray(feature_mapping[name], dtype=np.float32) for name in blob["feature_names"]],
        axis=-1,
    )
    prediction = predict_emulator(ordered, model_path=model_path)
    return {name: np.asarray(prediction[..., index], dtype=np.float32) for index, name in enumerate(blob["target_names"])}


class PythonBackendEmulator:
    def __init__(self, model_file: str | Path | None = None):
        self.model_file = model_file
        self._blob = _load_blob(model_file)
        self.input_names = tuple(self._blob["feature_names"])
        self.output_names = tuple(self._blob["target_names"])

    def infer(self, inputs: dict[str, np.ndarray], outputs: dict[str, np.ndarray]) -> None:
        ordered = np.stack([np.asarray(inputs[name], dtype=np.float32) for name in self.input_names], axis=-1)
        prediction = predict_emulator(ordered, model_path=self.model_file)
        for index, name in enumerate(self.output_names):
            outputs[name][...] = np.asarray(prediction[..., index], dtype=outputs[name].dtype)


def create_emulator(config: dict[str, Any]) -> PythonBackendEmulator:
    model_path = config.get("model_path") or None
    return PythonBackendEmulator(model_path)
