"""Declarative emulator YAML configuration. Requires Python 3.10+ and PyYAML.

Run directly to write the example: python emulator_config.py [output.yaml]
Relative paths are preserved and interpreted by the C++ reader relative to YAML.
Field order is preserved: LibTorch uses this order for arguments and results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import yaml


@dataclass(kw_only=True)
class Field:
    name: str
    dimensions: Sequence[str]

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "dimensions": list(self.dimensions),
        }


@dataclass(kw_only=True)
class Model:
    backend: str
    path: str | Path = ""
    input_layout: str = "FeatureMajor"
    output_layout: str = "FeatureMajor"
    options: Mapping[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict:
        backend = self.backend.lower()
        if backend not in {"python", "libtorch", "stub"}:
            raise ValueError(f"Unknown backend: {self.backend!r}")
        for layout in (self.input_layout, self.output_layout):
            if layout not in {"FeatureMajor", "SampleMajor"}:
                raise ValueError(f"Unknown layout: {layout!r}")
        if any(
            not isinstance(k, str) or not isinstance(v, str)
            for k, v in self.options.items()
        ):
            raise TypeError("Backend option keys and values must be strings")
        if backend == "python" and not self.options.get("python_module", "").strip():
            raise ValueError("Python backend requires options.python_module")
        path = str(self.path)
        if backend == "libtorch" and not path.strip():
            raise ValueError("LibTorch backend requires a model path")
        return {
            "path": path,
            "input_layout": self.input_layout,
            "output_layout": self.output_layout,
            "backend": backend,
            "options": dict(self.options),
        }


@dataclass(kw_only=True)
class EmulatorConfig:
    name: str
    model: Model
    dimensions: Sequence[str]
    inputs: Sequence[Field]
    outputs: Sequence[Field]
    type: str = "component"

    def to_dict(self) -> dict:
        if not self.name.strip():
            raise ValueError("Emulator name must not be empty")
        if isinstance(self.dimensions, str):
            raise TypeError("dimensions must be a sequence of dimension names")
        known = set(self.dimensions)
        if len(known) != len(self.dimensions):
            raise ValueError("Duplicate sample-space dimensions")
        if any(not isinstance(d, str) or not d.strip() for d in known):
            raise ValueError("Dimension names must be nonempty strings")

        def fields_to_list(fields: Sequence[Field]) -> list[dict]:
            result = []
            names = set()
            for item in fields:
                if not item.name.strip() or item.name in names:
                    raise ValueError(f"Empty or duplicate field name: {item.name!r}")
                if isinstance(item.dimensions, str):
                    raise TypeError(
                        f"Field {item.name!r}: dimensions must be a sequence"
                    )
                if not item.dimensions:
                    raise ValueError(f"Field {item.name!r} needs dimensions")
                unknown = set(item.dimensions) - known
                if unknown:
                    raise ValueError(
                        f"Field {item.name!r}: unknown dimensions {unknown}"
                    )
                names.add(item.name)
                result.append(item.to_dict())
            return result

        def fields_to_dict(fields: Sequence[Field]) -> dict:
            result = {}
            for item in fields:
                if item.name in result:
                    raise ValueError(f"Duplicate field name: {item.name!r}")
                result[item.name] = item.to_dict()
            return result

        return {
            "emulator": {
                "name": self.name,
                "type": self.type,
                "model": self.model.to_dict(),
            },
            "sample_space": {"dimensions": list(self.dimensions)},
            "input_order": [item.name for item in self.inputs],
            "output_order": [item.name for item in self.outputs],
            "inputs": fields_to_dict(self.inputs),
            "outputs": fields_to_dict(self.outputs),
        }

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.to_dict(), sort_keys=False, allow_unicode=True)

    def write(self, filename: str | Path) -> Path:
        """Write YAML, replacing filename if present; parent must already exist."""
        destination = Path(filename)
        destination.write_text(self.to_yaml(), encoding="utf-8")
        return destination


@dataclass(kw_only=True)
class CoupledField:
    """One coupling: sources export field, destinations import it."""

    field: Field
    kind: str
    sources: Sequence[EmulatorConfig]
    destinations: Sequence[EmulatorConfig]
    merge_type: str = "direct"

    def to_dict(self) -> dict:
        if not self.kind or not self.merge_type:
            raise ValueError("Coupling kind and merge_type must be nonempty")
        if not self.sources or not self.destinations:
            raise ValueError(
                f"{self.field.name}: at least one source and destination required"
            )
        for label, components in (
            ("sources", self.sources),
            ("destinations", self.destinations),
        ):
            names = [component.name for component in components]
            if any(not name.strip() for name in names) or len(set(names)) != len(names):
                raise ValueError(f"{self.field.name}: invalid or duplicate {label}")
        name = self.field.name
        for role, components in (
            ("outputs", self.sources),
            ("inputs", self.destinations),
        ):
            for component in components:
                matches = [f for f in getattr(component, role) if f.name == name]
                if len(matches) != 1:
                    raise ValueError(
                        f"{component.name}: expected exactly one {role} field named {name!r}"
                    )
                actual = matches[0]
        return {
            "kind": self.kind,
            "name": name,
            "merge_type": self.merge_type,
            "sources": [e.name for e in self.sources],
            "destinations": [e.name for e in self.destinations],
        }


@dataclass(kw_only=True)
class CouplingPlan:
    fields: Sequence[CoupledField]

    def to_dict(self) -> dict:
        fields = {}
        # Different objects with the same name can silently refer to different configs.
        components = {}
        for item in self.fields:
            for component in (*item.sources, *item.destinations):
                if (
                    component.name in components
                    and components[component.name] is not component
                ):
                    raise ValueError(
                        f"Multiple emulator objects named {component.name!r}"
                    )
                components[component.name] = component
            name = item.field.name
            if name in fields:
                raise ValueError(f"Duplicate coupled field: {name!r}")
            fields[name] = item.to_dict()
        return {"fields": fields}

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.to_dict(), sort_keys=False, allow_unicode=True)

    def write(self, filename: str | Path) -> Path:
        destination = Path(filename)
        destination.write_text(self.to_yaml(), encoding="utf-8")
        return destination


if __name__ == "__main__":
    temperature = Field(
        name="temperature",
        dimensions=["column"],
    )
    soil_moisture = Field(
        name="soil_moisture",
        dimensions=["column", "soil_levels"],
    )
    surface_flux = Field(
        name="surface_flux",
        dimensions=["column"],
    )

    forcing = EmulatorConfig(
        name="forcing-emulator",
        model=Model(
            backend="libtorch", path="forcing_emulator.pt", options={"dtype": "float64"}
        ),
        dimensions=["column", "soil_levels"],
        inputs=[surface_flux],
        outputs=[temperature, soil_moisture],
    )
    land = EmulatorConfig(
        name="lnd-emulator",
        model=Model(
            backend="libtorch", path="land_emulator.pt", options={"dtype": "float64"}
        ),
        dimensions=["column", "soil_levels"],
        inputs=[temperature, soil_moisture],
        outputs=[surface_flux],
    )

    print(forcing.write("forcing_emulator_config.yaml"))
    print(land.write("land_emulator_config.yaml"))

    plan = CouplingPlan(
        fields=[
            CoupledField(
                field=temperature, kind="state", sources=[forcing], destinations=[land]
            ),
            CoupledField(
                field=soil_moisture,
                kind="state",
                sources=[forcing],
                destinations=[land],
            ),
            CoupledField(
                field=surface_flux, kind="flux", sources=[land], destinations=[forcing]
            ),
        ]
    )
    plan.write("test_coupling_plan.yaml")
