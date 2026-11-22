import sys
from pathlib import Path

from emulator_config_api import (
    CoupledField,
    CouplingPlan,
    EmulatorConfig,
    Field,
    Model,
    PrognosticField,
    prognostic_fields,
)

filter_columns = Field(name="filter_columns", dimensions=["column"], dtype="int32")
temperature = Field(name="temperature", dimensions=["column"])
soil_moisture = Field(name="soil_moisture", dimensions=["column", "soil_levels"])
surface_flux = Field(name="surface_flux", dimensions=["column"])

# Output directory; the build passes its fixture dir. Defaults to the cwd.
fixture_dir = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()


def test_emulator_coupling():
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
        inputs=[temperature, soil_moisture, filter_columns],
        outputs=[surface_flux],
    )

    print(forcing.write(f"{fixture_dir}/forcing_emulator_config.yaml"))
    print(land.write(f"{fixture_dir}/land_emulator_config.yaml"))

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
    plan.write(f"{fixture_dir}/test_coupling_plan.yaml")


def test_harness():
    land = EmulatorConfig(
        name="lnd-emulator",
        model=Model(
            backend="libtorch", path="land_emulator.pt", options={"dtype": "float64"}
        ),
        dimensions=["column", "soil_levels"],
        inputs=[temperature, soil_moisture, filter_columns],
        outputs=[surface_flux],
    )
    no_model = Model(backend="", path="")
    external = EmulatorConfig(
        name="ext",
        model=no_model,
        dimensions=[],
        inputs=[surface_flux],
        outputs=[temperature, soil_moisture],
    )
    plan = CouplingPlan(
        fields=[
            CoupledField(
                field=temperature, kind="state", sources=[external], destinations=[land]
            ),
            CoupledField(
                field=soil_moisture,
                kind="state",
                sources=[external],
                destinations=[land],
            ),
            CoupledField(
                field=surface_flux, kind="flux", sources=[land], destinations=[external]
            ),
        ]
    )
    plan.write(f"{fixture_dir}/test_harness_plan.yaml")


def test_prognostic():
    state = Field(name="state", dimensions=["column"])
    emulator = EmulatorConfig(
        name="prognostic-emulator",
        model=Model(backend="stub"),
        dimensions=["column"],
        inputs=[state],
        outputs=[state],
    )
    external = EmulatorConfig(
        name="ext",
        model=Model(backend="stub"),
        dimensions=[],
        inputs=[state],
        outputs=[state],
    )
    emulator.write(f"{fixture_dir}/prognostic_emulator_config.yaml")
    plan = CouplingPlan(
        fields=[
            PrognosticField(field=item, component=emulator, external=external)
            for item in prognostic_fields(emulator)
        ]
    )
    plan.write(f"{fixture_dir}/test_prognostic_plan.yaml")


if __name__ == "__main__":
    test_emulator_coupling()
    test_harness()
    test_prognostic()
