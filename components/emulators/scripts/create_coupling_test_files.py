from create_emulator_config_files import (
    CoupledField,
    CouplingPlan,
    EmulatorConfig,
    Field,
    Model,
)

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


def test_harness():
    land = EmulatorConfig(
        name="lnd-emulator",
        model=Model(
            backend="libtorch", path="land_emulator.pt", options={"dtype": "float64"}
        ),
        dimensions=["column", "soil_levels"],
        inputs=[temperature, soil_moisture],
        outputs=[surface_flux],
    )
    no_model = Model(backend="", path="")
    cpl7 = EmulatorConfig(
        name="cpl7",
        model=no_model,
        dimensions=[],
        inputs=[surface_flux],
        outputs=[temperature,soil_moisture],
    )
    plan = CouplingPlan(
        fields=[
            CoupledField(
                field=temperature, kind="state", sources=[cpl7], destinations=[land]
            ),
            CoupledField(
                field=soil_moisture,
                kind="state",
                sources=[cpl7],
                destinations=[land],
            ),
            CoupledField(
                field=surface_flux, kind="flux", sources=[land], destinations=[cpl7]
            ),
        ]
    )
    plan.write("test_harness_plan.yaml")

if __name__ == "__main__":
    test_emulator_coupling()
    test_harness()
