(omega-design-frazil)=
# Frazil

## 1 Overview

The Frazil module computes the frazil tendencies in the ocean column from the local temperature, salinity, pressure, and layer-thickness state. It is used as a tendency term in the OMEGA model and accumulates completed-step totals for export to the coupler. The module is organized as a registry with a default instance and optional non-default instances, and it exposes per-column per-timestep totals for accumulated freshwater mass, salt mass, and energy.

The implementation is intended to be:
- local to a single ocean column,
- consistent with the thermodynamic assumptions in the EOS,
- conservative with respect to water mass, salt mass and energy
- compatible with the existing tendency infrastructure and enable flag.


## 2 Requirements

### 2.1 Requirement: Compute frazil tendencies from formation and melt

The module must compute the total pseudo-thickness (aka mass), salt and conservative temperature tendency contribution from frazil formation and melt as a function of the water column state.
The calculation must allow for tendency contributions from
- the formation of frazil ice from supercooled seawater (leading to the removal of the supercooled water)
- melt of existing frazil if frazil is formed below encounters a layer warm enough to melt it.
- salt and energy contributions from the above processes.


### 2.2 Requirement: Preserve conservation and physical consistency

The implementation must maintain column-level conservation: the sum of each Frazil tendency in a column should match the corresponding column-integrated Frazil term. it must reject or flag situations where the computed mass, salt or energy budgets cannot be reconciled. The implementation should also treat the 3 terms (mass, salt, and energy) in a self-consistent way, in line with the thermodynamic assumptions of each frazil option.

Incorrect mass or energy partitioning can corrupt the tracer and thickness budgets and lead to a mass, salt, or energy leak in the coupled model.

### 2.3 Requirement: Support configuration-based enable/disable

The frazil capability must be enabled by configuration. When disabled, the default frazil object must not be created and the tendency contribution must be skipped.

This requirement is needed because the model should be able to run without a frazil parameterization without introducing unresolved dependencies or unnecessary computation.

### 2.4 Requirement: Provide per-cell per-timestep accumulated diagnostics

The module must track column-accumulated totals for:
- frazil freshwater mass,
- salt mass,
- energy associated with the mass fluxes

These totals are used to compute completed-step contributions and are needed for coupling to the sea ice model.

### 2.5 Requirement: calculate frazil contributions consistent with the TEOS-10 thermodynamics

Omega uses ocean thermodynamics consistent with TEOS-10. We target to form (and melt) frazil in a way consistent with the choice of the TEOS-10 Equation Of State (EOS).

### 2.6 Requirement: Support IO field registration and output

The per-timestep column-wise frazil totals should be exposed as fields so they can be selected in output streams. The fields must be attached to the correct arrays and registered under a group.

This is needed because we desire to test the frazil contributions in standalone mode before coupling. This will also be useful for diagnostics of production simulations.

### 2.6 Requirement: Include a thickness limit
For stability reasons, the total mass tendency removed (or added) to a layer should be limited to a maximum fraction of the layer thickness. Ten percent is a threshold used in previous implementations (and the intended default value). The requirement is for there to be a parameter that caps the mass tendency.

### 2.7 Desired: Include a depth limit
In other implementations, deep formation of frazil (e.g. below 1000m) was an issue. Thus we prefer to have a parameter restricting the frazil formation (and melt) to depths shallower than this parameter.

### 2.7 Desired: have a simpler frazil scheme
To limit the risk of instability, we desire to implement a frazil scheme to provide a simpler alternative to the TEOS-10 option. As much as possible, this simpler frazil scheme should rely on assumptions that past frazil implementations have made to increase the chances of providing a stable frazil scheme in a coupled model.

### 2.8 Desired: limit computational expense
We desire to limit unnecessary calls and computations. Thus, we add a check on the water column state (is there super cooled water in this column?) before running through the frazil calculation. This is inspired by MPAS-O implementation.

## 3 Algorithmic Formulation

The frazil calculation is implemented as a local column process. The model state provides:
- conservative temperature $T$,
- absolute salinity $S$,
- pressure $P$,
- layer thickness $H$.

For a given active-layer cell, the module computes the formation and melt contributions to:
- ice mass,
- liquid mass,
- salt mass,
- total energy.

The calculation is organized around a set of vertically-accumulated quantities:
$$
M_{ice},\quad M_{liq},\quad M_{salt},\quad E_{ice},\quad E_{liq}
$$
and the tendency terms
$$
\frac{\partial H}{\partial t},\quad
\frac{\partial T}{\partial t},\quad
\frac{\partial S}{\partial t}.
$$

For each stage of a timestep, the calculation resets the (vertical) accumulator arrays, applies the local frazil algorithm, and stores the resulting stage-local contributions in the tendency arrays. At the end of the completed outer-step update, the accumulated totals are normalized by the full ocean timestep:
$$
\bar{F}_{mass} = \frac{\sum_k w_k M_{mass,k}}{\Delta t},
\quad
\bar{F}_{salt} = \frac{\sum_k w_k M_{salt,k}}{\Delta t},
\quad
\bar{F}_{energy} = \frac{\sum_k w_k E_k}{\Delta t},
$$
where $w_k$ are the timestepper-specific update weights. This preserves the correct outer-step mean flux while remaining compatible with Runge-Kutta stage logic.
We prefer that the public output arrays captures mean fluxes rather than stage or timestep totals to make the Polaris testing independent of the timestep used.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

The Frazil configuration is stored under the `Frazil` config group. Relevant entries include:
- `FrazilType`: fixed-property or TEOS-based formulation
- `LayerMassFracMax`: maximum total frazil mass changed allowed in an ocean layer
- `Phi`: target liquid mass fraction for the formed frazil
- `DepthLimit`: optional depth limit below which frazil is disabled
- `ConservationCheck`: whether conservation must be validated

These parameters are read during object creation in `Frazil.cpp`.

#### 4.1.2 Class/structs/data types

The relevant public structures in the Frazil implementation are:

```c++
enum class FrazilType {
   FixedPropertyFrazil,
   TeosFrazil
};

struct FrazilFormation { ... };
struct FrazilMelt { ... };
struct FixedPropertyFrazilFormation { ... };
struct FixedPropertyFrazilMelt { ... };

class Frazil {
 public:
   static void init();
   static Frazil *getDefault();
   static Frazil *get(const std::string &Name);
   static void clear();
   static void erase(std::string InName);

   Array2DReal FrazilTTend;
   Array2DReal FrazilSTend;
   Array2DReal FrazilHTend;

   Array1DReal AccMIce;
   Array1DReal AccEIce;
   Array1DReal AccMLiq;
   Array1DReal AccELiq;
   Array1DReal AccMSalt;

   Array1DReal FrazilMassFlux;
   Array1DReal FrazilSaltFlux;
   Array1DReal FrazilEnergyFlux;

   void computeFrazil(...);
   void resetOcnStepFluxes();
   void accumulateOcnStepFluxes(Real FinalUpdateWeight, R8 TimeStepSeconds);
   void registerFields();
   void unregisterFields();
};
```

The design keeps the arrays and accumulation logic local to the Frazil object, which makes the lifetime management explicit and avoids leaking the data model into time-stepping logic.

### 4.2 Methods

#### 4.2.1 Initialization

The static `Frazil::init()` method:
- ensures the mesh and vertical coordinate exist,
- initializes EOS dependencies,
- checks the `Tendencies` config for `FrazilTendencyEnable`,
- creates the default Frazil instance when enabled.

This ensures the object exists only when the feature is active.

#### 4.2.2 Creation and retrieval

The creation pattern mirrors the rest of the Omega module structure:
- a map of all Frazil objects,
- a pointer to the default instance,
- a name-based retrieval method,
- methods for erase and clear.

This design allows future non-default instances while keeping the default case fast and simple.

#### 4.2.3 Computation

The `computeFrazil()` method is the core thermodynamic routine. It:
- reads the local tracer state,
- evaluates the local freezing point,
- computes mass and energy changes associated with frazil formation or melt,
- writes the results into the layer-by-layer tendency arrays.

The method is later invoked from the tendency operator, which is the integration point between the frazil module and the model tendency machinery.

#### 4.2.4 Completed-step accumulation

The `accumulateOcnStepFluxes()` method updates the final per-cell ocean-step fluxes. It:
- reads the stage-local accumulation arrays,
- applies the update weight,
- divides by the full ocean timestep seconds,
- flips the sign, because the accumulators are counted positive into the frazil reservoir while the exported fluxes follow the coupler convention of being positive into the ocean,
- writes the mean flux into the public output arrays.

This is a key design choice because it keeps the output/coupling values independent of the internal stage ordering of the timestepper.

#### 4.2.5 Field registration and output

The `registerFields()` method creates the three default output fields:
- `FrazilMassFlux`
- `FrazilSaltFlux`
- `FrazilEnergyFlux`

All three are counted positive into the ocean, so frazil formation gives
negative values and net melt gives positive ones.

Each field is:
- 1D with `NCells` dimension,
- attached to the corresponding array,
- inserted into the `Frazil` field group.

The `unregisterFields()` method is guarded by a registration flag so only the object that successfully registered the fields destroys them. This avoids collisions when non-default instances are later introduced.

#### 4.2.6 Tendency hook

The `FrazilOnCell` object is the bridge between the Frazil module and the tendency infrastructure. It:
- checks whether frazil is enabled,
- resets the stage-local accumulators,
- calls `computeFrazil()`,
- accumulates the completed-step totals,
- writes the resultant tendency contributions into the thickness and tracer tendency arrays.

This is the integration point used by `TendencyTerms.cpp`.

## 5 Verification and Testing

### 5.1 Ctest - Test cold formation

The cold-case frazil formation test verifies that in a single layer set to super-cooled conditons, a call to FrazilFormation() produces:
- accumulated frazil ice mass is positive,
- accumulated frazil salt mass is positive,
- accumulated frazil energy is negative,
- layer temperature tendency is positive,
- layer salinity and thickness tendencies are negative.

This verifies the sign and directionality of the formation branch for a single layer.
The `testFrazilFormationCold` verifies the Teos-10 implementation, and `testFixedPropertyFrazilFormationCold` verifies the FixedProperty implementation.

### 5.2 Ctest - Test warm formation

The warm-case formation test verifies that the formation branch is zero when water is not supercooled. This checks that the algorithm respects the physical cutoff condition.
The `testFrazilFormationWarm` verifies the Teos-10 implementation, and `testFixedPropertyFrazilFormationWarm` verifies the FixedProperty implementation.

### 5.4 Ctest - Test column conservation and formation/melt switching

The column test verifies the full `Frazil::computeFrazil()` path across a single water column with multiple active layers.  The test constructs a vertical profile with cold layers separated by warm layers so that frazil can form in supercooled water, be carried upward through the column accumulators, and melt again where warmer water is encountered.

The test checks the branch switching through the signs of the layer tendencies:
- cold layers should form frazil, giving negative thickness and salinity tendencies and a positive temperature tendency;
- warm layers with accumulated frazil should produce melt, giving the opposite tendency signs.

The test also forcibly enables the internal frazil column conservation check, which verifies that the summed layer tendencies balance the accumulated frazil mass, salt, and energy reservoirs.  This provides coverage of the conservation relationship inside one call to `computeFrazil()`, while the separate timestep-level tests cover completed-step output and coupling diagnostics.

This case has also been used for manual checks by adding logging commands (which work on CPU). The checks on layer accummulation, left over after limited melt etc. are not part of the automated checks.

### 5.5 Test depth limit

The depth-limit test ensures that any layers below the configured depth threshold generate frazil tendencies of zero, which guards the parameterized cutoff behavior.

### 5.6 Test mass limit

The mass-limit tests verify that `LayerMassFracMax` caps the frazil mass change to a
fixed fraction of the layer thickness.  The tests are single layer, use the default value `LayerMassFracMax = 0.10` and exercise both available frazil formulations.

For formation, the initial conditions are set to produce a large frazil mass (i.e. small layer of 1m pseudothickness and very super-cooled water), which should be capped. `testFrazilFormationMassLimit()` verifies the TEOS-10 path and
`testFixedPropertyFrazilFormationMassLimit()` verifies the fixed-property path.
The TEOS-10 test checks that the produced solid ice mass, the `Phi`-dependent liquid
mass, and thickness tendency against the analytical capped values.  The fixed-property test checks that the produced ice thickness and thickness tendency are equal to the
expected capped value.

For melt, the initial ocean conditions are set to allow melting:  the small layer of 1m pseudothickness is warm water and a large amount of frazil terms is present for melting, which should be capped. `testFrazilMeltMassLimit()` verifies the TEOS-10 path and `testFixedPropertyFrazilMeltMassLimit()` verifies the fixed-property path.
These tests check that the computed fraction of frazil to be melted is bounded between zero and one,
that the 3 frazil tendencies match the expected values, and that the remaining frazil
reservoirs are reduced consistently.

### 5.7 Test Phi sensitivity

The `testFrazilFormationPhi()` test verifies that the TEOS-10 formation path
responds to the `Phi` parameter.  The test runs the same scalar state with
`Phi = 0.75` and `Phi = 0.85`.  With all other inputs unchanged, increasing
`Phi` should increase the total frazil mass, increase the salt content, and
make the total frazil energy more negative.  The test also checks that the
solid ice mass is unchanged so the comparison is not accidentally testing the
mass limiter.

### 5.8 Test per-timestep conservation

A future test should check that the ocean tendencies and the terms ready for coupling export are conserved. The above test 5.4 ensures conservation within the frazil call but does not check the full timestep totals, which will differ across time stepper choice.
