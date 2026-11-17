Frazil ice formation in MPAS-Ocean
=================================

date: 2026/11/16

Contributors: Alice Barthel, Co-pilot, others? 


Summary
-------

During each (inner) ocean timestep, the presence of super-cooled water in 
MPAS-Ocean triggers frazil formation (and possibly melting). 
MPAS-Ocean estimates freezing or melting energy
in ocean layers, forms frazil from supercooled water, permits some newly formed
frazil to melt in warmer layers above it, and accumulates the surviving
mass, salt, and energy for coupling. This document describes the behavior in
the current source, and the requirements. The aim of this document is to define 
the acceptable level of documentation. Ideally, no code changes would be necessary 
unless the physics cannot be reconciled under any assumption/documentation. 

The implementation is in
:file:`../../src/shared/mpas_ocn_frazil_forcing.F`; its namelist options are
defined in :file:`../../src/Registry.xml`.


Requirements
------------

Requirement 1: Frazil formation removes ocean supercooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Frazil formation adds heat to the ocean to remove supercooling in ocean
layers. If the implementation does not bring a layer to its freezing point,
the reason and resulting behavior must be documented, including any effects of
multiphase time steppers.

.. Current implementation: The code computes a layer freezing potential from
.. layer thickness, seawater heat capacity and density, and the difference
.. between layer temperature and freezing temperature. It converts negative
.. potential to freezing energy and calculates frazil thickness from latent heat
.. of fusion. The thickness is capped by a configured fraction of the layer, so
.. the implementation does not guarantee that the layer reaches its freezing
.. point when this cap is active. The temperature tendency uses latent energy
.. plus, by default, a sensible-heat adjustment.

Requirement 2: Consistent treatment of temperature, salinity, and volume
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assumptions linking temperature, salinity, and volume changes during frazil
formation and melting must be documented. Ideally, their tendencies must be mutually
consistent. 

.. Current implementation: Frazil formation decreases layer thickness using
.. frazil density divided by local ocean density, removes the computed
.. frazil-weighted salt, and warms the ocean using seawater density and heat
.. capacity. The source comments state that local density is used for the
.. thickness tendency to preserve buoyancy during formation. Melting applies the
.. opposite-sign thickness, salt, and temperature tendencies, with melted salt
.. and energy taken proportionally from the running frazil inventory. These
.. separate formulas are the current implementation; their joint conservation
.. properties are not established by the code alone.

Requirement 3: Column-integrated frazil fluxes match layer tendencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Column-integrated frazil fluxes must equal the sum of layer tendencies for
mass, salt, and energy.

.. Current implementation: The routine sums surviving frazil thickness, weighted
.. salt, and energy over the water-column pass and adds them to separate
.. accumulators for mass, salt, and energy. These accumulators are used to form
.. coupling fluxes. The implementation should conserve but this needs to be verified.

Requirement 4: Warm layers can melt frazil formed in layers below
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Frazil formed in a supercooled layer must be able to melt when a warmer layer
above it provides the energy for melting. Formation and melting must treat
mass, salt, and energy consistently; for example, mass must not be removed
without the corresponding salt exchange. Column-integrated quantities must
remain in balance.

.. Current implementation: Layers are traversed from the deepest eligible layer
.. toward the surface. During this pass, positive melting energy can melt the
.. frazil formed in deeper layers earlier in the same pass. Melting is limited by
.. the amount in this running, newly formed inventory and by the configured
.. fraction of the warm layer thickness. Previously accumulated frazil awaiting
.. coupling is not included in this meltable inventory.

Requirement 5: Compatibility with MPAS-SeaIce external frazil
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MPAS-Ocean implementation must remain compatible with the MPAS-SeaIce
``external`` frazil option (as well as the ``omega-fluxes`` option).

Current implementation: MPAS-Ocean exports :code:`Fioo_q`,
:code:`Fioo_frazilm`, :code:`Fioo_frazils`, and :code:`Fioo_frazilh`. The
current code does not define a :code:`Fioo_frazil` field. In the non-
:code:`omega-fluxes` MPAS-SeaIce coupling path, sea ice recalculates a frazil
mass from the energy potential and surface salinity, then returns the
difference from the received mass as :code:`frazilMassAdjust` in a subsequent
ocean flux exchange. The meaning, timing, and conservation behavior of this
correction should be described alongside the implementation.


Desired Outcomes
----------------

Desired 1: Frazil formation does not change buoyancy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The inline code comments explicitly notes the use of the local ocean density in the thickness tendency 
 as intended to preserve buoyancy during formation. However, given that other buoyancy tracers are modified in the timestep, I am not sure how essential this buoyancy conservation is to model stability. 
 This choice competes with the idea of treating all tracers consistently (R2). By construction, meeting D1 means D2 cannot be met to machine precision unless the definition of ocean salt content is modified. 

Desired 2: Salinity is unchanged when frazil salinity equals ocean salinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a conceptual model, removing a volume of ocean water to form frazil should not change the salinity if the frazil salinity equals the ocean salinity. 
Evaluate the expected salinity tendency when :math:`S_{\mathrm{frazil}} = S`
and document whether the implementation satisfies this expectation.

.. Current implementation: The salt removed with new frazil is computed from its
.. thickness, frazil density, assigned frazil salinity, and the salinity-unit
.. conversion factor. The layer-thickness tendency removes the corresponding
.. frazil mass using local ocean density. Whether these paired tendencies yield
.. the expected salinity behavior in the time integration should be verified.

Desired 3: Implement a depth limit for frazil formation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Document the purpose and behavior of the depth limiter, including its
configuration and interaction with the vertical coordinate.

Desired 4: Support an option for variable frazil salinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The current implementation includes the
:code:`config_use_frazil_ice_porosity` and
:code:`config_frazil_ice_porosity` options, which affect the salinity assigned
to internal frazil. This option only changes the salt content of the frazil (not its mass/thickness or energy content) and its limits need to be documented. 

Desired 5: Variable frazil salinity should be compatible with coupling choices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When variable frazil salinity is used, the salt removed from each ocean layer is impacted. 
However, the other coupling choices (eg. constant salinity exchange to sea ice) need to remain possible. 
This is why a check and correction for excess salinity is included in the frazil code. If turned on, this option allows a surface redistribution of the salt that does not meet the coupling constraint. 
This may impact ocean stability but is a stop-gap option to maintain cross-component compatibility. It can be thought of as an instantaneous surface brine rejection scheme. 


Known limitations of current implementation
------------------------------------------

Limit 1: frazil formation is based on "bulk" quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The frazil total mass is calculated based on the energy available divided by the latent heat of fusion of frazil. 
This user-configurable parameter uniquely sets the total mass of frazil removed. The default value (3.34e5) matches the specific latent heat of fusion for pure water.
It is converted to an ice-thickness accounting using the "bulk" frazil density. There is no check on the thermodynamical compatibility of these 2 parameters (frazil density and frazil latent heat).
These are not currently independent of other choices (such as porosity).

Limit 2: two limits to R1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Although the frazil formation calculation is based on energy and designed to remove exactly the super-cooled energy of the previous ocean layer state, there are two limits to Requirement 1. 
In multi-phase timesteppers, the total (outer timestep) tendency will not fully remove the super-cooled energy of the previous ocean state because it uses tendencies based on provisional intermediate states.
The second limit is that if the frazil variable of "ice-thickness" exceeds the 0.1 the ocean layer thickness, the energy removed will not be the total super-cooled energy of the previous layer state.


Limit 3: salt content definition is not consistent with Boussinesq equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Based on the definition of `saltflux_factor`, the salt content of an individual layer is based on the density of freshwater, not the reference density of seawater. 
Thus, the treatment of the salt tracer and the temperature tracer are not aligned since they use different density references. 
The current implementation is self-consistent for each individual tracer, but inconsistent with the Boussinesq approximation of using a single $\rho_0$ to switch from the mass-based theoretical tracer conservation to the volume-based tracer conservation of the implementation.

Limit 4: R2, D1 and D2 are in conflict
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Ideally, all tracer (thickness, salinity, temperature) should be treated in the same way based on the same density reference. 
However, the current choice to define salt and heat content based on different density references creates a conflict. Instead, salinity and temperature are each treated consistently with their definitions, but not consistent with each other (see Limit 3)..
In addition, giving priority to D1 means that D2 cannot be exactly met. 

Limit 5: use_porosity is not a full porosity-based calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In its current implementation, the `use_porosity` option can be misunderstood. Toggling this option does not change the amount of frazil formed, it only serves as an option to adjust the salinity based on a factor (`porosity` aka $\phi$) times the local layer salinity. 
It is used as a mass fraction, not a volume-based porosity. 
The current implementation meets D4 by adjusting the salinity based on this parameter, but it does not fully implement a porosity-based calculation for frazil formation.

Algorithm Design
----------------

Algorithm Design: Frazil formation and melt tendencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The forward-mode initialization sets the module's formation switch from
:code:`config_use_frazil_ice_formation`. When enabled, the forward timestep
builds frazil tendency arrays before the layer-thickness and tracer tendency
calculations apply them. If ecosystem tracers and
:code:`config_use_ecosys_frazil_tendencies` are enabled, a separate routine
computes their frazil tendencies.

For each cell, the routine uses land-ice and floating-land-ice masks together
with :code:`config_frazil_in_open_ocean` and
:code:`config_frazil_under_land_ice` to decide whether to process it. Grounded
land-ice cells are skipped. The eligible vertical range is bounded by the
layer satisfying
:code:`-zMid(k,iCell) < -ssh(iCell) + config_frazil_maximum_depth`. If the
minimum temperature in that range exceeds
:code:`config_frazil_maximum_freezing_temperature`, the column is skipped.

The remaining layers are visited from the deepest eligible layer to the
surface. At each layer, the freezing temperature is calculated from salinity
and pressure, with the land-ice-cavity flag passed to the equation-of-state
routine. The code defines the signed freezing/melting potential as

.. math::

   P_k = h_k c_{p,w} \rho_{sw} (T_k - T_{f,k}),

where :math:`h_k` is layer thickness, :math:`T_k` is layer temperature,
:math:`T_{f,k}` is its freezing temperature, :math:`c_{p,w}` is configured
seawater specific heat, and :math:`\rho_{sw}` is the seawater reference
density. Negative :math:`P_k` is treated as freezing energy
:math:`E_{f,k}=\max(0,-P_k)`; positive :math:`P_k` is melting energy
:math:`E_{m,k}=\max(0,P_k)`.

For freezing, the unbounded frazil thickness is
:math:`E_{f,k}/(L_f\rho_i)`, where :math:`L_f` is the configured heat of
fusion and :math:`\rho_i` is the configured frazil density. The actual
thickness is capped at :math:`\alpha h_k`, where :math:`\alpha` is
:code:`config_frazil_fractional_thickness_limit`. The code then adds
temperature, salinity, and layer-thickness tendencies using that capped
thickness. If sensible heat is enabled, the energy assigned to the frazil
includes a correction based on freezing temperature, seawater heat capacity,
reference density, and local ocean density; otherwise it uses latent heat
alone.

For a formed frazil thickness :math:`\delta h_f`, the current tendencies are
computed from

.. math::

   M_{s,f} &= \delta h_f \rho_i S_f c_s,\\
   \mathcal{T}_{h,k}^{frazil} &= -\frac{\delta h_f \rho_i}{\rho_k \Delta t},\\
   \mathcal{T}_{S,k}^{frazil} &= -\frac{M_{s,f}}{\Delta t},\\
   \mathcal{T}_{T,k}^{frazil} &= \frac{E_f^{ice}}{\rho_{sw}c_{p,w}\Delta t},

where :math:`S_f` is assigned frazil salinity, :math:`c_s` is the code's
salinity-unit conversion factor, :math:`\rho_k` is local ocean density, and
:math:`\Delta t` is the ocean timestep. With sensible heat enabled, the
energy assigned to frazil is

.. math::

   E_f^{ice} =
   \delta h_f \rho_i L_f -
   \rho_{sw}c_{p,w}T_{f,k}
   \frac{\delta h_f\rho_i}{\rho_k};

   E_f^{ice} = \delta h_f \rho_i L_f
   \quad \text{when sensible heat is disabled.}

For melting, the thickness derived from melting energy is capped both by the
surviving frazil formed earlier in the current bottom-to-surface pass and by
:math:`\alpha h_k`. The melted fraction removes the same fraction of the
running frazil salt and energy inventories. Melting increases layer thickness,
releases salt, and cools the ocean. The routine does not use frazil
accumulated during earlier coupling intervals as input to this melt
calculation.

For melted thickness :math:`\delta h_m`, the code applies the same
mass-to-volume conversion with the opposite sign, releases the proportional
salt inventory, and applies a cooling tendency from the proportional frazil
energy inventory. This proportional removal is based on the melted fraction
of the running newly formed frazil inventory, not on a separately stored
vertical frazil field. If :math:`f_m` is the melted fraction of that running
inventory, the corresponding tendency contributions are

.. math::

   \mathcal{T}_{h,k}^{melt} &= \frac{\delta h_m \rho_i}{\rho_k \Delta t},\\
   \mathcal{T}_{S,k}^{melt} &= \frac{f_m M_{s,f}^{running}}{\Delta t},\\
   \mathcal{T}_{T,k}^{melt} &= -\frac{f_m E_f^{ice,running}}{\rho_{sw}c_{p,w}\Delta t}.

These are the tendency-array contributions computed by the frazil routines;
they should not be read as a proof of column-integrated conservation after
the tracer and thickness time integration.

Algorithm Design: Frazil-ice porosity and passive tracers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Initialization sets an internal switch to zero unless
:code:`config_use_frazil_ice_porosity` is true, in which case it is set to
one. Denote this switch by :math:`u` and
:code:`config_frazil_ice_porosity` by :math:`\phi`. The salinity assigned to
new frazil is capped at the ocean salinity and is calculated as

.. math::

   S_f =
   \begin{cases}
      \min(S, u\phi S), & \text{under land ice},\\
      \min(S, u\phi S + (1-u)S_{ref}), & \text{in the open ocean},
   \end{cases}

where :math:`S` is the local ocean salinity and :math:`S_{ref}` is
:code:`config_frazil_sea_ice_reference_salinity`. Thus, when porosity is
disabled, open-ocean frazil uses the reference salinity and frazil under land
ice uses zero salinity. When enabled, internal frazil salinity is based on
the local ocean salinity and configured porosity.

When :code:`config_expel_excess_frazil_salt` is enabled, the code compares
the internal salt inventory with the salt implied by the fixed coupling
salinity. If the internal inventory is higher, the difference is added to the
surface-layer salinity tendency. If it is lower, the difference is recorded
in :code:`lowSalinityFrazilIce` and the accumulated frazil salt remains at the
computed internal value; the code does not add salt to make up the shortfall.
This correction is applied only when the porosity switch is enabled.

If ecosystem frazil tendencies are enabled, the separate ecosystem routine
removes a porosity-scaled part of each ecosystem tracer during formation,
returns a proportional part during melting, and puts the remaining
column-integrated tracer content into the surface layer. The source comments
state that sea ice calculates its own frazil ecosystem concentration.


Implementation
--------------

Implementation: Frazil formation and melt tendencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The implementation is in :file:`../../src/shared/mpas_ocn_frazil_forcing.F`.
Frazil formation is controlled by :code:`config_use_frazil_ice_formation`.
The current Registry defaults and controls include:

- :code:`config_use_frazil_ice_formation` (default: :code:`.false.`)
- :code:`config_frazil_in_open_ocean` (default: :code:`.true.`)
- :code:`config_frazil_under_land_ice` (default: :code:`.false.`)
- :code:`config_frazil_maximum_depth` (default: :code:`100.0` m)
- :code:`config_frazil_maximum_freezing_temperature` (default: :code:`0.0` C)
- :code:`config_frazil_fractional_thickness_limit` (default: :code:`0.1`)
- :code:`config_frazil_heat_of_fusion` (default: :code:`3.34e5` J kg\ :sup:`-1`)
- :code:`config_frazil_ice_density` (default: :code:`1000.0` kg m\ :sup:`-3`)
- :code:`config_specific_heat_sea_water` (default: :code:`3985.0` J kg\ :sup:`-1` C\ :sup:`-1`)
- :code:`config_frazil_sea_ice_reference_salinity` (default: :code:`4.0`, units: :code:`1.e-3`)
- :code:`config_include_sensible_heat_frazil` (default: :code:`.true.`)
- :code:`config_expel_excess_frazil_salt` (default: :code:`.true.`)
- :code:`config_frazil_use_surface_pressure` (default: :code:`.true.`)
- :code:`config_use_frazil_ice_porosity` (default: :code:`.false.`)
- :code:`config_frazil_ice_porosity` (default: :code:`1.0`)

When :code:`config_frazil_use_surface_pressure` is enabled, the routine
computes :code:`frazilSurfacePressure` as accumulated frazil mass times
gravity. Frazil formation itself is disabled by default.

The main routine accumulates surviving frazil mass, salt, and energy in state
arrays between coupling exchanges. When the ocean Sea-Ice coupling code
exports the accumulated inventory, it divides each quantity by the ocean
coupling interval and then resets the accumulators. In open-ocean cells,
:code:`Fioo_frazilm`, :code:`Fioo_frazils`, and :code:`Fioo_frazilh` carry the
mass, salt, and energy fluxes. When frazil mass is present,
:code:`Fioo_q` is set to zero; otherwise it carries the non-positive
near-surface melt potential. The field metadata defines :code:`Fioo_q` and
:code:`Fioo_frazilh` in W m\ :sup:`-2`, and the frazil mass and salt fields
in kg m\ :sup:`-2` s\ :sup:`-1`.

For cells under land ice, the formation routine accumulates land-ice frazil
mass and energy separately and sets :code:`frazilIceFreshwaterFlux`; this
freshwater flux is included in the land-ice freshwater flux. The Sea-Ice
coupling export suppresses the open-ocean frazil fields for land-ice cells.
In standalone or data land-ice-flux modes it instead exports the land-ice
frazil mass and energy fields.

Implementation: Frazil-ice porosity and passive tracers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When :code:`config_use_frazil_ice_porosity` is false, :code:`usePorosity` is
zero and internal open-ocean frazil uses the fixed reference salinity. When
it is true, :code:`usePorosity` is one and internal frazil salinity is the
configured porosity times local seawater salinity, capped at the seawater
salinity. The Registry allows :code:`config_frazil_ice_porosity` values from
0 to 1.

The porosity option does not change the default fixed salinity used by the
Sea-Ice coupling. With :code:`config_expel_excess_frazil_salt` enabled, the
routine returns internal salt in excess of the coupling assumption to the
surface layer; a shortfall is recorded in :code:`lowSalinityFrazilIce`.

Ecosystem-tracer tendencies are computed only when ecosystem tracers and
:code:`config_use_ecosys_frazil_tendencies` are enabled. The current routine
removes a porosity-scaled tracer amount with newly formed frazil, returns a
proportional amount if that frazil melts during the column pass, and releases
the remaining column inventory into the surface layer.

In the current MPAS-SeaIce Registry, :code:`config_frazil_coupling_type`
defaults to :code:`external`. For the non-:code:`omega-fluxes` branch, SeaIce
reconstructs frazil mass from the received energy potential and surface
salinity, then sets :code:`frazilMassAdjust` to the difference between
MPAS-Ocean's mass flux and the reconstructed mass flux. On the next coupling
exchange this adjustment is returned to the ocean through the
:code:`Fioi_meltw` freshwater and :code:`Fioi_salt` salt fluxes, using the
Sea-Ice reference salinity. The :code:`omega-fluxes` branch instead uses the
supplied frazil mass, salt, and heat fluxes directly and sets the mass
adjustment to zero.


Testing
-------

Testing and Validation: Frazil formation and melt tendencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Melting applies the opposite-sign thickness, salt, and temperature tendencies, with melted salt
and energy taken proportionally from the running frazil inventory. This ensures that melting applies equally to all 3 frazil quantities. 

TODO: No automated tests check the sign of the tendencies for formation and melting, but it has been manually checked in coupled simulations. 
TODO: No automated checks on the column-integrated tendencies vs final coupler fluxes exist. Nicole checked the energy closure. 

No checks on open-ocean and under-land-ice
behavior, or depth-limit behavior. No explicit checks for the external
Sea-Ice mass correction and for the low-salinity diagnostic.


