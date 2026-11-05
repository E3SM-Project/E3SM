Omega-type frazil coupling with variable salinity between MPAS-Ocean and MPAS-Seaice
====================================================================================

date: 2026/11/05

Contributors: Nicole Jeffery, Andrew Roberts, Darin Comeau, Erin Thomas,
Kara Peterson, Jon Wolfe

This design document describes the changes merged in E3SM PR #8957
(merge commit ``1946bb1d22``).


Summary
-------

In the default (original) E3SM frazil coupling, MPAS-Ocean computes frazil
ice formation and sends the sea ice model a frazil mass flux
(``Fioo_frazil``) and a frazil energy flux. The energy is carried in the same
coupler field (``Fioo_q``) as the ocean freezing/melting potential. The sea
ice model then recomputes the frazil mass from the energy, using its own
assumptions, and returns any mismatch to the ocean as a freshwater and salt
correction. Salt fluxes are assumed proportional to the ice mass through a
fixed reference salinity
(``config_frazil_sea_ice_reference_salinity`` in the ocean and
``seaiceReferenceSalinity`` in the sea ice). As a result, the ocean releases
any frazil salt in excess of the reference salinity back into its top layer.

This work adds an alternative, "omega-type" coupling. MPAS-Ocean is the single
authority for the mass, salt and energy of newly formed frazil ice:

* MPAS-Ocean accumulates frazil mass, salt mass and energy over the coupling
  interval and passes all three to the coupler as separate fields.
* MPAS-Seaice, when ``config_frazil_coupling_type = 'omega-fluxes'``, takes
  these fluxes directly. It uses them to set the volume, bulk salinity and
  enthalpy of new ice, rather than deriving new ice properties from the
  freezing potential and sea surface salinity.
* Because the sea ice uses the salinity it is given, the ocean no longer
  needs to expel "excess" frazil salt to match a fixed coupling salinity.
  The new ocean namelist option ``config_expel_excess_frazil_salt`` controls
  this behavior.
* MPAS-Seaice now allows ``config_salt_flux_coupling_type = 'prognostic'``,
  which turns on the existing Icepack option ``saltflux_option =
  'prognostic'``. Its main effect is on melt: the salt released to the ocean
  when ice melts is computed from the ice's prognostic salinity instead of a
  fixed reference salinity. The salt returned on melt then matches the salt
  that entered the ice as frazil.

The dual use of ``Fioo_q`` is removed. ``Fioo_q`` now carries only the melt
potential (negative), and a new field ``Fioo_frazilh`` carries the frazil
energy. New ice-to-ocean frazil fields are also added so that the coupler
budget diagnostics can track the frazil adjustment returned by the sea ice.

The capability is a stealth feature in MPAS-Ocean. It is off by default.
Success is defined as:

1. Default runs keep the original fixed-salinity frazil coupling. Answers
   differ from the previous code only at roundoff level, from the new frazil
   salt unit conversion.
2. Coupled runs with ``omega-fluxes`` run stably.
3. Frazil mass, salt and heat are conserved across the ocean, coupler and
   sea ice, as reported by the coupler budget tables and the MPAS-Ocean
   conservation check analysis member.


Requirements
------------

Requirement: Ocean is the single source of frazil mass, salt and energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ocean must provide the mass, salt mass and energy of frazil ice formed
during each coupling interval as independent quantities. The sea ice must be
able to form new ice that has exactly this mass, salt and energy, without
recomputing them under different assumptions.

Requirement: Variable frazil salinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The salinity of new frazil ice must be allowed to vary in space and time
according to the ocean's frazil formulation (for example, with fixed frazil reference salinity or the porosity-based
frazil salinity enabled by ``config_use_frazil_ice_porosity``).

Requirement: Prognostic salt flux from sea ice melt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When frazil enters the sea ice with a variable salinity, the salt flux
returned to the ocean as the ice melts must reflect the salt actually stored
in the ice, not a fixed reference salinity. Otherwise, salt is not conserved
between the ocean and the sea ice over a freeze-melt cycle.

Requirement: Distinct coupling fields for melt potential and frazil energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ocean freezing energy and the ocean melt potential must be passed in
separate coupler fields, so that each field has a single physical meaning
and sign convention.

Requirement: Conservation and diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Frazil mass, salt and heat must be conserved across components. The coupler
budget diagnostics (MCT and MOAB drivers) and the MPAS-Ocean conservation
check must report the frazil terms separately, for both the ocean-to-ice
and the ice-to-ocean exchanges.

Requirement: Backward compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default configuration must keep the original fixed-salinity frazil
coupling. Answers are not required to be bit-for-bit with the previous code:
roundoff-level differences from the frazil salt unit conversion are
acceptable, in addition to the baseline differences from the added coupler
fields.


Algorithm Design
----------------

Algorithm Design: Ocean is the single source of frazil mass, salt and energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In each ocean time step, frazil forms in layer :math:`k` when the layer
temperature falls below the freezing point. The thickness of new frazil ice is

.. math::

   \Delta h_k = \min\left( \frac{E_k}{L_f \rho_i},\;
                \alpha\, h_k \right),

where :math:`E_k` is the freezing energy, :math:`L_f` is
``config_frazil_heat_of_fusion``, :math:`\rho_i` is
``config_frazil_ice_density``, :math:`h_k` is the layer thickness and
:math:`\alpha` is ``config_frazil_fractional_thickness_limit``. The salt mass
per unit area carried by the new frazil is

.. math::

   \Delta s_k = \Delta h_k \, S_{f,k} \, \rho_i \, c_{ppt},

where :math:`S_{f,k}` is the frazil salinity (ppt) and :math:`c_{ppt}`
(``ppt_to_salt``, 10\ :sup:`-3`) converts ppt to kg of salt per kg. This is
the salt mass in kg m\ :sup:`-2`. It is also used directly as the
thickness-weighted salinity tendency (m ppt), because MPAS-Ocean treats
1 kg m\ :sup:`-2` of salt as 1 m ppt (``sflux_factor = 1``). Summing over the
column, the ocean accumulates three quantities over the coupling interval:

.. math::

   M_f = \rho_i \sum_k \Delta h_k, \qquad
   S_f = \sum_k \Delta s_k, \qquad
   Q_f = \sum_k \Delta e_k,

where :math:`\Delta e_k` is the frazil energy (latent heat, plus sensible
heat when ``config_include_sensible_heat_frazil`` is true). At the coupling
time, these are divided by the coupling interval :math:`\Delta t_c` to form
fluxes (:math:`M_f/\Delta t_c`, :math:`S_f/\Delta t_c` and
:math:`Q_f/\Delta t_c`).

On the sea ice side, with ``omega-fluxes`` and mushy thermodynamics, Icepack
sets the volume, bulk salinity and enthalpy of new ice directly from these
fluxes:

.. math::

   V_{new} = \frac{F_m \, \Delta t}{\rho_i}, \qquad
   S_{new} = \max\left( 1000\,\frac{F_s}{F_m},\; S_{min} \right), \qquad
   q_{new} = -\frac{F_h}{F_m}\,\rho_i,

where :math:`F_m`, :math:`F_s` and :math:`F_h` are the frazil mass, salt and
heat fluxes received from the ocean. When :math:`S_{new}` falls below the
Icepack minimum salinity, Icepack issues a non-fatal warning and limits the
salinity to that minimum. Because the sea ice does not recompute the frazil
mass, the frazil mass adjustment returned to the ocean is zero in this mode.

Algorithm Design: Variable frazil salinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With frazil porosity enabled (``config_use_frazil_ice_porosity = .true.``),
the frazil salinity in each layer is

.. math::

   S_{f,k} = \phi \, S_k,

where :math:`\phi` (``config_frazil_ice_porosity``, between 0 and 1) is the
frazil porosity and :math:`S_k` is the ocean salinity in layer :math:`k`.
With porosity disabled, the frazil salinity is the fixed reference salinity
:math:`S_{ref}` (``config_frazil_sea_ice_reference_salinity``), limited so
that frazil does not take more salt than the layer contains:

.. math::

   S_{f,k} = \min\left( S_{ref},\; S_k \right).

In the original fixed-salinity coupling, the sea ice assumes that new ice
has salinity :math:`S_{ref}`. With porosity enabled, the ocean therefore
returns salt in excess of :math:`S_{ref}` to the top ocean layer, to stay
consistent:

.. math::

   \Delta S_{expel} = \max\left(0,\; S_f - M_f \, S_{ref} \, c_{ppt}\right).

When :math:`S_f < M_f S_{ref} c_{ppt}`, the deficit is recorded in
``lowSalinityFrazilIce`` as a possible conservation error.

In the omega-type coupling, the sea ice uses :math:`S_f` directly, so this
step must be turned off. Setting ``config_expel_excess_frazil_salt = .false.``
skips the expulsion. All of the frazil salt is then passed to the sea ice
through ``Fioo_frazils``.

The combination of porosity on (``config_use_frazil_ice_porosity = .true.``)
and expulsion off (``config_expel_excess_frazil_salt = .false.``) conserves
salt in the coupled system only when it is used with
``config_frazil_coupling_type = 'omega-fluxes'`` and
``config_salt_flux_coupling_type = 'prognostic'`` in MPAS-Seaice. With any
other sea ice settings, the sea ice assumes a fixed frazil salinity (or
recomputes the frazil salinity itself). The salt the ocean passes in
``Fioo_frazils`` is then not the salt the ice stores, nor the salt it returns
on melt, so salt is not conserved.

Land-ice exception: frazil that forms under land ice (ice-shelf cavities) is
not passed to the sea ice, and its salinity is assumed to be 0. With porosity
on, however, the layer loop still removes salt :math:`\phi S_k` from the
ocean as this frazil forms. To keep salt in the ocean, the expulsion step is
always applied under land ice, regardless of
``config_expel_excess_frazil_salt``. Because the target salinity is 0 under
land ice, all of the extracted salt is returned to the top ocean layer of the
cavity column, and none is added to ``accumulatedFrazilIceSalinity``.

As part of this change, the salt content in the frazil routine was made
dimensionally consistent: salt mass is now computed as thickness times
salinity times ice density times ``ppt_to_salt`` (kg m\ :sup:`-2`), in both
the formation step and the expulsion step. The previous code used thickness
times salinity. With the default ice density of 1000 kg m\ :sup:`-3` the two
are mathematically equal, but multiplying by 1000 and then by 0.001 rounds
twice in floating point. Results are therefore not bit-for-bit with the
previous code; the differences are at roundoff level.

For the non-omega (default) sea ice path, MPAS-Seaice adds a spatially
varying forcing field ``oceanFrazilSalinityReduction`` (default 3 ppt). It is
passed to Icepack as ``dSin0_frazil_cell``, the bulk salinity reduction of
new frazil relative to the sea surface salinity.

Algorithm Design: Prognostic salt flux from sea ice melt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This capability uses an option that already exists in Icepack,
``saltflux_option``. MPAS-Seaice previously allowed only ``'constant'``.
This PR lets MPAS-Seaice select ``'prognostic'``.

With ``'constant'``, the salt flux to the ocean when ice melts is the melted
ice mass times a fixed reference salinity :math:`S_{ref}`
(``ice_ref_salinity``):

.. math::

   F_{salt} = \frac{\rho_i \, \Delta V_{melt} \, S_{ref} \times 10^{-3}}{\Delta t}.

With ``'prognostic'``, Icepack uses the prognostic salinity of the melting
ice, :math:`\bar{S}`, in place of :math:`S_{ref}`:

.. math::

   F_{salt} = \frac{\rho_i \, \Delta V_{melt} \, \bar{S} \times 10^{-3}}{\Delta t}.

The salt released as ice melts therefore tracks the salt that entered the
ice. With omega-type coupling, that salt is the variable frazil salt supplied
by the ocean.

The option requires mushy thermodynamics, the only Icepack scheme that
carries a prognostic ice salinity. In MPAS-Seaice, it is also required to be
used with ``omega-fluxes``.

Algorithm Design: Distinct coupling fields for melt potential and frazil energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Previously, ``Fioo_q`` was positive when frazil formed (frazil energy) and
negative when the ocean was above freezing (melt potential). The new scheme
splits these roles:

* ``Fioo_q``: melt potential only (:math:`\le 0`), set when no frazil has
  accumulated.
* ``Fioo_frazilh``: frazil energy flux (:math:`\ge 0`), set when frazil has
  accumulated.

The sea ice reconstructs the old combined freezing/melting potential from
these two fields: it uses ``Fioo_q`` when that field is negative and
``Fioo_frazilh`` otherwise. Existing sea ice code paths therefore see the
same input as before.

The frazil mass field is renamed from ``Fioo_frazil`` to ``Fioo_frazilm``,
and a frazil salt flux ``Fioo_frazils`` is added.

Algorithm Design: Conservation and diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three ice-to-ocean fields are added: ``Fioi_frazilm``, ``Fioi_frazils`` and
``Fioi_frazilh``. They carry the frazil mass adjustment the sea ice returns
to the ocean, and its associated salt and heat. These fields are already
included in ``Fioi_meltw`` and ``Fioi_salt``, so they are diagnostic: they let
the coupler budget separate the frazil terms from melting. The coupler budget
tables gain separate entries for frazil heat, water and salt, for both the
ocean side (``hfrz``, ``wfrz``, ``sfrazo``) and the ice side (``hfrazi``,
``wfrzi``, ``sfrazi``).

In the MPAS-Ocean conservation check, the frazil heat flux is now labeled
with ``o2x_Fioo_frazilh`` and the frazil mass flux with ``o2x_Fioo_frazilm``.
The frazil mass term is computed from the change in accumulated frazil mass,
instead of from the layer thickness tendency.


Implementation
--------------

Implementation: Ocean is the single source of frazil mass, salt and energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MPAS-Ocean:

* ``src/shared/mpas_ocn_frazil_forcing.F``: accumulates
  ``accumulatedFrazilIceMass``, ``accumulatedFrazilIceSalinity`` (salt mass,
  kg m\ :sup:`-2`) and ``accumulatedFrazilIceEnergy`` (J m\ :sup:`-2`) in
  every column, regardless of the porosity setting.
* ``driver/ocn_comp_mct.F``: at the coupling time, sets

  .. code-block:: fortran

     if (keepFrazil) then
        if (accumulatedFrazilIceMass(i) > 0) then
           Fioo_frazilh = accumulatedFrazilIceEnergy / ocn_cpl_dt
        else
           Fioo_q = min(melt potential of top layer, 0) / ocn_cpl_dt
        endif
        Fioo_frazilm = accumulatedFrazilIceMass     / ocn_cpl_dt
        Fioo_frazils = accumulatedFrazilIceSalinity / ocn_cpl_dt
     endif

  The accumulated values are saved before being reset, for use in the
  restart and conservation logic.
* ``driver/mpaso_cpl_indices.F``: adds the ``Fioo_frazilm``,
  ``Fioo_frazils`` and ``Fioo_frazilh`` indices.

MPAS-Seaice:

* ``driver/ice_comp_mct.F``: reads the new fields into the
  ``frazilMassFlux``, ``frazilSaltFlux`` and ``frazilHeatFlux`` registry
  arrays. With ``omega-fluxes``, it sets ``frazilMassAdjust = 0``. Otherwise,
  it calls ``frazil_mass`` as before, and the salt and heat fluxes are zero.
* ``src/shared/mpas_seaice_icepack.F``: passes ``frazilm``, ``frazils`` and
  ``frazilh`` (and ``dSin0_frazil_cell``) to Icepack's
  ``icepack_step_therm2``.
* Icepack (submodule updated to ``7eb2ef0``): in ``add_new_ice``, uses the
  frazil fluxes to set the volume, salinity and enthalpy of new ice when
  ``cpl_frazil == 'omega-fluxes'``.

Implementation: Variable frazil salinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* New MPAS-Ocean namelist option in the ``frazil_ice`` group (added to
  ``Registry.xml``, ``namelist_definition_mpaso.xml``,
  ``namelist_defaults_mpaso.xml``, ``build-namelist`` and
  ``build-namelist-section``):

  .. code-block:: fortran

     config_expel_excess_frazil_salt = .true.   ! default

  If true, the ocean releases frazil salt in the upper ocean layer in excess
  of the fixed salinity coupling assumption. Set it to false only together
  with ``config_frazil_coupling_type = 'omega-fluxes'`` and
  ``config_salt_flux_coupling_type = 'prognostic'`` in MPAS-Seaice;
  otherwise the coupled system does not conserve salt.
* In ``src/shared/mpas_ocn_frazil_forcing.F``, the expulsion step runs when

  .. code-block:: fortran

     if (sumNewFrazilIceThickness > 0.0_RKIND .and. usePorosity .ge. 1 .and. &
         (config_expel_excess_frazil_salt .or. underLandIce)) then

  so that frazil salt under land ice is always returned to the ocean.
* The MPAS-Seaice forcing field ``oceanFrazilSalinityReduction`` is added to
  the ``ocean_coupling`` pool and to the monthly and initial forcing groups.

Implementation: Prognostic salt flux from sea ice melt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* MPAS-Seaice namelist option (``Registry.xml``):

  .. code-block:: fortran

     config_salt_flux_coupling_type = 'constant'   ! default; or 'prognostic'

  Previously, any value other than ``'constant'`` caused a configuration
  error. ``'prognostic'`` is now accepted.
* ``src/shared/mpas_seaice_icepack.F`` (``check_column_package_configs``):
  stops the run with a critical error if ``'prognostic'`` is used without
  ``config_frazil_coupling_type = 'omega-fluxes'`` or without
  ``config_thermodynamics_type = 'mushy'``.
* The value is passed to Icepack as ``saltflux_option`` through
  ``icepack_init_parameters``. The prognostic salt flux calculation itself
  already existed in Icepack and was not changed by this PR.

Implementation: Distinct coupling fields for melt potential and frazil energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``driver-mct/shr/seq_flds_mod.F90`` and ``driver-moab/shr/seq_flds_mod.F90``:
  rename ``Fioo_frazil`` to ``Fioo_frazilm`` and add ``Fioo_frazils``
  (kg m\ :sup:`-2` s\ :sup:`-1`) and ``Fioo_frazilh`` (W m\ :sup:`-2`) to the
  ocean-to-ice field lists.
* ``components/data_comps/docn/src/docn_comp_mod.F90``: updated for the
  renamed field.

Implementation: Conservation and diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``seq_flds_mod.F90`` (MCT and MOAB): add ``Fioi_frazilm``,
  ``Fioi_frazils`` and ``Fioi_frazilh`` to the ice-to-ocean field lists.
* ``driver-mct/main/seq_diag_mct.F90`` and ``driver-moab/main/seq_diag_moab.F90``:
  add budget indices ``f_hfrazi``, ``f_wfrzi``, ``f_sfrazo`` and ``f_sfrazi``,
  and accumulate the frazil terms from the ocean and ice fields.
* ``src/analysis_members/mpas_ocn_conservation_check.F``: updates the labels
  for the frazil heat and mass fluxes, and computes the frazil mass term from
  the accumulated frazil mass.
* MPAS-Seaice ``cime_config/buildnml``: adds ``frazilMassFlux``,
  ``frazilSaltFlux`` and ``frazilHeatFlux`` to the sea ice output streams.

Implementation: Backward compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The defaults (``config_frazil_coupling_type = 'external'``,
``config_expel_excess_frazil_salt = .true.``,
``config_salt_flux_coupling_type = 'constant'``) keep the original behavior,
with roundoff-level differences from the frazil salt unit conversion.
To enable the omega-type coupling, set the following in ``user_nl_mpassi``:

.. code-block:: fortran

   config_frazil_coupling_type = 'omega-fluxes'
   config_ocean_surface_type = 'free'
   config_salt_flux_coupling_type = 'prognostic'

All three settings are required. ``omega-fluxes`` is checked against
``config_ocean_surface_type = 'free'`` at initialization in coupled runs.

To use variable (porosity-based) frazil salinity in the coupled system, also
set the following in ``user_nl_mpaso``:

.. code-block:: fortran

   config_use_frazil_ice_porosity = .true.
   config_expel_excess_frazil_salt = .false.

These ocean settings must only be used with the three MPAS-Seaice settings
above. The ocean does not check the sea ice configuration, so this is the
user's responsibility.


Testing
-------

Testing and Validation: Backward compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Existing E3SM tests that use frazil are not bit-for-bit with the previous
code. The frazil salt unit conversion introduces roundoff-level differences
in the ocean salinity, which grow over the run. Baseline comparisons also
show differences because new coupler fields were added and ``Fioo_q`` was
split into two fields. Tests are expected to run to completion, and the
default configuration should remain climatologically unchanged.

Testing and Validation: Omega-type coupling and conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two tests are in the ``e3sm_ocnice_stealth_features`` suite in
``cime_config/tests.py`` (also inherited by ``e3sm_ocnice_extra_coverage``):

.. code-block:: none

   SMS_PS.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.mpassi-frazil_omega_coupling
   SMS_PS.ne30pg2_r05_IcoswISC30E3r5.CRYO1850.mpaso-frazil_porosity_omega_coupling

The ``mpassi/frazil_omega_coupling`` testmod sets the three MPAS-Seaice
namelist options listed above, including
``config_salt_flux_coupling_type = 'prognostic'``, so the prognostic salt
melt flux is exercised together with the omega-type frazil coupling. It does
not change the MPAS-Ocean frazil options.

The ``mpaso/frazil_porosity_omega_coupling`` testmod tests the full
variable frazil salinity configuration. It sets the same three MPAS-Seaice
options, plus the following in ``user_nl_mpaso``:

.. code-block:: fortran

   config_use_frazil_ice_porosity = .true.
   config_frazil_ice_porosity = 0.85
   config_expel_excess_frazil_salt = .false.

A porosity below 1 is used so that the frazil salinity, :math:`0.85 \, S_k`,
is lower than the ocean salinity. This exercises a variable frazil salinity
that differs from both the ocean salinity and the fixed reference salinity.

The test uses the CRYO1850 compset, which has ice-shelf cavities with frazil
under land ice, so the land-ice exception is also exercised. The testmod
turns on
the conservation check analysis members in both MPAS-Ocean and MPAS-Seaice
(``config_AM_conservationCheck_enable = .true.``). Their energy, mass and
salt budgets are written to the component log files. Neither check stops the
run on a salt imbalance, so the logs must be inspected to confirm
conservation.

Conservation across components is also checked with the coupler budget
tables, which now list frazil heat, water and salt separately for the ocean
and ice.
