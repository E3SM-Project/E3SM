Omega-type frazil coupling with variable salinity between MPAS-Ocean and MPAS-Seaice
====================================================================================

date: 2026/11/04

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

1. Default runs reproduce the original answers apart from the additional
   coupler fields.
2. Coupled runs with ``omega-fluxes`` run stably.
3. Frazil mass, salt and heat are conserved across the ocean, coupler and
   sea ice, as reported by the coupler budget tables and the MPAS-Ocean
   conservation check analysis member.


Requirements
------------

Requirement: Ocean is the single source of frazil mass, salt and energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

The ocean must provide the mass, salt mass and energy of frazil ice formed
during each coupling interval as independent quantities. The sea ice must be
able to form new ice that has exactly this mass, salt and energy, without
recomputing them under different assumptions.

Requirement: Variable frazil salinity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

The salinity of new frazil ice must be allowed to vary in space and time
according to the ocean's frazil formulation (for example, with fixed frazil reference salinity or the porosity-based
frazil salinity enabled by ``config_use_frazil_ice_porosity``). It must not be
constrained to a fixed reference salinity at the coupling interface.

Requirement: Prognostic salt flux from sea ice melt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

When frazil enters the sea ice with a variable salinity, the salt flux
returned to the ocean as the ice melts must reflect the salt actually stored
in the ice, not a fixed reference salinity. Otherwise, salt is not conserved
between the ocean and the sea ice over a freeze-melt cycle.

Requirement: Distinct coupling fields for melt potential and frazil energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

The ocean freezing energy and the ocean melt potential must be passed in
separate coupler fields, so that each field has a single physical meaning
and sign convention.

Requirement: Conservation and diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

Frazil mass, salt and heat must be conserved across components. The coupler
budget diagnostics (MCT and MOAB drivers) and the MPAS-Ocean conservation
check must report the frazil terms separately, for both the ocean-to-ice
and the ice-to-ocean exchanges.

Requirement: Backward compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

The default configuration must keep the original fixed-salinity frazil
coupling and give answers that are bit-for-bit with the previous code. The
only expected baseline differences come from the added coupler fields.


Algorithm Design
----------------

Algorithm Design: Ocean is the single source of frazil mass, salt and energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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
(``ppt_to_salt``) converts ppt to kg of salt per kg. Summing over the
column, the ocean accumulates three quantities over the coupling interval:

.. math::

   M_f = \rho_i \sum_k \Delta h_k, \qquad
   S_f = \sum_k \Delta s_k, \qquad
   Q_f = \sum_k \Delta e_k,

where :math:`\Delta e_k` is the frazil energy (latent heat, plus sensible
heat when ``config_include_sensible_heat_frazil`` is true). At the coupling
time these are divided by the coupling interval :math:`\Delta t_c` to form
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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

With frazil porosity enabled (``config_use_frazil_ice_porosity = .true.``),
the frazil salinity in each layer is

.. math::

   S_{f,k} = \min\left( \phi \, S_k,\; S_k \right),

where :math:`\phi` is ``config_frazil_ice_porosity`` and :math:`S_k` is the
ocean salinity in layer :math:`k`. In the original fixed-salinity coupling,
the sea ice assumes that new ice has salinity :math:`S_{ref}`
(``config_frazil_sea_ice_reference_salinity``). To stay consistent, the
ocean returns salt in excess of :math:`S_{ref}` to the top ocean layer:

.. math::

   \Delta S_{expel} = \max\left(0,\; S_f - M_f \, S_{ref} \, c_{ppt}\right).

When :math:`S_f < M_f S_{ref} c_{ppt}`, the deficit is recorded in
``lowSalinityFrazilIce`` as a possible conservation error.

In the omega-type coupling, the sea ice uses :math:`S_f` directly, so this
step must be turned off. Setting ``config_expel_excess_frazil_salt = .false.``
skips the expulsion. All of the frazil salt is then passed to the sea ice
through ``Fioo_frazils``.

As part of this change, the salt content in the frazil routine was made
dimensionally consistent: salt mass is now computed as thickness times
salinity times ice density times ``ppt_to_salt`` (kg m\ :sup:`-2`), in both
the formation step and the expulsion step.

For the non-omega (default) sea ice path, MPAS-Seaice adds a spatially
varying forcing field ``oceanFrazilSalinityReduction`` (default 3 ppt). It is
passed to Icepack as ``dSin0_frazil_cell``, the bulk salinity reduction of
new frazil relative to the sea surface salinity.

Algorithm Design: Prognostic salt flux from sea ice melt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

* New MPAS-Ocean namelist option in the ``frazil_ice`` group (added to
  ``Registry.xml``, ``namelist_definition_mpaso.xml``,
  ``namelist_defaults_mpaso.xml``, ``build-namelist`` and
  ``build-namelist-section``):

  .. code-block:: fortran

     config_expel_excess_frazil_salt = .true.   ! default

  If true, the ocean releases frazil salt in the upper ocean layer in excess
  of the fixed salinity coupling assumption. Set it to false when using
  variable sea ice salinity coupling.
* The MPAS-Seaice forcing field ``oceanFrazilSalinityReduction`` is added to
  the ``ocean_coupling`` pool and to the monthly and initial forcing groups.

Implementation: Prognostic salt flux from sea ice melt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

* ``driver-mct/shr/seq_flds_mod.F90`` and ``driver-moab/shr/seq_flds_mod.F90``:
  rename ``Fioo_frazil`` to ``Fioo_frazilm`` and add ``Fioo_frazils``
  (kg m\ :sup:`-2` s\ :sup:`-1`) and ``Fioo_frazilh`` (W m\ :sup:`-2`) to the
  ocean-to-ice field lists.
* ``components/data_comps/docn/src/docn_comp_mod.F90``: updated for the
  renamed field.

Implementation: Conservation and diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

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

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

The defaults (``config_frazil_coupling_type = 'external'``,
``config_expel_excess_frazil_salt = .true.``,
``config_salt_flux_coupling_type = 'constant'``) keep the original behavior.
To enable the omega-type coupling, set the following in ``user_nl_mpassi``:

.. code-block:: fortran

   config_frazil_coupling_type = 'omega-fluxes'
   config_ocean_surface_type = 'free'
   config_salt_flux_coupling_type = 'prognostic'

All three settings are required. ``omega-fluxes`` is checked against
``config_ocean_surface_type = 'free'`` at initialization in coupled runs.


Testing
-------

Testing and Validation: Backward compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

Existing E3SM developer tests are expected to be bit-for-bit in model state.
They show baseline differences only because new coupler fields were added
and ``Fioo_q`` was split into two fields.

Testing and Validation: Omega-type coupling and conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2026/11/04

Contributors: Nicole Jeffery

A new test is added to the ``e3sm_integration`` suite in
``cime_config/tests.py``:

.. code-block:: none

   SMS_PS.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.mpassi-frazil_omega_coupling

The ``mpassi/frazil_omega_coupling`` testmod sets the three MPAS-Seaice
namelist options listed above, including
``config_salt_flux_coupling_type = 'prognostic'``, so the prognostic salt
melt flux is exercised together with the omega-type frazil coupling. Conservation is checked with the coupler
budget tables, which now list frazil heat, water and salt separately for the
ocean and ice, and with the MPAS-Ocean conservation check analysis member
(``mpaso-conservation_check`` testmod).
