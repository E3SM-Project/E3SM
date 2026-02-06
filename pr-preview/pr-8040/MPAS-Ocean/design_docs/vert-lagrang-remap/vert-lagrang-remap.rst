
Vertical Lagrangian-remap method
================================

date: 2020/12/21
Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis



Summary
-------

We propose to implement the vertical Lagrangian-remap method to support the
choice of hybrid vertical coordinates, which may be required for evolving ice-
shelf cavity geometries. It can be seen as an extension of the Arbitrary
Lagrangian-Eulerian method which is already implemented in MPAS-Ocean, a key
difference being that the solution of diasurface transport is coordinate-free.
This has the advantage of removing the vertical CFL constraint, which is of
particular concern for an evolving ice shelf geometry in terrain-following
coordinate schemes. This development may also reduce numerical diffusion in the
vertical, but we do not require this result for this implementation to be
successful.

Generalization of the vertical coordinate is not in the scope of this development.
Additionally, the infrastructure for specifying a target grid that differs from
the existing z-star and z-tilde options is reserved for future development.
The minimum criteria for success are not degrading the accuracy of the solution 
or computational performance.


Requirements
------------

Requirement: Ability to perform vertical Lagrangian-remapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

The core of the proposed development is the ability to perform vertical
Lagrangian-remapping (VLR). When VLR is performed, the momentum, volume and
tracer conservation equations are first solved in a Lagrangian sense; that is,
layer thicknesses evolve such that there is no vertical transport across layer
interfaces. A target grid is then specified. For this stage of development,
we only require that the algorithm be able to use the existing grid
specifications. However, the implementation should be general enough that
a different grid could be used for remapping. This regridding step is followed
by a conservative remapping of velocity and scalars to the target grid.

Requirement: Grid is sufficiently conditioned for accuracy and stability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

A severely deformed grid can lead to issues with the accuracy and stability of
the numerical methods used for the dynamics and thermodynamics. If layers are
allowed to become very thin they can also degrade the accuracy and stability of
the solution. Thus, a requirement of the VLR implementation is that there are
safeguards to prevent layers from becoming too thin or deformed.

Requirement: Ability to specify remapping method and its options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

There are several choices to make for the remapping operation to balance
accuracy, stability and computational cost. To the extent that it is feasible,
the user should have the ability to control the order of accuracy of the
remapping method, and other remapping options (e.g., via namelist options).

Requirement: Support for existing time stepping schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Vertical Lagrangian remapping should be supported (eventually) in both the RK4
and the split-explicit time stepping schemes currently implemented in MPAS-
Ocean. The scheme will also have a straightforward interface so it can be easily
added to new time stepping schemes.

Requirement: not climate changing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

For the same target grid, time-stepping scheme and problem, the solution for
layer thicknesses and key prognostic variables using VLR is within a non-
climate-changing threshold value of the solution when VLR is not performed.
Tests should span from idealized tests of simplified dynamics to full climate
simulations. Acceptable thresholds specified for each test could be similar to
those from changing time stepping schemes.

Requirement: tracer and volume conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2020/12/04

Contributors: Darren Engwirda, Xylar Asay-Davis, Carolyn Begeman, Mark Petersen

The column-wise tracer content should be conserved over the remapping step to a
specified threshold, and using a specified method to integrate the tracer. The
column-integrated volume must remain unchanged over the remapping step, as
expected.

Requirement: performance
^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2020/12/04

Contributors: Darren Engwirda, Xylar Asay-Davis, Carolyn Begeman, Mark Petersen

Performance is not degraded below an acceptable level. Compute time for
stand-alone ocean, without i/o, on full primitive equations, is expected to
remain very similar - within 5% of compute time when using vertical
advection method but likely less because remapping is a column-wise, in-cache
operation.

Requirement: no interference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Remapping does not introduce errors in the computation of prognostic and
diagnostic variables, such as through the use of an Eulerian vertical velocity
that is not consistent with the Lagrangian solution.

Requirement: modularity
^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2020/12/21

Contributors: Darren Engwirda, Xylar Asay-Davis, Carolyn Begeman

To the degree possible, the code for determining the target grid and performing 
remapping should be kept in its own Fortran module(s) for better readability.
These modules should be called by both timestepping routines to maximize code
reuse.


Algorithm Design (optional)
---------------------------

Algorithm Design: Ability to perform vertical Lagrangian-remapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2020/12/15

Contributors: Darren Engwirda, Xylar Asay-Davis, Carolyn Begeman

The conservation of momentum, volume, and tracer equations in MPAS-Ocean (
`Ringler et al. 2013 <https://www.sciencedirect.com/science/article/abs/pii/S1463500313000760>`_; 
`Petersen et al. 2014 <https://www.sciencedirect.com/science/article/abs/pii/S1463500314001796>`_)
are:

.. math::

   \frac{\partial u_k}{\partial t} + q_k h_k u_k^{normal} + \overline{w^t \delta z^t u} = -\frac{1}{\rho_0} \nabla p_k - \frac{\rho_k g}{\rho_0} \nabla z_k - \nabla K_k + [D_h^u]_k + [D_{\nu}^u]_k + F_k^u

   \frac{\partial h_k}{\partial t} + \nabla \cdot (h_k \mathbf{u}_k) + w_k^t - w_{k+1}^t = 0

   \frac{\partial (h_k \phi_k)}{\partial t} + \nabla \cdot (h_k \mathbf{u}_k \phi_k) + \overline{\phi}_k^t w_k^t - \overline{\phi}_{k+1}^t w_{k+1}^t = [D_h^{\phi}]_k + [D_v^{\phi}]_k + F_k^{\phi}

For the Lagrangian step, the vertical velocity through the top of the cell,
:math:`w_k^t`, is set to zero in all of the above equations. Thus, these
equations simplify to:

.. math::

   \frac{\partial u_k}{\partial t} + q_k h_k u_k^{\perp} = -\frac{1}{\rho_0} \nabla p_k - \frac{\rho_k g}{\rho_0} \nabla z_k - \nabla K_k + [D_h^u]_k + [D_v^u]_k + F_k^u

   \frac{\partial h_k}{\partial t} + \nabla \cdot (h_k \mathbf{u}_k) = 0

   \frac{\partial (h_k \phi_k)}{\partial t} + \nabla \cdot (h_k \mathbf{u}_k \phi_k) = [D_h^{\phi}]_k + [D_v^{\phi}]_k + F_k^{\phi}

The time-stepping algorithm (RK4 or split-explicit) advances the prognostic
variables and layer thickness from their values at time n
:math:`u_k^{n},\phi_k^{n},h_k^{n}`, to their values after the Lagrangian step,
designated by the superscript *lg*, :math:`u_k^{lg},h_k^{lg},\phi_k^{lg}`.

Note that the vertical mixing terms :math:`D_v^h, D_v^{\phi}`
are retained here. We opt to compute these terms prior to remapping as this
allow for future development in which the dynamics are subcycled relative to
the thermodynamics and remapping is scheduled on the thermodynamic timestep.
This computation of vertical mixing terms prior to remapping is similar to
both MOM6 and HYCOM. We anticipate that there could be a trade-off between (a)
loss of accuracy of vertical mixing terms when their computation precedes
remapping due to grid deformation and (b) loss of accuracy when their
computation follows remapping due to remapping errors in vertical gradients of
prognostic variables. We do not intend to test this at this time.

The target grid needs to be determined after the solution for prognostic
variables so that the vertical Lagrangian-remapping method is general enough to
be used with coordinate systems that depend on the ocean state (this includes
the z-star coordinate system in which SSH perturbations are vertically
distributed between layers). We do not present an algorithmic design for
regridding to coordinate systems not already supported in MPAS-Ocean, as this
will be the subject of future development. For now, the target grid is based on a
constant set of z-star levels that are specified at initialization.

For the grid selection step, the target grid, :math:`h_k^{target}`, is
determined, conserving volume:

.. math::

   \sum_{k=1}^{kmax}h_k^{target} = \sum_{k=1}^{kmax}h_k^{lg}


For scalar remapping, layer thicknesses at the next timestep,
:math:`h_k^{n+1}` are set according to the target grid and scalars are remapped
to the target grid using the remapping operations represented by the function
:math:`G`:

.. math::

   h_k^{n+1} = h_k^{target}

   hEdge_k^{n+1} = 0.5 (h_{k,cell1}^{n+1} + h_{k,cell2}^{n+1})

   u_k^{n+1} = G(u_k^{lg},hEdge_k^{lg},hEdge_k^{n+1})

   \phi_k^{n+1} = G(\phi_k^{lg},\phi_k^{lg},h_k^{n+1})

For velocity remapping, we solve for layer thicknesses at edges after the
lagrangian step and the regridded thickness. In this development, we only plan
to support centered edge layer thicknesses consistent with the centered
advection scheme. There does not appear to be a precedent among ocean models
(HYCOM, MOM6) at this time for using upwinded layer thickness in the remapping
operation. We touch on a few of the implementation challenges with using
upwinded layer thicknesses for remapping in the corresponding implementation
section. If VLR is run with an upwinded thickness flux, the horizontal momentum
flux will not be conserved as :math:`hEdge^{n+1}` will be reassigned to the
upwinded layer thickness (errors will likely increase as horizontal gradients
in layer thickness increase). Otherwise, the remapping operation :math:`G`
conserves volume flux and scalar concentration.

.. math::

   hEdge_k^{lg} = 0.5 (h_{k,cell1}^{lg} + h_{k,cell2}^{lg})

   hEdge_k^{n+1} = 0.5 (h_{k,cell1}^{n+1} + h_{k,cell2}^{n+1})

   u_k^{n+1} = G(u_k^{lg},hEdge_k^{lg},hEdge_k^{n+1})


.. math::

   \sum_{k=1}^{kmax} u_k^{n+1} h_k^{n+1} = \sum_{k=1}^{kmax} u_k^{lg} h_k^{lg}

   \sum_{k=1}^{kmax} \phi_k^{n+1} h_k^{n+1} = \sum_{k=1}^{kmax} \phi_k^{lg} h_k^{lg}

The vertical velocity across layer interfaces may be computed anytime after
regridding. It can be computed as

.. math::

   w = - \nabla \cdot (h_k \mathbf{u}_k) - (h_k^{t+1} - h_k^t)/dt

or

.. math::

   w = (h_k^{t+1} - h_k^{lg})/dt

The choice between the two is discussed in the Implementation section.


Implementation
--------------

Implementation: Ability to perform vertical Lagrangian-remapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Namelist options:

- To turn VLR on/off:
  :code:`advection, config_vert_advection_method = 'remap' or 'flux-form'`


Lagrangian step:

The solution for prognostic variables in RK4 and split-explicit remains
largely the same. The main difference is that the vertical velocity through
the top of layers is set to zero in the routine
:code:`ocn_vert_transport_velocity_top`. This is the same as what is done when
:code:`config_vert_coord_movement` is :code:`impermeable_interfaces`.


Grid selection step:

#. :math:`z_k^{target}`, the depth of the top of the layer, is determined based
   on an analytical expression for the grid.
   The simplest case is constant z-levels, :math:`z_k^{target} = z_k^{init}`.
   Since :math:`z_k^{target}` can be a function of the ocean state (e.g.,
   :math:`\rho` for isopycnal coordinates, grid selection doesn't occur until
   after the solution for prognostic variables.
#. Superimpose SSH perturbations according to one of the existing depth-
   dependent functions, :math:`z_k^{target} = z_k^{target} + c(z) \: \eta`. As
   in :code:`ocn_ALE_thickness`, layer thicknesses are adjusted from the
   seafloor upwards. Currently :code:`layerThicknessTarget` is computed in a
   new routine :code:`ocn_vert_regrid` following other instances in which a z-
   star grid is recomputed.
#. Apply conditioning steps outlined in the following section. (Not
   implemented until coordinates other than z-star are implemented.)
#. In preparation for remapping, compute :math:`z_k^{target}` from
   :code:`layerThicknessTarget`.

All of the grid selection steps will be performed from a separate module called
:code:`ocn_vert_regrid`. This topic is further addressed in section
Implementation: modularity.


Remapping step:

There is a new remapping routine :code:`ocn_remap_vert_state` which updates
all state variables with depth coordinates at :code:`timeLevel=2`. On input,
:code:`layerThickness(tlev=2)` reflects the Lagrangian layer thickness
determined by :code:`ocn_tend_thick`. On output, it is equal to
:code:`layerThicknessNew` as determined by the regridding routine
:code:`ocn_vert_regrid`.

Members of :code:`statePool` that will be remapped:

- All members of :code:`tracerPool` unless :code:`activeTracersOnly`, in which
  case only the :code:`activeTracers` will be remapped
- :code:`normalVelocity`
- :code:`highFreqThickness`
- :code:`lowFreqDivergence`
- :code:`normalBarotropicVelocity`, only for split-explicit time-stepping

In preparation for remapping, we compute the scratch variables
:code:`layerThickEdgeNew` as the mean of :code:`layerThicknessNew` at
neighboring cells. We do these locally rather than through a call to 
:code:`ocn_diagnostic_solve_layerThicknessEdge`. For the Lagrangian
:code:`layerThickEdgeMean` we use the existing solution from
:code:`ocn_diagnostic_solve_layerThicknessEdge`. This approach requires that
:code:`config_thickness_flux_type` is :code:`'centered'`. At initialization,
we throw an error but do not terminate the run if
:code:`config_thickness_flux_type` is not :code:`'centered'` and VLR is active.


A note about difficulties of implementing upwinded :code:`layerThicknessEdge`
for remapping:

Currently, the PPR library assumes that the total height is the same before and
after remapping (i.e., :code:`sum_k(layerThicknessOld)` equals
:code:`sum_k(layerThicknessTarget)`. Over the course of remapping, the upwind
cell could change for one or more layers and thus the total column height could
change. PPR would have to be carefully adapted to deal with this condition in
order to preserve the total volume flux as well as the vertical distribution of
momentum during remapping.

An alternative to modifying the remapping library is to use centered edge layer
thicknesses for remapping and correct :code:`uNormal` prior to remapping such
that :code:`uNormalCorr(k) * layerThicknessEdgeCntr(k) = uNormal(k) * layerThicknessEdgeUpwind(k)`.
When edge layer thicknesses are upwinded based on the remapped :code:`uNormal`,
:code:`uNormal` must be corrected again to preserve layerwise fluxes. There
will still be some error in the vertical distribution of volume flux with this
approach. Given the complexity of either of these implementation options, we
leave this issue for future development.


After determining the layer thicknesses to remap to, the vertical layer
interface locations are determined, :code:`heightCellNow`,
:code:`heightCellNew`, :code:`heightEdgeNow`, :code:`heightEdgeNew`.
These variables and the Lagrangian state variable are the inputs to the PPR
routine :code:`rmap1d`, and the output is the remapped state variable.


Some implementation considerations for PPR:
 
- Error-checking in PPR: make consistent with MPAS errors, consider additional
  error checks. (Not implemented)

After remapping, :code:`ocn_diagnostic_solve` is called. This is needed to
compute the density and pressure fields based on the remapped ocean state and
the diagnostic field :code:`vertVelocityTop` which is the vertical velocity
through the top of the layer. This is only used as a diagnostic variable for
computing the MOC streamfunction and tracer budgets. None of the mixing
parameterizations require a vertical velocity (Eulerian or diasurface velocity).

Note: if `vertVelocityTop` is computed between regridding and remapping then it
can be computed as

.. code::

   vertVelocityTop(k) = vertVelocityTop(k+1) - div_hu(k) - 
                        (layerThickness(k,tlev=2) - layerThickness(k,tlev=1))/dt

In our implementation, `vertVelocityTop` is computed after remapping and
:code:`div_hu` is no longer appropriate as it has been remapped. In this case,
the Lagrangian layer thickness is stored in the state variable
:code:`layerThicknessLag` and then the vertical velocity through the top of the
layer is computed as:

.. code::

   layerThicknessLag = layerThickness(tlev=2)

   layerThickness(tlev=2) = layerThicknessTarget

   vertVelocityTop = (layerThickness(tlev=2) - layerThicknessLag)/dt

The function that performs this computation is
:code:`ocn_diagnostic_solve_vertVel_remap`.

If :code:`normalGMBolusVelocity` is computed based on the remapped ocean state
then the computation of :code:`vertTransportVelocityTop` and
:code:`vertGMBolusVelocitytop` is unchanged as these fields represent Eulerian
velocities.

However, the current implementation will not compute
:code:`normalGMBolusVelocity` based on the remapped ocean state before
:code:`ocn_diagnostic_solve` is called. Thus, :code:`vertGMBolusVelocityTop` and
:code:`vertTransportVelocityTop` will be inaccurate. This will only pose an
issue when the number of vertical levels changes during regridding, a
capability which isn't included in this development scope.

Tracer tendencies that are computed as diagnostics will also be inaccurate
after regridding, as they will not be remapped. Remapping these variables does
not make physical sense without also computing vertical tracer fluxes, which
would be overly burdensome. Analysis members that currently use these diagnostic
variables are :code:`mpas_ocn_layer_volume_weighted_averages` and
:code:`mpas_ocn_mixed_layer_heat_budget`.

Implementation: Grid is sufficiently conditioned for accuracy and stability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Since the grid is remapped back to z-star in the current implementation, we
assume that the grid is sufficiently regular.

For future coordinate choices, consdider applying the following after
determining the target grid:

#. Assign :math:`h_k^{t+1}` to :math:`h_k^{lg}` if :math:`h_k^{t+1} - h_k^{lg}`
   is less than a minimum thickness alteration. This is motivated by accuracy
   considerations, as each remapping may introduce errors. (Not implemented)
#. Apply minimum layer thickness criterion. (Not implemented)
#. Smooth layers in space and time. (Not implemented)

Namelist options:

- Minimum layer thickness
- Minimum thickness change for remapping to occur


Implementation: Ability to specify remapping method and its options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Namelist options:

- Order of remapping: 
  :code:`advection, config_vert_remap_order = 1, 2, 3, or 5`
- Interval number of timesteps for remapping: 
  :code:`advection, config_vert_remap_interval = integer >= 0`
- Slope limiter: 
  :code:`advection, config_remap_limiter = 'none', 'monotonic', or 'weno'`

Implementation: Support for existing time stepping schemes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Currently, only the split-explicit time integration scheme is supported. We
explored the RK4 implementation but it was not found to be stable.

Implementation: performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

OMP calls were implemented for all cell loops in the remapping routine. Our
main strategy for improving performance is remapping an interval of time steps
rather than every time step. For EC60to30 global configurations, an interval of
3-4 timesteps was found to not significantly degrade the solution.

Options for improving performance (not implemented):

- Splitting the scalar and momentum timesteps
- Only remapping when the change in thickness exceeds given threshold
- Optimizing/parallelizing PPR

Implementation: no interference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2020/12/15

Contributors: Darren Engwirda, Xylar Asay-Davis, Carolyn Begeman

Ensure that no calculations are made outside of remapping itself using
prognostic (or diagnostic) variables both before and after remapping to avoid
introducing additional errors.

Look for places in the code where prognostic variables are used at the previous
timestep.

Implementation: modularity
^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Remapping operations (PPR) are performed in a separate routine,
:code:`ocn_vert_remap_state`.

Target grid levels are determined in a separate routine,
:code:`ocn_vert_regrid`.


Testing
-------

Testing and Validation: Ability to perform vertical Lagrangian-remapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Ability to handle strong vertical velocities: 

- Baroclinic channel test case

Evaluating spurious mixing due to remapping: Compare with and without VLR

- Internal wave test case
- Baroclinic channel test case


Testing and Validation: not climate changing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

Global ocean case :code:`GMPAS-JRA1p4.TL319_EC30to60E2r2` 50 year run comparison
between VLR branch with VLR on and master without remapping:
`Link to analysis <https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/ac.cbegeman/scratch/anvil/mpas_analysis_output/20211020.GMPAS-JRA1p4.TL319_EC30to60E2r2.vlr.anvil.control/yrs46-50/>`_

Testing and Validation: performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2022/05/21

Contributors: Carolyn Begeman, Darren Engwirda, Xylar Asay-Davis

For the same global ocean case above we show that remapping increases run time
by 20% for remapping every time step and that it scales almost perfectly with as
the remapping interval increases, i.e., 5% increase for remapping every 4 time
steps. We also show that 90% of the remapping time is spent in the PPR routine
:code:`rmap1d`, indicating that this would be the target for performance
improvements.
