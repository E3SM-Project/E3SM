
Prognostic non-cohesive sediment transport model
================================================

date: 2021/01/14

Contributors: Zhendong Cao



Summary
-------

The current MPAS-Ocean includes a sediment transport analysis member. But it lacks
the capabilities to simulate the morphological change and the potential interactions
between sediment transport and other system components (e.g., flow field, vegetation dynamics, etc.).
Therefore, we propose a prognostic sediment transport model to enable these missing capabilities.

The proposed model development will focus on the non-cohesive sediment transport.
The model requirements, algorithms and implementation mainly follow the sediment transport
development work in `Warner et al., 2005 <https://doi.org/10.1016/j.cageo.2008.02.012>`_ and
the `COAWST code (USGS version) <https://github.com/jcwarner-usgs/COAWST>`_.
The model components include bedload sediment transport and suspended sediment transport;
These are the two main transport forms considered in this development.

The primary implementation challenges lie in the suspended sediment transport component.
It requires solving the three-dimensional scalar equation for suspended sediment with
its own source/sink terms, which need to resolve the deposition and resuspension at each
time step.

The successful implementation will enable MPAS-Ocean to simulate the sediment transport (especially
the computation of suspended sediment concentration; SSC) with given initial and boundary conditions,
as well as the potential interactions between sediment transport, hydrodynamics and morphology.


Requirements
------------

Requirement: Bedload sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/14

Contributors: Zhendong Cao, Xylar Asay-Davis

The prognostic sediment transport capability in MPAS-Ocean must include bedload sediment transport.
Bedload sediment transport means the sediment particles have successive contacts with the bed during transport.
Bedload mass is exchanged horizontally at the bed layer.

Requirement: Suspended sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/17

Contributors: Zhendong Cao, Xylar Asay-Davis

The prognostic sediment transport capability in MPAS-Ocean must include suspended sediment transport.
If water flows significantly faster than the threshold of the sediment motion, sediment will be entrained off the bed and into suspension.
To remain in suspension, the upward velocity of turbulence eddies must be greater than the sediment settling velocity.
Suspended mass is exchanged vertically between the water column and the top bed layer.

Requirement: MPAS-Ocean unchanged if sediment transport is disabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/25

Contributors: Xylar Asay-Davis, Zhendong Cao

The model should run exactly as it currently does (results should be "bit-for-bit identical" in our regression tests)
when sediment transport is disabled.

Requirement: MPAS-Ocean performance essentially unchanged if sediment transport is disabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/25

Contributors: Xylar Asay-Davis, Zhendong Cao

The model should not take significantly more time than before the sediment transport capability was
added when the new capability is disabled. Timing will be required to increase by less than 1-2% in regressing timing
tests. The effect on timing should ideally be completely negligible.

When sediment transport is enabled, performance will necessarily be affected. For this implementation, we will not worry too much about this performance -- optimization will be left to follow-up work.

Algorithm Design
----------------

Algorithm Design: Bedload sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/16

Contributors: Zhendong Cao, Xylar Asay-Davis

We will implement three formulations to solve the bedload transport process, 
and users need to define which one to use (by setting :code:`config_sediment_bedload_formulation`)

- `Meyer-Peter-Mueller formulation (1948) <http://resolver.tudelft.nl/uuid:4fda9b61-be28-4703-ab06-43cdc2a21bd7>`_
  This formulation is for the unidirectional flow.

- `Soulsby-Damgaard formulation (2005) <https://doi.org/10.1016/j.coastaleng.2005.04.003>`_
  This formulation accounts for the combined effects of currents and waves.

- `Engelund-Hansen formulation (1967) <http://resolver.tudelft.nl/uuid:81101b08-04b5-4082-9121-861949c336c9>`_
  This formulation was derived for total transport but has been widely used for bedload transport in long-term sediment transport simulations.

In practice, a bedload coefficient (:math:`bedld_{coef}`) should be defined to adjust the bedlad transport rate;
besides, if multiple sediment classes are defined, the sediment class fraction (:math:`sed_{frac}`) should be also pre-defined.
(see `COAWST code here <https://github.com/jcwarner-usgs/COAWST/blob/master/ROMS/Nonlinear/Sediment/sed_bedload.F#L608-L612>`_).

The Meyer-Peter-Mueller and Soulsby-Dammgaard share the same expression of the bedload transport rate

.. math::

      q_b = \Phi \sqrt{(s-1)g D_{50}^3}\rho_s

in which :math:`q_b` is the bedload transport rate (:math:`kg\,m^{-1}\,s^{-1}` ), :math:`\Phi` is dimensionless bedload
transport rate, :math:`s=\rho_s-\rho_w` is the sediment specific density in water, :math:`\rho_s` and :math:`\rho_w`
are the sediment denisty and the (*ocean bottom layer*) water density, respectively (:math:`kg\,m^{-3}`), 
:math:`g` is the gravitational acceleration (:math:`m\,s^{-2}`), :math:`D_{50}` is the median grain diameter (:math:`m`).

The main difference of the two formulations is the computation of :math:`\Phi`, which is a function of another
dimensionless parameter called Shields parameter :math:`\theta`

.. math::

      \theta = \frac{\tau}{g(s-1)D_{50}}

in which :math:`\tau` is the magnitude of the bottom shear stress (:math:`m^2\,s^2`)

**Meyer-Peter-Mueller**

In Meyer-Peter-Mueller formulation, the dimensionless bedload transport rate is computed by

.. math::

      \Phi = \max \left( 8.0(\theta -\theta_{cr})^{1.5}, 0.0 \right)

with :math:`\theta_{cr}=0.047`.

**Soulsby-Damgaard**

In Soulsby-Damgaard formulation, the dimensionless bedload transport rate is computed by

.. math::

      \Phi = \max\left(12.0\theta ^{1/2}(\theta - \theta_{cr}) , 0.0\right)

in which :math:`\theta_{cr}=0.05`.

**Engelund-Hansen**

One form of the Engelund-Hansen formulation is

.. math::

      q_b = \frac{0.05\overline{U}^5\rho_s} {\sqrt{g}C^3(s-1)^2D_{50}}

in which :math:`\overline{U}` is the magnitude of the depth-averaged velocity (:math:`m\,s^{-1}`), :math:`C`
is the Chezy coefficient (:math:`m^{1/2}\,s^{-1}`).


Algorithm Design: Suspended sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/17

Contributors: Zhendong Cao, Xylar Asay-Davis

We will implement two empirical approaches and one process-based approach in MPAS-Ocean to solve 
the suspended sediment transport process.
The emprical approaches solve SSC profile by the classic Rouse profile assumption, which requires
the knowledge of the near-bottom reference SSC and the reference height.
The process-based approach solves the suspended sediment transport by the three-dimensional advection-diffusion
equation (i.e., the scaler transport equation),  with an additional source/sink term for sediment vertical
settling and exchange with the bed.

**Empirical approach**

When sediment is in suspension, the settling towards the bed is counterbalanced by the upward diffusion of the sediment due to
the turbulence eddies. Thus, SSC profile :math:`C(z)` can be obtained by solving the following equations

.. math::

      C(z) w_s + \epsilon_s \frac{dC(z)}{dz} = 0

in which :math:`C(z)` is the suspended sediment concentration (:math:`kg\,m^{-3}`),
:math:`w_s` is the sediment settling velocity (:math:`m\,s^{-1}`),
:math:`\epsilon_s` is the sediment diffusion coefficient (:math:`m^2\,s^{-1}`),
and :math:`z` is the vertical coordinate (:math:`m`) with positive upward.
Assuming the eddy diffusivity varies quadratically with height, the solution of
:math:`C(z)` can be written as

.. math::

      C(z) = C_{ref}\left[ \frac{z_{ref}}{z} \frac{H-z}{H-z_{ref}} \right]^b

in which :math:`C_{ref}` is the near-bottom reference SSC (:math:`kg\,m^{-3}`),
:math:`z_{ref}` is the reference height (:math:`m`) of :math:`C_{ref}`,
:math:`H` is the total water depth (:math:`m`),
and :math:`b` is the suspension parameter called Rouse number, which is calculated by

.. math::

      b = \frac{w_s}{\kappa u_*}

where :math:`\kappa` is the von Karman constant (0.4) and
:math:`u_*` is the shear friction velocity (:math:`m\,s^{-1}`).
The shear stress :math:`\tau = C_d \overline{U}^2={u_*}^2`, so :math:`u_*={C_d}^{0.5} \overline{U}`.

We provide two options to compute :math:`C_{ref}` and :math:`z_{ref}`:

- `Zyserman and Fredsoe, 1994 <https://ascelibrary.org/doi/pdf/10.1061/%28ASCE%290733-9429%281994%29120%3A9%281021%29>`_

  .. math::

      C_{ref} = \frac{0.331(\theta - 0.045)^{1.75}}{1+0.72(\theta-0.045)^{1.75}}

  and :math:`z_{ref}=2D_{50}`. Again, :math:`\theta` is the Shields parameter and
  :math:`D_{50}` is the sediment median grain diameter.

- `Goldstein et al., 2014 <https://esurf.copernicus.org/articles/2/67/2014/esurf-2-67-2014-discussion.html>`_

  .. math::

      C_{ref} = \left[ \frac{0.328U_b}{0.0688+1000D_{50}}\right]^2

  at height :math:`z_{ref}=0.01 m`, where :math:`U_b` is the magnitude of the bottom velocity (:math:`m\,s^{-1}`)

**Process-based approach**

The proces-based approach solves SSC by the three-dimensional advection-diffusion equations,
with an additional source/sink term for vertical settling and exchange with the bed

.. math::

      \frac{\partial{H(z)C(z)} }{\partial t} + \nabla \cdot \left( H(z)C(z) \bf{u}\right) = \nabla \cdot \left ( \epsilon_s \nabla H(z)C(z) \right) + C_{source}

with the addtional term written as

.. math::

      C_{source} = -w_s C(z) + E_s

in which :math:`H(z)` is the layer thickness (:math:`m`) in the ocean model,
:math:`t` is time (:math:`s`), :math:`\nabla` represents the advection process,
:math:`\bf{u}` is the velocity vector (:math:`m\,s^{-1}`),
:math:`E_s` is the erosion source (:math:`kg\,m^{-2}\,s^{-1}`), which is computed as

.. math::

      E_s = E_0(1-\phi)\frac{\tau_{sf}-\tau_{ce}}{\tau_{ce}}

where :math:`E_0` is the bed erodibility constant (:math:`kg\,m^{-2}\,s^{-1}`),
:math:`\phi` is the porosity of the bed layer,
:math:`\tau_{ce}` is the bed critical erosion stress (:math:`m^2\,s^{2}`)

**Sediment settling velocity**

Suspended sediment transport also requires the knowledge of the sediment settling velocity.
Sediment settling velocity is a function of the sediment geometry (e.g., particle size, shape and roundness),
the fluid characteristics (density and viscosity), and well as the fluid turbulence velocity.

Two sediment settling veloity formulations are selected to implement into MPAS-Ocean, including

- `Soulsby, 1997 <https://www.icevirtuallibrary.com/doi/abs/10.1680/doms.25844.fm>`_

  .. math::

      w_s = \frac{\nu}{D_{50}}\left[\sqrt{ (10.36^2+1.049D_*^3)}-10.36\right]

  where :math:`\nu` is the kinematic viscosity of water (:math:`m^2\,s^{-1}`).

- `Cheng, 1997 <https://doi.org/10.1061/(ASCE)0733-9429(1997)123:2(149)>`_

  .. math::

      w_s = \frac{\nu}{D_{50}} (\sqrt{25+1.2D_*^2}-5)^{1.5}

  in which :math:`D_* = D_{50}\left[\frac{(s-1)g}{\nu ^2}\right]^{1/3}` is dimensionless grain size.



Implementation
--------------

Implementation: Bedload sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/19

Contributors: Zhendong Cao, Xylar Asay-Davis

The sediment bedload transport requires the pre-definitions of the following namelist options:

- A logical variable to switch ON/OFF the bedload transport: :code:`config_sediment_bedload_enable`

- A character variable to define the bedload transport method: :code:`config_sediment_bedload_formulation`
  And there are three options for it:

  * :code:`Engelund-Hansen`

  * :code:`Meyer-Peter-Mueller`

  * :code:`Soulsby-Damgaard`

These can be defined in **Registry.xml** as follows:

.. code-block::

  <nml_record name="sediment_transport" mode="init;forward">
        <nml_option name="config_sediment_bedload_enable" type="logical" default_value=".false." units="unitless"
                    description="Controls if sediment bedload transport is computed."
                    possible_values=".true. or .false."
        />
        <nml_option name="config_sediment_bedload_formulation" type="chracter" default_value="Meyer-Peter-Mueller" units="unitless"
                    description="Select the sediment bedload transport formulation"
                    possible_values=" 'Engelund-Hansen', 'Meyer-Peter-Mueller', 'Soulsby-Damgaard'"
	/>
  </nml_record>

Bedload transport also requires the definitions of the following variables:

- Pre-defined sediment-related variables

  * sediment grain diameter (:math:`D_{50}`): :code:`sedimentGrainDiameter(nSedimentClasses, nCells)`

  * sediment grain density (:math:`\rho_s`): :code:`sedimentGrainDensity(nSedimentClasses, nCells)`

  * sediment class fraction on bed layer (:math:`sed_{frac}`): :code:`sedimentClassFraction(nSedimentClasses)`

  * bedload transport coefficient (:math:`bedld_{coef}`): :code:`sedimentBedloadCoefficient(nSedimentClasses)`

- Other variables

  * bedload transport rate (:math:`q_b`): :code:`sedimentBedloadTransportRate(nSedimentClasses,nCells)`


Implementation: Suspended sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/19

Contributors: Zhendong Cao, Xylar Asay-Davis

The following namelist options should be defined:

- A logical swith to turn ON/OFF suspended sediment transport: :code:`config_sediment_suspended_enable`

- A character-type variable to define the empirical suspended transport method: :code:`config_sediment_suspended_formulation`. There are three options for it:

  * :code:`None`

  * :code:`Zyserman-Fredsoe`

  * :code:`Goldstein`

  If this is set as :code:`None`, then the suspended transport is solved by the process-based approach.

- A character-type variable to define the method to compute the sediment settling velocity: :code:`config_sediment_settling_formulation`. There are three options for it:

  * :code:`None`

  * :code:`Soulsby`

  * :code:`Cheng`

  If this is set as :code:`None`, then the sediment settling velocity is pre-defined manually.

These can be defined in **Registry.xml** as follows:

.. code-block::

  <nml_record name="sediment_transport" mode="init;forward">
        <nml_option name="config_sediment_suspended_enable" type="logical" default_value=".false." units="unitless"
                    description="Controls if suspended sediment transport is computed."
                    possible_values=".true. or .false."
        />
        <nml_option name="config_sediment_suspended_formulation" type="chracter" default_value="None" units="unitless"
                    description="Select the reference SSC  formulation"
                    possible_values=" 'None', 'Zyserman-Fredsoe', 'Goldstein'"
	/>
        <nml_option name="config_sediment_settling_formulation" type="chracter" default_value="None" units="unitless"
                    description="Select the sediment settling velocity formulation"
                    possible_values=" 'None', 'Soulsby', 'Cheng'"
	/>
  </nml_record>

Suspended sediment transport also requires the definitions of the following variables:

- Pre-defined sediment-related variables

  * sediment settling velocity (:math:`w_s`): :code:`sedimentSettlingVelocity(nSedimentClasses)`

  * sediment critical erosion shear (:math:`\tau_{ce}`): :code:`sedimentCriticalErosionShear(nSedimentClasses)`

  * sediment erodibility constant (:math:`E_0`): :code:`sedimentErodibilityConstant(nSedimentClasses)`

  * sediment bed porosity (:math:`\psi`): :code:`sedimentBedPorosity(nSedimentClasses)`

- Other variables

  * sediment suspended concentration (:math:`C(z)`): :code:`sedimentSuspendedConcentration(nSedimentClasses, nCell)`

Testing
-------

Testing and Validation: name-of-topic-here (same as Requirement)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: YYYY/MM/DD

Contributors: (add your name to this list if it does not appear)

**To be continued ...**
