
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
The main model components include 
1) sediment classes; 
2) bedload sediment transport; 
3) suspended sediment transport; and 
4) morphological change. 
Sediment classes determine the properties of sediment-related variables. 
Bedload and suspended sediment transport are the two main transport forms considered in this development. 
Morphological change provides the trigger for interactions between different model components.

The primary implimentation challenges lie in the suspended sediment transport component.
It requires to solve the three-dimensional scalar equation for suspended sediment with 
its own source/sink terms, which need to resolve the deposition and resuspension at each
time step.

The sucessful implementation will enable MPAS-Ocean to simulate the sediment transport (especially
the computation of suspended sediment cocentration (SSC)) with given initial and boundary conditions,
as well as the mutual feedback between sediment transport, hydrodynamics and morphology.


Requirements
------------

Requirement: Sediment classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/14

Contributors: Zhendong Cao

We will enable the model the capability in representing multiple user-defined
sediment classes. Sediment classes determine the sediment properties such as
sediment grain size, grain density, settling veloity, crtical erosion shear
and erodibility constant. Apparently, the definitions of these sediment
properties require the assignment of the number of sediment classes.


Requirement: Bedload sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/14

Contributors: Zhendong Cao

Bedload sediment transport means the sediment particles have successive contacts with the bed during transport.
Bedload mass is exchanged horizontally at the bed layer.
The following three methods (formulations) are implemented to solve the bedload transport process:

- `Engelund-Hansen formulation (1967) <https://repository.tudelft.nl/islandora/object/uuid%3A81101b08-04b5-4082-9121-861949c336c9>`_ This formulation was derived for total transport but has been widely used for bedload transport in long-term sediment transport simulations.
- `Meyer-Peter-Mueller formulation (1948) <https://repository.tudelft.nl/islandora/object/uuid%3A4fda9b61-be28-4703-ab06-43cdc2a21bd7>`_ This formulation is for the unidirectional flow.
- `Soulsby-Damgaard formulation (2005) <https://doi.org/10.1016/j.coastaleng.2005.04.003>`_ This formulation accounts for the combined effects of currents and waves.

Requirement: Suspended sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/17

Contributors: Zhendong Cao

If water flows significantly beyond the threshold of the sediment motion, sediment will be entrained off the bed 
and into suspension. 
To remain in suspension, the upward turbulence eddies must be greater than the sediment settling velocity.
Suspended mass is exchanged vertically between the water column and the top bed layer.

Suspended sediment transport can be solved either by emprical approaches or by a process-based approach.
The emprical approaches solve SSC profile by the classic Rouse profile assumption. 
This approach requires the knowledge of the near-bottom reference SSC,
which can be computed by one of the following two methods:

- `Zyserman and Fredsoe, 1994 <https://ascelibrary.org/doi/pdf/10.1061/%28ASCE%290733-9429%281994%29120%3A9%281021%29>`_
- `Goldstein et al., 2014 <https://esurf.copernicus.org/articles/2/67/2014/esurf-2-67-2014-discussion.html>`_ 

The process-based approach solves the suspended sediment transport by the three-dimensional advection-diffusion equation 
(i.e., the scaler transport equation),  with an additional source/sink term for sediment vertical settling and exchange 
with the bed.

Suspended sediment transport requires the knowledge of the sediment settling velocity.
Sediment settling velocity is a function of the sediment geometry (e.g., particle size, shape and roundness), 
the fluid characteristics (density and viscosity), and well as the fluid turbulence velocity. 
Two sediment settling veloity formulations are selected to implement into MPAS-Ocean, including

- `Soulsby, 1997 <https://www.icevirtuallibrary.com/doi/abs/10.1680/doms.25844.fm>`_
- `Cheng, 1997 <https://doi.org/10.1061/(ASCE)0733-9429(1997)123:2(149)>`_


Requirement: Morphological change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/20

Contributors: Zhendong Cao

Morphological change here stands for the change of sea floor elevation (bathymetry) due to the convergence/divergence of net
sediment flux during a time period. The Exner equation (`Paola and Voller, 2005 <https://doi.org/10.1029/2004JF000274>`_) 
is implemented to compute the morphological change.


Algorithm Design
---------------------------

Algorithm Design: Sediment classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/16

Contributors: Zhendong Cao

The number of sediment classes will be defined manually and considered as a dimension of sediment-related variables.

Algorithm Design: Bedload sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/16

Contributors: Zhendong Cao

In practice, a bedload coefficient should be defined to adjust the bedlad transport rate (see `COAWST code here 
<https://github.com/jcwarner-usgs/COAWST/blob/master/ROMS/Nonlinear/Sediment/sed_bedload.F#L608-L612>`_). 

The Meyer-Peter-Mueller and Soulsby-Dammgaard share the same expression of the bedload transport rate

.. math::
    q_b = \Phi \sqrt{(s-1)g D_{50}^3}\rho_s

in which :math:`q_b` is the bedload transport rate (:math:`kg m^{-1}s^{-1}` ), :math:`\Phi` is dimensionless bedload 
transport rate, :math:`s=\rho_s-\rho_w` is the sediment specific density in water, :math:`\rho_s` and :math:`\rho_w` 
are the sediment denisty and the water density, respectively (:math:`kg m^{-3}`), :math:`g` is the gravitational 
acceleration (:math:`m s^{-2}`), :math:`D_{50}` is the median grain diameter (:math:`m`).

The main difference of the two formulations is the computation of :math:`\Phi`, which is a function of another
dimensionless parameter called Shields parameter :math:`\theta`

.. math::
	\theta = \frac{\tau}{g(s-1)D_{50}}

in which :math:`\tau` is the magnitude of the bottom shear stress (:math:`m^2 s^2`)

**Meyer-Peter-Mueller** 

In Meyer-Peter-Mueller formulation, the dimensionless bedload transport rate is computed by

.. math::
	\Phi = \max(8.0(\theta -\theta_{cr})^{1.5}, 0.0)

with :math:`\theta_{cr}=0.047`.

**Soulsby-Damgaard**

In Soulsby-Damgaard formulation, the dimensionless bedload transport rate is computed by

.. math::
	\Phi_{x1} = \max(12.0\theta ^{1/2}(\theta - \theta_{cr}) , 0.0)
.. math::

	\Phi_{x2} = 12.0(0.9534+0.1907\cos(2\psi))\theta_w^{0.5}\theta_m + 12(0.229\gamma_w\theta_w^{1.5}\cos\psi)

.. math::
	\Phi_x = \max[\Phi_{x1},\Phi_{x2}, 0.0]

.. math::
	\Phi_y = \frac{12(0.1907\theta_w^2)} {\theta_w^{1.5}+1.5 \theta_m^{1.5}} (\theta_m \sin (2\psi))+1.2\gamma_w\theta_w\sin\psi)

.. math::
	\tau_m = \tau_c \left(1+1.2\left(\frac{\tau_w}{\tau_c+\tau_w} \right) \right)

in which :math:`\theta_{cr}=0.05`, :math:`\gamma_w` is the asymmetry factor (constrained to be less than 0.2), 
:math:`\psi` is the angle between current direction and wave travelling direction, :math:`\theta_w` and :math:`\theta_m`
are wave-induced Shields parameter and the mean Shields paramter, respectively (i.e. computed by :math:`\tau_w` and 
:math:`\tau_m`, respectively, using the Shields parameter equation).
:math:`\tau_c` and :math:`\tau_w` are bottom shear stress from the currents and waves, respectively.
Notice that :math:`\Phi_x = \Phi_y =0` if :math:`\theta_m \leq \theta_{cr}`. 


**Engelund-Hansen**

One form of the Engelund-Hansen formulation is

.. math::
  q_b = \frac{0.05\overline{U}^5\rho_s} {\sqrt{g}C^3(s-1)^2D_{50}}

in which :math:`\overline{U}` is the magnitude of the depth-averaged velocity (:math:`m/s`), :math:`C` 
is the Chezy coefficient (:math:`m^{1/2} s^{-1}`).


Algorithm Design: Suspended sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/17

Contributors: Zhendong Cao

**Empirical approach**

When sediment is in suspension, the settling towards the bed is counterbalanced by the upward diffusion of the sediment due to
the turbulence eddies. Thus, SSC profile :math:`C(z)` can be obtained by solving the following equations

.. math::
	C(z) w_s + \epsilon_s \frac{dC(z)}{dz} = 0

in which :math:`C(z)` is the suspended sediment concentration (:math:`kg m^{-3}`), :math:`w_s` is the sediment settling velocity 
(:math:`ms^{-1}`), :math:`\epsilon_s` is the sediment diffusion coefficient (:math:`m^2s^{-1}`), and :math:`z` is the vertical 
coordinate (:math:`m`) with positive upward. Assuming the eddy diffusivity varying parabolically with height, the solution of
:math:`C(z)` can be written as

.. math::
	C(z) = C_{ref}\left[ \frac{z_{ref}}{z} \frac{H-z}{H-z_{ref}} \right]^b

in which :math:`C_{ref}` is the near-bottom reference SSC (:math:`kg m^{-3}`), 
:math:`z_{ref}` is the reference height (:math:`m`) of :math:`C_{ref}`, 
:math:`H` is the total water depth (:math:`m`),
and :math:`b` is the suspension parameter called Rouse number, which is calculated by

.. math::
	b = \frac{w_s}{\kappa u_*}

where :math:`\kappa` is the von Karman constant (0.4) and :math:`u_*` is the shear friction velocity (:math:`ms^{-1}`). 
The shear stress :math:`\tau = C_d \overline{U}^2={u_*}^2`, so :math:`u_*={C_d}^{0.5} \overline{U}`.

We provde two options to compute :math:`C_{ref}` and :math:`z_{ref}`:

- *Zyserman and Fredsoe, 1994*

.. math::
	C_{ref} = \frac{0.331(\theta - 0.045)^{1.75}}{1+0.72(\theta-0.045)^{1.75}}

and :math:`z_{ref}=2D_{50}`. Again, :math:`\theta` is the Shields parameter and 
:math:`D_{50}` is the sediment median grain diameter.

- *Goldstein et al., 2014*

.. math::
	C_{ref} = \left[ \frac{0.328U_b}{0.0688+1000D_{50}}\right]^2

at height :math:`z_{ref}=0.01 m`, where :math:`U_b` is the magnitude of the bottom velocity (:math:`ms^{-1}`)

**Process-based approach**

The proces-based approach solves SSC by the three-dimensional advection-diffusion equations, with an additional source/sink term 
for vertical settling and exchange with the bed

.. math::
	\frac{\partial{H(z)C(z)} }{\partial t} + \nabla \cdot \left( H(z)C(z) \bf{u}\right) = \nabla \cdot \left ( \epsilon_s \nabla H(z)C(z) \right) + C_{source}

with the addtional term written as

.. math::
	C_{source} = -w_s C(z) + E_s

in which :math:`H(z)` is the layer thickness (:math:`m`), :math:`t` is time (:math:`s`), :math:`\nabla` represents 
the advection process, :math:`\bf{u}` is the velocity vector (:math:`ms^{-1}`), 
:math:`E_s` is the erosion source (:math:`kg m^{-2} s^{-1}`), which is computed as

.. math::
	E_s = E_0(1-\phi)\frac{\tau_{sf}-\tau_{ce}}{\tau_{ce}}

where :math:`E_0` is the bed erodibility constant (:math:`kg m^{-2}s^{-1}`), :math:`\phi` is the porosity of the bed layer,
:math:`\tau_{ce}` is the bed critical erosion stress (:math:`m^2s^{2}`)

The methods to compute sediment settling velocity include

- *Soulsby 1997*

.. math::
	w_s = \frac{\nu}{D_{50}}\left[\sqrt{ (10.36^2+1.049D_*^3)}-10.36\right]

where :math:`\nu` is the kinematic viscosity of water (:math:`m^2s^{-1}`)

- *Cheng 1997*

.. math::
	w_s = \frac{\nu}{D_{50}} (\sqrt{25+1.2D_*^2}-5)^{1.5}

in which :math:`D_* = D_{50}\left[\frac{(s-1)g}{\nu ^2}\right]^{1/3}` is dimensionless grain size.


Algorithm Design: Morphological change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/20

Contributors: Zhendong Cao

Exner equation reads as

.. math::
	(1-\phi)\frac{\partial Z_b}{\partial t} + \left(\frac{\partial q_{bx}}{\partial x}+\frac{\partial q_{by}}{\partial y}\right) = - \frac{C_{source,b}}{\rho_s}

in which :math:`\phi` is bed porosity, :math:`Z_b` is the bathymetry (:math:`m`), :math:`q_{bx}` and :math:`q_{by}` are bottom
bedload transport flux (:math:`m^2s^{-1}`) in :math:`x` and :math:`y` direction, respectively; :math:`C_{source,b}` is the
bottom suspended transport rate (:math:`kg m^{-2}s^{-1}`).


Implementation
--------------

Implementation: Sediment classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/14

Contributors: Zhendong Cao

The number of sediment classes `nSedimentClasses` will be defined as a dimension before the 
definitions of the sediment properties. The definition will be in a `Registry.xml` file and the code is

.. code::

	dim name='nSedimentClasses' unit='unitless'
	    description='the number of sediment classes considered in the model'

Implementation: Bedload sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/19

Contributors: Zhendong Cao

The sediment bedload transport requires the pre-definitions of the following namelist options:

- A logical variable to switch ON/OFF the bedload transport: :code:`config_sediment_bedload`

- A character variable to define the bedload transport method: :code:`config_sediment_bedload_formulation` 
  And there are three options for it:
	* :code:`Engelund-Hansen`
	* :code:`Meyer-Peter-Mueller`
	* :code:`Soulsby-Damgaard`

These can be defined in **Registry.xml** as follows:

.. code::

  <nml_record name="sediment_transport" mode="init;forward">
        <nml_option name="config_sediment_bedload" type="logical" default_value=".false." units="unitless"
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
	* sediment grain diameter: :code:`sedimentGrainDiameter(nSedimentClasses, nCells)`
	* sediment grain density: :code:`sedimentGrainDensity(nSedimentClasses, nCells)`
	* sediment class fraction on bed layer: :code:`sedimentClassFraction(nSedimentClasses)`
	* bedload transport coefficient: :code:`sedimentBedloadCoefficient(nSedimentClasses)`
- Other variables
	* bedload transport rate: :code:`sedimentBedloadTransportRate(nCells)`


Implementation: Suspended sediment transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/19

Contributors: Zhendong Cao

The following namelist options should be defined:

- A logical swith to turn ON/OFF suspended sediment transport: :code:`config_sediment_suspended`.
- A character-type variable to define the empirical suspended transport method: 
  :code:`config_sediment_suspended_formulation`. There are three options for it:
	* :code:`None`
	* :code:`Zyserman-Fredsoe`
	* :code:`Goldstein`
  
  If this is set as :code:`None`, then the suspended transport is solved by the process-based approach.

- A character-type variable to define the method to compute the sediment settling velocity: 
  :code:`config_sediment_settling_formulation`. There are three options for it:
	* :code:`None`
	* :code:`Soulsby`
	* :code:`Cheng`

  If this is set as :code:`None`, then the sediment settling velocity is pre-defined manually.

These can be defined in **Registry.xml** as follows:

.. code::

  <nml_record name="sediment_transport" mode="init;forward">
        <nml_option name="config_sediment_suspended" type="logical" default_value=".false." units="unitless"
                    description="Controls if sediment bedload transport is computed."
                    possible_values=".true. or .false."
        />
        <nml_option name="config_sediment_suspended_formulation" type="chracter" default_value="None" units="unitless"
                    description="Select the sediment bedload transport formulation"
                    possible_values=" 'None', 'Zyserman-Fredsoe', 'Goldstein'"
	/>
  </nml_record>

Suspended sediment transport also requires the definitions of the following variables:

- Pre-defined sediment-related variables:
	* sediment settling velocity: :code:`sedimentSettlingVelocity(nSedimentClasses)`
	* sediment critical erosion shear: :code:`sedimentCriticalErosionShear(nSedimentClasses)`
	* sediment erodibility constant: :code:`sedimentErodibilityConstant(nSedimentClasses)`
	* sediment bed porosity: :code:`sedimentBedPorosity(nSedimentClasses)`
- Other variables:
	* sediment suspended concentration: :code:`sedimentSuspendedConcentration(nSedimentClasses, nCell)`

Implementation: Morphological change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: 2021/01/20

Contributors: Zhendong Cao

A logical variable need to be defined to swith ON/OFF the morphological change:
:code:`config_morphological_change`.

A user-specified bed layer thickness should be defined to account for the erodible bed sediment mass:
:code:`sedimentBedLayerThickness`

For long-term simulation, a morphological scale factor also need defined to accelerate the morphological change:
:code:`config_morphological_scale_factor`.
With this implementation, the bedload transport flux, erosion, and deposition rates will be multiplied by the
defined scale factor at each time step.

These can be defined in **Registry.xml** as follows:

.. code::

  <nml_record name="sediment_transport" mode="init;forward">
        <nml_option name="config_morphological_change" type="logical" default_value=".false." units="unitless"
                    description="Controls if morphological change is taken account."
                    possible_values=".true. or .false."
        />
        <nml_option name="config_morphological_scale_factor" type="real" default_value="1.0" units="unitless"
                    description="Acceleration rate of the morphological change"
                    possible_values=" Any value equal or larger than 1.0; 1.0 means no acceleration."
        />
  </nml_record>


Testing
-------

Testing and Validation: name-of-topic-here (same as Requirement)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Date last modified: YYYY/MM/DD

Contributors: (add your name to this list if it does not appear)

***To be continued ...**
