!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_types.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_types

  !> Holds type definitions for the derived types used by each 
  !> instance of the ice model. Originally, each of these types
  !> was a module containing variables, which were used as containers
  !> for global variables. However, the need to allow for multiple
  !> ice model instances meant that the nested derived types were instituted
  !> instead. However, there is probably one too many levels in this scheme. 
  !> It would be better if the different types here were contained in the 
  !> higher-level instance type (\texttt{glint\_instance}), rather than 
  !> the intermediate model type (\texttt{glide\_global\_type}). 
  !> 
  !> Note that this \emph{is} now where the defaults are defined for these
  !> variables.

!TODO - Clean up the glide_global type so it holds fewer subtypes?
!       For example, we could replace some work types (tempwk, velowk) with local arrays and parameters.

  use glimmer_sparse_type
  use glimmer_global, only: sp, dp
  use glimmer_ncdf
  use profile
  use glimmer_coordinates, only: coordsystem_type
  use glimmer_map_types
  use glimmer_physcon

  implicit none

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Constants that describe the options available
  ! We use these integer parameters elsewhere in the code to avoid
  !  hardwiring of option numbers

  ! basic Glimmer/Glide options

  integer, parameter :: GLOBAL_BC_PERIODIC = 0  ! doubly periodic
  integer, parameter :: GLOBAL_BC_OUTFLOW = 1   ! free outflow; scalars in global halo set to zero
  integer, parameter :: GLOBAL_BC_NO_PENETRATION = 2   ! no penetration; outflow set to zero at global boundaries
                                                       ! NOTE: The no-penetration option is currently supported for Glissade only

  integer, parameter :: DYCORE_GLIDE = 0     ! old shallow-ice dycore from Glimmer
  integer, parameter :: DYCORE_GLAM = 1      ! Payne-Price finite-difference solver
  integer, parameter :: DYCORE_GLISSADE = 2  ! prototype finite-element solver
  integer, parameter :: DYCORE_ALBANYFELIX = 3  ! External Albany-Felix finite-element solver
  integer, parameter :: DYCORE_BISICLES = 4     ! BISICLES-Chombo external FVM solver

  integer, parameter :: EVOL_PSEUDO_DIFF = 0    ! glide only
  integer, parameter :: EVOL_ADI = 1            ! glide only
  integer, parameter :: EVOL_DIFFUSION = 2      ! glide only
  integer, parameter :: EVOL_INC_REMAP = 3      ! glam/glissade only
  integer, parameter :: EVOL_UPWIND = 4         ! glam/glissade only
  integer, parameter :: EVOL_NO_THICKNESS = 5   ! glam/glissade only

  !NOTE: Use option 1 for prognostic temperature with any dycore
  !      Option 3 is under construction

  integer, parameter :: TEMP_SURFACE_AIR_TEMP = 0
  integer, parameter :: TEMP_PROGNOSTIC = 1
  integer, parameter :: TEMP_STEADY = 2
  integer, parameter :: TEMP_ENTHALPY = 3

  integer, parameter :: TEMP_INIT_ZERO = 0
  integer, parameter :: TEMP_INIT_ARTM = 1
  integer, parameter :: TEMP_INIT_LINEAR = 2

  integer, parameter :: FLWA_CONST_FLWA = 0
  integer, parameter :: FLWA_PATERSON_BUDD_CONST_TEMP = 1
  integer, parameter :: FLWA_PATERSON_BUDD = 2

  integer, parameter :: BTRC_ZERO = 0
  integer, parameter :: BTRC_CONSTANT = 1
  integer, parameter :: BTRC_CONSTANT_BWAT = 2
  integer, parameter :: BTRC_CONSTANT_TPMP = 3
  integer, parameter :: BTRC_LINEAR_BMLT = 4
  integer, parameter :: BTRC_TANH_BWAT = 5

  integer, parameter :: BWATER_NONE  = 0
  integer, parameter :: BWATER_LOCAL = 1
  integer, parameter :: BWATER_FLUX  = 2
  integer, parameter :: BWATER_CONST = 3
  integer, parameter :: BWATER_OCEAN_PENETRATION = 4   ! effective pressure calculation with pw=ocean pressure for grounding line parameterisation (Leguy, et al., TC, 2014)
  !integer, parameter :: BWATER_BASAL_PROC = 4  ! not currently supported

  integer, parameter :: BMLT_FLOAT_NONE = 0
  integer, parameter :: BMLT_FLOAT_CONSTANT = 1
  integer, parameter :: BMLT_FLOAT_MISMIP = 2

  integer, parameter :: BASAL_MBAL_NO_CONTINUITY = 0
  integer, parameter :: BASAL_MBAL_CONTINUITY = 1

  integer, parameter :: GTHF_UNIFORM = 0
  integer, parameter :: GTHF_PRESCRIBED_2D = 1
  integer, parameter :: GTHF_COMPUTE = 2

  integer, parameter :: RELAXED_TOPO_NONE = 0     ! Do nothing
  integer, parameter :: RELAXED_TOPO_INPUT = 1    ! Input topo is relaxed
  integer, parameter :: RELAXED_TOPO_COMPUTE = 2  ! Input topo in isostatic equilibrium
                                                  ! compute relaxed topo

  integer, parameter :: ISOSTASY_NONE = 0
  integer, parameter :: ISOSTASY_COMPUTE = 1

  integer, parameter :: LITHOSPHERE_LOCAL = 0
  integer, parameter :: LITHOSPHERE_ELASTIC = 1

  integer, parameter :: ASTHENOSPHERE_FLUID = 0
  integer, parameter :: ASTHENOSPHERE_RELAXING = 1

  integer, parameter :: CALVING_NONE = 0
  integer, parameter :: CALVING_FLOAT_ZERO = 1
  integer, parameter :: CALVING_FLOAT_FRACTION = 2
  integer, parameter :: CALVING_RELX_THRESHOLD = 3
  integer, parameter :: CALVING_TOPG_THRESHOLD = 4
  integer, parameter :: CALVING_HUYBRECHTS = 5
  integer, parameter :: CALVING_THCK_THRESHOLD = 6
  integer, parameter :: CALVING_DAMAGE = 7

  !WHL - added an option to determine whether calving occurs at initialization
  integer, parameter :: CALVING_INIT_OFF = 0
  integer, parameter :: CALVING_INIT_ON = 1

  !WHL - added an option to determine whether calving can occur everywhere the calving
  !      criterion is met, or only at the ocean edge
  integer, parameter :: CALVING_DOMAIN_OCEAN_EDGE = 0
  integer, parameter :: CALVING_DOMAIN_EVERYWHERE = 1
  integer, parameter :: CALVING_DOMAIN_OCEAN_CONNECT = 2

  integer, parameter :: VERTINT_STANDARD = 0
  integer, parameter :: VERTINT_KINEMATIC_BC = 1

  integer, parameter :: SIGMA_COMPUTE_GLIDE = 0
  integer, parameter :: SIGMA_EXTERNAL = 1
  integer, parameter :: SIGMA_CONFIG = 2
  integer, parameter :: SIGMA_COMPUTE_EVEN = 3
  integer, parameter :: SIGMA_COMPUTE_PATTYN = 4

  integer, parameter :: RESTART_FALSE = 0
  integer, parameter :: RESTART_TRUE = 1

  integer, parameter :: RESTART_EXTEND_VELO_FALSE = 0
  integer, parameter :: RESTART_EXTEND_VELO_TRUE = 1
  
  !basal proc option disabled for now
  integer, parameter :: BAS_PROC_DISABLED = 0
!!  integer, parameter :: BAS_PROC_FULLCALC = 1
!!  integer, parameter :: BAS_PROC_FASTCALC = 2

  ! higher-order options

  integer, parameter :: HO_EFVS_CONSTANT = 0
  integer, parameter :: HO_EFVS_FLOWFACT = 1
  integer, parameter :: HO_EFVS_NONLINEAR = 2

  integer, parameter :: HO_DISP_NONE = -1
  integer, parameter :: HO_DISP_SIA = 0
  integer, parameter :: HO_DISP_FIRSTORDER = 1

  integer, parameter :: HO_BABC_CONSTANT = 0
  integer, parameter :: HO_BABC_BETA_TPMP = 1
  integer, parameter :: HO_BABC_YIELD_PICARD = 2
  integer, parameter :: HO_BABC_BETA_BWAT = 3
  integer, parameter :: HO_BABC_LARGE_BETA = 4
  integer, parameter :: HO_BABC_EXTERNAL_BETA = 5
  integer, parameter :: HO_BABC_NO_SLIP = 6
  integer, parameter :: HO_BABC_YIELD_NEWTON = 7
  integer, parameter :: HO_BABC_ISHOMC = 8
  integer, parameter :: HO_BABC_POWERLAW = 9
  integer, parameter :: HO_BABC_COULOMB_FRICTION = 10
  integer, parameter :: HO_BABC_COULOMB_CONST_BASAL_FLWA = 11
  integer, parameter :: HO_BABC_COULOMB_POWERLAW_TSAI = 12
  integer, parameter :: HO_BABC_SIMPLE = 13

  integer, parameter :: HO_NONLIN_PICARD = 0
  integer, parameter :: HO_NONLIN_JFNK = 1

  integer, parameter :: HO_RESID_MAXU = 0
  integer, parameter :: HO_RESID_MAXU_NO_UBAS = 1
  integer, parameter :: HO_RESID_MEANU = 2
  integer, parameter :: HO_RESID_L2NORM = 3
  integer, parameter :: HO_RESID_L2NORM_RELATIVE = 4

  integer, parameter :: HO_SPARSE_PCG_INCH = -1
  integer, parameter :: HO_SPARSE_BICG = 0
  integer, parameter :: HO_SPARSE_GMRES = 1
  integer, parameter :: HO_SPARSE_PCG_STANDARD = 2
  integer, parameter :: HO_SPARSE_PCG_CHRONGEAR = 3
  integer, parameter :: HO_SPARSE_TRILINOS = 4

  integer, parameter :: HO_APPROX_LOCAL_SIA = -1
  integer, parameter :: HO_APPROX_SIA = 0
  integer, parameter :: HO_APPROX_SSA = 1
  integer, parameter :: HO_APPROX_BP = 2
  integer, parameter :: HO_APPROX_L1L2 = 3
  integer, parameter :: HO_APPROX_DIVA = 4

  integer, parameter :: HO_PRECOND_NONE = 0
  integer, parameter :: HO_PRECOND_DIAG = 1
  integer, parameter :: HO_PRECOND_SIA  = 2

  integer, parameter :: HO_GRADIENT_CENTERED = 0
  integer, parameter :: HO_GRADIENT_UPSTREAM = 1

  integer, parameter :: HO_GRADIENT_MARGIN_ALL = 0
  integer, parameter :: HO_GRADIENT_MARGIN_ICE_LAND = 1
  integer, parameter :: HO_GRADIENT_MARGIN_ICE_ONLY = 2

  integer, parameter :: HO_VERTICAL_REMAP_FIRST_ORDER = 0
  integer, parameter :: HO_VERTICAL_REMAP_SECOND_ORDER = 1

  integer, parameter :: HO_ASSEMBLE_BETA_STANDARD = 0
  integer, parameter :: HO_ASSEMBLE_BETA_LOCAL = 1

  integer, parameter :: HO_ASSEMBLE_TAUD_STANDARD = 0
  integer, parameter :: HO_ASSEMBLE_TAUD_LOCAL = 1

  integer, parameter :: HO_ASSEMBLE_BFRIC_STANDARD = 0
  integer, parameter :: HO_ASSEMBLE_BFRIC_LOCAL = 1

  integer, parameter :: HO_GROUND_NO_GLP = 0
  integer, parameter :: HO_GROUND_GLP = 1
  integer, parameter :: HO_GROUND_ALL = 2

  integer, parameter :: HO_FLOTATION_FUNCTION_PATTYN = 0
  integer, parameter :: HO_FLOTATION_FUNCTION_INVERSE_PATTYN = 1
  integer, parameter :: HO_FLOTATION_FUNCTION_LINEAR = 2

  integer, parameter :: HO_ICE_AGE_NONE = 0 
  integer, parameter :: HO_ICE_AGE_COMPUTE = 1 

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_general

    !> Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !> The number of grid-points in the E-W direction.
    integer :: nsn = 0  !> The number of grid-points in the N-S direction.
    integer :: upn = 1  !> The number of vertical levels in the model.

    type(coordsystem_type) :: ice_grid  !> coordinate system of the ice grid
    type(coordsystem_type) :: velo_grid !> coordinate system of the velocity grid

    real(dp), dimension(:),pointer :: x0 => null() !original x0 grid 
    real(dp), dimension(:),pointer :: y0 => null() !original y0 grid
    real(dp), dimension(:),pointer :: x1 => null() !original x1 grid
    real(dp), dimension(:),pointer :: y1 => null() !original y1 grid

    integer :: global_bc = 0     ! 0 for periodic, 1 for outflow, 2 for no-penetration

  end type glide_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_options

    !> Holds user options controlling the methods used in the ice-model integration.

    !-----------------------------------------------------------------------
    ! standard options
    !-----------------------------------------------------------------------

    integer :: whichdycore = 2

    ! Choice of two Glimmer dycores:
    !> \begin{description} 
    !> \item[0] Glide dycore (SIA, serial (SLAP) only)
    !> \item[1] SEACISM/Glam dycore (1st-order, FDM, serial (SLAP) or parallel (Trilinos))
    !> \item[2] Glissade dycore (1st-order, FEM, serial (SLAP) or parallel (F90 native PCG solver) )
    !> \item[3] FELIX-Albany dycore (1st-order, FEM, using Trilino/Albany, mesh information from Glissade)
    !> \item[4] BISICLES dycore (L1L2, FVM, parallel using Chombo AMR)
    !> \end{description}

    integer :: whichevol = 0

    !> Thickness evolution method:
    !> \begin{description}
    !> \item[0] Pseudo-diffusion 
    !> \item[1] Alternating direction implicit (ADI)
    !> \item[2] Diffusion (also calculates velocities) 
    !> \item[3] Incremental remapping
    !> \item[4] 1st-order upwind scheme
    !> \item[5] Temperature evolves but thickness does not
    !> \end{description}

    integer :: whichtemp = 1

    !> Method of ice temperature calculation:
    !> \begin{description} 
    !> \item[0] Set column to surface air temperature
    !> \item[1] Prognostic temperature solution 
    !> \item[2] Do NOTHING - hold temperatures steady at initial value  
    !> \item[3] Prognostic enthalpy solution
    !> \end{description}

    integer :: temp_init = 1

    ! Temperature initialization:
    !> \begin{description} 
    !> \item[0] Initialize temperature to 0 C
    !> \item[1] Initialize temperature to surface air temperature
    !> \item[2] Initialize temperature with a linear profile in each column
    !> \end{description}

    !> Method for calculating flow factor $A$:

    integer :: whichflwa = 2

    !> \begin{description} 
    !> \item[0] Set equal to $1\times 10^{-16}\,\mathrm{yr}^{-1}
    !> \item[1] \emph{Paterson and Budd} relationship, 
    !> with temperature set to $-5^{\circ}\mathrm{C}$ 
    !> \item[2] \emph{Paterson and Budd} relationship 
    !> \,\mathrm{Pa}^{-n}$
    !> \end{description}

    integer :: whichbtrc = 0

    !> Basal slip coefficient:
    !> \begin{description}
    !> \item[0] Set equal to zero everywhere
    !> \item[1] Set to (non--zero) constant
    !> \item[2] Set to (non--zero) constant where basal water is present, otherwise to zero
    !> \item[3] Set to (non--zero) constant where temperature is at pressure melting point of ice, otherwise to zero
    !> \item[4] linear function of basal melt rate
    !> \item[5] \texttt{tanh} function of basal water depth 
    !> \end{description}

    integer :: whichbwat = 0

    !> Basal water depth: 
    !> \begin{description} 
    !> \item[0] Set to zero everywhere 
    !> \item[1] Compute from local basal water balance 
    !> \item[2] Compute the basal water flux, then find depth via calculation
    !> \item[3] Set to constant (10 m) everywhere, to force T = Tpmp.
    !> \item[4] Calculated from till water content, in the basal processes module
    !> \end{description}

    integer :: whichbmlt_float = 0

    !> basal melt rate for floating ice:
    !> \begin{description}
    !> \item[0] Basal melt rate = 0 for floating ice
    !> \item[1] Basal melt rate = constant for floating ice (with option to selectively mask out melting)
    !> \item[2] Basal melt rate for floating ice as prescribed for MISMIP+
    !> \end{description}

    integer :: basal_mbal = 0

    !> basal mass balance:
    !> \begin{description}
    !> \item[0] Basal mass balance not included in continuity equation
    !> \item[1] Basal mass balance included in continuity equation
    !> \end{description}

    integer :: gthf = 0

    !> geothermal heat flux:
    !> \begin{description}
    !> \item[0] prescribed uniform geothermal heat flux
    !> \item[1] read 2D geothermal flux field from input file (if present)
    !> \item[2] calculate geothermal flux using 3d diffusion
    !> \end{description}

    ! This replaces model%isos%do_isos
    integer :: isostasy = 0

    !> isostasy:
    !> \begin{description}
    !> \item[0] no isostatic adjustment
    !> \item[1] compute isostatic adjustment using lithosphere/asthenosphere model
    !> \end{description}

    !TODO - Should whichrelaxed move from the options to the isostasy section?
    integer :: whichrelaxed = 0

    !> relaxed topography:
    !> \begin{description}
    !> \item[0] get relaxed topo from separate variable (in practice, do nothing)
    !> \item[1] first time slice of input topo is relaxed
    !> \item[2] first time slice of input topo is in isostatic equilibrium
    !> \end{description}

    integer :: whichcalving = 1

    !> Calving: 
    !> \begin{description} 
    !> \item[0] No calving 
    !> \item[1] Set thickness to zero if floating 
    !> \item[2] Lose a fraction of floating ice at marine margin
    !> \item[3] Set thickness to zero if relaxed bedrock is more than a
    !>          certain water depth (variable "marine_limit" in glide_types)  
    !> \item[4] Set thickness to zero if present bedrock topography lies below
    !>          a certain water depth (variable "marine_limit" in glide_types)  
    !> \item[5] Huybrechts grounding line scheme for Greenland initialization
    !> \item[6] Set thickness to zero if ice at marine margin is thinner than
    !>          a certain value (variable 'calving_minthck' in glide_types)
    !> \item[7] Calve ice that is sufficiently damaged
    !> \end{description}

    integer :: calving_init = 0
    !> \begin{description}
    !> \item[0] Do not calve at initialization
    !> \item[1] Calve at initialization
    !> \end{description}

    integer :: calving_domain = 0
    !> \begin{description}
    !> \item[0] Calve only at ocean edge
    !> \item[1] Calve wherever the calving criterion is met
    !> \item[2] Calve where the calving criterion is met, and there is a connected path
    !>          to the ocean through other cells where the criterion is met.
    !> \end{description}

    integer :: whichwvel = 0

    !> Vertical velocities: 
    !> \begin{description}
    !> \item[0] Usual vertical integration 
    !> \item[1] Vertical integration constrained so that 
    !> upper kinematic B.C. obeyed 
    !> \end{description}

    integer :: which_sigma = 0

    !> \begin{description}
    !> \item[0] compute standard Glimmer sigma coordinates
    !> \item[1] sigma coordinates are given in external file
    !> \item[2] sigma coordinates are given in configuration file
    !> \item[3] evenly spaced levels, as required for glam dycore
    !> \item[4] compute Pattyn sigma coordinates
    !> \end{description}

    !TODO - Make is_restart a logical variable?

    integer :: is_restart = 0
    !> if the run is a restart of a previous run
    !> \begin{description}
    !> \item[0] normal start-up (take init fields from .nc input file OR if absent, use default options)
    !> \item[1] restart model from previous run (do not calc. temp, rate factor, or vel)
    !> \end{description}

    integer :: restart_extend_velo = 0
    !> if velocity fields should be written on the extended staggered mesh
    !> \begin{description}
    !> \item[0] write uvel and vvel to restart file on standard staggered mesh
    !> \item[1] write uvel_extend and vvel_extend to restart file on extended staggered mesh
    !>          (required if restart velocities are nonzero on global boundaries)
    !> \end{description}

    ! This is a Glimmer serial option
    ! The parallel code enforces periodic EW and NS boundary conditions by default
    logical :: periodic_ew = .false.

    !> \begin{description}
    !> \item[0] no periodic EW boundary conditions
    !> \item[1] periodic EW boundary conditions
    !> \end{description}

    !-----------------------------------------------------------------------
    ! Higher-order options
    ! Associated with Payne-Price dycore (glam) and newer glissade dycore
    !-----------------------------------------------------------------------

    integer :: which_ho_efvs = 2

    !> Flag that indicates how effective viscosity is computed
    !> \begin{description}
    !> \item[0] constant value
    !> \item[1] multiple of flow factor
    !> \item[2] compute from effective strain rate

    integer :: which_ho_disp = 1

    !> Flag that indicates method for computing the dissipation during the temperature calc.
    !> \begin{description}
    !> \item[-1] for no dissipation
    !> \item[0] for 0-order SIA approx
    !> \item[1] for first-order dissipation (Blatter-Pattyn)
    !>      
    !> \end{description}

    integer :: which_ho_babc = 4

    !> Flag that describes basal boundary condition for HO dyn core: 
    !> \begin{description}
    !> \item[0] spatially uniform value (low value of 10 Pa/yr by default)
    !> \item[1] large value for frozen bed, lower value for bed at pressure melting point
    !> \item[2] treat beta value as a till yield stress (in Pa) using Picard iteration 
    !> \item[3] linear (inverse) function of bwat 
    !> \item[4] very large value for beta to enforce no slip everywhere 
    !> \item[5] beta field passed in from .nc input file as part of standard i/o
    !> \item[6] no slip everywhere (using Dirichlet BC rather than large beta)
    !> \item[7] treat beta value as till yield stress (in Pa) using Newton-type iteration (in development)
    !> \item[8] beta field as prescribed for ISMIP-HOM test C (serial only)
    !> \item[9] power law based using effective pressure
    !> \item[10] Coulomb friction law using effective pressure, with flwa from lowest ice layer
    !> \item[11] Coulomb friction law using effective pressure, with constant basal flwa
    !> \item[12] basal stress is the minimum of Coulomb and power-law values, as in Tsai et al. (2015)
    !> \item[13] simple hard-coded pattern (useful for debugging)
    !> \end{description}

    integer :: which_ho_nonlinear = 0
    !> Flag that indicates method for solving the nonlinear iteration when solving 
    !> the first-order momentum balance
    !> \item[0] use the standard Picard iteration
    !> \item[1] use Jacobian Free Newton Krylov (JFNK) method

    integer :: which_ho_resid = 3
    !> Flag that indicates method for computing residual in PP dyn core: 
    !> \begin{description}
    !> \item[0] maxval 
    !> \item[1] maxval ignoring basal velocity 
    !> \item[2] mean value
    !> \item[3] L2 norm of system residual, Ax-b=resid
    !> \item[4] L2 norm of system residual relative to rhs, |Ax-b|/|b|
    !> \begin{description}

    integer :: which_ho_sparse = 0
    !> Flag that indicates method for solving the sparse linear system
    !> that arises from the higher-order solver
    !> \begin{description}
    !> \item[-1] SLAP (serial): Preconditioned conjugate gradient, incomplete Cholesky preconditioner
    !> \item[0]  SLAP (serial): Biconjugate gradient, incomplete LU preconditioner
    !> \item[1]  SLAP (serial): GMRES, incomplete LU preconditioner
    !> \item[2] Native PCG, parallel-enabled, standard solver
    !> \item[3] Native PCG, parallel-enabled, Chronopoulos-Gear solver
    !> \item[4] standalone interface to Trilinos
    !> \end{description}

    ! parameters to store external dycore options/information -- Doug Ranken 04/20/12
    integer*4 :: external_dycore_type = 0
    integer*4 :: external_dycore_model_index = -1  
    !> Flag to select an external dynamic core.
    !> \begin{description}
    !> \item[0] Do not use an external dynamic core
    !> \item[1] Use the BISICLES external dynamic core
    !> \item[2] Use the ALBANY_FELIX external dynamic core
    !> \end{description}

    character(fname_length) :: dycore_input_file=''
    !FD Name of a file containing external dycore settings.

    integer :: which_ho_approx = 2 
    !> Flag that indicates which Stokes approximation to use with the glissade dycore.
    !> Not valid for other dycores 
    !> Compute Blatter-Pattyn HO momentum balance by default.
    !> Note: There are two SIA options:
    !>       Option -1 uses module glissade_velo_sia to compute local SIA velocities, similar to Glide
    !>       Option 0 uses module glissade_velo_higher to compute SIA velocities via an iterative solve
    !> \begin{description}
    !> \item[-1] Shallow-ice approximation, Glide-type calculation; uses glissade_velo_sia
    !> \item[0]  Shallow-ice approximation, vertical-shear stresses only; uses glissade_velo_higher
    !> \item[1]  Shallow-shelf approximation, horizontal-plane stresses only; uses glissade_velo_higher
    !> \item[2]  Blatter-Pattyn approximation with both vertical-shear and horizontal-plane stresses; uses glissade_velo_higher
    !> \item[3]  Vertically integrated 'L1L2' approximation with vertical-shear and horizontal-plane stresses; uses glissade_velo_higher
    !> \item[4]  Depth-integrated viscosity approximation based on Goldberg (2011); uses glissade_velo_higher 
    !> \end{description}

    integer :: which_ho_precond = 2    
    !> Flag that indicates which Stokes preconditioner to use in the glissade dycore.
    !> Not valid for other dycores 
    !> \begin{description}
    !> \item[0] No preconditioner
    !> \item[1] Diagonal preconditioner
    !> \item[2] Physics-based shallow-ice preconditioner
    !> \end{description}

    integer :: which_ho_gradient = 0    
    !> Flag that indicates which gradient operator to use in the glissade dycore.
    !> Not valid for other dycores
    !> NOTE: Option 1 may be better for ice evolution because it damps checkerboard noise.
    !> \begin{description}
    !> \item[0] Centered gradient
    !> \item[1] Upstream gradient

    integer :: which_ho_gradient_margin = 1
    !> Flag that indicates how to compute the gradient at the ice margin in the glissade dycore.
    !> Not valid for other dycores
    !> \begin{description}
    !> \item[0] Use info from all neighbor cells, ice-covered or ice-free
    !> \item[1] Use info from ice-covered and/or land cells, not ice-free ocean
    !> \item[2] Use info from ice-covered cells only

    !TODO: Change the default to 2nd order vertical remapping
    ! WHL: Keeping this 1st order for now so that standard tests are BFB
    integer :: which_ho_vertical_remap = 0
    !> Flag that indicates the order of accuracy for vertical remapping
    !> \begin{description}
    !> \item[0] first-order accurate in the vertical direction
    !> \item[1] second-order accurate in the vertical direction

    integer :: which_ho_assemble_beta = 0

    !> Flag that describes how beta terms are assembled in the glissade finite-element calculation
    !> \begin{description}
    !> \item[0] standard finite-element calculation (which effectively smooths beta at discontinuities)
    !> \item[1] apply local value of beta at each vertex

    integer :: which_ho_assemble_taud = 0

    !> Flag that describes how driving-stress terms are assembled in the glissade finite-element calculation
    !> \begin{description}
    !> \item[0] standard finite-element calculation (which effectively smooths the driving stress)
    !> \item[1] apply local value of driving stress at each vertex

    integer :: which_ho_assemble_bfric = 0

    !> Flag that describes how the basal friction heat flux is computed in the glissade finite-element calculation
    !> \begin{description}
    !> \item[0] standard finite-element calculation summing over quadrature points
    !> \item[1] apply local value of beta*(u^2 + v^2) at each vertex

    integer :: which_ho_ground = 0
    !> Flag that indicates how to compute the grounded fraction of each gridcell in the glissade dycore.
    !> Not valid for other dycores
    !> \begin{description}
    !> \item[0] fground = 0 in floating cells (based on flotation condition), else fground = 1 
    !> \item[1] 0 <= fground <= 1, based on a grounding line parameterization
    !> \item[2] fground = 1 in all cells

    !TODO - Change default to linear function 2?
    integer :: which_ho_flotation_function = 1
    !> Flag that indicates how to compute the flotation function at and near vertices in the glissade dycore
    !> Not valid for other dycores
    !> \begin{description}
    !> \item[0] f_flotation = (-rhow*b/rhoi*H) = f_pattyn; <=1 for grounded, > 1 for floating
    !> \item[1] f_flotation = (rhoi*H)/(-rhow*b) = 1/f_pattyn; >=1 for grounded, < 1 for floating
    !> \item[2] f_flotation = -rhow*b - rhoi*H = ocean cavity thickness; <=0 for grounded, > 0 for floating 

    integer :: which_ho_ice_age = 1    
    !> Flag that indicates whether to compute a 3d ice age tracer
    !> \item[0] ice age computation off
    !> \item[1] ice age computation on

    integer :: glissade_maxiter = 100    
    !> maximum number of nonlinear iterations to be used by the Glissade velocity solver

    ! The remaining options are not currently supported

    !integer :: which_bproc = 0
    !Options for the basal processes code
    !> \begin{description}
    !> \item[0] Disabled
    !> \item[1] Full calculation, with at least 3 nodes to represent the till layer
    !> \item[2] Fast calculation, using Tulaczyk empirical parametrization
    !> \end{description}

    !integer :: use_plume = 0   !! Option to be supported in future releases
    !> \begin{description}
    !> \item[0] standard bmlt calculation
    !> \item[1] use plume to calculate bmlt
    !> \end{description}

  end type glide_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geometry

    !> Holds fields and other information relating to the
    !> geometry of the ice sheet and bedrock.

    real(dp),dimension(:,:),pointer :: thck => null()
    !> The thickness of the ice, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: usrf => null()
    !> The elevation of the upper ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: lsrf => null() 
    !> The elevation of the lower ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: topg => null() 
    !> The elevation of the topography, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: f_flotation => null() 
    !> flotation function, (rhoi*thck) / (-rhoo*(topg-eus))
    !> previously was f_pattyn = -rhoo*(topg-eus)/(rhoi*thck)
    !    (computed by glissade dycore only)

    real(dp),dimension(:,:),pointer :: f_ground => null() 
    !> The fractional area at each vertex which is grounded 
    !    (computed by glissade dycore only)

    real(dp),dimension(:,:,:),pointer :: ice_age => null()
    !> The age of a given ice layer, divided by \texttt{tim0}.
    !> Used to be called 'age', but changed to 'ice_age' for easier grepping

    integer :: ntracers
    !> number of tracers to be transported

    real(dp),dimension(:,:,:,:), pointer :: tracers => null()
    !> all 3D tracers packed into one array

    real(dp),dimension(:,:,:), pointer :: tracers_usrf => null()
    !> value at upper surface for all tracers

    real(dp),dimension(:,:,:), pointer :: tracers_lsrf => null()
    !> value at lower surface for all tracers

    integer, dimension(:,:),pointer :: thkmask => null()
    !> see glide_mask.f90 for possible values

    integer, dimension(:,:),pointer :: stagmask => null()
    !> see glide_mask.f90 for possible values

    !TODO - Consider moving BISICLES variables to their own type at some point
    !* (DFM ----------------- added for BISICLES interface --------------)
    real(dp),dimension(:,:),pointer :: floating_mask => null()
    !*(DFM) Real-valued mask indicated where ice is grounded or floating

    !* (DFM ----------------- added for BISICLES interface --------------)
    real(dp),dimension(:,:),pointer :: ice_mask => null()
    !*(DFM) Real-valued mask indicating where ice is present or absent


    !* (DFM ----------------- added for BISICLES interface --------------)
    real(dp),dimension(:,:),pointer :: lower_cell_loc => null()
    !*(DFM) The z-location of the center of the lowest ice cell center

    !* (DFM ----------------- added for BISICLES interface --------------)
    real(dp),dimension(:,:),pointer :: lower_cell_temp => null()
    !*(DFM) The temperature in the cell located at lower_cell_loc

    integer, dimension(:,:),pointer :: thck_index => null()
    ! Set to nonzero integer for ice-covered cells (thck > 0), cells adjacent to ice-covered cells,
    !  and cells with acab > 0.  The non-zero points are numbered in sequence from the bottom left 
    !  to the top right, going along the rows.

    integer :: totpts = 0       ! total number of points with nonzero thck_index
    logical :: empty = .true.   ! true if totpts = 0

    real(dp) :: ivol, iarea,iareag, iareaf !> ice volume and ice area

  end type glide_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geomderv

    !> Holds the horizontal and temporal derivatives of the thickness and
    !> upper surface elevation, as well as the thickness on the staggered grid.

    !*tb* Added a bunch of stuff here to clean up the higher order code that
    !I've been writing.  Might be worth it to add a mechanism to conditionally
    !allocate these depending on whether they are needed by the SIA core or by
    !the higher-order extensions

    !First derivatives on a staggered grid
    real(dp),dimension(:,:),pointer :: dthckdew => null() !> E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !> E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !> N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !> N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dlsrfdew => null() !*tb* added
    real(dp),dimension(:,:),pointer :: dlsrfdns => null() !*tb* added

    !Second derivatives on a staggered grid
    !*tb* added all of these
    ! Used by glam_strs2
    real(dp),dimension(:,:),pointer :: d2usrfdew2 => null()
    real(dp),dimension(:,:),pointer :: d2usrfdns2 => null()
    real(dp),dimension(:,:),pointer :: d2thckdew2 => null()
    real(dp),dimension(:,:),pointer :: d2thckdns2 => null()

    !Time derivatives
    real(dp),dimension(:,:),pointer :: dthckdtm => null() !> Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm => null() !> Temporal derivative of upper surface elevation.

    !TODO - Move staggered variables from glide_geomderv type to glide_geometry?

    !Staggered grid versions of geometry variables
    real(dp),dimension(:,:),pointer :: stagthck => null() !> Thickness averaged onto the staggered grid.
    real(dp),dimension(:,:),pointer :: stagusrf => null() !> Upper surface averaged onto the staggered grid
    real(dp),dimension(:,:),pointer :: staglsrf => null() !> Lower surface averaged onto the staggered grid
    real(dp),dimension(:,:),pointer :: stagtopg => null() !> Bedrock topography averaged onto the staggered grid

  end type glide_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_tensor
    real(dp), dimension(:,:,:), pointer :: scalar => null()
    real(dp), dimension(:,:,:), pointer :: xz => null()
    real(dp), dimension(:,:,:), pointer :: yz => null()
    real(dp), dimension(:,:,:), pointer :: xx => null()
    real(dp), dimension(:,:,:), pointer :: yy => null()
    real(dp), dimension(:,:,:), pointer :: xy => null()
  end type glide_tensor
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velocity

    !> Holds the velocity fields in 2D and 3D. At least some of these fields
    real(dp),dimension(:,:,:),pointer :: uvel  => null()   !> 3D $x$-velocity.
    real(dp),dimension(:,:,:),pointer :: vvel  => null()   !> 3D $y$-velocity.
    real(dp),dimension(:,:,:),pointer :: velnorm => null() ! horizontal ice speed
    real(dp),dimension(:,:,:),pointer :: wvel  => null()   !> 3D $z$-velocity.
    real(dp),dimension(:,:,:),pointer :: wgrd  => null()   !> 3D grid vertical velocity.
    real(dp),dimension(:,:,:),pointer :: wvel_ho  => null()!> 3D $z$-velocity.from higher-order dycores
    real(dp),dimension(:,:)  ,pointer :: uflx  => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: vflx  => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: diffu => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: diffu_x => null() !*sfp* moved from velocity_hom deriv type
    real(dp),dimension(:,:)  ,pointer :: diffu_y => null() 
    real(dp),dimension(:,:)  ,pointer :: total_diffu => null() !> total diffusivity

    real(dp),dimension(:,:)  ,pointer :: uvel_2d  => null()   !> 2D vertically averaged $x$-velocity.
    real(dp),dimension(:,:)  ,pointer :: vvel_2d  => null()   !> 2D vertically averaged $y$-velocity.
    real(dp),dimension(:,:)  ,pointer :: ubas  => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: ubas_tavg  => null()
    real(dp),dimension(:,:)  ,pointer :: vbas  => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: vbas_tavg  => null() 

    !! next 3 used for output of residual fields (when relevant code in glam_strs2 is active)
!    real(dp),dimension(:,:,:),pointer :: ures => null() !> 3D $x$-residual.
!    real(dp),dimension(:,:,:),pointer :: vres  => null() !> 3D $y$-residual.
!    real(dp),dimension(:,:,:),pointer :: magres  => null() !> 3D $magnitude$-residual.

    ! Note: uvel_extend and vvel_extend can be used for input and output of uvel, vvel on a staggered grid 
    !       that is the same size as the unstaggered grid. This is required for exact restart if velocities
    !       are nonzero along the north and east boundaries of the global domain.
    real(dp),dimension(:,:,:),pointer :: uvel_extend => null()  !> 3D $x$-velocity on extended staggered grid
    real(dp),dimension(:,:,:),pointer :: vvel_extend => null()  !> 3D $y$-velocity on extended staggered grid
    real(dp),dimension(:,:)  ,pointer :: uvel_2d_extend => null()  !> 2D $x$-velocity on extended staggered grid
    real(dp),dimension(:,:)  ,pointer :: vvel_2d_extend => null()  !> 2D $y$-velocity on extended staggered grid

    real(dp),dimension(:,:)  ,pointer :: bed_softness => null() !> bed softness parameter
    real(dp),dimension(:,:)  ,pointer :: btrc  => null()        !>  basal traction (scaler field)
    real(dp),dimension(:,:,:),pointer :: btraction => null()    !> x(1,:,:) and y(2,:,:) "consistent" basal traction fields 
    real(dp),dimension(:,:)  ,pointer :: beta  => null()        !> basal shear coefficient on velo grid (Pa yr/m)
    real(dp),dimension(:,:)  ,pointer :: beta_internal => null()!> beta weighted by f_ground or otherwise adjusted (glissade only)
    real(dp),dimension(:,:)  ,pointer :: unstagbeta  => null()  !> basal shear coefficient on ice grid (Pa yr/m)
    real(dp),dimension(:,:)  ,pointer :: tau_x => null()        !> SIA basal shear stress, x-dir
    real(dp),dimension(:,:)  ,pointer :: tau_y => null()        !> SIA basal shear stress, y-dir

    !WHL - A reasonable value of beta_grounded_min might be 10 Pa yr/m.  
    !      However, this choice is not BFB for the confined-shelf test case, so I am choosing a default value of 0 for now.
    !      The default can be overridden in the config file.
    !TODO: Set beta_grounded_min = 10?
    real(dp) :: beta_grounded_min = 0.d0     !> minimum value of beta for grounded ice, Pa yr/m (glissade only; scaled during init)
    real(dp) :: ho_beta_const = 10.d0        !> spatially uniform beta for HO dycores, Pa yr/m (scaled during init)

    !> mask that specifies where the velocity being read in should be held constant as a dirichlet condition
    integer, dimension(:,:), pointer  :: kinbcmask => null()    

    !> masks that specify where the outflow velocities on the global boundary should be set to zero
    integer, dimension(:,:), pointer  :: umask_no_penetration => null()
    integer, dimension(:,:), pointer  :: vmask_no_penetration => null()

    !*sfp* mask on vel grid showing which dyn bc is applied at each grid cell (mainly for debugging)
    integer, dimension(:,:), pointer    :: dynbcmask => null()    

    ! for viewing the spatial pattern of residuals
    real(dp),dimension(:,:,:),pointer :: resid_u => null()     ! u component of residual Ax - b where x is the velocity
    real(dp),dimension(:,:,:),pointer :: resid_v => null()     ! v component of residual Ax - b where x is the velocity

    ! for viewing the driving stress on the RHS
    real(dp),dimension(:,:,:),pointer :: rhs_u => null()     ! u component of b in Ax = b
    real(dp),dimension(:,:,:),pointer :: rhs_v => null()     ! v component of b in Ax = b

  end type glide_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_stress_t      

    type(glide_tensor) :: tau    ! HO only
    real(dp),dimension(:,:,:),pointer :: efvs => null()    !> effective viscosity
    real(dp),dimension(:,:),  pointer :: btractx => null() !> basal traction (Pa), x comp
    real(dp),dimension(:,:),  pointer :: btracty => null() !> basal traction (Pa), y comp
    !WHL - The extended versions are needed for exact restart if using DIVA solver for a problem with nonzero traction at global boundaries
    real(dp),dimension(:,:),  pointer :: btractx_extend => null() !> basal traction (Pa), x comp, on extended staggered grid
    real(dp),dimension(:,:),  pointer :: btracty_extend => null() !> basal traction (Pa), y comp, on extended staggered grid
    real(dp),dimension(:,:),  pointer :: taudx => null()   !> driving stress (Pa), x comp
    real(dp),dimension(:,:),  pointer :: taudy => null()   !> driving stress (Pa), y comp

  end type glide_stress_t      

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!TODO - Make eus a config file parameter.
!TODO - Rename acab in glide_climate type to avoid confusion over units? (e.g., acab_ice?)
!       Here, acab has units of m/y ice, whereas in Glint, acab has units of m/y water equiv.

  type glide_climate
     !> Holds fields used to drive the model
     real(dp),dimension(:,:),pointer :: acab      => null() !> Annual mass balance (m/y ice)
     real(dp),dimension(:,:),pointer :: acab_tavg => null() !> Annual mass balance (time average).
     real(dp),dimension(:,:),pointer :: artm      => null() !> Annual mean air temperature (degC)

     real(dp) :: eus = 0.d0                                !> eustatic sea level

  end type glide_climate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_calving
     !> holds fields and paramters related to calving
     !> Note: The 3D damage field is prognostic; the 2D damage_column field is diagnosed from the 3D damage field.
     real(dp),dimension(:,:),  pointer :: calving_thck => null()   !> thickness loss in grid cell due to calving
                                                                   !< scaled by thk0 like mass balance, thickness, etc.
     real(dp),dimension(:,:,:),pointer :: damage => null()         !> 3D damage tracer, 0 > damage < 1
     real(dp),dimension(:,:),  pointer :: damage_column => null()  !> 2D vertically integrated damage tracer, 0 > damage_column < 1
  
     real(dp) :: marine_limit =    -200.d0  !> minimum value of topg/relx before floating ice calves
                                            !> (whichcalving = CALVING_RELX_THRESHOLD, CALVING_TOPG_THRESHOLD)
     real(dp) :: calving_fraction = 0.2d0   !> fractional thickness of floating ice that calves
                                            !> (whichcalving = CALVING_FLOAT_FRACTION)
                                            !> WHL - previously defined as the fraction of floating ice that does not calve
     real(dp) :: calving_timescale = 0.0d0  !> time scale (yr) for calving (Glissade only); calving_thck = thck * max(dt/calving_timescale, 1)
                                            !> if calving_timescale = 0, then calving_thck = thck
     real(dp) :: calving_minthck = 100.d0   !> minimum thickness (m) of floating ice at marine edge before it calves
                                            !> (whichcalving = CALVING_THCK_THRESHOLD)
     real(dp) :: damage_threshold = 1.0d0   !> threshold at which ice column is deemed sufficiently damaged to calve
                                            !> assuming that 0 = no damage, 1 = total damage

  end type glide_calving

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type eismint_climate_type

     ! holds parameters for the eismint climate

     ! For EISMINT2:
     ! airt(1) = Tmin = summit surface temperature (K)
     ! airt(2) = S_T = horizontal temperature gradient (K/m)
     ! nmsb(1) = M_max = max accumulation (m/yr)
     ! nmsb(2) = S_b = horizontal smb gradient (m/yr/m)
     ! nmsb(3) = R_el = radial distance from summit where mass balance = 0 (m)   
     !

     integer :: eismint_type = 0
     !> select EISMINT experiment
     !> \begin{description}
     !> \item[{\bf 1}] EISMINT-1 fixed margin
     !> \item[{\bf 2}] EISMINT-1 moving margin
     !> \item[{\bf 3}] EISMINT-2
     !> \item[{\bf 4}] MISMIP-1 (not EISMINT but has similar climate parameters)
     !> \item[{\bf 5}] Exact verification (not EISMINT but has similar climate parameters)
     !> \end{description}

     ! NOTE: The initial nmsb values in the declarations below are appropriate
     !       for EISMINT-2, but the initial airt values are not.
     ! TODO: Change default airt values in eismint_type to be consistent with EISMINT-2?

     !> air temperature parameterisation K, K km$^{-3}$
     real(dp), dimension(2) :: airt = (/ -3.15d0, 1.d-2 /)  

     !> mass balance parameterisation:
     real(dp), dimension(3) :: nmsb = (/ 0.5d0, 1.05d-5, 450.0d3 /)

     !> EISMINT time-dep climate forcing period, switched off when set to 0
     real(dp) :: period = 0.d0

     !> EISMINT amplitude of mass balance time-dep climate forcing
     real(dp) :: mb_amplitude = 0.2d0

  end type eismint_climate_type


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_temper

    !> Holds fields relating to temperature.

    !Note: In the Glide dycore, temp, flwa and dissip live on the unstaggered vertical grid
    !       at layer interfaces and have vertical dimension (1:upn).
    !      In the Glam/Glissade dycore, with remapping advection of temperature, 
    !       temp, flwa and dissip live on the staggered vertical grid at layer midpoints.  
    !       The vertical dimensions are (0:upn) for temp and (1:upn-1) for flwa and dissip.
    !
    !      bheatflx, ucondflx, and lcondflx are defined as positive down,
    !       so they will often be < 0.  
    !      However, bfricflx and dissipcol are defined to be >= 0.
    !
    !      If bheatflx is read from a data file, be careful about the sign!
    !      In input data, the geothermal heat flux is likely to be defined as positive upward.

    real(dp),dimension(:,:,:),pointer :: temp => null()      !> 3D temperature field.
    real(dp),dimension(:,:),  pointer :: bheatflx => null()  !> basal heat flux (W/m^2) (geothermal, positive down)
    real(dp),dimension(:,:,:),pointer :: flwa => null()      !> Glen's flow factor $A$.
    real(dp),dimension(:,:,:),pointer :: dissip => null()    !> interior heat dissipation rate, divided by rhoi*Ci (deg/s)
    real(dp),dimension(:,:),  pointer :: bwat => null()      !> Basal water depth
    real(dp),dimension(:,:),  pointer :: bwatflx => null()   !> Basal water flux 
    real(dp),dimension(:,:),  pointer :: stagbwat => null()  !> Basal water depth on velo grid
    real(dp),dimension(:,:),  pointer :: bmlt_ground =>null()!> Basal melt-rate for grounding ice (> 0 for melt, < 0 for freeze-on)
    real(dp),dimension(:,:),  pointer :: bmlt_float => null()!> Basal melt rate for floating ice (> 0 for melt, < 0 for freeze-on) 
    real(dp),dimension(:,:),  pointer :: stagbtemp => null() !> Basal temperature on velo grid
    real(dp),dimension(:,:),  pointer :: bpmp => null()      !> Basal pressure melting point
    real(dp),dimension(:,:),  pointer :: stagbpmp => null()  !> Basal pressure melting point on velo grid
    real(dp),dimension(:,:),  pointer :: bfricflx => null()  !> basal heat flux (W/m^2) from friction (>= 0)
    real(dp),dimension(:,:,:),pointer :: waterfrac => null() !> fractional water content in layer (0 <= waterfrac <= 1)
    real(dp),dimension(:,:,:),pointer :: enthalpy => null()  !> specific enthalpy in layer (J m-3)
                                                             !> = rhoi * Ci * T for cold ice
    !TODO - Remove ucondflx, lcondflx, dissipcol and make these local to glissade_therm?
    !       Probably cannot remove ucondflx because it may be needed for coupling.
    real(dp),dimension(:,:),  pointer :: ucondflx => null()  !> conductive heat flux (W/m^2) at upper sfc (positive down)
    real(dp),dimension(:,:),  pointer :: lcondflx => null()  !> conductive heat flux (W/m^2) at lower sfc (positive down)
    real(dp),dimension(:,:),  pointer :: dissipcol => null() !> total heat dissipation rate (W/m^2) in column (>= 0)
    integer  :: niter   = 0   
    real(dp) :: perturb = 0.d0
    real(dp) :: grid    = 0.d0 
    integer  :: tpt     = 0      !> Pointer to time series data
    logical  :: first1  = .true. !>
    logical  :: newtemps = .false. !> new temperatures

    ! parameters and fields for MISMIP+ experiments with basal melting
    ! Note: Parameters with units yr^{-1} are scaled to s^{-1} in subroutine glide_scale_params
    real(dp) :: bmlt_float_omega = 0.2d0           !> time scale for basal melting (yr-1)
                                                   !> default value = 0.2 yr^{-1} for MISIMP+
    real(dp) :: bmlt_float_h0 = 75.d0              !> scale for sub-shelf cavity thickness (m)
                                                   !> default value = 75 m for MISMIP+
    real(dp) :: bmlt_float_z0 = -100.d0            !> scale for ice draft, relative to sea level (m)
                                                   !> default value = -100 m for MISMIP+
    real(dp) :: bmlt_float_rate = 100.d0           !> constant melt rate (m/yr)
                                                   !> default value = 100 m/yr for MISMIP+ experiment Ice2r
    integer, dimension(:,:), pointer :: bmlt_float_mask => null()   !> switch off melt where mask = 1
                                                                    !> mask = 1 where x < 480 km for MISMIP+ experiment Ice2r

  end type glide_temper

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_basal_physics
      !< Holds variables related to basal physics associated with ice dynamics

      ! see glissade_basal_traction.F90 for usage details
      ! Note: It may make sense to move effecpress to a hydrology model when one is available.
      real(dp), dimension(:,:), pointer :: effecpress => null()          !< effective pressure  
      real(dp), dimension(:,:), pointer :: effecpress_stag => null()     !< effective pressure on staggered grid
      real(dp), dimension(:,:), pointer :: C_space_factor => null()      !< spatial factor for basal shear stress (no dimension)
      real(dp), dimension(:,:), pointer :: C_space_factor_stag => null() !< spatial factor for basal shear stress on staggered grid (no dimension)
      real(dp) :: friction_powerlaw_k = 8.4d-9    !< the friction coefficient for the power-law friction law (m y^-1 Pa^-2).  
                                                  !< The default value is from Bindschadler (1983) based on fits to observations, converted to CISM units.

      ! Parameters for Coulomb friction sliding law (default values from Pimentel et al. 2010)
      real(dp) :: Coulomb_C = 0.42d0              !< basal stress constant (no dimension)
                                                  !< Pimentel et al. have Coulomb_C = 0.84*m_max, where m_max = Coulomb_Bump_max_slope
      real(dp) :: Coulomb_bump_wavelength = 2.0d0 !< bed rock wavelength at subgrid scale precision (m)
      real(dp) :: Coulomb_bump_max_slope = 0.5d0  !< maximum bed bump slope at subgrid scale precision (no dimension) 
      real(dp) :: flwa_basal = 1.0d-16            !< Glen's A at the bed for the Schoof (2005) Coulomb friction law, in units Pa^{-n} yr^{-1} 
                                                  !< = 3.1688d-24 Pa{-n} s{-1}, the value used by Leguy et al. (2014)

      ! parameters for power law, taub_b = C * u_b^(1/m); used for HO_BABC_COULOMB_POWERLAW_TSAI
      ! The default values are from Asay-Davis et al. (2015).
      ! The value of powerlaw_C suggested by Tsai et al. (2015) is 7.624d6 Pa m^(-1/3) s^(1/3).
      ! This value can be converted to CISM units by dividing by scyr^(1/3), to obtain 2.413d4 Pa m^(-1/3) yr^(1/3).
      ! Note: The Tsai et al. Coulomb friction law uses Coulomb_C above, with
      !       effective pressure N as in Leguy et al. (2014) with p_ocean_penetration = 1.
      ! 
      real(dp) :: powerlaw_C = 1.0d4              !< friction coefficient in power law, units of Pa m^(-1/3) yr^(1/3)
      real(dp) :: powerlaw_m = 3.d0               !< exponent in power law (unitless)
      
  end type glide_basal_physics

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_lithot_type
     !> holds variables for temperature calculations in the lithosphere

     real(dp),dimension(:,:,:),pointer :: temp => null()    !> Three-dimensional temperature field.
     logical, dimension(:,:), pointer :: mask => null()     !> whether the point has been ice covered at some time

     integer :: num_dim = 1                                 !> either 1 or 3 for 1D/3D calculations

     ! The sparse matrix and linearised arrays
     type(sparse_matrix_type) :: fd_coeff, fd_coeff_slap
     integer :: all_bar_top
     real(dp), dimension(:), pointer :: rhs
     real(dp), dimension(:), pointer :: answer
     real(dp), dimension(:), pointer :: supd,diag,subd

     ! work arrays for solver
     real(dp), dimension(:), pointer :: rwork
     integer, dimension(:), pointer :: iwork
     integer mxnelt

     real(dp), dimension(:), pointer :: deltaz => null()    !> array holding grid spacing in z
     real(dp), dimension(:,:), pointer :: zfactors => null()!> array holding factors for finite differences of vertical diffu
     real(dp) :: xfactor,yfactor !> factors for finite differences of horizontal diffu


     real(dp) :: surft = 2.d0          !> surface temperature, used for calculating initial temperature distribution
     real(dp) :: mart  = 2.d0          !> sea floor temperature 
     integer  :: nlayer = 20           !> number of layers in lithosphere
     real(dp) :: rock_base = -5000.d0  !> depth below sea-level at which geothermal heat gradient is applied
     
     integer :: numt = 0        !> number time steps for spinning up GTHF calculations

     real(dp) :: rho_r = 3300.0d0 !> The density of lithosphere (kg m$^{-3}$)
     real(dp) :: shc_r = 1000.0d0 !> specific heat capcity of lithosphere (J kg$^{-1}$ K$^{-1}$)
     real(dp) :: con_r = 3.3d0    !> thermal conductivity of lithosphere (W m$^{-1}$ K$^{-1}$)

     real(dp) :: diffu = 0. !> diffusion coefficient

  end type glide_lithot_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type isos_elastic
     !> Holds data used by isostatic adjustment calculations

     real(dp) :: d = 0.24d25                !> flexural rigidity  !TODO - What are units of d?
     real(dp) :: lr                         !> radius of relative stiffness
     real(dp) :: a                          !> radius of disk
     real(dp) :: c1,c2,cd3,cd4              !> coefficients
     real(dp), dimension(:,:), pointer :: w !> matrix operator for lithosphere deformation
     integer :: wsize                       !> size of operator (0:rbel_wsize, 0:rbel_wsize), operator is axis symmetric
  end type isos_elastic

  type isostasy_type
     !> contains isostasy configuration

     integer :: lithosphere = 0
     !> method for calculating equilibrium bedrock depression
     !> \begin{description}
     !> \item[0] local lithosphere, equilibrium bedrock depression is found using Archimedes' principle
     !> \item[1] elastic lithosphere, flexural rigidity is taken into account
     !> \end{description}

     integer :: asthenosphere = 0
     !> method for approximating the mantle
     !> \begin{description}
     !> \item[0] fluid mantle, isostatic adjustment happens instantaneously
     !> \item[1] relaxing mantle, mantle is approximated by a half-space
     !> \end{description}

     real(dp) :: relaxed_tau = 4000.d0    ! characteristic time constant of relaxing mantle (yr)
     real(dp) :: period = 500.d0          ! lithosphere update period (yr)
     real(dp) :: next_calc                ! when to update lithosphere
     logical :: new_load = .false.        ! set to true if there is a new surface load
     type(isos_elastic) :: rbel           ! structure holding elastic lithosphere setup

     real(dp),dimension(:,:),pointer :: relx => null()  ! elevation of relaxed topography, by \texttt{thck0}.
     real(dp),dimension(:,:),pointer :: load => null()  ! load imposed on lithosphere
     real(dp),dimension(:,:),pointer :: load_factors => null() ! temporary used for load calculation

  end type isostasy_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_funits
    character(fname_length) :: sigfile=''                      !> sigma coordinates file
    character(fname_length) :: ncfile=''                       !> configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !> first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !> first element of linked list defining netCDF inputs
    type(glimmer_nc_input), pointer :: frc_first=>NULL()       !> first element of linked list defining netCDF forcings
    ! Note: forcing files are of the same type as input files since they share a lot in common.
  end type glide_funits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_numerics

    !> Parameters relating to the model numerics
    real(dp) :: tstart =     0.d0 !> starting time
    real(dp) :: tend   =  1000.d0 !> end time
    real(dp) :: time   =     0.d0 !> main loop counter in years
    real(dp) :: tinc   =     1.d0 !> time step of main loop in years 
    real(dp) :: ntem   =     1.d0 !> multiplier of main time step; allows longer temperature time step
    real(dp) :: alpha  =    0.5d0 !> richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas =    0.5d0 !> was a parameter in the original
    real(dp) :: thklim =   100.d0 ! min thickness for computing ice dynamics (m) 
    real(dp) :: thklim_temp =   1.d0 ! min thickness for computing vertical temperature (m) (higher-order only)
    real(dp) :: dew    =    20.d3
    real(dp) :: dns    =    20.d3
    real(dp) :: dt     =     0.d0     ! ice dynamics timestep
    real(dp) :: dttem  =     0.d0     ! temperature timestep
    real(dp) :: dt_transport = 0.d0   ! timestep for subcycling transport within the dynamics timestep dt
    real(dp) :: nshlf  =     0.d0          !TODO - not currently used; remove?
    integer  :: subcyc =     1
    real(dp) :: periodic_offset_ew = 0.d0 ! optional periodic_offsets for ismip-hom and similar tests
    real(dp) :: periodic_offset_ns = 0.d0 ! These may be needed to ensure continuous ice geometry at
                                          !  the edges of the global domain.

    integer  :: timecounter = 0   !> count time steps
    
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),pointer :: sigma => null() !> Sigma values for vertical spacing of 
                                                     !> model levels
    real(dp),dimension(:),pointer :: stagsigma => null() !> Staggered values of sigma (layer midpts)
    real(dp),dimension(:),pointer :: stagwbndsigma => null() !> Staggered values of sigma (layer midpts) with boundaries

    integer :: profile_period = 100            ! profile frequency

    real(dp) :: dt_diag = 0.d0     ! diagnostic interval (write diagnostic output every dt_diag years)
                                   ! dt_diag = 0 => never write diagnostic output
    integer  :: ndiag = 0          ! diagnostic interval (write diagnostic output every ndiag timesteps)
                                   ! ndiag = 0 => never write diagnostic output
    integer  :: idiag = 1          ! global grid indices for diagnostic point
    integer  :: jdiag = 1          ! 
    integer  :: idiag_local = 1    ! local grid indices for diagnostic point
    integer  :: jdiag_local = 1
    integer  :: rdiag_local = 0    ! task number for diagnostic point

    real(dp) :: adv_cfl_dt = 0.0d0  ! maximum allowable dt (yrs) based on advective CFL (calculated by model for each time step)
    real(dp) :: diff_cfl_dt = 0.0d0 ! maximum allowable dt (yrs) based on diffusive CFL (calculated by model for each time step)
  end type glide_numerics

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !TODO - Is the glide_grnd type still needed?
  type glide_grnd
    ! variables for tracking the grounding line    
    real(dp),dimension(:,:),pointer :: gl_ew => null()
    real(dp),dimension(:,:),pointer :: gl_ns => null()
    real(dp),dimension(:,:),pointer :: gline_flux => null() !> flux at the
                                                            !grounding line
  end type glide_grnd

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velowk
    real(dp),dimension(:),  pointer :: depth    => null()
    real(dp),dimension(:),  pointer :: dupsw    => null()
    real(dp),dimension(:),  pointer :: depthw   => null()
    real(dp),dimension(:),  pointer :: suvel    => null()
    real(dp),dimension(:),  pointer :: svvel    => null()
    real(dp),dimension(:,:),pointer :: fslip    => null()
    real(dp),dimension(:,:),pointer :: dintflwa => null()
    real(dp),dimension(:),  pointer :: dups     => null()
    real(dp),dimension(4) :: fact
    real(dp),dimension(4) :: c = 0.d0
    real(dp) :: watwd  = 3.0d0
    real(dp) :: watct  = 10.0d0
    real(dp) :: trc0   = 0.d0
    real(dp) :: trcmin = 0.0d0
    real(dp) :: marine = 1.0d0
    real(dp) :: trcmax = 10.0d0
    real(dp) :: btrac_const = 0.0d0  !TODO - Remove from glide_velowk type; already in glide_paramets type.
    real(dp) :: btrac_slope = 0.0d0
    real(dp) :: btrac_max = 0.d0
  end type glide_velowk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_thckwk
     real(dp),dimension(:,:),  pointer :: oldthck   => null()
     real(dp),dimension(:,:),  pointer :: oldthck2  => null()
     real(dp),dimension(:,:,:),pointer :: olds      => null()
     integer  :: nwhich  = 2
     real(dp) :: oldtime = 0.d0
     
     ! next four are for ADI evolution only
     real(dp), dimension(:), pointer :: alpha => null()
     real(dp), dimension(:), pointer :: beta  => null()
     real(dp), dimension(:), pointer :: gamma => null()
     real(dp), dimension(:), pointer :: delta => null()

  end type glide_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !WHL - Moved dissip to glide_temper
  type glide_tempwk
    real(dp),dimension(:,:,:),pointer :: inittemp => null()
    real(dp),dimension(:,:,:),pointer :: compheat => null()
    real(dp),dimension(:,:,:),pointer :: initadvt => null()
    real(dp),dimension(:),    pointer :: dupa     => null()
    real(dp),dimension(:),    pointer :: dupb     => null()
    real(dp),dimension(:),    pointer :: dupc     => null()
    real(dp),dimension(:),    pointer :: c1       => null()
    real(dp),dimension(:,:),  pointer :: dups     => null()
    real(dp),dimension(:,:),  pointer :: wphi     => null()
    real(dp),dimension(:,:),  pointer :: smth     => null()
    real(dp),dimension(:,:,:),pointer :: hadv_u   => null()
    real(dp),dimension(:,:,:),pointer :: hadv_v   => null()

    !*sfp** added space to the next 2 (cons, f) for use w/ HO and SSA dissip. calc. 
    real(dp),dimension(5)             :: cons     = 0.d0
    real(dp),dimension(5)             :: f        = 0.d0
    real(dp),dimension(8)             :: c        = 0.d0
    real(dp),dimension(2)             :: slide_f
    real(dp) :: noflow      = -1
    real(dp),dimension(2) :: advconst = 0.d0
    real(dp) :: zbed        = 0.d0
    real(dp) :: dupn        = 0.d0
    real(dp) :: wmax        = 0.d0
    real(dp) :: dt_wat      = 0.d0
    real(dp) :: watvel      = 0.d0
    integer  :: nwat        = 0
  end type glide_tempwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_paramets
    real(dp),dimension(5) :: bpar = (/ 0.2d0, 0.5d0, 0.0d0 ,1.0d-2, 1.0d0/)
    real(dp) :: btrac_const = 0.d0     ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_slope = 0.0d0    ! Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_max = 0.d0       ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: geot   = -5.0d-2       ! W m^{-2}, positive down
    real(dp) :: flow_enhancement_factor = 1.0d0   ! flow enhancement parameter for the Arrhenius relationship;
                                                  ! typically used in SIA model to speed up the ice
                                       ! (NOTE change relative to prev. versions of code - used to be 3)
    real(dp) :: slip_ratio = 1.0d0     ! Slip ratio, used only in higher order code when the slip ratio beta computation is requested
    real(dp) :: hydtim = 1000.0d0      ! years, converted to s^{-1} and scaled
                                       ! 0 if no drainage
    real(dp) :: bwat_smooth = 0.01d0   ! basal water field smoothing strength
    real(dp) :: default_flwa = 1.0d-16 ! Glen's A to use in isothermal case, in units Pa^{-n} yr^{-1} 
                                       ! (would change to e.g. 4.6e-18 in EISMINT-ROSS case)
    real(dp) :: efvs_constant = 2336041.d0  ! value of efvs to use in constant efvs case, in units Pa yr
                                       ! = 0.5*A^(-1), where A = 2.140373 Pa^(-1) yr^(1) is the value used in ISMIP-HOM Test F
    real(dp) :: p_ocean_penetration = 0.0d0  ! p-exponent parameter for ocean penetration parameterization
    real(dp) :: max_slope = 1.0d0      ! maximum surface slope allowed in Glissade dycore (unitless)
                                       ! Note: It may be necessary to reduce max_slope to ~0.1 to prevent huge velocities
                                       !       in regions of rough coastal topography
            
  end type glide_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Should the glide_basalproc type be removed?
  !       Keeping it for now because glam_strs2 uses mintauf (but this could be moved to another type).
  type glide_basalproc
    !Tuneables, set in the config file 
!    real(dp):: fric=0.45d0                   ! Till coeff of internal friction: ND
!    real(dp):: etillo=0.7d0                  ! Till void ratio at No
!    real(dp):: No=1000.d0                    ! Reference value of till effective stress
!    real(dp):: Comp=0.12d0                   ! Till coeff of compressibility: ND
!    real(dp):: Cv = 1.0d-8                   ! Till hydraulic diffusivity: m2/s
!    real(dp):: Kh = 1.0d-10                  !Till hydraulic conductivity: m/s
!    real(dp):: Zs = 3.0d0                    ! Solid till thickness: m
!    real(dp):: aconst=994000000d0            ! Constant in till strength eq. (Pa)
!    real(dp):: bconst=21.7d0                 ! Constant in till strength eq. (ND)
!    integer:: till_hot = 0
!    integer:: tnodes = 5

    real(dp), dimension (:) , pointer :: till_dz => null()  !holds inital till layer spacing - 
    
    !Model variables that will be passed to other subroutines
    real(dp),dimension(:,:)  ,pointer :: mintauf => null() !Bed strength calculated with basal proc. mod.
!    real(dp),dimension(:,:)  ,pointer :: Hwater  => null() !Water available from till layer (m)
    !Model variables necessary for restart
!    real(dp),dimension(:,:,:)  ,pointer :: u => null()     !Till excess pore pressure (Pa)
!    real(dp),dimension(:,:,:)  ,pointer :: etill  => null()  !Till void ratio (ND)  
    
  end type glide_basalproc

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_prof_type
     integer :: geomderv
     integer :: hvelos
     integer :: ice_mask1
     integer :: temperature
     integer :: ice_evo
     integer :: ice_mask2
     integer :: isos_water
     integer :: isos
  end type glide_prof_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!TODO - Remove the glide_phaml type?  Commented out for now
!!  type glide_phaml
!!    real(dp),dimension(:,:),pointer :: uphaml => null()
!!    real(dp),dimension(:,:),pointer :: init_phaml => null()
!!    real(dp),dimension(:,:),pointer :: rs_phaml => null()
!!    !maybe put the x/y vectors here too just for simplicity
!!  end type glide_phaml
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! for JFNK, NOX Trilinos solver
  type, public :: glissade_solver

    integer ,dimension(:,:) ,allocatable :: ui 
    integer ,dimension(:,:) ,allocatable :: um 
    real(dp) ,dimension(:,:) ,allocatable :: d2thckcross
    real(dp) ,dimension(:,:) ,allocatable :: d2usrfcross
    integer ,dimension(2) :: pcgsize
    integer ,dimension(:) ,allocatable :: gxf
    real(dp)  :: L2norm
    type(sparse_matrix_type)  :: matrix
    type(sparse_matrix_type)  :: matrixA
    type(sparse_matrix_type)  :: matrixC
    real(dp),dimension(:),pointer :: rhsd    => null()
    real(dp),dimension(:),pointer :: answ    => null()
    integer :: ct     = 0

    !TODO - KJE - Remove ewn/nsn from glissade_solver type once new glide_global_type is working and we can use those ewn/nsn
    integer :: ewn
    integer :: nsn

  end type glissade_solver

       
  type glide_global_type    ! type containing all of the above for an ice sheet model instance
    integer              :: model_id !> Used in the global model list for error handling purposes
    type(glide_general)  :: general
    type(glide_options)  :: options
    type(glide_geometry) :: geometry
    type(glide_geomderv) :: geomderv
    type(glide_velocity) :: velocity
    type(glide_stress_t) :: stress   
    type(glide_climate)  :: climate
    type(eismint_climate_type) :: eismint_climate
    type(glide_calving)  :: calving
    type(glide_temper)   :: temper
    type(glide_basal_physics)  :: basal_physics
    type(glide_lithot_type) :: lithot
    type(glide_funits)   :: funits
    type(glide_numerics) :: numerics
    type(glide_velowk)   :: velowk
    type(glide_thckwk)   :: thckwk
    type(glide_tempwk)   :: tempwk
    type(glide_paramets) :: paramets
    type(glimmap_proj)   :: projection
    type(glide_basalproc):: basalproc
    type(profile_type)   :: profile
    type(glide_prof_type):: glide_prof
    type(isostasy_type)  :: isostasy
!!    type(glide_phaml)    :: phaml
    type(glide_grnd)     :: ground
    type(glissade_solver):: solver_data

  end type glide_global_type

contains

  subroutine glide_allocarr(model)    

    !> Allocates the model arrays, and initialises some of them to zero.
    !> These are the arrays allocated, and their dimensions:

    !TODO - Make sure the itemized lists in subroutine glide_allocarr are complete.

    !> In \texttt{model\%temper}:
    !> \begin{itemize}
    !> \item \texttt{temp(upn,0:ewn+1,0:nsn+1))}   !WHL - 2 choices
    !> \item \texttt{bheatflx(ewn,nsn))}
    !> \item \texttt{flwa(upn,ewn,nsn))}           !WHL - 2 choices
    !> \item \texttt{dissip(upn,ewn,nsn))}         !WHL - 2 choices
    !> \item \texttt{bwat(ewn,nsn))}
    !> \item \texttt{bmlt_ground(ewn,nsn))}
    !> \item \texttt{bmlt_float(ewn,nsn))}
    !> \item \texttt{bmlt_float_mask(ewn,nsn))}
    !> \item \texttt{bfricflx(ewn,nsn))}
    !> \item \texttt{ucondflx(ewn,nsn))}
    !> \item \texttt{lcondflx(ewn,nsn))}
    !> \item \texttt{dissipcol(ewn,nsn))}
    !> \item \texttt{waterfrac(upn-1,ewn,nsn))}
    !> \item \texttt{enthalpy(0:upn,ewn,nsn))}
    !> \end{itemize}

    !> In \texttt{model\%velocity}:
    !> \begin{itemize}
    !> \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !> \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !> \item \texttt{wvel(upn,ewn,nsn))}
    !> \item \texttt{wgrd(upn,ewn,nsn))}
    !> \item \texttt{uflx(ewn-1,nsn-1))}
    !> \item \texttt{vflx(ewn-1,nsn-1))}
    !> \item \texttt{diffu(ewn,nsn))}
    !> \item \texttt{btrc(ewn,nsn))}
    !> \item \texttt{ubas(ewn,nsn))}
    !> \item \texttt{vbas(ewn,nsn))}
    !> \end{itemize}

    !> In \texttt{model\%climate}:
    !> \begin{itemize}
    !> \item \texttt{acab(ewn,nsn))}
    !> \item \texttt{artm(ewn,nsn))}
    !> \end{itemize}

    !> In \texttt{model\%geomderv}:
    !> \begin{itemize}
    !> \item \texttt{dthckdew(ewn,nsn))}
    !> \item \texttt{dusrfdew(ewn,nsn))}
    !> \item \texttt{dthckdns(ewn,nsn))}
    !> \item \texttt{dusrfdns(ewn,nsn))}
    !> \item \texttt{dthckdtm(ewn,nsn))}
    !> \item \texttt{dusrfdtm(ewn,nsn))}
    !> \item \texttt{stagthck(ewn-1,nsn-1))}
    !> \end{itemize}
  
    !> In \texttt{model\%geometry}:
    !> \begin{itemize}
    !> \item \texttt{thck(ewn,nsn))}
    !> \item \texttt{usrf(ewn,nsn))}
    !> \item \texttt{lsrf(ewn,nsn))}
    !> \item \texttt{topg(ewn,nsn))}
    !> \item \texttt{mask(ewn,nsn))}
    !> \item \texttt{age(upn-1,ewn,nsn))}
    !> \item \texttt{tracers(ewn,nsn,ntracers,upn-1)}
    !> \item \texttt{f_flotation(ewn,nsn)}
    !> \item \texttt{f_ground(ewn-1,nsn-1)}
    !* (DFM) added floating_mask, ice_mask, lower_cell_loc, and lower_cell_temp
    !> \item \texttt{floating_mask(ewn,nsn))}
    !> \item \texttt{ice_mask(ewn,nsn))}
    !> \item \texttt{lower_cell_loc(ewn,nsn))}
    !> \item \texttt{lower_cell_temp(ewn,nsn))}
    !> \end{itemize}

    !> In \texttt{model\%thckwk}:
    !> \begin{itemize}
    !> \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !> \end{itemize}

    !> In \texttt{model\%numerics}:
    !> \begin{itemize}
    !> \item \texttt{sigma(upn))}
    !> \end{itemize}

    !> In \texttt{model\%numerics}:
    !> \begin{itemize}
    !> \item \texttt{stagsigma(upn-1))}
    !> \end{itemize}

    use glimmer_log
    use glimmer_coordinates, only: coordsystem_allocate
    use glimmer_paramets, only: unphys_val

    implicit none

    type(glide_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn

    ! for simplicity, copy these values...

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
    
    ! horizontal coordinates

    allocate(model%general%x0(ewn-1))!; model%general%x0 = 0.d0  ! velocity grid
    allocate(model%general%y0(nsn-1))!; model%general%y0 = 0.d0
    allocate(model%general%x1(ewn))!; model%general%x1 = 0.d0    ! ice grid (for scalars)
    allocate(model%general%y1(nsn))!; model%general%y1 = 0.d0

    ! vertical sigma coordinates
    ! If we already have sigma, don't reallocate

    if (associated(model%numerics%sigma)) then
       if (size(model%numerics%sigma) /= upn) then
          call write_log('Wrong number of sigma levels given',GM_FATAL)
       end if
    else
       allocate(model%numerics%sigma(upn))
    endif

    allocate(model%numerics%stagsigma(upn-1))
    allocate(model%numerics%stagwbndsigma(0:upn))  !MJH added (0:upn) as separate variable

    ! temperature arrays

    !NOTE: In the glide dycore (whichdycore = DYCORE_GLIDE), the temperature and 
    !       flow factor live on the unstaggered vertical grid, and extra rows and columns 
    !       (with indices 0:ewn+1, 0:nsn+1) are needed.
    !      In the glam/glissade dycore, the temperature and flow factor live on
    !       the staggered vertical grid, with temp and flwa defined at the
    !       center of each layer k = 1:upn-1.  The temperature (but not flwa)
    !       is defined at the upper surface (k = 0) and lower surface (k = upn).

    if (model%options%whichdycore == DYCORE_GLIDE) then
       allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1))
       call coordsystem_allocate(model%general%ice_grid, upn, model%temper%flwa)
       call coordsystem_allocate(model%general%ice_grid, upn, model%temper%dissip)
    else    ! glam/glissade dycore
       allocate(model%temper%temp(0:upn,1:ewn,1:nsn))
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%temper%flwa)
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%temper%dissip)
    endif

    ! MJH - Set temp and flwa to physically unrealistic values so we can tell later if 
    !       arrays were initialized correctly
    model%temper%temp(:,:,:) = unphys_val  ! unphys_val = -999.d0
    model%temper%flwa(:,:,:) = unphys_val
    model%temper%dissip(:,:,:) = 0.d0

    call coordsystem_allocate(model%general%ice_grid,  model%temper%bheatflx)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bwat)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bwatflx)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbwat)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bmlt_ground)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bmlt_float)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bmlt_float_mask)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bpmp)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbpmp)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbtemp)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%ucondflx)

    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glam/glissade only
       call coordsystem_allocate(model%general%ice_grid, model%temper%bfricflx)
       call coordsystem_allocate(model%general%ice_grid, model%temper%lcondflx)
       call coordsystem_allocate(model%general%ice_grid, model%temper%dissipcol)
       ! water fraction and enthalpy live at the midpoint of each layer (with temp and flwa)
       ! enthalpy (like temp) is defined at the upper and lower surfaces as well
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%temper%waterfrac)
       allocate(model%temper%enthalpy(0:upn,1:ewn,1:nsn))
       model%temper%enthalpy(:,:,:) = 0.d0
    endif

    ! velocity arrays

    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%uvel)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%vvel)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%velnorm)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%uflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%bed_softness)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%btrc)
    call coordsystem_allocate(model%general%velo_grid, 2, model%velocity%btraction)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%resid_u)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%resid_v)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%rhs_u)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%rhs_v)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%uvel_2d)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vvel_2d)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas_tavg)

    ! The following are on the extended staggered grid, which is the same size as the ice grid.
    call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%uvel_extend)
    call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%vvel_extend)
    call coordsystem_allocate(model%general%ice_grid,  model%velocity%uvel_2d_extend)
    call coordsystem_allocate(model%general%ice_grid,  model%velocity%vvel_2d_extend)

    if (model%options%whichdycore == DYCORE_GLIDE) then
       call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%wvel)
       call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%wgrd)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu_x)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu_y)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%total_diffu)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%tau_x)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%tau_y)
    else   ! glam/glissade dycore
       call coordsystem_allocate(model%general%velo_grid, model%velocity%beta)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%beta_internal)
       call coordsystem_allocate(model%general%ice_grid, model%velocity%unstagbeta)
       ! Set beta and unstagbeta to physically unrealistic values so we can tell later 
       ! if these fields were read correctly from an input file
       model%velocity%beta(:,:) = unphys_val   ! unphys_val = -999.0d0
       model%velocity%unstagbeta(:,:) = unphys_val

       call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%wvel_ho)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%kinbcmask)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%dynbcmask)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%umask_no_penetration)
       call coordsystem_allocate(model%general%velo_grid, model%velocity%vmask_no_penetration)
         ! next 3 used for output of residual fields (when relevant code in glam_strs2 is active)
!       call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%ures)
!       call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%vres)
!       call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%magres)
    endif

    ! higher-order stress arrays

    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glam/glissade dycore
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%efvs)
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%scalar) 
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%xz)
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%yz)
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%xx)
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%yy)
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%xy)
       call coordsystem_allocate(model%general%velo_grid, model%stress%btractx)
       call coordsystem_allocate(model%general%velo_grid, model%stress%btracty)
       call coordsystem_allocate(model%general%ice_grid, model%stress%btractx_extend)
       call coordsystem_allocate(model%general%ice_grid, model%stress%btracty_extend)
       call coordsystem_allocate(model%general%velo_grid, model%stress%taudx)
       call coordsystem_allocate(model%general%velo_grid, model%stress%taudy)
    endif

    ! geometry arrays
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thck)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%usrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%lsrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%topg)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thkmask)
    call coordsystem_allocate(model%general%velo_grid, model%geometry%stagmask)

    call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagthck)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdns)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdns)

    !* (DFM) -- added floating_mask, ice_mask, lower_cell_loc, and lower_cell_temp here
    call coordsystem_allocate(model%general%ice_grid, model%geometry%floating_mask)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%ice_mask)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%lower_cell_loc)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%lower_cell_temp)

    if (model%options%whichdycore == DYCORE_GLIDE) then
       call coordsystem_allocate(model%general%ice_grid, model%geometry%thck_index)
       call coordsystem_allocate(model%general%ice_grid,  model%geomderv%dthckdtm)
       call coordsystem_allocate(model%general%ice_grid,  model%geomderv%dusrfdtm)
       allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
       model%thckwk%olds = 0.0d0
       call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck)
       call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck2)
    else   ! glam/glissade dycore
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%geometry%ice_age)
       call coordsystem_allocate(model%general%ice_grid,  model%geometry%f_flotation)
       call coordsystem_allocate(model%general%velo_grid, model%geometry%f_ground)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%dlsrfdew)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%dlsrfdns)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%staglsrf)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagusrf)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagtopg)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2usrfdew2)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2usrfdns2)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2thckdew2)
       call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2thckdns2)
       !Note: model%geometry%tracers and related arrays are allocated later, in glissade_transport_setup

       ! Basal Physics
       !WHL - Since the number of basal BC options is proliferating, simplify the logic by allocating the following arrays
       !      whenever running glam/glissade
!!       if ( (model%options%which_ho_babc == HO_BABC_POWERLAW) .or. &
!!            (model%options%which_ho_babc == HO_BABC_COULOMB_FRICTION) .or. &
!!            (model%options%which_ho_babc == HO_BABC_COULOMB_CONST_BASAL_FLWA) .or. &
!!            (model%options%whichbwat == BWATER_OCEAN_PENETRATION)     ) then
       call coordsystem_allocate(model%general%ice_grid, model%basal_physics%effecpress)
       call coordsystem_allocate(model%general%velo_grid, model%basal_physics%effecpress_stag)
       call coordsystem_allocate(model%general%ice_grid, model%basal_physics%C_space_factor)
       call coordsystem_allocate(model%general%velo_grid, model%basal_physics%C_space_factor_stag)
!!       endif

    endif  ! glam/glissade

    ! climate arrays
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab)
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab_tavg)
    call coordsystem_allocate(model%general%ice_grid, model%climate%artm)

    ! calving arrays
    call coordsystem_allocate(model%general%ice_grid, model%calving%calving_thck)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%calving%damage)
    call coordsystem_allocate(model%general%ice_grid, model%calving%damage_column)

    ! matrix solver arrays

    allocate (model%solver_data%rhsd(ewn*nsn))
    allocate (model%solver_data%answ(ewn*nsn))

    call new_sparse_matrix(ewn*nsn, 5*ewn*nsn, model%solver_data%matrix)

    !TODO - In model%lithot%temp, put the vertical index 3rd as in model%temper%temp?

    ! lithosphere arrays

    if (model%options%gthf == GTHF_COMPUTE) then
       allocate(model%lithot%temp(1:ewn,1:nsn,model%lithot%nlayer)); model%lithot%temp = 0.d0
       call coordsystem_allocate(model%general%ice_grid, model%lithot%mask)
    endif

    ! isostasy arrays

    call coordsystem_allocate(model%general%ice_grid, model%isostasy%relx)  ! MJH: relx needs to be allocated always.
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call coordsystem_allocate(model%general%ice_grid, model%isostasy%load)
       call coordsystem_allocate(model%general%ice_grid, model%isostasy%load_factors)
    endif

    ! The remaining arrays are not currently used (except mintauf)
    ! phaml arrays
!!    call coordsystem_allocate(model%general%ice_grid, model%phaml%init_phaml)
!!    call coordsystem_allocate(model%general%ice_grid, model%phaml%rs_phaml)
!!    call coordsystem_allocate(model%general%ice_grid, model%phaml%uphaml)

    ! grounding line arrays (not currently supported)

!!    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glam/glissade dycore
!!       allocate (model%ground%gl_ew(ewn-1,nsn))
!!       allocate (model%ground%gl_ns(ewn,nsn-1))
!!       allocate (model%ground%gline_flux(ewn,nsn)) 
!!    endif

    ! basal process arrays
    ! not currently supported, except that glam_strs2 uses mintauf

    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glam/glissade dycore
!!       call coordsystem_allocate(model%general%ice_grid, model%basalproc%Hwater)
       call coordsystem_allocate(model%general%velo_grid, model%basalproc%mintauf)
!!       allocate(model%basalproc%u (ewn-1,nsn-1,model%basalproc%tnodes)); model%basalproc%u=41.0d3
!!       allocate(model%basalproc%etill (ewn-1,nsn-1,model%basalproc%tnodes));model%basalproc%etill=0.5d0
    endif

  end subroutine glide_allocarr


  subroutine glide_deallocarr(model)

    !> deallocate model arrays
    !TODO - Check that all arrays allocated above are deallocated here.

    implicit none
    type(glide_global_type),intent(inout) :: model

    ! horizontal coordinates

    if (associated(model%general%x0)) &
        deallocate(model%general%x0) 
    if (associated(model%general%y0)) &
        deallocate(model%general%y0) 
    if (associated(model%general%x1)) &
        deallocate(model%general%x1) 
    if (associated(model%general%y1)) &
        deallocate(model%general%y1) 

    ! vertical sigma coordinates

    if (associated(model%numerics%sigma)) &
        deallocate(model%numerics%sigma)
    if (associated(model%numerics%stagsigma)) &
        deallocate(model%numerics%stagsigma)
    if (associated(model%numerics%stagwbndsigma)) &
        deallocate(model%numerics%stagwbndsigma)

    ! temperature arrays

    if (associated(model%temper%temp)) &
        deallocate(model%temper%temp)
    if (associated(model%temper%bheatflx)) &
        deallocate(model%temper%bheatflx)
    if (associated(model%temper%bwat)) &
        deallocate(model%temper%bwat)
    if (associated(model%temper%bwatflx)) &
        deallocate(model%temper%bwatflx)
    if (associated(model%temper%stagbwat)) &
        deallocate(model%temper%stagbwat)
    if (associated(model%temper%bmlt_ground)) &
        deallocate(model%temper%bmlt_ground)
    if (associated(model%temper%bmlt_float)) &
        deallocate(model%temper%bmlt_float)
    if (associated(model%temper%bmlt_float_mask)) &
        deallocate(model%temper%bmlt_float_mask)
    if (associated(model%temper%bpmp)) &
        deallocate(model%temper%bpmp)
    if (associated(model%temper%stagbpmp)) &
        deallocate(model%temper%stagbpmp)
    if (associated(model%temper%stagbtemp)) &
        deallocate(model%temper%stagbtemp)
    if (associated(model%temper%bfricflx)) &
        deallocate(model%temper%bfricflx)
    if (associated(model%temper%ucondflx)) &
        deallocate(model%temper%ucondflx)
    if (associated(model%temper%lcondflx)) &
        deallocate(model%temper%lcondflx)
    if (associated(model%temper%dissipcol)) &
        deallocate(model%temper%dissipcol)
    if (associated(model%temper%waterfrac)) &
        deallocate(model%temper%waterfrac)
    if (associated(model%temper%enthalpy)) &
        deallocate(model%temper%enthalpy)
    if (associated(model%temper%flwa)) &
        deallocate(model%temper%flwa)
    if (associated(model%temper%dissip)) &
        deallocate(model%temper%dissip)

    ! velocity arrays

    if (associated(model%velocity%uvel)) &
        deallocate(model%velocity%uvel)
    if (associated(model%velocity%vvel)) &
        deallocate(model%velocity%vvel)
    if (associated(model%velocity%uvel_2d)) &
        deallocate(model%velocity%uvel_2d)
    if (associated(model%velocity%vvel_2d)) &
        deallocate(model%velocity%vvel_2d)
    if (associated(model%velocity%velnorm)) &
        deallocate(model%velocity%velnorm)
    if (associated(model%velocity%wvel)) &
        deallocate(model%velocity%wvel)
    if (associated(model%velocity%uflx)) &
        deallocate(model%velocity%uflx)
    if (associated(model%velocity%vflx)) &
        deallocate(model%velocity%vflx)
    if (associated(model%velocity%bed_softness)) &
        deallocate(model%velocity%bed_softness)
    if (associated(model%velocity%btrc)) &
        deallocate(model%velocity%btrc)
    if (associated(model%velocity%btraction)) &
        deallocate(model%velocity%btraction)
    if (associated(model%velocity%resid_u)) &
        deallocate(model%velocity%resid_u)
    if (associated(model%velocity%resid_v)) &
        deallocate(model%velocity%resid_v)
    if (associated(model%velocity%rhs_u)) &
        deallocate(model%velocity%rhs_u)
    if (associated(model%velocity%rhs_v)) &
        deallocate(model%velocity%rhs_v)
    if (associated(model%velocity%uvel_extend)) &
        deallocate(model%velocity%uvel_extend)
    if (associated(model%velocity%vvel_extend)) &
        deallocate(model%velocity%vvel_extend)
    if (associated(model%velocity%uvel_2d_extend)) &
        deallocate(model%velocity%uvel_2d_extend)
    if (associated(model%velocity%vvel_2d_extend)) &
        deallocate(model%velocity%vvel_2d_extend)
    if (associated(model%velocity%ubas)) &
        deallocate(model%velocity%ubas)
    if (associated(model%velocity%ubas_tavg)) &
        deallocate(model%velocity%ubas_tavg)
    if (associated(model%velocity%vbas)) &
        deallocate(model%velocity%vbas)
    if (associated(model%velocity%vbas_tavg)) &
        deallocate(model%velocity%vbas_tavg)

    if (associated(model%velocity%wgrd)) &
        deallocate(model%velocity%wgrd)
    if (associated(model%velocity%diffu)) &
        deallocate(model%velocity%diffu)
    if (associated(model%velocity%diffu_x)) &
        deallocate(model%velocity%diffu_x)
    if (associated(model%velocity%diffu_y)) &
        deallocate(model%velocity%diffu_y)
    if (associated(model%velocity%total_diffu)) &
        deallocate(model%velocity%total_diffu)
    if (associated(model%velocity%tau_x)) &
        deallocate(model%velocity%tau_x)
    if (associated(model%velocity%tau_y)) &
        deallocate(model%velocity%tau_y)

!!    if (associated(model%velocity%velmask)) &   ! no longer used
!!       deallocate(model%velocity%velmask) 

    if (associated(model%velocity%beta)) &
        deallocate(model%velocity%beta)
    if (associated(model%velocity%unstagbeta)) &
        deallocate(model%velocity%unstagbeta)
    if (associated(model%velocity%beta_internal)) &
        deallocate(model%velocity%beta_internal)
    if (associated(model%velocity%wvel_ho)) &
        deallocate(model%velocity%wvel_ho)
    if (associated(model%velocity%kinbcmask)) &
        deallocate(model%velocity%kinbcmask)
    if (associated(model%velocity%dynbcmask)) &
        deallocate(model%velocity%dynbcmask)
    if (associated(model%velocity%umask_no_penetration)) &
        deallocate(model%velocity%umask_no_penetration)
    if (associated(model%velocity%vmask_no_penetration)) &
        deallocate(model%velocity%vmask_no_penetration)

    !! next 3 used for output of residual fields (when relevant code in glam_strs2 is active)
!    if (associated(model%velocity%ures)) & 
!        deallocate(model%velocity%ures) 
!    if (associated(model%velocity%vres)) & 
!        deallocate(model%velocity%vres) 
!    if (associated(model%velocity%magres)) & 
!        deallocate(model%velocity%magres) 

    ! higher-order stress arrays

    if (associated(model%stress%efvs)) &
        deallocate(model%stress%efvs)
    if (associated(model%stress%tau%scalar)) &
        deallocate(model%stress%tau%scalar)
    if (associated(model%stress%tau%xz)) &
        deallocate(model%stress%tau%xz)
    if (associated(model%stress%tau%yz)) &
        deallocate(model%stress%tau%yz)
    if (associated(model%stress%tau%xx)) &
        deallocate(model%stress%tau%xx)
    if (associated(model%stress%tau%yy)) &
        deallocate(model%stress%tau%yy)
    if (associated(model%stress%tau%xy)) &
        deallocate(model%stress%tau%xy)
    if (associated(model%stress%btractx)) &
        deallocate(model%stress%btractx)
    if (associated(model%stress%btracty)) &
        deallocate(model%stress%btracty)
    if (associated(model%stress%btractx_extend)) &
        deallocate(model%stress%btractx_extend)
    if (associated(model%stress%btracty_extend)) &
        deallocate(model%stress%btracty_extend)
    if (associated(model%stress%taudx)) &
        deallocate(model%stress%taudx)
    if (associated(model%stress%taudy)) &
        deallocate(model%stress%taudy)

    ! basal physics arrays
    if (associated(model%basal_physics%effecpress)) &
        deallocate(model%basal_physics%effecpress)
    if (associated(model%basal_physics%effecpress_stag)) &
        deallocate(model%basal_physics%effecpress_stag)
    if (associated(model%basal_physics%C_space_factor)) &
        deallocate(model%basal_physics%C_space_factor)
    if (associated(model%basal_physics%C_space_factor_stag)) &
        deallocate(model%basal_physics%C_space_factor_stag)

    ! geometry arrays

    if (associated(model%geometry%thck)) &
        deallocate(model%geometry%thck)
    if (associated(model%geometry%usrf)) &
        deallocate(model%geometry%usrf)
    if (associated(model%geometry%lsrf)) &
        deallocate(model%geometry%lsrf)
    if (associated(model%geometry%topg)) &
        deallocate(model%geometry%topg)
    if (associated(model%geometry%thkmask)) &
        deallocate(model%geometry%thkmask)
    if (associated(model%geometry%stagmask)) &
        deallocate(model%geometry%stagmask)
    if (associated(model%geomderv%stagthck)) &
        deallocate(model%geomderv%stagthck)
    if (associated(model%geomderv%dthckdew)) &
        deallocate(model%geomderv%dthckdew)
    if (associated(model%geomderv%dthckdns)) &
        deallocate(model%geomderv%dthckdns)
    if (associated(model%geomderv%dusrfdew)) &
        deallocate(model%geomderv%dusrfdew)
    if (associated(model%geomderv%dusrfdns)) &
        deallocate(model%geomderv%dusrfdns)
!!    if (associated(model%geometry%marine_bc_normal)) &
!!       deallocate(model%geometry%marine_bc_normal)

    !*SFP: fields that need to be passed to POP for ice ocean coupling
    !* (DFM -- deallocate floating_mask, ice_mask, lower_cell_loc, and lower_cell_temp)
    if (associated(model%geometry%floating_mask)) &
       deallocate(model%geometry%floating_mask)
    if (associated(model%geometry%ice_mask)) &
       deallocate(model%geometry%ice_mask)
    if (associated(model%geometry%lower_cell_loc)) &
       deallocate(model%geometry%lower_cell_loc)
    if (associated(model%geometry%lower_cell_temp)) &
       deallocate(model%geometry%lower_cell_temp)

    if (associated(model%geometry%thck_index)) &
        deallocate(model%geometry%thck_index)
    if (associated(model%geomderv%dthckdtm)) &
        deallocate(model%geomderv%dthckdtm)
    if (associated(model%geomderv%dusrfdtm)) &
        deallocate(model%geomderv%dusrfdtm)
    if (associated(model%thckwk%olds)) &
        deallocate(model%thckwk%olds)
    if (associated(model%thckwk%oldthck)) &
        deallocate(model%thckwk%oldthck)
    if (associated(model%thckwk%oldthck2)) &
        deallocate(model%thckwk%oldthck2)

    if (associated(model%geometry%ice_age)) &
        deallocate(model%geometry%ice_age)
    if (associated(model%geometry%tracers)) &
        deallocate(model%geometry%tracers)
    if (associated(model%geometry%f_flotation)) &
        deallocate(model%geometry%f_flotation)
    if (associated(model%geometry%f_ground)) &
        deallocate(model%geometry%f_ground)
    if (associated(model%geomderv%dlsrfdew)) &
        deallocate(model%geomderv%dlsrfdew)
    if (associated(model%geomderv%dlsrfdns)) &
        deallocate(model%geomderv%dlsrfdns)
    if (associated(model%geomderv%staglsrf)) &
        deallocate(model%geomderv%staglsrf)
    if (associated(model%geomderv%stagusrf)) &
        deallocate(model%geomderv%stagusrf)
    if (associated(model%geomderv%stagtopg)) &
        deallocate(model%geomderv%stagtopg)
    if (associated(model%geomderv%d2usrfdew2)) &
        deallocate(model%geomderv%d2usrfdew2)
    if (associated(model%geomderv%d2usrfdns2)) &
        deallocate(model%geomderv%d2usrfdns2)
    if (associated(model%geomderv%d2thckdew2)) &
        deallocate(model%geomderv%d2thckdew2)
    if (associated(model%geomderv%d2thckdns2)) &
        deallocate(model%geomderv%d2thckdns2)

    ! climate arrays

    if (associated(model%climate%acab)) &
        deallocate(model%climate%acab)
    if (associated(model%climate%acab_tavg)) &
        deallocate(model%climate%acab_tavg)
    if (associated(model%climate%artm)) &
        deallocate(model%climate%artm)

    ! calving arrays
    if (associated(model%calving%calving_thck)) &
        deallocate(model%calving%calving_thck)
    if (associated(model%calving%damage)) &
        deallocate(model%calving%damage)
    if (associated(model%calving%damage_column)) &
        deallocate(model%calving%damage_column)

    ! matrix solver arrays

    if (associated(model%solver_data%rhsd))  &  
        deallocate(model%solver_data%rhsd)
    if (associated(model%solver_data%answ))  &
        deallocate(model%solver_data%answ)

    !KJE do we need this here? The parts within are allocated in glam_strs2
    call del_sparse_matrix(model%solver_data%matrix)

    ! lithosphere arrays

    if (associated(model%lithot%temp)) &
        deallocate(model%lithot%temp)
    if (associated(model%lithot%mask)) &
        deallocate(model%lithot%mask)

    ! isostasy arrays

    if (associated(model%isostasy%relx)) &
        deallocate(model%isostasy%relx)
    if (associated(model%isostasy%load)) &
        deallocate(model%isostasy%load)
    if (associated(model%isostasy%load_factors)) &
        deallocate(model%isostasy%load_factors)

    ! The remaining arrays are not currently used (except mintauf)
    ! phaml arrays

!!    if (associated(model%phaml%init_phaml)) &
!!       deallocate(model%phaml%init_phaml)
!!    if (associated(model%phaml%rs_phaml)) &
!!       deallocate(model%phaml%rs_phaml)    
!!    if (associated(model%phaml%uphaml)) &
!!       deallocate(model%phaml%uphaml)

    ! grounding line arrays (not currently supported)

!!    if (associated(model%ground%gl_ns)) &
!!        deallocate(model%ground%gl_ns)
!!    if (associated(model%ground%gl_ew)) &
!!        deallocate(model%ground%gl_ew)
!!    if (associated(model%ground%gline_flux)) &
!!        deallocate(model%ground%gline_flux)

    ! basal process arrays
    ! not currently supported, except that glam_strs2 uses mintauf

!!    if (associated(model%basalproc%Hwater)) &
!!       deallocate(model%basalproc%Hwater)
    if (associated(model%basalproc%mintauf)) &
       deallocate(model%basalproc%mintauf)
!!    if (associated(model%basalproc%u)) &
!!       deallocate(model%basalproc%u)
!!    if (associated(model%basalproc%etill)) &
!!       deallocate(model%basalproc%etill)

  end subroutine glide_deallocarr


  ! some accessor functions
  function get_dew(model)
    !> return scaled x node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dew
    type(glide_global_type) :: model

    get_dew = model%numerics%dew * len0
  end function get_dew

  function get_dns(model)
    !> return scaled y node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dns
    type(glide_global_type) :: model

    get_dns = model%numerics%dns * len0
  end function get_dns

  function get_tstart(model)
    !> return start time
    implicit none
    real(dp) :: get_tstart
    type(glide_global_type) :: model
    
    get_tstart = model%numerics%tstart
  end function get_tstart

  function get_tend(model)
    !> return end time
    implicit none
    real(dp) :: get_tend
    type(glide_global_type) :: model
    
    get_tend = model%numerics%tend
  end function get_tend

  function get_tinc(model)
    !> return time increment
    implicit none
    real(dp) :: get_tinc
    type(glide_global_type) :: model
    
    get_tinc = model%numerics%tinc
  end function get_tinc

  function get_ewn(model)
    !> get number of nodes in x dir
    implicit none
    integer get_ewn
    type(glide_global_type) :: model

    get_ewn = model%general%ewn
  end function get_ewn

  function get_nsn(model)
    !> get number of nodes in y dir
    implicit none
    integer get_nsn
    type(glide_global_type) :: model

    get_nsn = model%general%nsn
  end function get_nsn
  
  subroutine set_time(model,time)
    !> Set the model time counter --- useful for
    !> fractional year output
    implicit none
    type(glide_global_type) :: model
    real(dp) :: time

    model%numerics%time = time
  end subroutine set_time

end module glide_types
