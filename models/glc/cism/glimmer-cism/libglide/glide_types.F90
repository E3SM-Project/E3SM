!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_types.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_types

  !*FD Holds type definitions for the derived types used by each 
  !*FD instance of the ice model. Originally, each of these types
  !*FD was a module containing variables, which were used as containers
  !*FD for global variables. However, the need to allow for multiple
  !*FD ice model instances meant that the nested derived types were instituted
  !*FD instead. However, there is probably one too many levels in this scheme. 
  !*FD It would be better if the different types here were contained in the 
  !*FD higher-level instance type (\texttt{glint\_instance}), rather than 
  !*FD the intermediate model type (\texttt{glide\_global\_type}). 
  !*FD 
  !*FD Note that this \emph{is} now where the defaults are defined for these
  !*FD variables.

!TODO - Reorganize the types, as suggested above?
!       The higher-level glint_instance type is defined in glint_type.F90.
!       The various instances are contained in the highest-level type, glint_params.

!TODO - We might consider cleaning up the glide_global type so that it doesn't
!        include as many subtypes.  For example, we might be able to remove
!        replace some work types (tempwk, velowk) with local arrays and parameters.

  use glimmer_sparse_type
  use glimmer_global
  use glimmer_ncdf
  use profile
  use glimmer_coordinates
  use glimmer_map_types

  implicit none

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Constants that describe the options available
  ! We use these integer parameters elsewhere in the code to avoid
  !  hardwiring of option numbers

  ! basic Glimmer/Glide options

  integer, parameter :: DYCORE_GLIDE = 0     ! old shallow-ice dycore from Glimmer
  integer, parameter :: DYCORE_GLAM = 1      ! Payne-Price finite-difference solver
  integer, parameter :: DYCORE_GLISSADE = 2  ! prototype finite-element solver

  !WHL - Removed -1 option (replaced by new glint option: evolve_ice)
  integer, parameter :: EVOL_PSEUDO_DIFF = 0    ! glide only
  integer, parameter :: EVOL_ADI = 1            ! glide only
  integer, parameter :: EVOL_DIFFUSION = 2      ! glide only
  integer, parameter :: EVOL_INC_REMAP = 3      ! glam/glissade only
  integer, parameter :: EVOL_UPWIND = 4         ! glam/glissade only
  integer, parameter :: EVOL_NO_THICKNESS = 5   ! glam/glissade only

  !NOTE: Option 3 is now deprecated.
  !      Use option 1 for prognostic temperature with any dycore
  !TODO: Remove option 3
  integer, parameter :: TEMP_SURFACE_AIR_TEMP = 0
  integer, parameter :: TEMP_PROGNOSTIC = 1
  integer, parameter :: TEMP_STEADY = 2
  integer, parameter :: TEMP_REMAP_ADV = 3

  integer, parameter :: TEMP_INIT_ZERO = 0
  integer, parameter :: TEMP_INIT_ARTM = 1
  integer, parameter :: TEMP_INIT_LINEAR = 2

  !WHL - Swapped 0 and 2
!  integer, parameter :: FLWA_PATERSON_BUDD = 0
!  integer, parameter :: FLWA_PATERSON_BUDD_CONST_TEMP = 1
!  integer, parameter :: FLWA_CONST_FLWA = 2
  integer, parameter :: FLWA_CONST_FLWA = 0
  integer, parameter :: FLWA_PATERSON_BUDD_CONST_TEMP = 1
  integer, parameter :: FLWA_PATERSON_BUDD = 2

  !WHL: Swapped 3 and 5
!!  integer, parameter :: BTRC_ZERO = 0
!!  integer, parameter :: BTRC_CONSTANT = 1
!!  integer, parameter :: BTRC_CONSTANT_BWAT = 2
!!  integer, parameter :: BTRC_TANH_BWAT = 3
!!  integer, parameter :: BTRC_LINEAR_BMLT = 4
!!  integer, parameter :: BTRC_CONSTANT_TPMP = 5
  integer, parameter :: BTRC_ZERO = 0
  integer, parameter :: BTRC_CONSTANT = 1
  integer, parameter :: BTRC_CONSTANT_BWAT = 2
  integer, parameter :: BTRC_CONSTANT_TPMP = 3
  integer, parameter :: BTRC_LINEAR_BMLT = 4
  integer, parameter :: BTRC_TANH_BWAT = 5

  !WHL - Set NONE = 0, LOCAL = 1, FLUX = 2
!!  integer, parameter :: BWATER_LOCAL = 0
!!  integer, parameter :: BWATER_FLUX  = 1
!!  integer, parameter :: BWATER_NONE  = 2
!!  integer, parameter :: BWATER_CONST = 3
!!  !integer, parameter :: BWATER_BASAL_PROC = 4  ! not currently supported
  integer, parameter :: BWATER_NONE  = 0
  integer, parameter :: BWATER_LOCAL = 1
  integer, parameter :: BWATER_FLUX  = 2
  integer, parameter :: BWATER_CONST = 3
  !integer, parameter :: BWATER_BASAL_PROC = 4  ! not currently supported

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

  !WHL - Swapped 2 and 3
!!  integer, parameter :: MARINE_NONE = 0
!!  integer, parameter :: MARINE_FLOAT_ZERO = 1
!!  integer, parameter :: MARINE_RELX_THRESHOLD = 2
!!  integer, parameter :: MARINE_FLOAT_FRACTION = 3
!!  integer, parameter :: MARINE_TOPG_THRESHOLD = 4
!!  integer, parameter :: MARINE_HUYBRECHTS = 5
  integer, parameter :: MARINE_NONE = 0
  integer, parameter :: MARINE_FLOAT_ZERO = 1
  integer, parameter :: MARINE_FLOAT_FRACTION = 2
  integer, parameter :: MARINE_RELX_THRESHOLD = 3
  integer, parameter :: MARINE_TOPG_THRESHOLD = 4
  integer, parameter :: MARINE_HUYBRECHTS = 5

  integer, parameter :: VERTINT_STANDARD = 0
  integer, parameter :: VERTINT_KINEMATIC_BC = 1

  integer, parameter :: SIGMA_COMPUTE_GLIDE = 0
  integer, parameter :: SIGMA_EXTERNAL = 1
  integer, parameter :: SIGMA_CONFIG = 2
  integer, parameter :: SIGMA_COMPUTE_EVEN = 3
  integer, parameter :: SIGMA_COMPUTE_PATTYN = 4

  !TODO - Make this a logical variable?
  integer, parameter :: RESTART_FALSE = 0
  integer, parameter :: RESTART_TRUE = 1

  !basal proc option disabled for now
  integer, parameter :: BAS_PROC_DISABLED = 0
!!  integer, parameter :: BAS_PROC_FULLCALC = 1
!!  integer, parameter :: BAS_PROC_FASTCALC = 2


  ! higher-order options

  !WHL: Swapped 0 and 2
!  integer, parameter :: HO_EFVS_NONLINEAR = 0
!  integer, parameter :: HO_EFVS_FLOWFACT = 1
!  integer, parameter :: HO_EFVS_CONSTANT = 2
  integer, parameter :: HO_EFVS_CONSTANT = 0
  integer, parameter :: HO_EFVS_FLOWFACT = 1
  integer, parameter :: HO_EFVS_NONLINEAR = 2

  integer, parameter :: SIA_DISP = 0
  integer, parameter :: FIRSTORDER_DISP = 1
!!  integer, parameter :: SSA_DISP = 2  ! not supported

  integer, parameter :: HO_BABC_CONSTANT = 0
  integer, parameter :: HO_BABC_SIMPLE = 1
  integer, parameter :: HO_BABC_YIELD_PICARD = 2
  integer, parameter :: HO_BABC_CIRCULAR_SHELF = 3
  integer, parameter :: HO_BABC_LARGE_BETA = 4
  integer, parameter :: HO_BABC_EXTERNAL_BETA = 5
  integer, parameter :: HO_BABC_NO_SLIP = 6
  integer, parameter :: HO_BABC_YIELD_NEWTON = 7

  integer, parameter :: HO_NONLIN_PICARD = 0
  integer, parameter :: HO_NONLIN_JFNK = 1

  integer, parameter :: HO_RESID_MAXU = 0
  integer, parameter :: HO_RESID_MAXU_NO_UBAS = 1
  integer, parameter :: HO_RESID_MEANU = 2
  integer, parameter :: HO_RESID_L2NORM = 3

  integer, parameter :: HO_SPARSE_BICG = 0
  integer, parameter :: HO_SPARSE_GMRES = 1
  integer, parameter :: HO_SPARSE_PCG_INCH = 2
  integer, parameter :: HO_SPARSE_PCG_STRUC = 3
  integer, parameter :: HO_SPARSE_TRILINOS = 4

!WHL - added options for different Stokes approximations
!      (for glissade dycore only)
!      commented out for now
!!  integer, parameter :: HO_APPROX_SIA = 0
!!  integer, parameter :: HO_APPROX_SSA = 1
!!  integer, parameter :: HO_APPROX_BP = 2

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_general

    !*FD Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !*FD The number of grid-points in the E-W direction.
    integer :: nsn = 0  !*FD The number of grid-points in the N-S direction.
    integer :: upn = 1  !*FD The number of vertical levels in the model.

    type(coordsystem_type) :: ice_grid  !*FD coordinate system of the ice grid
    type(coordsystem_type) :: velo_grid !*FD coordinate system of the velocity grid

    real(sp), dimension(:),pointer :: x0 => null() !original x0 grid 
    real(sp), dimension(:),pointer :: y0 => null() !original y0 grid
    real(sp), dimension(:),pointer :: x1 => null() !original x1 grid
    real(sp), dimension(:),pointer :: y1 => null() !original y1 grid

  end type glide_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_options

    !*FD Holds user options controlling the methods used in the ice-model integration.

    !-----------------------------------------------------------------------
    ! standard options
    !-----------------------------------------------------------------------

    integer :: whichdycore = 1

    ! Choice of two Glimmer dycores:
    !*FD \begin{description} 
    !*FD \item[0] Glide dycore (SIA, serial only)
    !*FD \item[1] Glissade dycore (HO, serial or parallel)
    !*FD \end{description}

    integer :: whichevol = 0

    !*FD Thickness evolution method:
    !*FD \begin{description}
    !*FD \item[0] Pseudo-diffusion 
    !*FD \item[1] Alternating direction implicit (ADI)
    !*FD \item[2] Diffusion (also calculates velocities) 
    !*FD \item[3] Incremental remapping
    !*FD \item[4] 1st-order upwind scheme
    !*FD \item[5] Temperature evolves but thickness does not
    !*FD \end{description}

    integer :: whichtemp = 1

    !TODO: Remove option 3 (after cleaning up config files)

    !*FD Method of ice temperature calculation:
    !*FD \begin{description} 
    !*FD \item[0] Set column to surface air temperature
    !*FD \item[1] Full prognostic temperature solution 
    !*FD \item[2] Do NOTHING - hold temperatures steady at initial value  
    !*FD \item[3] Use remapping to advect temperature (deprecated; now combined with [1])
    !*FD \end{description}

    !TODO: At some point, change the default to temp_init = 2.
    !      Setting default = 1 for now so that existing config files will get the same results.

    integer :: temp_init = 1

    ! Temperature initialization:
    !*FD \begin{description} 
    !*FD \item[0] Initialize temperature to 0 C
    !*FD \item[1] Initialize temperature to surface air temperature
    !*FD \item[2] Initialize temperature with a linear profile in each column
    !*FD \end{description}

    !*FD Method for calculating flow factor $A$:

    integer :: whichflwa = 2

    !*FD \begin{description} 
    !*FD \item[0] Set equal to $1\times 10^{-16}\,\mathrm{yr}^{-1}
    !*FD \item[1] \emph{Paterson and Budd} relationship, 
    !*FD with temperature set to $-10^{\circ}\mathrm{C}$ 
    !*FD \item[2] \emph{Paterson and Budd} relationship 
    !*FD \,\mathrm{Pa}^{-n}$
    !*FD \end{description}

    integer :: whichbtrc = 0

    !*FD Basal slip coefficient:
    !*FD \begin{description}
    !*FD \item[0] Set equal to zero everywhere
    !*FD \item[1] Set to (non--zero) constant
    !*FD \item[2] Set to (non--zero) constant where basal water is present, otherwise to zero
    !*FD \item[3] Set to (non--zero) constant where temperature is at pressure melting point of ice, otherwise to zero
    !*FD \item[4] linear function of basal melt rate
    !*FD \item[5] \texttt{tanh} function of basal water depth 
    !*FD \end{description}

    integer :: whichbwat = 0

    !*FD Basal water depth: 
    !*FD \begin{description} 
    !*FD \item[0] Set to zero everywhere 
    !*FD \item[1] Compute from local basal water balance 
    !*FD \item[2] Compute the basal water flux, then find depth via calculation
    !*FD \item[3] Set to constant (10 m) everywhere, to force T = Tpmp.
    !*FD \item[4] Calculated from till water content, in the basal processes module
    !*FD \end{description}

    integer :: basal_mbal = 0

    !*FD basal melt rate:
    !*FD \begin{description}
    !*FD \item[0] Basal melt rate not included in continuity equation
    !*FD \item[1] Basal melt rate included in continuity equation
    !*FD \end{description}

    integer :: gthf = 0

    !*FD geothermal heat flux:
    !*FD \begin{description}
    !*FD \item[0] prescribed uniform geothermal heat flux
    !*FD \item[1] read 2D geothermal flux field from input file (if present)
    !*FD \item[2] calculate geothermal flux using 3d diffusion
    !*FD \end{description}

    !WHL - new isostasy option; replaces model%isos%do_isos
    integer :: isostasy = 0

    !*FD isostasy:
    !*FD \begin{description}
    !*FD \item[0] no isostatic adjustment
    !*FD \item[1] compute isostatic adjustment using lithosphere/asthenosphere model
    !*FD \end{description}

    !TODO - Should this move from the options to the isostasy section?
    !TODO - Should the default be = 1?  Nothing happens for case 0.
    integer :: whichrelaxed = 0

    !*FD relaxed topography:
    !*FD \begin{description}
    !*FD \item[0] get relaxed topo from separate variable (in practice, do nothing)
    !*FD \item[1] first time slice of input topo is relaxed
    !*FD \item[2] first time slice of input topo is in isostatic equilibrium
    !*FD \end{description}

    integer :: whichmarn = 1

    !*FD Marine limit: 
    !*FD \begin{description} 
    !*FD \item[0] No action 
    !*FD \item[1] Set thickness to zero if floating 
    !*FD \item[2] Lose fraction of ice when edge cell
    !*FD \item[3] Set thickness to zero if relaxed bedrock is more than
    !*FD          certain water depth (variable "mlimit" in glide_types)  
    !*FD \item[4] Set thickness to zero if present bedrock is more than
    !*FD          certain water depth (variable "mlimit" in glide_types)  
    !*FD \item[5] Huybrechts grounding line scheme for Greenland initialization
    !*FD \end{description}

    integer :: whichwvel = 0

    !*FD Vertical velocities: 
    !*FD \begin{description}
    !*FD \item[0] Usual vertical integration 
    !*FD \item[1] Vertical integration constrained so that 
    !*FD upper kinematic B.C. obeyed 
    !*FD \end{description}

    integer :: which_sigma = 0

    !*FD \begin{description}
    !*FD \item[0] compute standard Glimmer sigma coordinates
    !*FD \item[1] sigma coordinates are given in external file
    !*FD \item[2] sigma coordinates are given in configuration file
    !*FD \item[3] evenly spaced levels, as required for glam dycore
    !*FD \item[2] compute Pattyn sigma coordinates
    !*FD \end{description}

    !TODO - Make is_restart a logical variable?

    integer :: is_restart = 0
    !*FD if the run is a restart of a previous run
    !*FD \begin{description}
    !*FD \item[0] normal start-up
    !*FD \item[1] restart model from previous run
    !*FD \end{description}

    ! This is a Glimmer serial option
    ! The parallel code enforces periodic EW and NS boundary conditions by default
    logical :: periodic_ew = .false.

    !*FD \begin{description}
    !*FD \item[0] no periodic EW boundary conditions
    !*FD \item[1] periodic EW boundary conditions
    !*FD \end{description}

    !-----------------------------------------------------------------------
    ! Higher-order options
    ! Associated with Payne-Price dycore (glam) and newer glissade dycore
    !-----------------------------------------------------------------------

    integer :: which_ho_efvs = 2

    !*FD Flag that indicates how effective viscosity is computed
    !*FD \begin{description}
    !*FD \item[0] constant value
    !*FD \item[1] multiple of flow factor
    !*FD \item[2] compute from effective strain rate

    integer :: which_disp = 0

    !*FD Flag that indicates method for computing the dissipation during the temperature calc.
    !*FD \begin{description}
    !*FD \item[0] for 0-order SIA approx
    !*FD \item[1] for 1-st order solution (e.g. Blatter-Pattyn)
    !*FD \item[2] for 1-st order depth-integrated solution (SSA)
    !*FD \end{description}

    integer :: which_ho_babc = 4

    !*FD Flag that describes basal boundary condition for PP dyn core: 
    !*FD \begin{description}
    !*FD \item[0] constant value (hardcoded in, useful for debugging)
    !*FD \item[1] simple pattern ("     ")
    !*FD \item[2] use till yield stress (Picard-type iteration)
    !*FD \item[3] circular ice shelf
    !*FD \item[4] no slip everywhere (using stress basal bc and large value for B^2)
    !*FD \item[5] beta^2 field passed in from CISM
    !*FD \item[6] no slip everywhere (using Dirichlet, no slip basal bc)
    !*FD \item[7] use till yield stress (Newton-type iteration)
    !*FD \end{description}

    integer :: which_ho_nonlinear = 0
    !*FD Flag that indicates method for solving the nonlinear iteration when solving 
    !*FD the first-order momentum balance
    !*FD \item[0] use the standard Picard iteration
    !*FD \item[1] use Jacobian Free Newton Krylov (JFNK) method

    integer :: which_ho_resid = 3
    !*FD Flag that indicates method for computing residual in PP dyn core: 
    !*FD \begin{description}
    !*FD \item[0] maxval 
    !*FD \item[1] maxval ignoring basal velocity 
    !*FD \item[2] mean value
    !*FD \item[3] L2 norm of system residual, Ax-b=resid
    !*FD \begin{description}

    integer :: which_ho_sparse = 0
    !*FD Flag that indicates method for solving the sparse linear system
    !*FD that arises from the higher-order solver
    !*FD \begin{description}
    !*FD \item[0] Biconjugate gradient, incomplete LU preconditioner
    !*FD \item[1] GMRES, incomplete LU preconditioner
    !*FD \item[2] Conjugate gradient, incomplete LU preconditioner
    !*FD \item[3] Conjugate gradient, structured grid, parallel-enabled
    !*FD \item[4] standalone interface to Trilinos
    !*FD \end{description}

    ! parameters to store external dycore options/information -- Doug Ranken 04/20/12
    integer :: external_dycore_type = 0  
    !*FD Flag to select an external dynamic core.
    !*FD \begin{description}
    !*FD \item[0] Do not use an external dynamic core
    !*FD \item[1] Use the BISICLES external dynamic core
    !*FD \end{description}

    character(fname_length) :: dycore_input_file=''
    !FD Name of a file containing external dycore settings.

!WHL - Added a glissade option to choose which Stokes approximation (SIA, SSA or Blatter-Pattyn HO)
!      Commented out for now

    ! Blatter-Pattyn HO by default
!!    integer :: which_ho_approx = 2    
    !*FD Flag that indicates which Stokes approximation to use in the glissade dycore.
    !*FD Not valid for other dycores 
    !*FD \begin{description}
    !*FD \item[0] Shallow-ice approximation, vertical shear stress only
    !*FD \item[1] Shallow-shelf approximation, horizontal-plane stresses only
    !*FD \item[2] Blatter-Pattyn with both vertical-shear and horizontal-plane stresses
    !*FD \end{description}

    ! The remaining options are not currently supported

    !integer :: which_bproc = 0
    !Options for the basal processes code
    !*FD \begin{description}
    !*FD \item[0] Disabled
    !*FD \item[1] Full calculation, with at least 3 nodes to represent the till layer
    !*FD \item[2] Fast calculation, using Tulaczyk empirical parametrization
    !*FD \end{description}

    !integer :: use_plume = 0   !! Option to be supported in future releases
    !*FD \begin{description}
    !*FD \item[0] standard bmlt calculation
    !*FD \item[1] use plume to calculate bmlt
    !*FD \end{description}

  end type glide_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geometry

    !*FD Holds fields and other information relating to the
    !*FD geometry of the ice sheet and bedrock.

    real, dimension(:,:), pointer :: temporary0 => null()
    !*FD temporary array used for masking velocity grid
    real, dimension(:,:), pointer :: temporary1 => null()
    !*FD temporary array used for masking temperature grid

    real(dp),dimension(:,:),pointer :: thck => null()
    !*FD The thickness of the ice, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: usrf => null()
    !*FD The elevation of the upper ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: lsrf => null() 
    !*FD The elevation of the lower ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: topg => null() 
    !*FD The elevation of the topography, divided by \texttt{thk0}.

    real(dp),dimension(:,:,:),pointer :: age => null()
    !*FD The age of a given ice layer, divided by \texttt{tim0}.

    integer, dimension(:,:),pointer :: thkmask => null()
    !*FD see glide_mask.f90 for possible values

    real(dp),dimension(:,:),pointer :: marine_bc_normal => null()
    !*FD NaN for all points except those that occur on the marine
    !*FD margin of an ice shelf, in which case contains the angle
    !*FD of the normal to the ice front. 

    integer, dimension(:,:),pointer :: thck_index => null()
    ! Set to nonzero integer for ice-covered cells (thck > 0), cells adjacent to ice-covered cells,
    !  and cells with acab > 0.  The non-zero points are numbered in sequence from the bottom left 
    !  to the top right, going along the rows.

    integer :: totpts = 0       ! total number of points with nonzero thck_index
    logical :: empty = .true.   ! true if totpts = 0

    real(dp) :: ivol, iarea,iareag, iareaf !*FD ice volume and ice area

  end type glide_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Some of these are PBJ variables that can be removed.

  type glide_geomderv

    !*FD Holds the horizontal and temporal derivatives of the thickness and
    !*FD upper surface elevation, as well as the thickness on the staggered grid.

    !*tb* Added a bunch of stuff here to clean up the higher order code that
    !I've been writing.  Might be worth it to add a mechanism to conditionally
    !allocate these depending on whether they are needed by the SIA core or by
    !the higher-order extensions

    !First derivatives on a staggered grid
    real(dp),dimension(:,:),pointer :: dthckdew => null() !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dlsrfdew => null() !*tb* added
    real(dp),dimension(:,:),pointer :: dlsrfdns => null() !*tb* added

    !Second derivatives on a staggered grid
    !*tb* added all of these
    ! Used by glam_strs2
    real(dp),dimension(:,:),pointer :: d2usrfdew2 => null()
    real(dp),dimension(:,:),pointer :: d2usrfdns2 => null()
    real(dp),dimension(:,:),pointer :: d2thckdew2 => null()
    real(dp),dimension(:,:),pointer :: d2thckdns2 => null()

    !First derivatives on a nonstaggered grid
    !*tb* added all of these
    !TODO - I think these can be removed
    real(dp),dimension(:,:),pointer :: dthckdew_unstag => null() !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew_unstag => null() !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns_unstag => null() !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns_unstag => null() !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dlsrfdew_unstag => null()
    real(dp),dimension(:,:),pointer :: dlsrfdns_unstag => null()

    !Second derivatives on a nonstaggered grid
    !*tb* added all of these
    !TODO - I think these can be removed
    real(dp),dimension(:,:),pointer :: d2usrfdew2_unstag => null()
    real(dp),dimension(:,:),pointer :: d2usrfdns2_unstag => null()
    real(dp),dimension(:,:),pointer :: d2thckdew2_unstag => null()
    real(dp),dimension(:,:),pointer :: d2thckdns2_unstag => null()

    !Time derivatives
    real(dp),dimension(:,:),pointer :: dthckdtm => null() !*FD Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm => null() !*FD Temporal derivative of upper surface elevation.

    !TODO - These staggered variables are not derivatives.
    !       Shoud they be part of the glide_geometry type?

    !Staggered grid versions of geometry variables
    real(dp),dimension(:,:),pointer :: stagthck => null() !*FD Thickness averaged onto the staggered grid.

    !*tb* added everything below
    real(dp),dimension(:,:),pointer :: stagusrf => null() !*FD Upper surface averaged onto the staggered grid
    real(dp),dimension(:,:),pointer :: staglsrf => null() !*FD Lower surface averaged onto the staggered grid
    real(dp),dimension(:,:),pointer :: stagtopg => null() !*FD Bedrock topography averaged onto the staggered grid
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

    !*FD Holds the velocity fields in 2D and 3D. At least some of these fields
    real(dp),dimension(:,:,:),pointer :: uvel  => null()   !*FD 3D $x$-velocity.
    real(dp),dimension(:,:,:),pointer :: vvel  => null()   !*FD 3D $y$-velocity.
    real(dp),dimension(:,:,:),pointer :: velnorm => null() ! horizontal ice speed
    real(dp),dimension(:,:,:),pointer :: wvel  => null()   !*FD 3D $z$-velocity.
    real(dp),dimension(:,:,:),pointer :: wgrd  => null()   !*FD 3D grid vertical velocity.
    real(dp),dimension(:,:)  ,pointer :: uflx  => null()   !*FD 
    real(dp),dimension(:,:)  ,pointer :: vflx  => null()   !*FD 
    real(dp),dimension(:,:)  ,pointer :: diffu => null()   !*FD 
    real(dp),dimension(:,:)  ,pointer :: diffu_x => null() !*sfp* moved from velocity_hom deriv type
    real(dp),dimension(:,:)  ,pointer :: diffu_y => null() 
    real(dp),dimension(:,:)  ,pointer :: total_diffu => null() !*FD total diffusivity
    real(dp),dimension(:,:)  ,pointer :: ubas  => null()   !*FD 
    real(dp),dimension(:,:)  ,pointer :: ubas_tavg  => null()
    real(dp),dimension(:,:)  ,pointer :: vbas  => null()   !*FD 
    real(dp),dimension(:,:)  ,pointer :: vbas_tavg  => null() 

    !! next 3 used for output of residual fields (when relevant code in glam_strs2 is active)
!    real(dp),dimension(:,:,:),pointer :: ures => null() !*FD 3D $x$-residual.
!    real(dp),dimension(:,:,:),pointer :: vres  => null() !*FD 3D $y$-residual.
!    real(dp),dimension(:,:,:),pointer :: magres  => null() !*FD 3D $magnitude$-residual.

    !! WHL - next 2 used for output of uvel, vvel on ice grid 
    !! (e.g., for problems with periodic BC, where the number of velocity points is
    !!  equal to the number of grid cells)
    real(dp),dimension(:,:,:),pointer :: uvel_icegrid => null() !*FD 3D $x$-velocity
    real(dp),dimension(:,:,:),pointer :: vvel_icegrid => null() !*FD 3D $x$-velocity

    real(dp),dimension(:,:)  ,pointer :: bed_softness => null() !*FD bed softness parameter
    real(dp),dimension(:,:)  ,pointer :: btrc  => null()        !*FD  basal traction (scaler field)
    real(dp),dimension(:,:,:),pointer :: btraction => null()    !*FD x(1,:,:) and y(2,:,:) "consistent" basal traction fields 
    real(dp),dimension(:,:)  ,pointer :: beta  => null()        !*FD basal shear coefficient
                                                                !*FD calculated from matrix coeffs in PP dyn core

    !*FD A mask similar to glide_geometry%thck_index, but on the velocity grid instead of the
    !*FD ice grid.  This is to aid in converging higher-order velocities
    integer, dimension(:,:), pointer    :: velmask => null()

    !*FD A mask that specifies where the velocity being read in should be held constant as a dirichlet condition
    integer, dimension(:,:), pointer    :: kinbcmask => null()

    !*sfp* mask on vel grid showing which dyn bc is applied at each grid cell (mainly for debugging)
    integer, dimension(:,:), pointer    :: dynbcmask => null()

  end type glide_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_stress_t      

    type(glide_tensor) :: tau
    real(dp),dimension(:,:,:),pointer :: efvs => null()
    real(dp),dimension(:,:)  ,pointer :: tau_x => null() !*FD SIA basal shear stress, x-dir
    real(dp),dimension(:,:)  ,pointer :: tau_y => null() !*FD SIA basal shear stress, y-dir
!*sfp* neither of the next two are currently used anywhere in the code (moved here from velocity_hom)
!    real(dp),dimension(:,:,:)  ,pointer :: gdsx => null() !*FD basal shear stress, x-dir
!    real(dp),dimension(:,:,:)  ,pointer :: gdsy => null() !*FD basal shear stress, y-dir

  end type glide_stress_t      

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!TODO - Should calving and eus be part of some other type?
!      Commented out slidconst, used by unsupported Huybrechts basal traction option
!TODO - Change to dp

  type glide_climate
     !*FD Holds fields used to drive the model
     real(sp),dimension(:,:),pointer :: acab     => null() !*FD Annual mass balance.
     real(sp),dimension(:,:),pointer :: acab_tavg     => null() !*FD Annual mass balance (time average).
     real(sp),dimension(:,:),pointer :: artm     => null() !*FD Annual mean air temperature
     real(sp),dimension(:,:),pointer :: lati     => null() !*FD Latitudes of model grid points
     real(sp),dimension(:,:),pointer :: loni     => null() !*FD Longitudes of model grid points
     real(sp),dimension(:,:),pointer :: calving  => null() !*FD Calving flux (scaled as mass balance, thickness, etc)
     real(sp) :: eus = 0.0                                 !*FD eustatic sea level
!!     real(sp) :: slidconst = 0.0     ! not currently used
  end type glide_climate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_temper

    !*FD Holds fields relating to temperature.

    !Note: In the glide dycore, temp and flwa live on the unstaggered vertical grid
    !       at layer interfaces and have vertical dimension (1:upn).
    !      In the glam/glissade dycore, with remapping advection of temperature, 
    !       temp and flwa live on the staggered vertical grid at layer midpoints.  
    !       The vertical dimensions are temp(0:upn) and flwa(1:upn-1).
    !
    !      bheatflx, ucondflx, and lcondflx are defined as positive down,
    !       so they will often be < 0.  
    !      However, bfricflx and dissipcol are defined to be >= 0.
    !
    !      If bheatflx is read from a data file, be careful about the sign!
    !      In input data, the geothermal heat flux is likely to be defined as positive upward.
    !
    !TODO: Create separate fields for basal melt beneath grounded and floating ice.


    real(dp),dimension(:,:,:),pointer :: temp => null()      !*FD 3D temperature field.
    real(dp),dimension(:,:),  pointer :: bheatflx => null()  !*FD basal heat flux (geothermal, positive down)
    real(dp),dimension(:,:,:),pointer :: flwa => null()      !*FD Glen's flow factor $A$.
    real(dp),dimension(:,:),  pointer :: bwat => null()      !*FD Basal water depth
    real(dp),dimension(:,:),  pointer :: bwatflx => null()   !*FD Basal water flux 
    real(dp),dimension(:,:),  pointer :: stagbwat => null()  !*FD Basal water depth on velo grid
    real(dp),dimension(:,:),  pointer :: bmlt => null()      !*FD Basal melt-rate (> 0 for melt, < 0 for freeze-on)
    real(dp),dimension(:,:),  pointer :: bmlt_tavg => null() !*FD Basal melt-rate
    real(dp),dimension(:,:),  pointer :: stagbtemp => null() !*FD Basal temperature on velo grid
    real(dp),dimension(:,:),  pointer :: bpmp => null()      !*FD Basal pressure melting point
    real(dp),dimension(:,:),  pointer :: stagbpmp => null()  !*FD Basal pressure melting point on velo grid
    real(dp),dimension(:,:),  pointer :: bfricflx => null()  !*FD basal heat flux from friction (>= 0)
    real(dp),dimension(:,:),  pointer :: ucondflx => null()  !*FD conductive heat flux at upper sfc (positive down)
    real(dp),dimension(:,:),  pointer :: lcondflx => null()  !*FD conductive heat flux at lower sfc (positive down)
    real(dp),dimension(:,:),  pointer :: dissipcol => null() !*FD total heat dissipation in column (>= 0)
    
    integer  :: niter   = 0      !*FD
    real(sp) :: perturb = 0.0    !*FD
    real(sp) :: grid    = 0.0    !*FD
    integer  :: tpt     = 0      !*FD Pointer to time series data
    logical  :: first1  = .true. !*FD
    logical  :: newtemps = .false. !*FD new temperatures
  end type glide_temper

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_lithot_type
     !*FD holds variables for temperature calculations in the lithosphere

     real(dp),dimension(:,:,:),pointer :: temp => null()    !*FD Three-dimensional temperature field.
     logical, dimension(:,:), pointer :: mask => null()     !*FD whether the point has been ice covered at some time

     integer :: num_dim = 1                                 !*FD either 1 or 3 for 1D/3D calculations

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

     real(dp), dimension(:), pointer :: deltaz => null()    !*FD array holding grid spacing in z
     real(dp), dimension(:,:), pointer :: zfactors => null()!*FD array holding factors for finite differences of vertical diffu
     real(dp) :: xfactor,yfactor !*FD factors for finite differences of horizontal diffu


     real :: surft = 2.         !*FD surface temperature, used for calculating initial temperature distribution
     real :: mart  = 2.         !*FD sea floor temperature 
     integer :: nlayer = 20     !*FD number of layers in lithosphere
     real :: rock_base = -5000. !*FD depth below sea-level at which geothermal heat gradient is applied
     
     integer :: numt = 0        !*FD number time steps for spinning up GTHF calculations

     real(dp) :: rho_r = 3300.0d0 !*FD The density of lithosphere (kg m$^{-3}$)
     real(dp) :: shc_r = 1000.0d0 !*FD specific heat capcity of lithosphere (J kg$^{-1}$ K$^{-1}$)
     real(dp) :: con_r = 3.3d0    !*FD thermal conductivity of lithosphere (W m$^{-1}$ K$^{-1}$)

     real(dp) :: diffu = 0. !*FD diffusion coefficient

  end type glide_lithot_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type isos_elastic
     !*FD Holds data used by isostatic adjustment calculations

     real(dp) :: d = 0.24e25                !*FD flexural rigidity  !TODO - Units?
     real(dp) :: lr                         !*FD radius of relative stiffness
     real(dp) :: a                          !*FD radius of disk
     real(dp) :: c1,c2,cd3,cd4              !*FD coefficients
     real(dp), dimension(:,:), pointer :: w !*FD matrix operator for lithosphere deformation
     integer :: wsize                       !*FD size of operator (0:rbel_wsize, 0:rbel_wsize), operator is axis symmetric
  end type isos_elastic

  type isostasy_type
     !*FD contains isostasy configuration

     ! do_isos has been replaced by model%options%isostasy
!!     logical :: do_isos = .false.    ! set to .true. if isostatic adjustment should be handled

     integer :: lithosphere = 0
     !*FD method for calculating equilibrium bedrock depression
     !*FD \begin{description}
     !*FD \item[0] local lithosphere, equilibrium bedrock depression is found using Archimedes' principle
     !*FD \item[1] elastic lithosphere, flexural rigidity is taken into account
     !*FD \end{description}

     integer :: asthenosphere = 0
     !*FD method for approximating the mantle
     !*FD \begin{description}
     !*FD \item[0] fluid mantle, isostatic adjustment happens instantaneously
     !*FD \item[1] relaxing mantle, mantle is approximated by a half-space
     !*FD \end{description}

     !TODO - Make these dp
     real :: relaxed_tau = 4000.    ! characteristic time constant of relaxing mantle (yr)
     real :: period = 500.          ! lithosphere update period (yr)
     real :: next_calc              ! when to update lithosphere
     logical :: new_load = .false.  ! set to true if there is a new surface load
     type(isos_elastic) :: rbel     ! structure holding elastic lithosphere setup

     real(dp),dimension(:,:),pointer :: relx => null()  ! elevation of relaxed topography, by \texttt{thck0}.
     real(dp),dimension(:,:),pointer :: load => null()  ! load imposed on lithosphere
     real(dp),dimension(:,:),pointer :: load_factors => null() ! temporary used for load calculation

  end type isostasy_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_funits
    character(fname_length) :: sigfile=''                      !*FD sigma coordinates file
    character(fname_length) :: ncfile=''                       !*FD configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !*FD first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !*FD first element of linked list defining netCDF inputs
  end type glide_funits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_numerics

    !*FD Parameters relating to the model numerics
    real(dp) :: tstart =    0.d0 !*FD starting time
    real(dp) :: tend   = 1000.d0 !*FD end time
    real(dp) :: time   =    0.d0 !*FD main loop counter in years
    real(dp) :: tinc   =    1.d0 !*FD time step of main loop in years 
    real(dp) :: ntem   =    1.d0 !*FD temperature time step (multiplier of main time step)
    real(dp) :: nvel   =    1.d0 !*FD velocity time step (multiplier of main time step)
    real(dp) :: alpha  =    0.5d0 !*FD richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas =    0.5d0 !*FD was a parameter in the original
    real(dp) :: thklim =  100.0   
    real(dp) :: mlimit = -200.0d0
    real(dp) :: calving_fraction = 0.8d0
    real(dp) :: dew    =   20.0d3
    real(dp) :: dns    =   20.0d3
    real(dp) :: dt     =    0.0
    real(dp) :: dttem  =    0.0
    real(sp) :: nshlf  =    0.0   !TODO - Change to dp
    integer  :: subcyc =    1
    real(dp) :: periodic_offset_ew = 0.d0 ! optional periodic_offsets for ismip-hom and similar tests
    real(dp) :: periodic_offset_ns = 0.d0 ! These may be needed to ensure continuous ice geometry at
                                          !  the edges of the global domain.

    integer  :: timecounter = 0   !*FD count time steps
    
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),pointer :: sigma => null() !*FD Sigma values for vertical spacing of 
                                                     !*FD model levels
    real(dp),dimension(:),pointer :: stagsigma => null() !*FD Staggered values of sigma (layer midpts)
    real(dp),dimension(:),pointer :: stagwbndsigma => null() !*FD Staggered values of sigma (layer midpts) with boundaries

    integer :: profile_period = 100            ! profile frequency

    !TODO - Compute ndiag as a function of dt_diag and pass to glide_diagnostics?
    !       This is more robust than computing mods of real numbers. 
    !TODO - Change names of idiag_global and jdiag_global?
    !       These are indices for the full ice sheet grid (before decomposition), but not a true global grid.
    real(dp) :: dt_diag = 0.d0            ! diagnostic time interval (write diagnostics every dt_diag years)
    integer  :: ndiag = -999              ! diagnostic period (write output every ndiag steps)
    integer  :: idiag_global = 1          ! grid indices for diagnostic point
    integer  :: jdiag_global = 1
  end type glide_numerics

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! variables for tracking the grounding line    
  type glide_grnd
    real(dp),dimension(:,:),pointer :: gl_ew => null()
    real(dp),dimension(:,:),pointer :: gl_ns => null()
    real(dp),dimension(:,:),pointer :: gline_flux => null() !*FD flux at the
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
    real(dp) :: trc0   = 0.0
    real(dp) :: trcmin = 0.0d0
    real(dp) :: marine = 1.0d0
    real(dp) :: trcmax = 10.0d0
    real(dp) :: btrac_const = 0.0d0  !TODO - Do we need two of these?  Also in glide_paramets below.
    real(dp) :: btrac_slope = 0.0d0
    real(dp) :: btrac_max = 0.d0
  end type glide_velowk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_thckwk
     real(dp),dimension(:,:),  pointer :: oldthck   => null()
     real(dp),dimension(:,:),  pointer :: oldthck2  => null()
     real(dp),dimension(:,:),pointer :: float => null()
     real(dp),dimension(:,:,:),pointer :: olds      => null()
     integer  :: nwhich  = 2
     real(dp) :: oldtime = 0.d0
     
     real(dp), dimension(:), pointer :: alpha => null()
     real(dp), dimension(:), pointer :: beta  => null()
     real(dp), dimension(:), pointer :: gamma => null()
     real(dp), dimension(:), pointer :: delta => null()

  end type glide_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_tempwk
    real(dp),dimension(:,:,:),pointer :: inittemp => null()
    real(dp),dimension(:,:,:),pointer :: dissip   => null()
    real(dp),dimension(:,:,:),pointer :: compheat => null()
    real(dp),dimension(:,:,:),pointer :: initadvt => null()
    real(dp),dimension(:),    pointer :: dupa     => null()
    real(dp),dimension(:),    pointer :: dupb     => null()
    real(dp),dimension(:),    pointer :: dupc     => null()
    real(dp),dimension(:),    pointer :: c1       => null()
    real(dp),dimension(:,:),  pointer :: dups     => null()
    real(dp),dimension(:,:),  pointer :: wphi     => null()
    real(dp),dimension(:,:),  pointer :: bwatu    => null()
    real(dp),dimension(:,:),  pointer :: bwatv    => null()
    real(dp),dimension(:,:),  pointer :: fluxew   => null()
    real(dp),dimension(:,:),  pointer :: fluxns   => null()
    real(dp),dimension(:,:),  pointer :: bint     => null()
    real(dp),dimension(:,:),  pointer :: smth     => null()
    real(dp),dimension(:,:,:),pointer :: hadv_u   => null()
    real(dp),dimension(:,:,:),pointer :: hadv_v   => null()

    !TODO - Do we need these?
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

  !TODO - Set btrac_const to nonzero value?

  type glide_paramets
    real(dp),dimension(5) :: bpar = (/ 0.2d0, 0.5d0, 0.0d0 ,1.0d-2, 1.0d0/)
    real(dp) :: btrac_const = 0.d0     ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_slope = 0.0d0    ! Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_max = 0.d0       ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: geot   = -5.0d-2       ! W m^{-2}, positive down
    real(dp) :: flow_factor = 3.0d0    ! "fiddle" parameter for the Arrhenius relationship
    real(dp) :: slip_ratio = 1.0d0     ! Slip ratio, used only in higher order code when the slip ratio beta computation is requested
    real(dp) :: hydtim = 1000.0d0      ! years, converted to s^{-1} and scaled
                                       ! 0 if no drainage
    real(dp) :: bwat_smooth = 0.01d0   ! basal water field smoothing strength
    real(dp) :: default_flwa = 1.0d-16 ! Glen's A to use in isothermal case, in units Pa^{-n} yr^{-1} 
                                       ! (would change to e.g. 4.6e-18 in EISMINT-ROSS case)
  end type glide_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - This type is not currently used.  Should it be removed?
  type glide_basalproc
    !Tuneables, set in the config file 
    real(dp):: fric=0.45d0                   ! Till coeff of internal friction: ND
    real(dp):: etillo=0.7d0                  ! Till void ratio at No
    real(dp):: No=1000.d0                    ! Reference value of till effective stress
    real(dp):: Comp=0.12d0                   ! Till coeff of compressibility: ND
    real(dp):: Cv = 1.0d-8                   ! Till hydraulic diffusivity: m2/s
    real(dp):: Kh = 1.0d-10                  !Till hydraulic conductivity: m/s
    real(dp):: Zs = 3.0d0                    ! Solid till thickness: m
    real(dp):: aconst=994000000d0            ! Constant in till strength eq. (Pa)
    real(dp):: bconst=21.7                   ! Constant in till strength eq. (ND)
    integer:: till_hot = 0
    integer:: tnodes = 5

    real(dp), dimension (:) , pointer :: till_dz => null()  !holds inital till layer spacing - 
    
    !Model variables that will be passed to other subroutines
    real(dp),dimension(:,:)  ,pointer :: minTauf => null() !Bed strength calculated with basal proc. mod.
    real(dp),dimension(:,:)  ,pointer :: Hwater  => null() !Water available from till layer (m)
    !Model variabled necessary for restart
    real(dp),dimension(:,:,:)  ,pointer :: u => null()     !Till excess pore pressure (Pa)
    real(dp),dimension(:,:,:)  ,pointer :: etill  => null()  !Till void ratio (ND)  
    
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

!TODO - Is this type used?
  type glide_phaml
    real(dp),dimension(:,:),pointer :: uphaml => null()
    real(dp),dimension(:,:),pointer :: init_phaml => null()
    real(dp),dimension(:,:),pointer :: rs_phaml => null()
    !maybe put the x/y vectors here too just for simplicity
  end type glide_phaml
  
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

!TODO - KJE - remove once new glide_global_type is working and we can use those ewn and nsn
    integer :: ewn
    integer :: nsn

  end type glissade_solver

       
  type glide_global_type    ! type containing all of the above for an ice sheet model instance
    integer              :: model_id !*FD Used in the global model list for error handling purposes
    type(glide_general)  :: general
    type(glide_options)  :: options
    type(glide_geometry) :: geometry
    type(glide_geomderv) :: geomderv
    type(glide_velocity) :: velocity
    type(glide_stress_t) :: stress   
    type(glide_climate)  :: climate
    type(glide_temper)   :: temper
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
    type(glide_phaml)    :: phaml
    type(glide_grnd)     :: ground
    type(glissade_solver):: solver_data

  end type glide_global_type

contains

  !TODO - Make sure these itemized lists are complete.

  subroutine glide_allocarr(model)    
    !*FD Allocates the model arrays, and initialises some of them to zero.
    !*FD These are the arrays allocated, and their dimensions:
    !*FD
    !*FD In \texttt{model\%temper}:
    !*FD \begin{itemize}
    !*FD \item \texttt{temp(upn,0:ewn+1,0:nsn+1))}   !WHL - 2 choices
    !*FD \item \texttt{bheatflx(ewn,nsn))}
    !*FD \item \texttt{flwa(upn,ewn,nsn))}           !WHL - 2 choices
    !*FD \item \texttt{bwat(ewn,nsn))}
    !*FD \item \texttt{bmlt(ewn,nsn))}
    !*FD \item \texttt{bfricflx(ewn,nsn))}
    !*FD \item \texttt{ucondflx(ewn,nsn))}
    !*FD \item \texttt{lcondflx(ewn,nsn))}
    !*FD \item \texttt{dissipcol(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%velocity}:
    !*FD \begin{itemize}
    !*FD \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{wvel(upn,ewn,nsn))}
    !*FD \item \texttt{wgrd(upn,ewn,nsn))}
    !*FD \item \texttt{uflx(ewn-1,nsn-1))}
    !*FD \item \texttt{vflx(ewn-1,nsn-1))}
    !*FD \item \texttt{diffu(ewn,nsn))}
    !*FD \item \texttt{btrc(ewn,nsn))}
    !*FD \item \texttt{ubas(ewn,nsn))}
    !*FD \item \texttt{vbas(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%climate}:
    !*FD \begin{itemize}
    !*FD \item \texttt{acab(ewn,nsn))}
    !*FD \item \texttt{artm(ewn,nsn))}
    !*FD \item \texttt{lati(ewn,nsn))}
    !*FD \item \texttt{loni(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%geomderv}:
    !*FD \begin{itemize}
    !*FD \item \texttt{dthckdew(ewn,nsn))}
    !*FD \item \texttt{dusrfdew(ewn,nsn))}
    !*FD \item \texttt{dthckdns(ewn,nsn))}
    !*FD \item \texttt{dusrfdns(ewn,nsn))}
    !*FD \item \texttt{dthckdtm(ewn,nsn))}
    !*FD \item \texttt{dusrfdtm(ewn,nsn))}
    !*FD \item \texttt{stagthck(ewn-1,nsn-1))}
    !*FD \end{itemize}
  
    !*FD In \texttt{model\%geometry}:
    !*FD \begin{itemize}
    !*FD \item \texttt{thck(ewn,nsn))}
    !*FD \item \texttt{usrf(ewn,nsn))}
    !*FD \item \texttt{lsrf(ewn,nsn))}
    !*FD \item \texttt{topg(ewn,nsn))}
    !*FD \item \texttt{mask(ewn,nsn))}
    !*FD \item \texttt{age(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%thckwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%numerics}:
    !*FD \begin{itemize}
    !*FD \item \texttt{sigma(upn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%numerics}:
    !*FD \begin{itemize}
    !*FD \item \texttt{stagsigma(upn-1))}
    !*FD \end{itemize}

!KJE add the new ones, once working and complete

    use glimmer_log
    use parallel

    implicit none

    type(glide_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn

    ! for simplicity, copy these values...

    ewn=model%general%ewn
    nsn=model%general%nsn
    upn=model%general%upn
    
    ! Allocate appropriately

    allocate(model%general%x0(ewn-1))!; model%general%x0 = 0.0
    allocate(model%general%y0(nsn-1))!; model%general%y0 = 0.0
    allocate(model%general%x1(ewn))!; model%general%x1 = 0.0
    allocate(model%general%y1(nsn))!; model%general%y1 = 0.0
    call coordsystem_allocate(model%general%ice_grid, model%temper%bheatflx)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bwat)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bwatflx)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbwat)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bmlt)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bmlt_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbtemp)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bpmp)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bwatflx)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbpmp)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bfricflx)
    call coordsystem_allocate(model%general%ice_grid, model%temper%ucondflx)
    call coordsystem_allocate(model%general%ice_grid, model%temper%lcondflx)
    call coordsystem_allocate(model%general%ice_grid, model%temper%dissipcol)

!NOTE: - In the glide dycore (whichdycore = DYCORE_GLIDE), temperature and 
!        flow factor live on the unstaggered vertical grid, and extra rows and columns 
!        (with indices 0:ewn+1, 0:nsn+1) are needed.
!      - In the glam/glissade dycore, temperature and flow factor live on the 
!        staggered vertical grid, with temperature and flwa defined at the
!        center of each layer k = 1:upn-1.  The temperature (but not flwa)
!        is defined at the upper surface (k = 0) and lower surface (k = upn).

    if (model%options%whichdycore == DYCORE_GLIDE) then
       allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1))
       call coordsystem_allocate(model%general%ice_grid, upn, model%temper%flwa)
    else   ! glam/glissade dycore
       allocate(model%temper%temp(0:upn,1:ewn,1:nsn))
       call coordsystem_allocate(model%general%ice_grid, upn-1, model%temper%flwa)
    endif

    ! MJH set these to physically unrealistic values so we can tell later if 
    !  arrays were initialized correctly
    model%temper%temp(:,:,:) = -999.0d0
    model%temper%flwa(:,:,:) = -999.0d0
 
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%uvel)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%vvel)

    !! next 3 used for output of residual fields (when relevant code in glam_strs2 is active)
!    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%ures)
!    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%vres)
!    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%magres)

    call coordsystem_allocate(model%general%ice_grid, upn, model%velocity%uvel_icegrid)
    call coordsystem_allocate(model%general%ice_grid, upn, model%velocity%vvel_icegrid)

    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%velnorm)
    call coordsystem_allocate(model%general%ice_grid, upn, model%velocity%wvel)
    call coordsystem_allocate(model%general%ice_grid, upn, model%velocity%wgrd)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%uflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu_x)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu_y)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%total_diffu)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%bed_softness)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%btrc)
    call coordsystem_allocate(model%general%velo_grid, 2, model%velocity%btraction)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%beta)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%velmask)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%kinbcmask)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%dynbcmask)

    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%scalar)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%xz)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%yz)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%xx)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%yy)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%tau%xy)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%stress%efvs)
    call coordsystem_allocate(model%general%velo_grid, model%stress%tau_x)
    call coordsystem_allocate(model%general%velo_grid, model%stress%tau_y)

    call coordsystem_allocate(model%general%ice_grid, model%climate%acab)
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab_tavg)
    call coordsystem_allocate(model%general%ice_grid, model%climate%artm)
    call coordsystem_allocate(model%general%ice_grid, model%climate%lati)
    call coordsystem_allocate(model%general%ice_grid, model%climate%loni)
    call coordsystem_allocate(model%general%ice_grid, model%climate%calving)

    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dlsrfdew)    
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdns)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdns)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dlsrfdns)
    
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2usrfdew2)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2usrfdns2)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2thckdew2)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%d2thckdns2)
    
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dthckdew_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dusrfdew_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dlsrfdew_unstag)    
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dthckdns_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dusrfdns_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dlsrfdns_unstag)
    
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%d2usrfdew2_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%d2usrfdns2_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%d2thckdew2_unstag)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%d2thckdns2_unstag)

    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dthckdtm)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dusrfdtm)

    call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagthck)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%staglsrf)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagusrf)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagtopg)
  
    call coordsystem_allocate(model%general%velo_grid, model%geometry%temporary0)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%temporary1)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thck)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%usrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%lsrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%topg)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thck_index)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thkmask)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%marine_bc_normal)
    call coordsystem_allocate(model%general%ice_grid, upn-1, model%geometry%age)

    allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
    model%thckwk%olds = 0.0d0
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck)
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck2)
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%float)

    ! If we already have sigma, don't reallocate
    if (associated(model%numerics%sigma)) then
       if (size(model%numerics%sigma)/=upn) then
          call write_log('Wrong number of sigma levels given',GM_FATAL)
       end if
    else
       allocate(model%numerics%sigma(upn))
    endif

    allocate(model%numerics%stagsigma(upn-1))
    allocate(model%numerics%stagwbndsigma(0:upn))  !MJH added (0:upn) as separate variable

    ! allocate memory for grounding line
    allocate (model%ground%gl_ew(ewn-1,nsn))
    allocate (model%ground%gl_ns(ewn,nsn-1))
    allocate (model%ground%gline_flux(ewn,nsn)) 
    allocate (model%solver_data%rhsd(ewn*nsn))
    allocate (model%solver_data%answ(ewn*nsn))

    call new_sparse_matrix(ewn*nsn, 5*ewn*nsn, model%solver_data%matrix)

    !TODO - The lithosphere temperature has the vertical layer as the 3rd index,
    !        whereas the ice temperature has the vertical layer as the 1st index.
    !       Should we switch to the ice temperature convention?
    allocate(model%lithot%temp(1:ewn,1:nsn,model%lithot%nlayer)); model%lithot%temp = 0.d0
    call coordsystem_allocate(model%general%ice_grid, model%lithot%mask)

    ! allocate isostasy grids

    call coordsystem_allocate(model%general%ice_grid, model%isostasy%relx)
    call coordsystem_allocate(model%general%ice_grid, model%isostasy%load)
    call coordsystem_allocate(model%general%ice_grid, model%isostasy%load_factors)

    !TODO - Are these needed?
    !allocate phaml variables
    call coordsystem_allocate(model%general%ice_grid, model%phaml%init_phaml)
    call coordsystem_allocate(model%general%ice_grid, model%phaml%rs_phaml)
    call coordsystem_allocate(model%general%ice_grid, model%phaml%uphaml)

    !allocate basal processes variables
    call coordsystem_allocate(model%general%ice_grid, model%basalproc%Hwater)
    call coordsystem_allocate(model%general%velo_grid, model%basalproc%minTauf)
    allocate(model%basalproc%u (ewn-1,nsn-1,model%basalproc%tnodes)); model%basalproc%u=41.0d3
    allocate(model%basalproc%etill (ewn-1,nsn-1,model%basalproc%tnodes));model%basalproc%etill=0.5d0

  end subroutine glide_allocarr

  subroutine glide_deallocarr(model)

    !*FD deallocate model arrays
    !TODO - Verify that all arrays allocated above are deallocated here.

    implicit none
    type(glide_global_type),intent(inout) :: model

    if (associated(model%general%x0)) &
       deallocate(model%general%x0) 
    if (associated(model%general%y0)) &
       deallocate(model%general%y0) 
    if (associated(model%general%x1)) &
       deallocate(model%general%x1) 
    if (associated(model%general%y1)) &
       deallocate(model%general%y1) 

    if (associated(model%temper%temp)) &
       deallocate(model%temper%temp)
    if (associated(model%temper%flwa)) &
       deallocate(model%temper%flwa)
    if (associated(model%temper%bheatflx)) &
       deallocate(model%temper%bheatflx)
    if (associated(model%temper%bwat)) &
       deallocate(model%temper%bwat)
    if (associated(model%temper%bwatflx)) &
       deallocate(model%temper%bwatflx)
    if (associated(model%temper%stagbwat)) &
       deallocate(model%temper%stagbwat)
    if (associated(model%temper%bmlt)) &
       deallocate(model%temper%bmlt)
    if (associated(model%temper%bmlt_tavg)) &
       deallocate(model%temper%bmlt_tavg)
    if (associated(model%temper%bfricflx)) &
       deallocate(model%temper%bfricflx)
    if (associated(model%temper%ucondflx)) &
       deallocate(model%temper%ucondflx)
    if (associated(model%temper%lcondflx)) &
       deallocate(model%temper%lcondflx)
    if (associated(model%temper%dissipcol)) &
       deallocate(model%temper%dissipcol)
    if (associated(model%temper%stagbtemp)) &
       deallocate(model%temper%stagbtemp)
    if (associated(model%temper%bpmp)) &
       deallocate(model%temper%bpmp)
    if (associated(model%temper%stagbpmp)) &
       deallocate(model%temper%stagbpmp)
    if (associated(model%ground%gl_ns)) &
       deallocate(model%ground%gl_ns)
    if (associated(model%ground%gl_ew)) &
       deallocate(model%ground%gl_ew)
    if (associated(model%ground%gline_flux)) &
       deallocate(model%ground%gline_flux)

    if (associated(model%lithot%temp)) &
       deallocate(model%lithot%temp)
    if (associated(model%lithot%mask)) &
       deallocate(model%lithot%mask)

    if (associated(model%velocity%uvel)) &
       deallocate(model%velocity%uvel)
    if (associated(model%velocity%vvel)) &
       deallocate(model%velocity%vvel)

    !! next 3 used for output of residual fields (when relevant code in glam_strs2 is active)
!    deallocate(model%velocity%ures) 
!    deallocate(model%velocity%vres)
!    deallocate(model%velocity%magres)

    if (associated(model%velocity%uvel_icegrid)) &
       deallocate(model%velocity%uvel_icegrid)
    if (associated(model%velocity%vvel_icegrid)) &
       deallocate(model%velocity%vvel_icegrid)

    if (associated(model%velocity%velnorm)) &
       deallocate(model%velocity%velnorm)
    if (associated(model%velocity%wvel)) &
       deallocate(model%velocity%wvel)
    if (associated(model%velocity%wgrd)) &
       deallocate(model%velocity%wgrd)

    if (associated(model%velocity%uflx)) &
       deallocate(model%velocity%uflx)
    if (associated(model%velocity%vflx)) &
       deallocate(model%velocity%vflx)
    if (associated(model%velocity%diffu)) &
       deallocate(model%velocity%diffu)
    if (associated(model%velocity%diffu_x)) &
       deallocate(model%velocity%diffu_x)
    if (associated(model%velocity%diffu_y)) &
       deallocate(model%velocity%diffu_y)
    if (associated(model%velocity%total_diffu)) &
       deallocate(model%velocity%total_diffu)
    if (associated(model%velocity%bed_softness)) &
       deallocate(model%velocity%bed_softness)
    if (associated(model%velocity%btrc)) &
       deallocate(model%velocity%btrc)
    if (associated(model%velocity%btraction)) &
       deallocate(model%velocity%btraction)
    if (associated(model%velocity%beta)) &
       deallocate(model%velocity%beta)
    if (associated(model%velocity%ubas)) &
       deallocate(model%velocity%ubas)
    if (associated(model%velocity%ubas_tavg)) &
       deallocate(model%velocity%ubas_tavg)
    if (associated(model%velocity%vbas)) &
       deallocate(model%velocity%vbas)
    if (associated(model%velocity%vbas_tavg)) &
       deallocate(model%velocity%vbas_tavg)
    if (associated(model%velocity%velmask)) &
       deallocate(model%velocity%velmask)
    if (associated(model%velocity%kinbcmask)) &
       deallocate(model%velocity%kinbcmask)
    if (associated(model%velocity%dynbcmask)) &
       deallocate(model%velocity%dynbcmask)

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
    if (associated(model%stress%efvs)) &
       deallocate(model%stress%efvs)
    if (associated(model%stress%tau_x)) &
       deallocate(model%stress%tau_x)
    if (associated(model%stress%tau_y)) &
       deallocate(model%stress%tau_y)

    if (associated(model%climate%acab)) &
       deallocate(model%climate%acab)
    if (associated(model%climate%acab_tavg)) &
       deallocate(model%climate%acab_tavg)
    if (associated(model%climate%artm)) &
       deallocate(model%climate%artm)
    if (associated(model%climate%lati)) &
       deallocate(model%climate%lati)
    if (associated(model%climate%loni)) &
       deallocate(model%climate%loni)

    if (associated(model%geomderv%dthckdew)) &
       deallocate(model%geomderv%dthckdew)
    if (associated(model%geomderv%dusrfdew)) &
       deallocate(model%geomderv%dusrfdew)
    if (associated(model%geomderv%dlsrfdew)) &
       deallocate(model%geomderv%dlsrfdew)
    if (associated(model%geomderv%dthckdns)) &
       deallocate(model%geomderv%dthckdns)
    if (associated(model%geomderv%dusrfdns)) &
       deallocate(model%geomderv%dusrfdns)
    if (associated(model%geomderv%dlsrfdns)) &
       deallocate(model%geomderv%dlsrfdns)

    if (associated(model%geomderv%d2usrfdew2)) &
       deallocate(model%geomderv%d2usrfdew2)
    if (associated(model%geomderv%d2thckdew2)) &
       deallocate(model%geomderv%d2thckdew2)
    if (associated(model%geomderv%d2usrfdns2)) &
       deallocate(model%geomderv%d2usrfdns2)
    if (associated(model%geomderv%d2thckdns2)) &
       deallocate(model%geomderv%d2thckdns2)

    if (associated(model%geomderv%dthckdew_unstag)) &
       deallocate(model%geomderv%dthckdew_unstag)
    if (associated(model%geomderv%dusrfdew_unstag)) &
       deallocate(model%geomderv%dusrfdew_unstag)
    if (associated(model%geomderv%dlsrfdew_unstag)) &
       deallocate(model%geomderv%dlsrfdew_unstag)
    if (associated(model%geomderv%dthckdns_unstag)) &
       deallocate(model%geomderv%dthckdns_unstag)
    if (associated(model%geomderv%dusrfdns_unstag)) &
       deallocate(model%geomderv%dusrfdns_unstag)
    if (associated(model%geomderv%dlsrfdns_unstag)) &
       deallocate(model%geomderv%dlsrfdns_unstag)

    if (associated(model%geomderv%d2usrfdew2_unstag)) &
       deallocate(model%geomderv%d2usrfdew2_unstag)
    if (associated(model%geomderv%d2thckdew2_unstag)) &
       deallocate(model%geomderv%d2thckdew2_unstag)
    if (associated(model%geomderv%d2usrfdns2_unstag)) &
       deallocate(model%geomderv%d2usrfdns2_unstag)
    if (associated(model%geomderv%d2thckdns2_unstag)) &
       deallocate(model%geomderv%d2thckdns2_unstag)

    if (associated(model%geomderv%dthckdtm)) &
       deallocate(model%geomderv%dthckdtm)
    if (associated(model%geomderv%dusrfdtm)) &
       deallocate(model%geomderv%dusrfdtm)
    if (associated(model%geomderv%stagthck)) &
       deallocate(model%geomderv%stagthck)
    if (associated(model%geomderv%stagusrf)) &
       deallocate(model%geomderv%stagusrf)
    if (associated(model%geomderv%staglsrf)) &
       deallocate(model%geomderv%staglsrf)
    if (associated(model%geomderv%stagtopg)) &
       deallocate(model%geomderv%stagtopg)

    if (associated(model%geometry%temporary0)) &
       deallocate(model%geometry%temporary0)
    if (associated(model%geometry%temporary1)) &
       deallocate(model%geometry%temporary1)
    if (associated(model%geometry%thck)) &
       deallocate(model%geometry%thck)
    if (associated(model%geometry%usrf)) &
       deallocate(model%geometry%usrf)
    if (associated(model%geometry%lsrf)) &
       deallocate(model%geometry%lsrf)
    if (associated(model%geometry%topg)) &
       deallocate(model%geometry%topg)
    if (associated(model%geometry%age)) &
       deallocate(model%geometry%age)
    if (associated(model%geometry%thck_index)) &
       deallocate(model%geometry%thck_index)
    if (associated(model%geometry%thkmask)) &
       deallocate(model%geometry%thkmask)
    if (associated(model%geometry%marine_bc_normal)) &
       deallocate(model%geometry%marine_bc_normal)

    if (associated(model%thckwk%olds)) &
       deallocate(model%thckwk%olds)
    if (associated(model%thckwk%oldthck)) &
       deallocate(model%thckwk%oldthck)
    if (associated(model%thckwk%oldthck2)) &
       deallocate(model%thckwk%oldthck2)
    if (associated(model%thckwk%float)) &
       deallocate(model%thckwk%float)
    if (associated(model%numerics%sigma)) &
       deallocate(model%numerics%sigma)
    if (associated(model%numerics%stagsigma)) &
       deallocate(model%numerics%stagsigma)
    if (associated(model%numerics%stagwbndsigma)) &
       deallocate(model%numerics%stagwbndsigma)
    if (associated(model%solver_data%rhsd,model%solver_data%answ)) &
       deallocate(model%solver_data%rhsd,model%solver_data%answ)

!KJE do we need this at all here, the parts within are allocated in glam_strs2
    call del_sparse_matrix(model%solver_data%matrix)

    ! deallocate isostasy grids

    ! new isostasy
    if (associated(model%isostasy%relx)) &
       deallocate(model%isostasy%relx)
    if (associated(model%isostasy%load)) &
       deallocate(model%isostasy%load)
    if (associated(model%isostasy%load_factors)) &
       deallocate(model%isostasy%load_factors)

    !deallocate phaml variables
    if (associated(model%phaml%init_phaml)) &
       deallocate(model%phaml%init_phaml)
    if (associated(model%phaml%rs_phaml)) &
       deallocate(model%phaml%rs_phaml)    
    if (associated(model%phaml%uphaml)) &
       deallocate(model%phaml%uphaml)

    ! deallocate till variables
    if (associated(model%basalproc%Hwater)) &
       deallocate(model%basalproc%Hwater)
    if (associated(model%basalproc%minTauf)) &
       deallocate(model%basalproc%minTauf)
    if (associated(model%basalproc%u)) &
       deallocate(model%basalproc%u)
    if (associated(model%basalproc%etill)) &
       deallocate(model%basalproc%etill)

  end subroutine glide_deallocarr

  ! some accessor functions
  function get_dew(model)
    !*FD return scaled x node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dew
    type(glide_global_type) :: model

    get_dew = model%numerics%dew * len0
  end function get_dew

  function get_dns(model)
    !*FD return scaled y node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dns
    type(glide_global_type) :: model

    get_dns = model%numerics%dns * len0
  end function get_dns

  function get_tstart(model)
    !*FD return start time
    implicit none
    real(dp) :: get_tstart
    type(glide_global_type) :: model
    
    get_tstart = model%numerics%tstart
  end function get_tstart

  function get_tend(model)
    !*FD return end time
    implicit none
    real(dp) :: get_tend
    type(glide_global_type) :: model
    
    get_tend = model%numerics%tend
  end function get_tend

  function get_tinc(model)
    !*FD return time increment
    implicit none
    real(dp) :: get_tinc
    type(glide_global_type) :: model
    
    get_tinc = model%numerics%tinc
  end function get_tinc

  function get_ewn(model)
    !*FD get number of nodes in x dir
    implicit none
    integer get_ewn
    type(glide_global_type) :: model

    get_ewn = model%general%ewn
  end function get_ewn

  function get_nsn(model)
    !*FD get number of nodes in y dir
    implicit none
    integer get_nsn
    type(glide_global_type) :: model

    get_nsn = model%general%nsn
  end function get_nsn
  
  subroutine set_time(model,time)
    !*FD Set the model time counter --- useful for
    !*FD fractional year output
    implicit none
    type(glide_global_type) :: model
    real(dp) :: time

    model%numerics%time = time
  end subroutine set_time

end module glide_types
