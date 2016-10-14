module SoilcolType
  
  ! -------------------------------------------------------- 
  ! ALM sub-grid hierarchy:
  ! Define Soil column unit data types, with Init and Clean for each
  ! -------------------------------------------------------- 
  ! 10 Oct 2016, PET
  
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval, spval
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak

  implicit none
  save
  private

  ! sub-grid geospatial and physical properties defined at the soil_column level
  ! migrate variables list from ColumnType.F90

  type, public :: soilcol_properties
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: landunit             (:)   ! index into landunit level quantities
     real(r8), pointer :: wtlunit              (:)   ! weight (relative to landunit)
     integer , pointer :: gridcell             (:)   ! index into gridcell level quantities
     real(r8), pointer :: wtgcell              (:)   ! weight (relative to gridcell)
     integer , pointer :: pfti                 (:)   ! beginning pft index for each column
     integer , pointer :: pftf                 (:)   ! ending pft index for each column
     integer , pointer :: npfts                (:)   ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype                (:)   ! column type
     logical , pointer :: active               (:)   ! true=>do computations on this column 

     ! topography
     real(r8), pointer :: glc_topo             (:)   ! surface elevation (m)
     real(r8), pointer :: micro_sigma          (:)   ! microtopography pdf sigma (m)
     real(r8), pointer :: n_melt               (:)   ! SCA shape parameter
     real(r8), pointer :: topo_slope           (:)   ! gridcell topographic slope
     real(r8), pointer :: topo_std             (:)   ! gridcell elevation standard deviation

     ! vertical levels
     integer , pointer :: snl                  (:)   ! number of snow layers
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi                   (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii                  (:)   ! convective boundary height [m]
     real(r8), pointer :: dz_lake              (:,:) ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake               (:,:) ! layer depth for lake (m)
     real(r8), pointer :: lakedepth            (:)   ! variable lake depth (m)                             

   contains

     procedure, public :: Init => init_col_pp
     procedure, public :: Clean => clean_col_pp

  end type soilcol_properties


  type, public :: soilcol_energy_state

  contains
     procedure, public :: Init => init_col_es
     procedure, public :: Clean => clean_col_es
  end type soilcol_energy_state


  type, public :: soilcol_water_state

  contains
     procedure, public :: Init => init_col_ws
     procedure, public :: Clean => clean_col_ws
  end type soilcol_water_state


  type, public :: soilcol_carbon_state

  contains
     procedure, public :: Init => init_col_cs
     procedure, public :: Clean => clean_col_cs
  end type soilcol_carbon_state

  type, public :: soilcol_nitrogen_state

  contains
     procedure, public :: Init => init_col_ns
     procedure, public :: Clean => clean_col_ns
  end type soilcol_nitrogen_state

  type, public :: soilcol_phosphorus_state

  contains
     procedure, public :: Init => init_col_ps
     procedure, public :: Clean => clean_col_ps
  end type soilcol_phosphorus_state


  ! declare the public instances of soilcolumn types
  type(soilcol_properties) , public, target :: col_pp
  type(soilcol_energy_state), public, target :: col_es
  type(soilcol_water_state), public, target :: col_ws

  contains 

  subroutine init_col_pp(this, begc, endc)
    class(soilcol_properties) :: this
    integer, intent(in) :: begc   ! beginning soil column index
    integer, intent(in) :: endc   ! ending soil column index

    allocate(this%lnd         (begc:endc)) ; this%lnd         (:) = ispval
   
  end subroutine init_col_pp

  
  subroutine clean_col_pp(this)
    class(soilcol_properties) :: this
  
    deallocate(this%lnd)
  end subroutine clean_col_pp

end module Soilcol_type

