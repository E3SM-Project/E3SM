module SoilorderConType

  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use abortutils     , only : endrun
  !
  implicit none
  save
  public
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: soilorderconInit
  !
  ! !PUBLIC TYPES:
  type, public :: soilordercon_type
     real(r8), allocatable :: smax(:)
     real(r8), allocatable :: ks_sorption(:)
     real(r8), allocatable :: r_weather(:)
     real(r8), allocatable :: r_adsorp(:)
     real(r8), allocatable :: r_desorp(:)
     real(r8), allocatable :: r_occlude(:)
     real(r8), allocatable :: k_s1_biochem(:)
     real(r8), allocatable :: k_s2_biochem(:)
     real(r8), allocatable :: k_s3_biochem(:)
     real(r8), allocatable :: k_s4_biochem(:)


  end type soilordercon_type

  type(soilordercon_type), public :: soilordercon ! soil order constants

contains

  !-----------------------------------------------------------------------
  subroutine soilorderconInit()
    !
    ! !USES:
    use elm_varpar, only : nsoilorder
    use soilorder_varcon, only : smax
    use soilorder_varcon, only : ks_sorption
    use soilorder_varcon, only : r_weather,r_adsorp,r_desorp,r_occlude
    use soilorder_varcon, only : k_s1_biochem,k_s2_biochem,k_s3_biochem,k_s4_biochem
    !
    ! !LOCAL VARIABLES:
    integer :: m, ib
    !------------------------------------------------------------------------


    allocate(soilordercon%smax           (0:nsoilorder))        ; soilordercon%smax(:)        =nan
    allocate(soilordercon%ks_sorption    (0:nsoilorder))        ; soilordercon%ks_sorption(:) =nan
    allocate(soilordercon%r_weather      (0:nsoilorder))        ; soilordercon%r_weather(:)   =nan
    allocate(soilordercon%r_adsorp       (0:nsoilorder))        ; soilordercon%r_adsorp(:)    =nan
    allocate(soilordercon%r_desorp       (0:nsoilorder))        ; soilordercon%r_desorp(:)    =nan
    allocate(soilordercon%r_occlude      (0:nsoilorder))        ; soilordercon%r_occlude(:)    =nan

    allocate(soilordercon%k_s1_biochem      (0:nsoilorder))        ; soilordercon%k_s1_biochem(:)    =nan
    allocate(soilordercon%k_s2_biochem      (0:nsoilorder))        ; soilordercon%k_s2_biochem(:)    =nan
    allocate(soilordercon%k_s3_biochem      (0:nsoilorder))        ; soilordercon%k_s3_biochem(:)    =nan
    allocate(soilordercon%k_s4_biochem      (0:nsoilorder))        ; soilordercon%k_s4_biochem(:)    =nan


    do m = 0,nsoilorder


       soilordercon%smax(m)         = smax(m)
       soilordercon%ks_sorption(m)         = ks_sorption(m)
       soilordercon%r_weather(m)         = r_weather(m)
       soilordercon%r_adsorp(m)         = r_adsorp(m)
       soilordercon%r_desorp(m)         = r_desorp(m)
       soilordercon%r_occlude(m)         = r_occlude(m)
       soilordercon%k_s1_biochem(m)         = k_s1_biochem(m)
       soilordercon%k_s2_biochem(m)         = k_s2_biochem(m)
       soilordercon%k_s3_biochem(m)         = k_s3_biochem(m)
       soilordercon%k_s4_biochem(m)         = k_s4_biochem(m)


    end do
  end subroutine soilorderconInit

end module SoilorderConType
                               
