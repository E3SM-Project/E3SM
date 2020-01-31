!--------------------------------------------------------------------------------
! Manages the absorber concentrations in the layers RRTMG operates 
! including an extra layer over the model if needed.
!
! Creator: Francis Vitt 
! 9 May 2011
!--------------------------------------------------------------------------------
module rrtmg_state

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, pverp

  implicit none
  private
  save
  
  public :: rrtmg_state_t
  public :: rrtmg_state_init
  public :: rrtmg_state_create
  public :: rrtmg_state_update
  public :: rrtmg_state_destroy
  public :: num_rrtmg_levs

  type rrtmg_state_t

     real(r8), allocatable :: h2ovmr(:,:)   ! h2o volume mixing ratio
     real(r8), allocatable :: o3vmr(:,:)    ! o3 volume mixing ratio
     real(r8), allocatable :: co2vmr(:,:)   ! co2 volume mixing ratio 
     real(r8), allocatable :: ch4vmr(:,:)   ! ch4 volume mixing ratio 
     real(r8), allocatable :: o2vmr(:,:)    ! o2  volume mixing ratio 
     real(r8), allocatable :: n2ovmr(:,:)   ! n2o volume mixing ratio 
     real(r8), allocatable :: cfc11vmr(:,:) ! cfc11 volume mixing ratio
     real(r8), allocatable :: cfc12vmr(:,:) ! cfc12 volume mixing ratio
     real(r8), allocatable :: cfc22vmr(:,:) ! cfc22 volume mixing ratio
     real(r8), allocatable :: ccl4vmr(:,:)  ! ccl4 volume mixing ratio

     real(r8), allocatable :: pmidmb(:,:)   ! Level pressure (hPa)
     real(r8), allocatable :: pintmb(:,:)   ! Model interface pressure (hPa)
     real(r8), allocatable :: tlay(:,:)     ! mid point temperature
     real(r8), allocatable :: tlev(:,:)     ! interface temperature

  end type rrtmg_state_t

  integer :: num_rrtmg_levs ! number of pressure levels greate than 1.e-4_r8 mbar

  real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
  real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
  real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
  real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
  real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
  real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
  real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
  real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

contains

!--------------------------------------------------------------------------------
! sets the number of model levels RRTMG operates
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_init

    use ref_pres,only : pref_edge
    implicit none

    ! The following cuts off RRTMG at roughly the point where it becomes
    ! invalid due to low pressure.
    num_rrtmg_levs = count( pref_edge(:) > 1._r8 ) ! pascals (1.e-2 mbar)

  end subroutine rrtmg_state_init
  
!--------------------------------------------------------------------------------
! creates (alloacates) an rrtmg_state object
!--------------------------------------------------------------------------------

  function rrtmg_state_create( pstate, cam_in ) result( rstate )
    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t
    use physconst,        only: stebol
#ifdef MAML
    use seq_comm_mct,       only : num_inst_atm
#endif
    implicit none

    type(physics_state), intent(in) :: pstate
    type(cam_in_t),      intent(in) :: cam_in

    type(rrtmg_state_t), pointer  :: rstate

    real(r8) dy                   ! Temporary layer pressure thickness
    real(r8) :: tint(pcols,pverp)    ! Model interface temperature
    integer  :: ncol, i, kk, k
#ifdef MAML
    real(r8) :: lwupavg_in(pcols)
    real(r8) :: factor_xy
    integer :: ii
#endif
    allocate( rstate )

    allocate( rstate%h2ovmr(pcols,num_rrtmg_levs) )
    allocate( rstate%o3vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%co2vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%ch4vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%o2vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%n2ovmr(pcols,num_rrtmg_levs) )
    allocate( rstate%cfc11vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%cfc12vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%cfc22vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%ccl4vmr(pcols,num_rrtmg_levs) )

    allocate( rstate%pmidmb(pcols,num_rrtmg_levs) )
    allocate( rstate%pintmb(pcols,num_rrtmg_levs+1) )
    allocate( rstate%tlay(pcols,num_rrtmg_levs) )
    allocate( rstate%tlev(pcols,num_rrtmg_levs+1) )

    ncol = pstate%ncol

    ! Calculate interface temperatures (following method
    ! used in radtpl for the longwave), using surface upward flux and
    ! stebol constant in mks units
    do i = 1,ncol
       tint(i,1) = pstate%t(i,1)
#ifdef MAML
       lwupavg_in(i) =0._r8
       factor_xy = 1._r8 / dble(num_inst_atm)
       do ii=1,num_inst_atm
          lwupavg_in(i) = lwupavg_in(i)+cam_in%lwup(i,ii)*factor_xy
       enddo
       tint(i,pverp) = sqrt(sqrt(lwupavg_in(i)/stebol))
#else
       tint(i,pverp) = sqrt(sqrt(cam_in%lwup(i)/stebol))
#endif
       do k = 2,pver
          dy = (pstate%lnpint(i,k) - pstate%lnpmid(i,k)) / (pstate%lnpmid(i,k-1) - pstate%lnpmid(i,k))
          tint(i,k) = pstate%t(i,k) - dy * (pstate%t(i,k) - pstate%t(i,k-1))
       end do
    end do

    do k = 1, num_rrtmg_levs

       kk = max(k + (pverp-num_rrtmg_levs)-1,1)

       rstate%pmidmb(:ncol,k) = pstate%pmid(:ncol,kk) * 1.e-2_r8
       rstate%pintmb(:ncol,k) = pstate%pint(:ncol,kk) * 1.e-2_r8

       rstate%tlay(:ncol,k) = pstate%t(:ncol,kk)
       rstate%tlev(:ncol,k) = tint(:ncol,kk)

    enddo

    ! bottom interface
    rstate%pintmb(:ncol,num_rrtmg_levs+1) = pstate%pint(:ncol,pverp) * 1.e-2_r8 ! mbar
    rstate%tlev(:ncol,num_rrtmg_levs+1) = tint(:ncol,pverp)

    ! top layer thickness
    if (num_rrtmg_levs==pverp) then
       rstate%pmidmb(:ncol,1) = 0.5_r8 * rstate%pintmb(:ncol,2) 
       rstate%pintmb(:ncol,1) = 1.e-4_r8 ! mbar
    endif

  endfunction rrtmg_state_create

!--------------------------------------------------------------------------------
! updates the concentration fields
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_update(pstate,pbuf,icall,rstate)
    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc
    use rad_constituents, only: rad_cnst_get_gas

    implicit none

    type(physics_state), intent(in), target :: pstate
    type(physics_buffer_desc),  pointer :: pbuf(:)
    integer,             intent(in) :: icall                     ! index through climate/diagnostic radiation calls
    type(rrtmg_state_t), pointer    :: rstate

    real(r8), pointer, dimension(:,:) :: sp_hum ! specific humidity
    real(r8), pointer, dimension(:,:) :: n2o    ! nitrous oxide mass mixing ratio
    real(r8), pointer, dimension(:,:) :: ch4    ! methane mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o2     ! O2 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc11  ! cfc11 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc12  ! cfc12 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o3     ! Ozone mass mixing ratio
    real(r8), pointer, dimension(:,:) :: co2    ! co2   mass mixing ratio
    
    integer  :: ncol, i, kk, k

    ncol = pstate%ncol

    ! Get specific humidity
    call rad_cnst_get_gas(icall,'H2O', pstate, pbuf, sp_hum)
    ! Get oxygen mass mixing ratio.
    call rad_cnst_get_gas(icall,'O2',  pstate, pbuf, o2)
    ! Get ozone mass mixing ratio.
    call rad_cnst_get_gas(icall,'O3',  pstate, pbuf, o3)
    ! Get CO2 mass mixing ratio
    call rad_cnst_get_gas(icall,'CO2', pstate, pbuf, co2)
    ! Get N2O mass mixing ratio
    call rad_cnst_get_gas(icall,'N2O', pstate, pbuf, n2o)
    ! Get CH4 mass mixing ratio
    call rad_cnst_get_gas(icall,'CH4', pstate, pbuf, ch4)
    ! Get CFC mass mixing ratios
    call rad_cnst_get_gas(icall,'CFC11', pstate, pbuf, cfc11)
    call rad_cnst_get_gas(icall,'CFC12', pstate, pbuf, cfc12)

    do k = 1, num_rrtmg_levs

       kk = max(k + (pverp-num_rrtmg_levs)-1,1)

       rstate%ch4vmr(:ncol,k)   = ch4(:ncol,kk) * amdm
       rstate%h2ovmr(:ncol,k)   = (sp_hum(:ncol,kk) / (1._r8 - sp_hum(:ncol,kk))) * amdw
       rstate%o3vmr(:ncol,k)    = o3(:ncol,kk) * amdo
       rstate%co2vmr(:ncol,k)   = co2(:ncol,kk) * amdc
       rstate%ch4vmr(:ncol,k)   = ch4(:ncol,kk) * amdm
       rstate%o2vmr(:ncol,k)    = o2(:ncol,kk) * amdo2
       rstate%n2ovmr(:ncol,k)   = n2o(:ncol,kk) * amdn
       rstate%cfc11vmr(:ncol,k) = cfc11(:ncol,kk) * amdc1
       rstate%cfc12vmr(:ncol,k) = cfc12(:ncol,kk) * amdc2
       rstate%cfc22vmr(:ncol,k) = 0._r8
       rstate%ccl4vmr(:ncol,k)  = 0._r8

    enddo

  end subroutine rrtmg_state_update

!--------------------------------------------------------------------------------
! de-allocates an rrtmg_state object
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_destroy(rstate)

    implicit none

    type(rrtmg_state_t), pointer   :: rstate

    deallocate(rstate%h2ovmr)
    deallocate(rstate%o3vmr)
    deallocate(rstate%co2vmr)
    deallocate(rstate%ch4vmr)
    deallocate(rstate%o2vmr)
    deallocate(rstate%n2ovmr)
    deallocate(rstate%cfc11vmr)
    deallocate(rstate%cfc12vmr)
    deallocate(rstate%cfc22vmr)
    deallocate(rstate%ccl4vmr)

    deallocate(rstate%pmidmb)
    deallocate(rstate%pintmb)
    deallocate(rstate%tlay)
    deallocate(rstate%tlev)

    deallocate( rstate )
    nullify(rstate)

  endsubroutine rrtmg_state_destroy

end module rrtmg_state
