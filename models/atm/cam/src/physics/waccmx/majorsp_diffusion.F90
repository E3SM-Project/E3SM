
module majorsp_diffusion

!--------------------------------------------------------------------------
! This module computes the diffusion of major species (O2 and O) mass mixing
! ratio. This routine computes both the molecular and eddy diffusivity. This
! is adapted from the major species diffusion calculation of TIME-GCM.
! 
! Calling sequence:
!   initialization:
!      init
!         call mspd_init
!
!   interfacing:
!      tphysac
!         (after vertical_diffusion_tend)
!         call mspd_intr
!            call mspdiff
!
!---------------------------Code history--------------------------------
! Adapted from TIME-GCM (comp.F): H.-L. Liu, Nov 2003
!--------------------------------------------------------------------------


  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver, pverp
  use constituents, only: pcnst, cnst_name, cnst_get_ind, cnst_mw
  use spmd_utils,   only: masterproc

  implicit none

  private          ! Make default type private to the module
  save
!-----------------------
! Public interfaces
!-----------------------
  public mspd_init   ! Initialization
  public mspd_intr   ! Full routine
!-----------------------
! Private data
!-----------------------

  real(r8) :: rmass_o2, rmass_o1, rmass_n2               ! molecular weight kg/kmol
  real(r8) :: rmassinv_o2, rmassinv_o1, rmassinv_n2      ! 1/rmass_o2...
  real(r8) :: phi(2,3)                                   ! mutual diffusion constants of
                                                         ! major constituents
  real(r8) :: delta(2,2)                                 ! unit matrix

  real(r8), parameter :: t00=273._r8                     ! reference temperature
  real(r8), parameter :: ptref=5.e-5_r8                  ! thermosphere reference pressure (Pa)
  real(r8), parameter :: tau=1.86e3_r8                   ! diffusive time constant (sec).
  real(r8), parameter :: protonmass=1.6726e-27_r8        ! Proton mass (kg)
  real(r8), parameter :: mmrMin=1.e-20_r8                ! lower limit of o2 and o mixing ratio
  real(r8), parameter :: N2mmrMin=1.e-6_r8               ! lower limit of o2, o, and h mixing ratios

  integer :: indx_O2                                     ! cnst index for o2
  integer :: indx_O                                      ! cnst index for o
  integer :: indx_H                                      ! cnst index for h
  integer, parameter :: io2=1, io1=2                     ! local indices to o2 , o respectively
  logical :: fixed_ubc(2)                                ! flag for fixed upper boundary condition

  real(r8) :: o2mmr_ubc(pcols)                           ! MMR of O2 at top boundary (specified)
  real(r8) :: ommr_ubc(pcols)                            ! MMR of O at top boundary
  real(r8) :: t_ubc(pcols)                               ! T at top boundary


  character(len=8), private :: mjdiffnam(2)              ! names of v-diff tendencies

contains

!===============================================================================
  subroutine mspd_init()

    !-------------------------------------------------------------------------------
    ! Define constants and coeficient matrices, phi and delta, in the initialization.
    !-------------------------------------------------------------------------------
    use constituents, only: cnst_mw, cnst_fixed_ubc
    use cam_history,  only: addfld, add_default
    use phys_control, only: phys_getopts

    !------------------------------Arguments--------------------------------

    !---------------------------Local storage-------------------------------
    integer :: k, m

    !-----------------------------------------------------------
    ! Get required molecular weights
    !-----------------------------------------------------------
    call cnst_get_ind('O2', indx_O2, abort=.false.)
    call cnst_get_ind('O',  indx_O, abort=.false.)
    call cnst_get_ind('H',  indx_H, abort=.false.)

    rmass_o2 = cnst_mw(indx_O2)
    rmass_o1 = cnst_mw(indx_O)
    rmass_n2 = 28._r8

    rmassinv_o2 = 1._r8/rmass_o2
    rmassinv_o1 = 1._r8/rmass_o1
    rmassinv_n2 = 1._r8/rmass_n2

    !--------------------------------------------------------------------
    ! Get fixed upper boundary flags and set vertical range for diffusion
    !--------------------------------------------------------------------
    fixed_ubc(io2) = cnst_fixed_ubc(indx_O2)
    fixed_ubc(io1) = cnst_fixed_ubc(indx_O)

    !------------------------------------------------
    ! Set diffusion constants and setup matrix
    !------------------------------------------------
    phi(:,1)=(/0._r8  ,0.673_r8/)
    phi(:,2)=(/1.35_r8,0._r8   /)
    phi(:,3)=(/1.11_r8,0.769_r8/)
    delta(:,1)=(/1._r8,0._r8/)
    delta(:,2)=(/0._r8,1._r8/)

   ! Set names of major diffusion tendencies and declare them as history variables
    mjdiffnam(1) = 'MD'//cnst_name(indx_O2)
    call addfld (mjdiffnam(1),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_O2))
    mjdiffnam(2) = 'MD'//cnst_name(indx_O)
    call addfld (mjdiffnam(2),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_O))
    call add_default (mjdiffnam(1), 1, ' ')
    call add_default (mjdiffnam(2), 1, ' ')

  return
  end subroutine mspd_init

!===============================================================================
  subroutine mspd_intr(ztodt    ,state    ,ptend)

!-------------------------------------------------------------------------------
! interface routine. output tendency.
!-------------------------------------------------------------------------------
    use physics_types,  only: physics_state, physics_ptend
    use cam_history,    only: outfld
    use upper_bc,       only: ubc_get_vals
    use physconst,      only: rairv, mbarv
    use ref_pres,       only: ntop_molec

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    type(physics_state), intent(in)     :: state   ! Physics state variables
    type(physics_ptend), intent(inout)  :: ptend   ! indivdual parameterization tendencies
!---------------------------Local storage-------------------------------
    real(r8) :: rztodt                             ! 1/ztodt
    real(r8) :: tendo2o(pcols,pver,2)              ! temporary array for o2 and o tendency
    real(r8) :: ubc_mmr(pcols,pcnst)               ! upper bndy mixing ratios (kg/kg)
    real(r8) :: ubc_t(pcols)                       ! upper bndy temperature (K)
    real(r8) :: ubc_flux(pcnst)                    ! upper bndy flux (kg/s/m^2)
    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
    integer :: i, k, m                             ! indexing integers

    !--------------------------------------------------------------------------------------------
    ! local constants
    !--------------------------------------------------------------------------------------------
    rztodt = 1._r8/ztodt
    lchnk = state%lchnk
    ncol  = state%ncol

    !----------------------------------------------------------------------------------------------
    ! Store the o2 and o tendencies calculated from vertical_diffusion (due to eddy diffusion only)
    !----------------------------------------------------------------------------------------------
    tendo2o(:ncol,:,io2) = ptend%q(:ncol,:,indx_O2)
    tendo2o(:ncol,:,io1) = ptend%q(:ncol,:,indx_O)
    
    !----------------------------------------------------------------------
    ! Operate on copies of the input states, convert to tendencies at end.
    !----------------------------------------------------------------------
    ptend%q(:ncol,:,indx_O2) = state%q(:ncol,:,indx_O2)
    ptend%q(:ncol,:,indx_O) = state%q(:ncol,:,indx_O)

    if (fixed_ubc(io2) .or. fixed_ubc(io1)) then
       !-------------------------------------------
       ! set upper boundary values of O2 and O MMR.
       !-------------------------------------------
       call ubc_get_vals (lchnk, ncol, ntop_molec, state%pint, state%zi, ubc_t, ubc_mmr, ubc_flux)
       o2mmr_ubc(:ncol) = ubc_mmr(:ncol,indx_O2)
       ommr_ubc(:ncol) = ubc_mmr(:ncol,indx_O)
       t_ubc(:ncol) = ubc_t(:ncol)
    endif

    ! Since this is a combined tendency, retain the old name for output
    ! and debugging purposes.
    ptend%name  = trim(ptend%name)//"+mspd"
    ptend%lq(indx_O2) = .TRUE.
    ptend%lq(indx_O) = .TRUE.
    !---------------------------------------------
    ! Call the major species diffusion subroutine.
    !---------------------------------------------
    call mspdiff (lchnk      ,ncol       ,                                     &
                  state%t    ,ptend%q    ,state%pmid ,state%pint ,             &
                  state%pdel ,state%rpdel,ztodt      ,rairv(:,:,lchnk),  mbarv(:,:,lchnk))

    !---------------------------------------------
    ! Update O2 and O tendencies and output
    !---------------------------------------------
    do k=1,pver
       do i=1,ncol
          ptend%q(i,k,indx_O2) = (ptend%q(i,k,indx_O2)-state%q(i,k,indx_O2))*rztodt  &
                                 +tendo2o(i,k,io2)
          ptend%q(i,k,indx_O) = (ptend%q(i,k,indx_O)-state%q(i,k,indx_O))*rztodt     &
                                 +tendo2o(i,k,io1)
       enddo
    enddo

    call outfld(mjdiffnam(1),ptend%q(1,1,indx_O2),pcols,lchnk)
    call outfld(mjdiffnam(2),ptend%q(1,1,indx_O),pcols,lchnk)

  end subroutine mspd_intr


!===============================================================================
  subroutine mspdiff (lchnk      ,ncol       ,                                     &
                      t          ,q          ,pmid       ,pint       ,             &
                      pdel       ,rpdel      ,ztodt      ,rairv      ,mbarv)
!-----------------------------------------------------------------------
! Driver routine to compute major species diffusion (O2 and O).

! Turbulent diffusivities and boundary layer nonlocal transport terms are 
! obtained from the turbulence module.
!---------------------------Arguments------------------------------------
    use ref_pres,     only: lev1 => ntop_molec, lev0 => nbot_molec
    use physconst,    only: gravit

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    real(r8), intent(in) :: t(pcols,pver)          ! temperature input
    real(r8), intent(in) :: pmid(pcols,pver)       ! midpoint pressures
    real(r8), intent(in) :: pint(pcols,pverp)      ! interface pressures
    real(r8), intent(in) :: pdel(pcols,pver)       ! thickness between interfaces
    real(r8), intent(in) :: rpdel(pcols,pver)      ! 1./pdel
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    real(r8), intent(in) :: rairv(pcols,pver)                  ! composition dependent gas "constant"
    real(r8), intent(in) :: mbarv(pcols,pver)                  ! composition dependent mean mass

    real(r8), intent(inout) :: q(pcols,pver,pcnst) ! constituents 

!---------------------------Local storage-------------------------------
    real(r8) :: o2(pcols,pver), o1(pcols,pver)     ! o2, o1 mixing ratio (kg/kg moist air)
    real(r8) :: h_atom(pcols,pver)                 ! H mixing ratio
    real(r8) :: rztodt                             ! 1/ztodt
    real(r8) :: dz(pcols,pver)                     ! log-pressure interval between interfaces
    real(r8) :: rdz(pcols,pver)                    ! 1./dz, defined on midpoints
    real(r8) :: dzmid(pcols,pverp)                 ! log-pressure interval between midpoints
    real(r8) :: rdzmid(pcols,pverp)                ! 1./dzmid, defined on interfaces
    real(r8) :: ak(pcols,2,2,2)                    ! coefficient matrix "Alfa"
    real(r8) :: ep(pcols,2,2)                      ! coefficient matrix
    real(r8) :: difk(pcols,pverp)                  ! eddy diffusion normalized by scale height (1/sec)
    real(r8) :: expzm(pcols,pver)                  ! exp(-z)=pmid/ptref
    real(r8) :: expzi(pcols,pverp)                 ! exp(-z)=pint/ptref
    real(r8) :: wks1(pcols),wks2(pcols)
    real(r8) :: wks3(pcols),wks4(pcols)            ! temporary working arrays
    real(r8) :: psclht(pcols,pverp)                ! pressure scale height
    real(r8) :: p_ubc(pcols)                       ! extrapolated pressure at upper boundary level, pmid(1)^2=pmid(2)*p_ubc
    real(r8) :: rair_ubc(pcols)                    ! extrapolated rair at upper boundary level
    real(r8) :: mbar_ubc(pcols)                    ! extrapolated mbar at upper boundary level
    real(r8) :: flb(pcols,2)                       ! lower boundary condition for o2 and o, now
                                                   ! calculated locally.
    real(r8) :: fub(pcols,2)                       ! upper boundary condition for o2 and o, now
                                                   ! calculated locally.
    real(r8) :: fk(pcols,2)                        ! temporary working array for rhs
    real(r8) :: pk(pcols,2,2)                      ! temp array for coefficients on lower diagonal
    real(r8) :: rk(pcols,2,2)                      ! temp array for coefficients on upper diagonal
    real(r8) :: qk(pcols,2,2)                      ! temp array for coefficients on diagonal
    real(r8) :: apk(2,2,pcols,pver)                ! coefficients on lower diagonal
    real(r8) :: ark(2,2,pcols,pver)                ! coefficients on upper diagonal
    real(r8) :: aqk(2,2,pcols,pver)                ! coefficients on diagonal
    real(r8) :: rfk(2,pcols,pver)                  ! rhs of the array equation
    real(r8) :: betawk(2,2,pcols,pver)             ! working arrays for blktri solver.
    real(r8) :: gammawk(2,2,pcols,pver)            ! working arrays for blktri solver.
    real(r8) :: ywk(2,pcols,pver)                  ! working arrays for blktri solver.
    real(r8) :: xwk(2,pcols,pver)                  ! working arrays for blktri solver.
    integer  :: lev1p1                             ! ntop_molec+1
    integer  :: nlevs                              ! lev0-lev1p1+1
    integer  :: i, k, km, kp, m, ktmp, isp, kk, kr

    !---------------------------------------------------
    ! Set vertical grid and get time step for diffusion
    !---------------------------------------------------
    lev1p1 = lev1+1
    nlevs = lev0-lev1+1
    rztodt = 1._r8/ztodt

    !------------------------------------------------------
    ! Get species to diffuse and set upper/lower boundaries
    !------------------------------------------------------
    o2(:ncol,:) = q(:ncol,:,indx_O2)
    o1(:ncol,:) = q(:ncol,:,indx_O)
    h_atom(:ncol,:) = q(:ncol,:,indx_H)

    flb(:ncol,1) = o2(:ncol,lev0+1)      ! fixed lower boundary condition
    flb(:ncol,2) = o1(:ncol,lev0+1)
    if(fixed_ubc(io2).or.fixed_ubc(io1)) then
       fub(:ncol,1) = o2mmr_ubc(:ncol)      ! fixed upper boundary condition
       fub(:ncol,2) = ommr_ubc(:ncol)
    endif

    !------------------------------------------------------------------
    ! Get log-pressure intervals between midpoints and interface points
    !------------------------------------------------------------------
    dz(:ncol,:) = pdel(:ncol,:)/pmid(:ncol,:)
    rdz(:ncol,:) = 1._r8/dz(:ncol,:)
    do k=2,pver
       do i=1,ncol
          dzmid(i,k) = (pmid(i,k)-pmid(i,k-1))/pint(i,k)
          rdzmid(i,k) = 1._r8/dzmid(i,k)
       enddo
    enddo
    do i=1,ncol
       p_ubc(i) = pmid(i,1)*pmid(i,1)/pmid(i,2)
       dzmid(i,1) = (pmid(i,1)-p_ubc(i))/pint(i,1)
       rdzmid(i,1) = 1._r8/dzmid(i,1)
    enddo
    do i=1,ncol
       dzmid(i,pverp) = dzmid(i,pver)
       rdzmid(i,pverp) = 1._r8/dzmid(i,pverp)
    enddo

    !------------------------------------------------------------------
    ! Get log-pressure intervals between midpoints and interface points
    !------------------------------------------------------------------
    expzi(:ncol,:) = pint(:ncol,:)/ptref
    expzm(:ncol,:) = pmid(:ncol,:)/ptref

    !------------------------------------------------------------------
    ! Get pressure scale height
    !------------------------------------------------------------------
    do k=2,pver
       do i=1,ncol
          psclht(i,k) = .5_r8*(rairv(i,k)*t(i,k)+rairv(i,k-1)*t(i,k-1))/gravit
       enddo
    enddo
    do i=1,ncol
       rair_ubc(i) = 1.5_r8*rairv(i,1)-.5_r8*rairv(i,2)
       psclht(i,1) = .5_r8*(rairv(i,1)*t(i,1)+rair_ubc(i)*t_ubc(i))/gravit
       psclht(i,pverp) = psclht(i,pver)
    enddo

    !------------------------------------------------------------------
    ! Initialize scale height normalized eddy diffusion
    !------------------------------------------------------------------
    do k=1,pverp
       do i=1,ncol
          difk(i,k) = 0._r8          ! eddy diffusion already calculated in vertical_diffusion
       enddo
    enddo


    !------------------------------------------------------------------
    ! Set up mean mass working array
    !------------------------------------------------------------------
    ! ep, ak at the interface level immediately below midpoint level nbot_molec
    
    ! WKS4 = .5*(DMBAR/DZ)/MBAR
    do i=1,ncol
       wks4(i) = (mbarv(i,lev0)-mbarv(i,lev0+1))/                               &
                 (dzmid(i,lev0+1)*(mbarv(i,lev0+1)+mbarv(i,lev0)))
    enddo

    !-----------------------------------
    ! Calculate coefficient matrices
    !-----------------------------------
    km = 1 
    kp = 2
    do i=1, ncol
       ep(i,io2,kp) = 1._r8-(2._r8/(mbarv(i,lev0+1)+mbarv(i,lev0)))*                  &
                   (rmass_o2+(mbarv(i,lev0)-mbarv(i,lev0+1))*rdzmid(i,lev0+1))
       ep(i,io1,kp) = 1._r8-(2._r8/(mbarv(i,lev0+1)+mbarv(i,lev0)))*                  &
                   (rmass_o1+(mbarv(i,lev0)-mbarv(i,lev0+1))*rdzmid(i,lev0+1))
    enddo


    do m=1,2
      do i=1,ncol
         ak(i,io2,m,kp) =                                           &
            -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*    &
            .5_r8*(o2(i,lev0+1)+o2(i,lev0)))-(1._r8-delta(io2,m))*        &
            (phi(io2,m)-phi(io2,3))*.5_r8*(o2(i,lev0+1)+o2(i,lev0))
         ak(i,io1,m,kp) =                                           &
            -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*    &
            .5_r8*(o1(i,lev0+1)+o1(i,lev0)))-(1._r8-delta(io1,m))*        &
            (phi(io1,m)-phi(io1,3))*.5_r8*(o1(i,lev0+1)+o1(i,lev0))
      enddo
    enddo
!
! WKS1=MBAR/M3*(T00/(T0+T))*0.25/(TAU*DET(ak)) ak at the interface level 
! immediately below midpoint level nbot_molec. 
    do i=1,ncol
      wks1(i) = 0.5_r8*(mbarv(i,lev0+1)+mbarv(i,lev0))*rmassinv_n2*       &
        (2._r8*t00/(t(i,lev0+1)+t(i,lev0)))**0.25_r8/                      &
        (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*ak(i,2,1,kp)))
    enddo
!
! Complete calculation of ak at the interface level immediately below midpoint
! level nbot_molec.
    do m=1,2
      do i=1,ncol
        ak(i,io2,m,kp) = ak(i,io2,m,kp)*wks1(i)
        ak(i,io1,m,kp) = ak(i,io1,m,kp)*wks1(i)
      enddo
    enddo

    km = 1
    kp = 2
    do k=lev0,lev1p1,-1
       ktmp = km
       km = kp
       kp = ktmp
       do i=1,ncol
          ep(i,io2,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbarv(i,k-1)))*(rmass_o2+      &
                         (mbarv(i,k-1)-mbarv(i,k))*rdzmid(i,k))
          ep(i,io1,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbarv(i,k-1)))*(rmass_o1+      &
                         (mbarv(i,k-1)-mbarv(i,k))*rdzmid(i,k))
       enddo

       do m=1,2
          do i=1,ncol
             ak(i,io2,m,kp) =                                               &
                  -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*      &
                  .5_r8*(o2(i,k)+o2(i,k-1)))-                                  &
                  (1._r8-delta(io2,m))*(phi(io2,m)-phi(io2,3))*                &
                  .5_r8*(o2(i,k)+o2(i,k-1))

             ak(i,io1,m,kp) =                                               &
                  -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*      &
                  .5_r8*(o1(i,k)+o1(i,k-1)))-                                  &
                  (1._r8-delta(io1,m))*(phi(io1,m)-phi(io1,3))*                &
                  .5_r8*(o1(i,k)+o1(i,k-1))

          enddo
       enddo
       
    !---------------------------------------------
    ! Calculate coefficients for diagonals and rhs
    !---------------------------------------------
!
! WKS1=MBAR/M3*(T00/(T0+T))**0.25/(TAU*DET(ALFA))
       do i=1,ncol
          wks1(i) = 0.5_r8*(mbarv(i,k)+mbarv(i,k-1))*rmassinv_n2*              &
               (2._r8*t00/(t(i,k)+t(i,k-1)))**0.25_r8/                          &
               (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*ak(i,2,1,kp)))
          wks3(i) = wks4(i)
          wks4(i) = (mbarv(i,k-1)-mbarv(i,k))/                              &
               (dzmid(i,k)*(mbarv(i,k)+mbarv(i,k-1)))
       enddo

!
! FINISH CALCULATING AK(K+1/2) AND GENERATE PK, QK, RK
       do m=1,2
          do isp=io2,io1
             do i=1,ncol
                ak(i,isp,m,kp) = ak(i,isp,m,kp)*wks1(i)

                pk(i,isp,m) = (ak(i,isp,m,km)*(rdzmid(i,k+1)+ep(i,m,km)/2._r8)-   &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)-                  &
                     wks3(i))*delta(isp,m))*rdz(i,k)
                
                rk(i,isp,m) = (ak(i,isp,m,kp)*(rdzmid(i,k)-ep(i,m,kp)/2._r8)-     &
                     expzi(i,k)*difk(i,k)*(rdzmid(i,k)+                        &
                     wks4(i))*delta(isp,m))*rdz(i,k)
   
                qk(i,isp,m) = -(ak(i,isp,m,km)*(rdzmid(i,k+1)-ep(i,m,km)/2._r8)+  &
                     ak(i,isp,m,kp)*(rdzmid(i,k)+ep(i,m,kp)/2._r8))*rdz(i,k)+     &
                     ((expzi(i,k)*difk(i,k)*(rdzmid(i,k)-wks4(i))+             &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)+wks3(i)))*        &
                     rdz(i,k)+expzm(i,k)*rztodt)*delta(isp,m)

             enddo
          enddo
       enddo

       do i=1,ncol
          fk(i,io2) = expzm(i,k)*o2(i,k)*rztodt
          fk(i,io1) = expzm(i,k)*o1(i,k)*rztodt
       enddo
 
       !----------------------------
       ! Lower boundary
       !----------------------------
       if (k==lev0) then
          do m=1,2
            do i=1,ncol
              fk(i,io2) = fk(i,io2)-pk(i,io2,m)*flb(i,m)
              fk(i,io1) = fk(i,io1)-pk(i,io1,m)*flb(i,m)
              pk(i,:,m) = 0._r8
            enddo
          enddo
       endif

       kr = lev0-k+1

       do i=1,ncol
          do m=1,2
             do kk=1,2
                apk(kk,m,i,kr) = pk(i,kk,m)
                aqk(kk,m,i,kr) = qk(i,kk,m)
                ark(kk,m,i,kr) = rk(i,kk,m)
             enddo
          enddo
       enddo

       do i=1,ncol
          rfk(io2,i,kr) = fk(i,io2)
          rfk(io1,i,kr) = fk(i,io1)
       enddo

    enddo

    !----------------------------
    ! Upper boundary
    !----------------------------
    k=lev1
    ktmp = km
    km = kp
    kp = ktmp
    if(fixed_ubc(io2).or.fixed_ubc(io1)) then
       do i=1,ncol
          mbar_ubc(i) = 1._r8/(o2mmr_ubc(i)*rmassinv_o2+ommr_ubc(i)*rmassinv_o1+ &
               (1._r8-o2mmr_ubc(i)-ommr_ubc(i))*rmassinv_n2)
          ep(i,io2,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbar_ubc(i)))*(rmass_o2+      &
               (mbar_ubc(i)-mbarv(i,k))*rdzmid(i,k))
          ep(i,io1,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbar_ubc(i)))*(rmass_o1+      &
               (mbar_ubc(i)-mbarv(i,k))*rdzmid(i,k))
       enddo

       do m=1,2
          do i=1,ncol
             ak(i,io2,m,kp) =                                               &
                  -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*      &
                  .5_r8*(o2(i,k)+o2mmr_ubc(i)))-                                  &
                  (1._r8-delta(io2,m))*(phi(io2,m)-phi(io2,3))*                &
                  .5_r8*(o2(i,k)+o2mmr_ubc(i))

             ak(i,io1,m,kp) =                                               &
                  -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*      &
                  .5_r8*(o1(i,k)+ommr_ubc(i)))-                                  &
                  (1._r8-delta(io1,m))*(phi(io1,m)-phi(io1,3))*                &
                  .5_r8*(o1(i,k)+ommr_ubc(i))
             
          enddo
       enddo
       
!
! WKS1=MBAR/M3*(T00/(T0+T))**0.25/(TAU*DET(ALFA))
       do i=1,ncol
          wks1(i) = 0.5_r8*(mbarv(i,k)+mbar_ubc(i))*rmassinv_n2*              &
               (2._r8*t00/(t(i,k)+t_ubc(i)))**0.25_r8/                          &
               (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*ak(i,2,1,kp)))
          wks3(i) = wks4(i)
          wks4(i) = (mbar_ubc(i)-mbarv(i,k))/                              &
               (dzmid(i,k)*(mbarv(i,k)+mbar_ubc(i)))
       enddo

!
! FINISH CALCULATING AK(K+1/2) AND GENERATE PK, QK, RK
       do m=1,2
          do isp=io2,io1
             do i=1,ncol
                ak(i,isp,m,kp) = ak(i,isp,m,kp)*wks1(i)
                
                pk(i,isp,m) = (ak(i,isp,m,km)*(rdzmid(i,k+1)+ep(i,m,km)/2._r8)-   &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)-                  &
                     wks3(i))*delta(isp,m))*rdz(i,k)

                rk(i,isp,m) = (ak(i,isp,m,kp)*(rdzmid(i,k)-ep(i,m,kp)/2._r8)-     &
                     expzi(i,k)*difk(i,k)*(rdzmid(i,k)+                        &
                     wks4(i))*delta(isp,m))*rdz(i,k)

                qk(i,isp,m) = -(ak(i,isp,m,km)*(rdzmid(i,k+1)-ep(i,m,km)/2._r8)+  &
                     ak(i,isp,m,kp)*(rdzmid(i,k)+ep(i,m,kp)/2._r8))*rdz(i,k)+     &
                     ((expzi(i,k)*difk(i,k)*(rdzmid(i,k)-wks4(i))+             &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)+wks3(i)))*        &
                     rdz(i,k)+expzm(i,k)*rztodt)*delta(isp,m)

             enddo
          enddo
       enddo

       do i=1,ncol
          fk(i,io2) = expzm(i,k)*o2(i,k)*rztodt
          fk(i,io1) = expzm(i,k)*o1(i,k)*rztodt
       enddo
       do m=1,2
          do i=1,ncol
             fk(i,io2) = fk(i,io2)-rk(i,io2,m)*fub(i,m)
             fk(i,io1) = fk(i,io1)-rk(i,io1,m)*fub(i,m)
             rk(i,:,m) = 0._r8
          enddo
       enddo

    else
       
       do i=1,ncol
          wks3(i) = wks4(i)
       enddo
       do m=1,2
          do isp=io2,io1
             do i=1,ncol
                pk(i,isp,m) = (ak(i,isp,m,km)*(rdzmid(i,k+1)+ep(i,m,km)/2._r8)-   &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)-                  &
                     wks3(i))*delta(isp,m))*rdz(i,k)

                qk(i,isp,m) = -(ak(i,isp,m,km)*(rdzmid(i,k+1)-ep(i,m,km)/2._r8))  &
                     *rdz(i,k)+     &
                     (expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)+wks3(i))*        &
                     rdz(i,k)+expzm(i,k)*rztodt)*delta(isp,m)

             enddo
          enddo
       enddo

       do i=1,ncol
          fk(i,io2) = expzm(i,k)*o2(i,k)*rztodt
          fk(i,io1) = expzm(i,k)*o1(i,k)*rztodt
       enddo

    endif

    kr = lev0-k+1

    do i=1,ncol
       do m=1,2
          do kk=1,2
             apk(kk,m,i,kr) = pk(i,kk,m)
             aqk(kk,m,i,kr) = qk(i,kk,m)
             ark(kk,m,i,kr) = rk(i,kk,m)
          enddo
       enddo
    enddo

    do i=1,ncol
       rfk(io2,i,kr) = fk(i,io2)
       rfk(io1,i,kr) = fk(i,io1)
    enddo

    !------------------------------------
    ! Call solver to get diffused species
    !------------------------------------
    call blktri(apk,aqk,ark,rfk,pcols,1,ncol,pver,1,nlevs,    &
                betawk, gammawk, ywk, xwk)

    do k=lev0,lev1,-1
       kr = lev0-k+1
       do i=1,ncol
          o2(i,k) = xwk(1,i,kr)
          o1(i,k) = xwk(2,i,kr)
       enddo
    enddo

    !---------------------------------------------------------------
    ! Ensure non-negative O2 and O and check for N2 greater than one
    !---------------------------------------------------------------
    do i=1,ncol
       do k=lev1,lev0

          if (o2(i,k) < mmrMin) o2(i,k) = mmrMin
          if (o1(i,k) < mmrMin) o1(i,k) = mmrMin

          if(1._r8-mmrMin-o2(i,k)-o1(i,k)-h_atom(i,k) < 0._r8) then
             o2(i,k) = o2(i,k)*((1._r8-N2mmrMin-h_atom(i,k))/(o2(i,k)+o1(i,k)))
             o1(i,k) = o1(i,k)*((1._r8-N2mmrMin-h_atom(i,k))/(o2(i,k)+o1(i,k)))
          endif
       enddo
    enddo

    q(:ncol,:,indx_O2) = o2(:ncol,:)
    q(:ncol,:,indx_O)  = o1(:ncol,:)

  end subroutine mspdiff


!===============================================================================
      SUBROUTINE BLKTRI(A,B,C,F,IF,I1,I2,KF,K1,K2,BETA,GAMMA,Y,X)
      implicit none
!     ****
!     ****     This procedure solves (I2-I1+1) tridiagonal block matrix
!     ****     systems in which all blocks are 2 x 2 matrices.
!     ****
!     ****     Each system may be written:
!     ****
!     ****      A(K) * X(K-1) + B(K) * X(K) + Z(K) * X(K+1) = F(K)
!     ****
!     ****      where:
!     ****
!     ****       K = K1,K2,1
!     ****
!     ****       A(K), B(K), C(K) are given (2 x 2) matrices.
!     ****
!     ****       The F(k) are given two componente vectors.
!     ****
!     ****       The system is to be solved for the two component
!     ****       vectors, X(K).
!     ****
!     ****       A(K1) = C(K2) = 0.
!     ****
!     ****      BETA(K), GAMMA(K), (K = K1,K2,1), are work space for
!     ****      (2 x 2) matrices.
!     ****
!     ****      Y(K), (K = K1,K2,1), is work space for two component
!     ****      vectors.
!     ****
!     ****     Algorithm: (See Isaacson and Keller p55)
!     ****
!     ****      Forward sweep from K = K1 to K = K2:
!     ****
!     ****       BETA(K1) = B(K1)**(-1)
!     ****
!     ****       Y(K1) = BETA(K1)*F(K1)
!     ****
!     ****       GAMMA(K) = BETA(K)*C(K),        K = K1,(K2-1),1
!     ****
!     ****       BETA(K) = (B(K) - A(K)*GAMMA(K-1))**(-1),
!     ****                                       K = K1+1,K2,1
!     ****
!     ****       Y(K) = BETA(K)*(F(K) - A(K)*Y(K-1)),
!     ****                                       K = K1+1,K2,1
!     ****
!     ****      Backward sweep, K = K2,K1,-1
!     ****
!     ****       X(K2) = Y(K2)
!     ****
!     ****       X(K) = Y(K) - GAMMA(K)*X(K+1),  K = K2-1,K1,-1
!     ****
!     ****     Dimension statements:
!     ****
!     ****      Block matrices are dimensioned thus:
!     ****
!     ****       MATRIX(2,2,IF,KF)
!     ****
!     ****      Two component vectors are similarly treated:
!     ****
!     ****       VECTOR(2,IF,KF)
!     ****
!     ****     Our block matrix scheme spans the range, (K = K1,K2,1),
!     ****     where (1 .LE. K1 .LT. K2 .LE. KF)
!     ****
!     ****     Similarly, we are solving (I1-I2+1) systems
!     ****     simultaneously as the index, I, spans the range,
!     ****     (I = I1,I2,1), where (1 .LE. I1 .LT. I2 .LE. IF)
!     ****
!     ****
!     ****     Dimension statements:
!     SUBROUTINE BLKTRI(A,B,C,F,IF,I1,I2,KF,K1,K2,BETA,GAMMA,Y,X)
!     ****
! Args:
      integer,intent(in)   :: if,i1,i2,kf,k1,k2
      real(r8),intent(in)  :: a(2,2,if,kf), b(2,2,if,kf), c(2,2,if,kf),      &
                              f(2,if,kf)
      real(r8),intent(out) :: beta(2,2,if,kf), gamma(2,2,if,kf),             &
                              y(2,if,kf), x(2,if,kf)
!
!     DIMENSION A(2,2,IF,KF), B(2,2,IF,KF), C(2,2,IF,KF), F(2,IF,KF),
!    1  BETA(2,2,IF,KF), GAMMA(2,2,IF,KF), Y(2,IF,KF), X(2,IF,KF)
!
! Local:
      integer :: i,k
!     ****
!     ****     Lower boundary at K = K1
!     ****
      DO I = I1,I2
!     ****
!     ****     Y(1,I,K1) = determinant(B(K1))
!     ****
        Y(1,I,K1) = B(1,1,I,K1)*B(2,2,I,K1) - B(1,2,I,K1)*B(2,1,I,K1)
!     ****
!     ****     BETA(K1) = B(K1)**(-1)
!     ****
        BETA(1,1,I,K1) = B(2,2,I,K1)/Y(1,I,K1)
        BETA(1,2,I,K1) = -B(1,2,I,K1)/Y(1,I,K1)
        BETA(2,1,I,K1) = -B(2,1,I,K1)/Y(1,I,K1)
        BETA(2,2,I,K1) = B(1,1,I,K1)/Y(1,I,K1)
!     ****
!     ****     Y(K1) = BETA(K1)*F(K1)
!     ****
        Y(1,I,K1) = BETA(1,1,I,K1)*F(1,I,K1) + BETA(1,2,I,K1)*F(2,I,K1)
        Y(2,I,K1) = BETA(2,1,I,K1)*F(1,I,K1) + BETA(2,2,I,K1)*F(2,I,K1)
      ENDDO
!     ****
!     ****     Now deal with levels (K1+1),K2,1
!     ****
      DO K = K1+1,K2
        DO I = I1,I2
!         ****
!         ****     GAMMA(K-1) = BETA(K-1)*C(K-1)
!         ****
          GAMMA(1,1,I,K-1) = BETA(1,1,I,K-1)*C(1,1,I,K-1) +              &
            BETA(1,2,I,K-1)*C(2,1,I,K-1)
          GAMMA(1,2,I,K-1) = BETA(1,1,I,K-1)*C(1,2,I,K-1) +              &
            BETA(1,2,I,K-1)*C(2,2,I,K-1)
          GAMMA(2,1,I,K-1) = BETA(2,1,I,K-1)*C(1,1,I,K-1) +              &
            BETA(2,2,I,K-1)*C(2,1,I,K-1)
          GAMMA(2,2,I,K-1) = BETA(2,1,I,K-1)*C(1,2,I,K-1) +              &
            BETA(2,2,I,K-1)*C(2,2,I,K-1)
!         ****
!         ****     GAMMA(K) = B(K) - A(K)*GAMMA(K-1)
!         ****
          GAMMA(1,1,I,K) = B(1,1,I,K) - A(1,1,I,K)*GAMMA(1,1,I,K-1) -    &
            A(1,2,I,K)*GAMMA(2,1,I,K-1)
          GAMMA(1,2,I,K) = B(1,2,I,K) - A(1,1,I,K)*GAMMA(1,2,I,K-1) -    &
            A(1,2,I,K)*GAMMA(2,2,I,K-1)
          GAMMA(2,1,I,K) = B(2,1,I,K) - A(2,1,I,K)*GAMMA(1,1,I,K-1) -    &
            A(2,2,I,K)*GAMMA(2,1,I,K-1)
          GAMMA(2,2,I,K) = B(2,2,I,K) - A(2,1,I,K)*GAMMA(1,2,I,K-1) -    &
            A(2,2,I,K)*GAMMA(2,2,I,K-1)
!         ****
!         ****     Y(1,I,K) = determinant(GAMMA(K))
!         ****
          Y(1,I,K) = GAMMA(1,1,I,K)*GAMMA(2,2,I,K) -                     &
            GAMMA(1,2,I,K)*GAMMA(2,1,I,K)
!         ****
!         ****     BETA(K) = GAMMA(K)**(-1)
!         ****
          BETA(1,1,I,K) = GAMMA(2,2,I,K)/Y(1,I,K)
          BETA(1,2,I,K) = -GAMMA(1,2,I,K)/Y(1,I,K)
          BETA(2,1,I,K) = -GAMMA(2,1,I,K)/Y(1,I,K)
          BETA(2,2,I,K) = GAMMA(1,1,I,K)/Y(1,I,K)
!         ****
!         ****     X(K) = F(K) - A(K)*Y(K-1)
!         ****
          X(1,I,K) = F(1,I,K) - A(1,1,I,K)*Y(1,I,K-1) -                  &
            A(1,2,I,K)*Y(2,I,K-1)
          X(2,I,K) = F(2,I,K) - A(2,1,I,K)*Y(1,I,K-1) -                  &
            A(2,2,I,K)*Y(2,I,K-1)
!         ****
!         ****     Y(K) = BETA(K)*X(K)
!         ****
          Y(1,I,K) = BETA(1,1,I,K)*X(1,I,K) + BETA(1,2,I,K)*X(2,I,K)
          Y(2,I,K) = BETA(2,1,I,K)*X(1,I,K) + BETA(2,2,I,K)*X(2,I,K)
        ENDDO
      ENDDO
!     ****
!     ****     Backward sweep to determine final solution, X(K) for
!     ****     K = K2,K1,-1
!     ****
!     ****      X(K2) = Y(K2)
!     ****
      DO I = I1,I2
        X(1,I,K2) = Y(1,I,K2)
        X(2,I,K2) = Y(2,I,K2)
      ENDDO
!     ****
!     ****      X(K) = Y(K) - GAMMA(K)*X(K+1)
!     ****
      DO K = K2-1,K1,-1
        DO I = I1,I2
          X(1,I,K) = Y(1,I,K) - GAMMA(1,1,I,K)*X(1,I,K+1) -             &
            GAMMA(1,2,I,K)*X(2,I,K+1)
          X(2,I,K) = Y(2,I,K) - GAMMA(2,1,I,K)*X(1,I,K+1) -             &
            GAMMA(2,2,I,K)*X(2,I,K+1)

        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE BLKTRI

end module majorsp_diffusion
