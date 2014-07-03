module molec_diff

  !------------------------------------------------------------------------------------------------- !
  ! Module to compute molecular diffusivity for various constituents                                 !
  !                                                                                                  !    
  ! Public interfaces :                                                                              !
  !                                                                                                  !
  !    init_molec_diff           Initializes time independent coefficients                           !
  !    init_timestep_molec_diff  Time-step initialization for molecular diffusivity                  ! 
  !    compute_molec_diff        Computes constituent-independent terms for moleculuar diffusivity   !
  !    vd_lu_qdecomp             Computes constituent-dependent terms for moleculuar diffusivity and !
  !                              updates terms in the triadiagonal matrix used for the implicit      !
  !                              solution of the diffusion equation                                  !
  !                                                                                                  !
  !---------------------------Code history---------------------------------------------------------- !
  ! Modularized     :  J. McCaa, September 2004                                                      !
  ! Lastly Arranged :  S. Park,  January.  2010                                                      !
  !                    M. Mills, November  2011
  !------------------------------------------------------------------------------------------------- !

  use perf_mod
  use physconst,    only : mbarv
  use constituents, only : pcnst
  use phys_control, only : waccmx_is             !WACCM-X runtime switch
  use ref_pres,     only : nbot_molec, ntop_molec

  implicit none
  private       
  save

  public init_molec_diff  
  public init_timestep_molec_diff
  public compute_molec_diff 
  public vd_lu_qdecomp

  ! ---------- !
  ! Parameters ! 
  ! ---------- !

  integer,  parameter   :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8), parameter   :: km_fac = 3.55E-7_r8         ! Molecular viscosity constant [ unit ? ]
  real(r8), parameter   :: pr_num = 1._r8              ! Prandtl number [ no unit ]
  real(r8), parameter   :: pwr    = 2._r8/3._r8        ! Exponentiation factor [ unit ? ]
  real(r8), parameter   :: d0     = 1.52E20_r8         ! Diffusion factor [ m-1 s-1 ] molec sqrt(kg/kmol/K) [ unit ? ]
                                                       ! Aerononmy, Part B, Banks and Kockarts (1973), p39
                                                       ! Note text cites 1.52E18 cm-1 ...

  real(r8)              :: rair                        ! Gas constant for dry air
  real(r8)              :: mw_dry                      ! Molecular weight of dry air
  real(r8)              :: n_avog                      ! Avogadro's number [ molec/kmol ]
  real(r8)              :: gravit     
  real(r8)              :: cpair
  real(r8)              :: kbtz                        ! Boltzman constant

  real(r8), allocatable :: mw_fac(:)                   ! sqrt(1/M_q + 1/M_d) in constituent diffusivity [  unit ? ]
  real(r8), allocatable :: alphath(:)                  ! Thermal diffusion factor, -0.38 for H, 0 for others
  
contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_molec_diff( kind, ncnst, rair_in, mw_dry_in, n_avog_in, gravit_in, &
                              cpair_in, kbtz_in, errstring)
    
    use constituents,     only : cnst_mw, cnst_get_ind
    use upper_bc,         only : ubc_init
    use physics_buffer,   only : physics_buffer_desc

    integer,  intent(in)  :: kind           ! Kind of reals being passed in
    integer,  intent(in)  :: ncnst          ! Number of constituents
    real(r8), intent(in)  :: rair_in
    real(r8), intent(in)  :: mw_dry_in      ! Molecular weight of dry air
    real(r8), intent(in)  :: n_avog_in      ! Avogadro's number [ molec/kmol ]
    real(r8), intent(in)  :: gravit_in
    real(r8), intent(in)  :: cpair_in
    real(r8), intent(in)  :: kbtz_in        ! Boltzman constant

    character(len=*), intent(out) :: errstring
    
    ! Local
    
    integer               :: k              ! Level index
    integer               :: m              ! Constituent index
    integer               :: indx_H         ! Constituent index for H
    integer               :: ierr           ! Allocate error check

    errstring = ' '
    
    rair       = rair_in
    mw_dry     = mw_dry_in
    n_avog     = n_avog_in
    gravit     = gravit_in
    cpair      = cpair_in
    kbtz       = kbtz_in

    if( kind /= r8 ) then
       errstring = 'inconsistent KIND of reals passed to init_molec_diff'
       return
    end if

  ! Initialize upper boundary condition variables

    call ubc_init()

  ! Molecular weight factor in constitutent diffusivity
  ! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****
 
    allocate(mw_fac(ncnst))
    do m = 1, ncnst
       mw_fac(m) = d0 * mw_dry * sqrt(1._r8/mw_dry + 1._r8/cnst_mw(m)) / n_avog
    end do

    !--------------------------------------------------------------------------------------------
    ! For WACCM-X, get H data index and initialize thermal diffusion coefficient
    !--------------------------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
    
      call cnst_get_ind('H',  indx_H)
 
      allocate(alphath(ncnst), stat=ierr)
      if ( ierr /= 0 ) then
         errstring = 'allocate failed in init_molec_diff'
         return
      end if
      alphath(:ncnst) = 0._r8
      alphath(indx_H) = -0.38_r8
         
    endif

  end subroutine init_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_timestep_molec_diff(pbuf2d, state)
    !--------------------------- !
    ! Timestep dependent setting ! 
    !--------------------------- !
    use upper_bc,     only : ubc_timestep_init
    use physics_types,only: physics_state
    use ppgrid,       only: begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call ubc_timestep_init( pbuf2d, state)
    
  end subroutine init_timestep_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  integer function compute_molec_diff( lchnk             ,                                          &
       pcols             , pver                , ncnst     , ncol     , t      , pmid   , pint   ,  &
       zi                , ztodt               , kvm       , kvt      , tint   , rhoi   , tmpi2  ,  &
       kq_scal           , ubc_t               , ubc_mmr   , ubc_flux , dse_top, cc_top , cd_top ,  &
       cnst_mw_out       , cnst_fixed_ubc_out  , cnst_fixed_ubflx_out , mw_fac_out      ,           &
       ntop_molec_out    , nbot_molec_out      , kvt_returned )

    use upper_bc,        only : ubc_get_vals
    use constituents,    only : cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx
    use physconst,       only : cpairv, rairv, kmvis, kmcnd

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    integer,  intent(in)    :: lchnk                     ! Chunk identifier
    real(r8), intent(in)    :: t(pcols,pver)             ! Temperature input
    real(r8), intent(in)    :: pmid(pcols,pver)          ! Midpoint pressures
    real(r8), intent(in)    :: pint(pcols,pver+1)        ! Interface pressures
    real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t
    
    real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Viscosity ( diffusivity for momentum )
    real(r8), intent(out)   :: kvt(pcols,pver+1)         ! Kinematic molecular conductivity
    real(r8), intent(inout) :: tint(pcols,pver+1)        ! Interface temperature
    real(r8), intent(inout) :: rhoi(pcols,pver+1)        ! Density ( rho ) at interfaces
    real(r8), intent(inout) :: tmpi2(pcols,pver+1)       ! dt*(g*rho)**2/dp at interfaces

    real(r8), intent(out)   :: kq_scal(pcols,pver+1)     ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(out)   :: ubc_mmr(pcols,ncnst)      ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(out)   :: ubc_flux(ncnst)           ! Upper boundary flux [ kg/s/m^2 ]
    real(r8), intent(out)   :: cnst_mw_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubc_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubflx_out(ncnst)
    real(r8), intent(out)   :: mw_fac_out(pcols,pver+1,ncnst) ! composition dependent mw_fac on interface level
    real(r8), intent(out)   :: dse_top(pcols)            ! dse on top boundary
    real(r8), intent(out)   :: cc_top(pcols)             ! Lower diagonal at top interface
    real(r8), intent(out)   :: cd_top(pcols)             ! cc_top * dse ubc value
    integer,  intent(out)   :: ntop_molec_out
    integer,  intent(out)   :: nbot_molec_out
    logical,  intent(out)   :: kvt_returned              ! Whether we actually returned kvt (vs. kvh).

    ! --------------- !
    ! Local variables !
    ! --------------- !

    integer                 :: m                          ! Constituent index
    integer                 :: i                          ! Column index
    integer                 :: k                          ! Level index

    real(r8)                :: mbarvi                     ! mbarv on interface level
    real(r8)                :: km_top(pcols)              ! molecular conductivity at the top

    real(r8)                :: mkvisc                     ! Molecular kinematic viscosity c*tint**(2/3)/rho
    real(r8)                :: ubc_t(pcols)               ! Upper boundary temperature (K)

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! We don't apply cpairv to kvt if WACCM-X is on.
    kvt_returned = ( waccmx_is('ionosphere') .or. waccmx_is('neutral') )

  ! Get upper boundary values

    call ubc_get_vals( lchnk, ncol, ntop_molec, pint, zi, ubc_t, ubc_mmr, ubc_flux )

  ! Below are already computed, just need to be copied for output

    cnst_mw_out(:ncnst)          = cnst_mw(:ncnst)
    cnst_fixed_ubc_out(:ncnst)   = cnst_fixed_ubc(:ncnst)
    cnst_fixed_ubflx_out(:ncnst) = cnst_fixed_ubflx(:ncnst)
    ntop_molec_out               = ntop_molec
    nbot_molec_out               = nbot_molec
    !
    !  Need variable mw_fac for kvt and constant otherwise
    !
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
      do m = 1, ncnst
        do k = ntop_molec+1, nbot_molec-1
          do i = 1, ncol
             mbarvi = 0.5_r8 * (mbarv(i,k-1,lchnk)+mbarv(i,k,lchnk))
             mw_fac_out(i,k,m) = d0 * mbarvi * sqrt(1._r8/mbarvi + 1._r8/cnst_mw(m)) / n_avog
          enddo
        enddo
        mw_fac_out(:ncol,ntop_molec,m) = 1.5_r8*mw_fac_out(:ncol,ntop_molec+1,m)-.5_r8*mw_fac_out(:ncol,ntop_molec+2,m)
        do k = nbot_molec, pver+1
          mw_fac_out(:ncol,k,m) = mw_fac_out(:ncol,nbot_molec-1,m)
        enddo
      end do
    else
      do k = 1, pver+1
        do i = 1, ncol
          mw_fac_out(i,k,:ncnst) = mw_fac(:ncnst)
        enddo
      enddo
    endif
    
  ! Density and related factors for molecular diffusion and ubc.
  ! Always have a fixed upper boundary T if molecular diffusion is active. Why ?
  ! For kvt, set ubc temperature to average of next two lower interface level temperatures

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
      tint(:ncol,ntop_molec) = 1.5_r8*tint(:ncol,ntop_molec+1)-.5_r8*tint(:ncol,ntop_molec+2)
    else
      tint (:ncol,ntop_molec) = ubc_t(:ncol)    
    endif
    
    rhoi (:ncol,ntop_molec) = pint(:ncol,ntop_molec) / ( rairv(:ncol,ntop_molec,lchnk) * tint(:ncol,ntop_molec) )
    tmpi2(:ncol,ntop_molec) = ztodt * ( gravit * rhoi(:ncol,ntop_molec))**2 &
                                    / ( pmid(:ncol,ntop_molec) - pint(:ncol,ntop_molec) )
    
  ! Compute molecular kinematic viscosity, heat diffusivity and factor for constituent diffusivity
  ! This is a key part of the code.  For WACCM-X, use constituent dependent molecular viscosity and conductivity

    kvt     = 0._r8
    kq_scal = 0._r8
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
      do k = ntop_molec, nbot_molec
         do i = 1, ncol
           mkvisc  = kmvis(i,k,lchnk) / rhoi(i,k)
           kvm(i,k) = kvm(i,k) + mkvisc
           mkvisc  = kmcnd(i,k,lchnk) / rhoi(i,k)
           kvt(i,k) = mkvisc
           kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
         end do
      end do
    else
      do k = ntop_molec, nbot_molec
        do i = 1, ncol
          mkvisc   = km_fac * tint(i,k)**pwr / rhoi(i,k)
          kvm(i,k) = kvm(i,k) + mkvisc
          kvt(i,k) = mkvisc * pr_num * cpairv(i,k,lchnk)  
          kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
        end do
      end do
    endif
    
  ! Top boundary condition for dry static energy

    dse_top(:ncol) = cpairv(:ncol,ntop_molec,lchnk) * tint(:ncol,ntop_molec) + gravit * zi(:ncol,ntop_molec)

  ! Top value of cc for dry static energy

    do i = 1, ncol
      cc_top(i) = ztodt * gravit**2 * rhoi(i,ntop_molec) * km_fac * ubc_t(i)**pwr / &
    		 ( ( pint(i,2) - pint(i,1) ) * ( pmid(i,1) - pint(i,1) ) )
    enddo

    cd_top(:ncol) = cc_top(:ncol) * dse_top(:ncol)
    
    compute_molec_diff = 1
    return
  end function compute_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  integer function vd_lu_qdecomp( pcols , pver   , ncol       , fixed_ubc  , mw     , ubc_mmr , &
                                  kv    , kq_scal, mw_facm    , tmpi       , rpdel  ,           &
                                  decomp, rhoi   ,           &
                                  tint  , ztodt  , ntop_molec , nbot_molec , cd_top ,           &
                                  lchnk , pmid   , pint       , t          , m      )

    use infnan, only: nan, assignment(=)
    use vdiff_lu_solver, only: lu_decomp

    !------------------------------------------------------------------------------ !
    ! Add the molecular diffusivity to the turbulent diffusivity for a consitutent. !
    ! Update the superdiagonal (ca(k)), diagonal (cb(k)) and subdiagonal (cc(k))    !
    ! coefficients of the tridiagonal diffusion matrix, also ze and denominator.    !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncol                  ! Number of atmospheric columns

    integer,  intent(in)    :: ntop_molec
    integer,  intent(in)    :: nbot_molec

    logical,  intent(in)    :: fixed_ubc             ! Fixed upper boundary condition flag
    real(r8), intent(in)    :: kv(pcols,pver+1)      ! Eddy diffusivity
    real(r8), intent(in)    :: kq_scal(pcols,pver+1) ! Molecular diffusivity ( kq_fac*sqrt(T)*m_d/rho )
    real(r8), intent(in)    :: mw                    ! Molecular weight for this constituent
    real(r8), intent(in)    :: ubc_mmr(pcols)        ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(in)    :: mw_facm(pcols,pver+1) ! composition dependent sqrt(1/M_q + 1/M_d) for this constituent
    real(r8), intent(in)    :: tmpi(pcols,pver+1)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in)    :: rpdel(pcols,pver)     ! 1./pdel ( thickness bet interfaces )
    real(r8), intent(in)    :: rhoi(pcols,pver+1)    ! Density at interfaces [ kg/m3 ]
    real(r8), intent(in)    :: tint(pcols,pver+1)    ! Interface temperature [ K ]
    real(r8), intent(in)    :: ztodt                 ! 2 delta-t [ s ]

    ! LU decomposition information.
    type(lu_decomp), intent(inout) :: decomp

    real(r8), intent(out)   :: cd_top(pcols)         ! Term for updating top level with ubc

    integer,  intent(in)    :: lchnk		    ! Chunk number
    real(r8), intent(in)    :: pmid(pcols,pver)     ! midpoint pressures
    real(r8), intent(in)    :: pint(pcols,pver+1)   ! interface pressures
    real(r8), intent(in)    :: t(pcols,pver)	    ! temperature
    integer,  intent(in)    :: m 		    ! cnst index 

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer                 :: i                     ! Longitude index
    integer                 :: k, kp1                ! Vertical indicies

    real(r8)                :: rghd(pcols,pver+1)    ! (1/H_i - 1/H) * (rho*g)^(-1)
    real(r8)                :: kmq(ncol,pver+1)      ! Molecular diffusivity for constituent
    real(r8)                :: wrk0(ncol)            ! Work variable
    real(r8)                :: wrk1(ncol)            ! Work variable

    real(r8)                :: cb(pcols,pver)        ! - Diagonal

    real(r8)                :: gradm(pcols,pver+1)   ! 1/mbar * d(mbar)/dp *(rho*g)^2 * dt
    real(r8)                :: gradt(pcols,pver+1)   ! alphaTh*(rho*g)^2 1/T * dT/dp * dt, for now alphaTh is non-zero only for H.
    real(r8)                :: mbarvi                ! mbarv at interface

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !   

    ! --------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are !
    ! a combination of ca and cc; they are not required by the solver.      !
    !---------------------------------------------------------------------- !

    call t_startf('vd_lu_qdecomp')

    kmq(:,:)  = 0._r8
    cd_top(:) = 0._r8
    rghd(:,:) = nan
    cb(:,:) = nan

  ! Compute difference between scale heights of constituent and dry air

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     

      do k = ntop_molec+1, nbot_molec-1
    	 do i = 1, ncol
    	    mbarvi = 0.5_r8 * (mbarv(i,k-1,lchnk)+mbarv(i,k,lchnk))
    	    rghd(i,k)  = gravit  / (kbtz * n_avog * tint(i,k)) * (mw - mbarvi)
    	    rghd(i,k)  = ztodt * gravit * rhoi(i,k) * rghd(i,k)
    	    gradm(i,k) = (mbarv(i,k,lchnk)-mbarv(i,k-1,lchnk))/(pmid(i,k)-pmid(i,k-1))/mbarvi
    	    gradm(i,k) = ztodt * (rhoi(i,k) * gravit)**2 * gradm(i,k)
    	    gradt(i,k) = (t(i,k)-t(i,k-1))/(pmid(i,k)-pmid(i,k-1))/tint(i,k)
    	    gradt(i,k) = ztodt * alphath(m) * (rhoi(i,k) * gravit)**2 * gradt(i,k)
    	 enddo
      enddo
      do k = ntop_molec,ntop_molec
    	 do i = 1, ncol
    	    mbarvi = .75_r8*mbarv(i,k,lchnk)+0.5_r8*mbarv(i,k+1,lchnk)-.25_r8*mbarv(i,k+2,lchnk)
    	    rghd(i,k)  = gravit  / (kbtz * n_avog * tint(i,k)) * (mw - mbarvi)
    	    rghd(i,k)  = ztodt * gravit * rhoi(i,k) * rghd(i,k)
    	    gradm(i,k) = (mbarv(i,k,lchnk)-mbarvi)/(pmid(i,k)-pint(i,k))/(mbarv(i,k,lchnk)+mbarvi)*2._r8
    	    gradm(i,k) = ztodt * (rhoi(i,k) * gravit)**2 * gradm(i,k)
    	    gradt(i,k) = (t(i,k)-tint(i,k))/(pmid(i,k)-pint(i,k))/(t(i,k)+tint(i,k))*2._r8
    	    gradt(i,k) = ztodt * alphath(m) * (rhoi(i,k) * gravit)**2 * gradt(i,k)
    	 enddo
      enddo
      do k = nbot_molec,nbot_molec
    	 do i = 1, ncol
    	    mbarvi = mbarv(i,k-1,lchnk)
    	    rghd(i,k)  = gravit  / (kbtz * n_avog * tint(i,k)) * (mw - mbarvi)
    	    rghd(i,k)  = ztodt * gravit * rhoi(i,k) * rghd(i,k)
    	    gradm(i,k) = 0._r8
    	    gradt(i,k) = 0._r8						     ! set to zero because molecular diffusion is small at the lower boundary
    	 enddo
      enddo

    else

      do k = ntop_molec, nbot_molec
    	 do i = 1, ncol
    	    rghd(i,k) = gravit / ( kbtz * n_avog * tint(i,k) ) * ( mw - mw_dry )
    	    rghd(i,k) = ztodt * gravit * rhoi(i,k) * rghd(i,k)
    	 enddo
      enddo
      
    endif

    !-------------------- !
    ! Molecular diffusion !
    !-------------------- !

    do k = nbot_molec - 1, ntop_molec, -1
       kp1 = k + 1
       kmq(:ncol,kp1)  = kq_scal(:ncol,kp1) * mw_facm(:ncol,kp1)
       wrk0(:ncol) = ( kv(:ncol,kp1) + kmq(:ncol,kp1) ) * tmpi(:ncol,kp1)
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
         wrk1(:ncol) = .5_r8 * (kmq(:ncol,kp1) * rghd(:ncol,kp1) &
                            - (kv(:ncol,kp1) + kmq(:ncol,kp1)) * gradm(:ncol,kp1) &
                            - kmq(:ncol,kp1) * gradt(:ncol,kp1))
       else
         wrk1(:ncol) = kmq(:ncol,kp1) * 0.5_r8 * rghd(:ncol,kp1)
       endif
     ! Add species separation term
       decomp%ca(:ncol,k  )  = ( wrk0(:ncol) - wrk1(:ncol) ) * rpdel(:ncol,k)
       decomp%cc(:ncol,kp1)  = ( wrk0(:ncol) + wrk1(:ncol) ) * rpdel(:ncol,kp1)
    end do
    
    if( fixed_ubc ) then
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
          kmq(:ncol,ntop_molec)  = kq_scal(:ncol,ntop_molec) * mw_facm(:ncol,ntop_molec)
          wrk0(:ncol) = (kv(:ncol,ntop_molec) + kmq(:ncol,ntop_molec)) * tmpi(:ncol,ntop_molec)/2._r8
          ! /2. is to extrapolate/(pmid(1)-pint(1)) to /(pmid(1)-pmid(0))
          wrk1(:ncol) = .5_r8 * (kmq(:ncol,ntop_molec) * rghd(:ncol,ntop_molec) &
                      - (kv(:ncol,ntop_molec) + kmq(:ncol,ntop_molec)) * gradm(:ncol,ntop_molec) &
                      - kmq(:ncol,ntop_molec) * gradt(:ncol,ntop_molec))
          decomp%cc(:ncol,ntop_molec)  = (wrk0(:ncol) + wrk1(:ncol)) * rpdel(:ncol,ntop_molec)
       else
          decomp%cc(:ncol,ntop_molec) = kq_scal(:ncol,ntop_molec) * mw_facm(:ncol,ntop_molec) &
                               * ( tmpi(:ncol,ntop_molec) + rghd(:ncol,ntop_molec) )   &
                               * rpdel(:ncol,ntop_molec)
       endif
    end if

  ! Calculate diagonal elements

    do k = nbot_molec - 1, ntop_molec + 1, -1
       kp1 = k + 1
       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
          cb(:ncol,k) = 1._r8 + decomp%ca(:ncol,k) + decomp%cc(:ncol,k) &
                      + rpdel(:ncol,k) * (kmq(:ncol,kp1)*rghd(:ncol,kp1) &
                      - kmq(:ncol,k)*rghd(:ncol,k) &
                      -(kv(:ncol,kp1)+kmq(:ncol,kp1)) * gradm(:ncol,kp1)  &
                      +(kv(:ncol,k)+kmq(:ncol,k)) * gradm(:ncol,k) &
                      -kmq(:ncol,kp1) *gradt(:ncol,kp1) &
                      +kmq(:ncol,k) *gradt(:ncol,k))
       else
          cb(:ncol,k) = 1._r8 + decomp%ca(:ncol,k) + decomp%cc(:ncol,k)     &
                      + rpdel(:ncol,k) * ( kmq(:ncol,kp1) * rghd(:ncol,kp1) &
                      - kmq(:ncol,k) * rghd(:ncol,k) )
       endif
    end do

    k   = ntop_molec
    kp1 = k + 1
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
       if( fixed_ubc ) then
          cb(:ncol,k) = 1._r8 + decomp%ca(:ncol,k) + decomp%cc(:ncol,k) &
                      + rpdel(:ncol,k) * (kmq(:ncol,kp1)*rghd(:ncol,kp1) &
                      - kmq(:ncol,k)*rghd(:ncol,k) &
                      -(kv(:ncol,kp1)+kmq(:ncol,kp1)) * gradm(:ncol,kp1)  &
                      +(kv(:ncol,k)+kmq(:ncol,k)) * gradm(:ncol,k)  &
                      -kmq(:ncol,kp1) * gradt(:ncol,kp1)  &
                      +kmq(:ncol,k) * gradt(:ncol,k))
       else
          cb(:ncol,k) = 1._r8 + decomp%ca(:ncol,k) &
                      + rpdel(:ncol,k) * (kmq(:ncol,kp1)*rghd(:ncol,kp1) &
                      - (kv(:ncol,kp1)+kmq(:ncol,kp1)) * gradm(:ncol,kp1) &
                      -kmq(:ncol,kp1) * gradt(:ncol,kp1))
       end if
    else
       if( fixed_ubc ) then
          cb(:ncol,k) = 1._r8 + decomp%ca(:ncol,k)                            &
                      + rpdel(:ncol,k) * kmq(:ncol,kp1) * rghd(:ncol,kp1)     &
                      + kq_scal(:ncol,ntop_molec) * mw_facm(:ncol,ntop_molec) &
                      * ( tmpi(:ncol,ntop_molec) - rghd(:ncol,ntop_molec) )   &
                      * rpdel(:ncol,ntop_molec)
       else
          cb(:ncol,k) = 1._r8 + decomp%ca(:ncol,k) &
                      + rpdel(:ncol,k) * kmq(:ncol,kp1) * rghd(:ncol,kp1)
       end if
    endif

    k   = nbot_molec
    cb(:ncol,k) = 1._r8 + decomp%cc(:ncol,k) + decomp%ca(:ncol,k) &
                - rpdel(:ncol,k) * kmq(:ncol,k)*rghd(:ncol,k)

  ! Compute term for updating top level mixing ratio for ubc

    if( fixed_ubc ) then
        cd_top(:ncol) = decomp%cc(:ncol,ntop_molec) * ubc_mmr(:ncol)
    end if

    !-------------------------------------------------------- !
    ! Calculate e(k).                                         !
    ! This term is required in solution of tridiagonal matrix ! 
    ! defined by implicit diffusion equation.                 !
    !-------------------------------------------------------- !

    do k = nbot_molec, ntop_molec + 1, -1
       decomp%dnom(:ncol,k) = 1._r8 / &
            ( cb(:ncol,k) - decomp%ca(:ncol,k) * decomp%ze(:ncol,k+1) )
       decomp%ze(:ncol,k)   = decomp%cc(:ncol,k) * decomp%dnom(:ncol,k)
    end do
    k = ntop_molec
    decomp%dnom(:ncol,k) = 1._r8 / &
         ( cb(:ncol,k) - decomp%ca(:ncol,k) * decomp%ze(:ncol,k+1) )

    vd_lu_qdecomp = 1
    call t_stopf('vd_lu_qdecomp')
    return

  end function vd_lu_qdecomp

end module molec_diff
