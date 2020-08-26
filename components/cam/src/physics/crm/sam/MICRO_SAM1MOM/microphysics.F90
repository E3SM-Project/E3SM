module microphysics
  use cloud_mod
  use precip_init_mod
  use precip_proc_mod

  ! module for original SAM bulk microphysics
  ! Marat Khairoutdinov, 2006

  use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information
  use params, only: doprecip, docloud, crm_rknd, asyncid
  use micro_params
  implicit none

  !----------------------------------------------------------------------
  !!! required definitions:

  integer, parameter :: nmicro_fields = 2   ! total number of prognostic water vars

  !!! microphysics prognostic variables are storred in this array:


  integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
  integer, parameter :: index_cloud_ice = 1   ! index for cloud ice (sedimentation)

  ! both variables correspond to mass, not number
  ! SAM1MOM 3D microphysical fields are output by default.
  integer, allocatable :: flag_precip    (:)


  !!! these arrays are needed for output statistics:


  !======================================================================
  ! UW ADDITIONS

  !bloss: arrays with names/units for microphysical outputs in statistics.

  ! END UW ADDITIONS
  !======================================================================

  !------------------------------------------------------------------
  ! Optional (internal) definitions)

  ! make aliases for prognostic variables:
  ! note that the aliases should be local to microphysics

  real(crm_rknd) vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity

  real(crm_rknd), allocatable, target :: micro_field(:,:,:,:,:)
  real(crm_rknd), allocatable :: fluxbmk (:,:,:,:) ! surface flux of tracers
  real(crm_rknd), allocatable :: fluxtmk (:,:,:,:) ! top boundary flux of tracers
  real(crm_rknd), allocatable :: mkwle  (:,:,:)  ! resolved vertical flux
  real(crm_rknd), allocatable :: mkwsb  (:,:,:)  ! SGS vertical flux
  real(crm_rknd), allocatable :: mkadv  (:,:,:)  ! tendency due to vertical advection
  real(crm_rknd), allocatable :: mkdiff (:,:,:)  ! tendency due to vertical diffusion
  character*3   , allocatable :: mkname       (:)
  character*80  , allocatable :: mklongname   (:)
  character*10  , allocatable :: mkunits      (:)
  real(crm_rknd), allocatable :: mkoutputscale(:)
  real(crm_rknd), allocatable :: qn(:,:,:,:)  ! cloud condensate (liquid + ice)
  real(crm_rknd), allocatable :: qpsrc(:,:)  ! source of precipitation microphysical processes
  real(crm_rknd), allocatable :: qpevp(:,:)  ! sink of precipitating water due to evaporation


CONTAINS


  subroutine allocate_micro(ncrms)
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    integer :: icrm
    real(crm_rknd) :: zero
    allocate( micro_field(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields))
    allocate( fluxbmk(ncrms,nx, ny, 1:nmicro_fields) )
    allocate( fluxtmk(ncrms,nx, ny, 1:nmicro_fields) )
    allocate( mkwle(ncrms,nz,1:nmicro_fields)  )
    allocate( mkwsb(ncrms,nz,1:nmicro_fields)  )
    allocate( mkadv(ncrms,nz,1:nmicro_fields)  )
    allocate( mkdiff(ncrms,nz,1:nmicro_fields)  )
    allocate( mkname       (nmicro_fields))
    allocate( mklongname   (nmicro_fields))
    allocate( mkunits      (nmicro_fields))
    allocate( mkoutputscale(nmicro_fields))
    allocate( qn(ncrms,nx,ny,nzm)  )
    allocate( qpsrc(ncrms,nz)  )
    allocate( qpevp(ncrms,nz)  )
    allocate( flag_precip    (nmicro_fields) )

    call prefetch(micro_field  )
    call prefetch(fluxbmk   )
    call prefetch(fluxtmk   )
    call prefetch(mkwle    )
    call prefetch(mkwsb    )
    call prefetch(mkadv    )
    call prefetch(mkdiff   )
    call prefetch(mkoutputscale  )
    call prefetch(qn  )
    call prefetch(qpsrc  )
    call prefetch(qpevp  )
    call prefetch(flag_precip    )

    zero = 0

    micro_field = zero
    fluxbmk  = zero
    fluxtmk  = zero
    mkwle   = zero
    mkwsb   = zero
    mkadv   = zero
    mkdiff  = zero
    mkname        = ''
    mklongname    = ''
    mkunits       = ''
    mkoutputscale = zero
    qn = zero
    qpsrc = zero
    qpevp = zero
    flag_precip    (:)  = (/0,1/)
  end subroutine allocate_micro


  subroutine deallocate_micro()
    implicit none
    deallocate(micro_field  )
    deallocate(fluxbmk   )
    deallocate(fluxtmk   )
    deallocate(mkwle    )
    deallocate(mkwsb    )
    deallocate(mkadv    )
    deallocate(mkdiff   )
    deallocate(mkname         )
    deallocate(mklongname     )
    deallocate(mkunits        )
    deallocate(mkoutputscale  )
    deallocate(qn  )
    deallocate(qpsrc  )
    deallocate(qpevp  )
    deallocate(flag_precip    )
  end subroutine deallocate_micro


  ! required microphysics subroutines and function:
  !----------------------------------------------------------------------
  !!! Read microphysics options from prm file

  subroutine micro_setparm()
    ! no user-definable options in SAM1MOM microphysics.
  end subroutine micro_setparm

  !----------------------------------------------------------------------
  !!! Initialize microphysics:


  subroutine micro_init(ncrms)

    use grid, only: nrestart
    use vars
    implicit none
    integer, intent(in) :: ncrms
    integer k, n,icrm, i, j, l

    a_bg = 1./(tbgmax-tbgmin)
    a_pr = 1./(tprmax-tprmin)
    a_gr = 1./(tgrmax-tgrmin)

    if(nrestart.eq.0) then

    !$acc parallel loop collapse(4) async(asyncid)
    do l=1,nmicro_fields
      do j=1,ny
        do i=1,nx
          do icrm=1,ncrms
            fluxbmk(icrm,i,j,l) = 0.
            fluxtmk(icrm,i,j,l) = 0.
          enddo
        enddo
      enddo
    enddo

      if(docloud) then
        call micro_diagnose(ncrms)
      end if
    end if

    !$acc parallel loop collapse(3) async(asyncid)
    do l=1,nmicro_fields
      do k=1,nz
        do icrm = 1 , ncrms
          mkwle (icrm,k,l) = 0.
          mkwsb (icrm,k,l) = 0.
          mkadv (icrm,k,l) = 0.
          mkdiff(icrm,k,l) = 0.
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nz
      do icrm=1,ncrms
        qpsrc(icrm,k) = 0.
        qpevp(icrm,k) = 0.
      enddo
    enddo

    mkname(1) = 'QT'
    mklongname(1) = 'TOTAL WATER (VAPOR + CONDENSATE)'
    mkunits(1) = 'g/kg'
    mkoutputscale(1) = 1.e3

    mkname(2) = 'QP'
    mklongname(2) = 'PRECIPITATING WATER'
    mkunits(2) = 'g/kg'
    mkoutputscale(2) = 1.e3

  end subroutine micro_init

  !----------------------------------------------------------------------
  !!! fill-in surface and top boundary fluxes:
  subroutine micro_flux(ncrms)
    use vars, only: fluxbq, fluxtq
    implicit none
    integer, intent(in) :: ncrms
    integer :: icrm, i, j

    !$acc parallel loop collapse(3) async(asyncid)
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          fluxbmk(icrm,i,j,index_water_vapor) = fluxbq(icrm,i,j)
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) async(asyncid)
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          fluxtmk(icrm,i,j,index_water_vapor) = fluxtq(icrm,i,j)
        enddo
      enddo
    enddo
  end subroutine micro_flux

  !----------------------------------------------------------------------
  !!! compute local microphysics processes (bayond advection and SGS diffusion):
  !
  subroutine micro_proc(ncrms)
    use grid, only: nstep,dt,icycle
    use cloud_mod
    use precip_init_mod
    use precip_proc_mod

    implicit none
    integer, intent(in) :: ncrms
    integer :: icrm

    ! Update bulk coefficient
    if(doprecip.and.icycle.eq.1) call precip_init(ncrms)

    if(docloud) then
      call cloud(ncrms, micro_field(:,:,:,:,1), micro_field(:,:,:,:,2), qn)
      if(doprecip) call precip_proc(ncrms, qpsrc, qpevp, micro_field(:,:,:,:,1), micro_field(:,:,:,:,2), qn)
      call micro_diagnose(ncrms)
    end if
  end subroutine micro_proc

  !----------------------------------------------------------------------
  !!! Diagnose arrays nessesary for dynamical core and statistics:
  !
  subroutine micro_diagnose(ncrms)
    use vars
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) omn, omp
    integer i,j,k,icrm

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            qv(icrm,i,j,k) = micro_field(icrm,i,j,k,1) - qn(icrm,i,j,k)
            omn = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tbgmin)*a_bg))
            qcl(icrm,i,j,k) = qn(icrm,i,j,k)*omn
            qci(icrm,i,j,k) = qn(icrm,i,j,k)*(1.-omn)
            omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
            qpl(icrm,i,j,k) = micro_field(icrm,i,j,k,2)*omp
            qpi(icrm,i,j,k) = micro_field(icrm,i,j,k,2)*(1.-omp)
          end do
        end do
      end do
    end do
  end subroutine micro_diagnose

  !----------------------------------------------------------------------
  !!! function to compute terminal velocity for precipitating variables:
  ! In this particular case there is only one precipitating variable.

  subroutine term_vel_qp(ncrms,icrm,i,j,k,ind,qploc,rho,tabs,qp_threshold,tprmin,&
                                      a_pr,vrain,crain,tgrmin,a_gr,vgrau,cgrau,vsnow,csnow,term_vel)
    !$acc routine seq
    implicit none
    integer, intent(in) :: ncrms,icrm
    integer, intent(in) :: i,j,k,ind
    real(crm_rknd), intent(in) :: qploc
    real(crm_rknd), intent(in) :: rho(ncrms,nzm), tabs(ncrms,nx, ny, nzm)
    real(crm_rknd), intent(in) :: qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,a_gr,vgrau,cgrau,vsnow,csnow
    real(crm_rknd), intent(out) :: term_vel
    real(crm_rknd) wmax, omp, omg, qrr, qss, qgg

    term_vel = 0.
    if(qploc.gt.qp_threshold) then
      omp = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
      if(omp.eq.1.) then
        term_vel = vrain*(rho(icrm,k)*qploc)**crain
      elseif(omp.eq.0.) then
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
        qgg=omg*qploc
        qss=qploc-qgg
        term_vel = (omg*vgrau*(rho(icrm,k)*qgg)**cgrau &
        +(1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow)
      else
        omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))
        qrr=omp*qploc
        qss=qploc-qrr
        qgg=omg*qss
        qss=qss-qgg
        term_vel = (omp*vrain*(rho(icrm,k)*qrr)**crain + (1.-omp)*(omg*vgrau*(rho(icrm,k)*qgg)**cgrau + &
                      (1.-omg)*vsnow*(rho(icrm,k)*qss)**csnow))
      endif
    end if
  end subroutine term_vel_qp

  !----------------------------------------------------------------------
  !!! compute sedimentation
  !
  subroutine micro_precip_fall(ncrms)
    use vars
    use params, only : pi
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), allocatable :: omega(:,:,:,:)
    integer ind
    integer i,j,k,icrm

    allocate(omega(ncrms,nx,ny,nzm))
    call prefetch( omega )

    crain = b_rain / 4.
    csnow = b_snow / 4.
    cgrau = b_grau / 4.
    vrain = a_rain * gamr3 / 6. / (pi * rhor * nzeror) ** crain
    vsnow = a_snow * gams3 / 6. / (pi * rhos * nzeros) ** csnow
    vgrau = a_grau * gamg3 / 6. / (pi * rhog * nzerog) ** cgrau

    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            omega(icrm,i,j,k) = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tprmin)*a_pr))
          end do
        end do
      end do
    end do

    call precip_fall(ncrms, 2, omega, ind)

    deallocate(omega)

  end subroutine micro_precip_fall


  subroutine precip_fall(ncrms,hydro_type, omega, ind)
    !     positively definite monotonic advection with non-oscillatory option
    !     and gravitational sedimentation
    use vars
    use openacc_utils
    use params
    implicit none
    integer, intent(in) :: ncrms
    integer :: hydro_type   ! 0 - all liquid, 1 - all ice, 2 - mixed
    real(crm_rknd) :: omega(ncrms,nx,ny,nzm)   !  = 1: liquid, = 0: ice;  = 0-1: mixed : used only when hydro_type=2
    integer :: ind
    ! Terminal velocity fnction
    ! Local:
    real(crm_rknd), allocatable :: mx     (:,:,:,:)
    real(crm_rknd), allocatable :: mn     (:,:,:,:)
    real(crm_rknd), allocatable :: lfac   (:,:,:,:)
    real(crm_rknd), allocatable :: www    (:,:,:,:)
    real(crm_rknd), allocatable :: fz     (:,:,:,:)
    real(crm_rknd), allocatable :: wp     (:,:,:,:)
    real(crm_rknd), allocatable :: tmp_qp (:,:,:,:)
    real(crm_rknd), allocatable :: irhoadz(:,:)
    real(crm_rknd), allocatable :: iwmax  (:,:)
    real(crm_rknd), allocatable :: rhofac (:,:)
    real(crm_rknd) :: prec_cfl
    real(crm_rknd) :: eps
    integer :: i,j,k,kc,kb,icrm
    logical :: nonos
    real(crm_rknd) :: y,pp,pn
    real(crm_rknd) :: lat_heat, wmax
    integer nprec, iprec
    real(crm_rknd) :: flagstat, tmp

    !Statement functions
    pp(y)= max(real(0.,crm_rknd),y)
    pn(y)=-min(real(0.,crm_rknd),y)

    eps = 1.e-10
    nonos = .true.

    allocate( mx     (ncrms,nx,ny,nzm) )
    allocate( mn     (ncrms,nx,ny,nzm) )
    allocate( lfac   (ncrms,nx,ny,nz ) )
    allocate( www    (ncrms,nx,ny,nz ) )
    allocate( fz     (ncrms,nx,ny,nz ) )
    allocate( wp     (ncrms,nx,ny,nzm) )
    allocate( tmp_qp (ncrms,nx,ny,nzm) )
    allocate( irhoadz(ncrms,nzm) )
    allocate( iwmax  (ncrms,nzm) )
    allocate( rhofac (ncrms,nzm) )
    
    call prefetch( mx      )
    call prefetch( mn      )
    call prefetch( lfac    )
    call prefetch( www     )
    call prefetch( fz      )
    call prefetch( wp      )
    call prefetch( tmp_qp  )
    call prefetch( irhoadz )
    call prefetch( iwmax   )
    call prefetch( rhofac  )

    !$acc parallel loop gang vector collapse(2) async(asyncid)
    do k = 1,nzm
      do icrm = 1 , ncrms
        rhofac(icrm,k) = sqrt(1.29/rho(icrm,k))
        irhoadz(icrm,k) = 1./(rho(icrm,k)*adz(icrm,k)) ! Useful factor
        kb = max(1,k-1)
        wmax       = dz(icrm)*adz(icrm,kb)/dtn   ! Velocity equivalent to a cfl of 1.0.
        iwmax(icrm,k)   = 1./wmax
      enddo
    enddo

    ! 	Add sedimentation of precipitation field to the vert. vel.
    prec_cfl = 0.
    !$acc parallel loop gang vector collapse(4) reduction(max:prec_cfl) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            select case (hydro_type)
            case(0)
              lfac(icrm,i,j,k) = fac_cond
              flagstat = 1.
            case(1)
              lfac(icrm,i,j,k) = fac_sub
              flagstat = 1.
            case(2)
              lfac(icrm,i,j,k) = fac_cond + (1-omega(icrm,i,j,k))*fac_fus
              flagstat = 1.
            case(3)
              lfac(icrm,i,j,k) = 0.
              flagstat = 0.
            case default
              if(masterproc) then
                !print*, 'unknown hydro_type in precip_fall. exitting ...'
                !call task_abort
              endif
            end select
            call  term_vel_qp(ncrms,icrm,i,j,k,ind,micro_field(icrm,i,j,k,2),rho(1,1),&
                                                      tabs(1,1,1,1),qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,&
                                                      a_gr,vgrau,cgrau,vsnow,csnow,tmp)
            wp(icrm,i,j,k)=rhofac(icrm,k)*tmp
            tmp = wp(icrm,i,j,k)*iwmax(icrm,k)
            prec_cfl = max(prec_cfl,tmp) ! Keep column maximum CFL
            wp(icrm,i,j,k) = -wp(icrm,i,j,k)*rhow(icrm,k)*dtn/dz(icrm)
            if (k == 1) then
              fz(icrm,i,j,nz)=0.
              www(icrm,i,j,nz)=0.
              lfac(icrm,i,j,nz)=0
            endif
          enddo  ! k
        enddo
      enddo
    enddo

    ! If maximum CFL due to precipitation velocity is greater than 0.9,
    ! take more than one advection step to maintain stability.
    if (prec_cfl.gt.0.9) then
      nprec = CEILING(prec_cfl/0.9)
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k = 1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              ! wp already includes factor of dt, so reduce it by a
              ! factor equal to the number of precipitation steps.
              wp(icrm,i,j,k) = wp(icrm,i,j,k)/real(nprec,crm_rknd)
            enddo
          enddo
        enddo
      enddo
    else
      nprec = 1
    endif

    !  loop over iterations
    do iprec = 1,nprec
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k = 1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              tmp_qp(icrm,i,j,k) = micro_field(icrm,i,j,k,2) ! Temporary array for qp in this column
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              if(nonos) then
                kc=min(nzm,k+1)
                kb=max(1,k-1)
                mx(icrm,i,j,k)=max(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k))
                mn(icrm,i,j,k)=min(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k))
              endif  ! nonos
              ! Define upwind precipitation flux
              fz(icrm,i,j,k)=tmp_qp(icrm,i,j,k)*wp(icrm,i,j,k)
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc=k+1
              tmp_qp(icrm,i,j,k)=tmp_qp(icrm,i,j,k)-(fz(icrm,i,j,kc)-fz(icrm,i,j,k))*irhoadz(icrm,k) !Update temporary qp
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              ! Also, compute anti-diffusive correction to previous
              ! (upwind) approximation to the flux
              kb=max(1,k-1)
              ! The precipitation velocity is a cell-centered quantity,
              ! since it is computed from the cell-centered
              ! precipitation mass fraction.  Therefore, a reformulated
              ! anti-diffusive flux is used here which accounts for
              ! this and results in reduced numerical diffusion.
              www(icrm,i,j,k) = 0.5*(1.+wp(icrm,i,j,k)*irhoadz(icrm,k))*(tmp_qp(icrm,i,j,kb)*wp(icrm,i,j,kb) - &
                                     tmp_qp(icrm,i,j,k)*wp(icrm,i,j,k)) ! works for wp(k)<0
            enddo
          enddo
        enddo
      enddo

      !---------- non-osscilatory option ---------------
      if(nonos) then
        !$acc parallel loop gang vector collapse(4) async(asyncid)
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              do icrm = 1 , ncrms
                kc=min(nzm,k+1)
                kb=max(1,k-1)
                mx(icrm,i,j,k)=max(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k),mx(icrm,i,j,k))
                mn(icrm,i,j,k)=min(tmp_qp(icrm,i,j,kb),tmp_qp(icrm,i,j,kc),tmp_qp(icrm,i,j,k),mn(icrm,i,j,k))
                kc=min(nzm,k+1)
                mx(icrm,i,j,k)=rho(icrm,k)*adz(icrm,k)*(mx(icrm,i,j,k)-tmp_qp(icrm,i,j,k))/(pn(www(icrm,i,j,kc)) + pp(www(icrm,i,j,k))+eps)
                mn(icrm,i,j,k)=rho(icrm,k)*adz(icrm,k)*(tmp_qp(icrm,i,j,k)-mn(icrm,i,j,k))/(pp(www(icrm,i,j,kc)) + pn(www(icrm,i,j,k))+eps)
              enddo
            enddo
          enddo
        enddo
        !$acc parallel loop gang vector collapse(4) async(asyncid)
        do k=1,nzm
          do j=1,ny
            do i=1,nx
              do icrm = 1 , ncrms
                kb=max(1,k-1)
                ! Add limited flux correction to fz(k).
                fz(icrm,i,j,k) = fz(icrm,i,j,k) + pp(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,k), mn(icrm,i,j,kb)) - &
                                                  pn(www(icrm,i,j,k))*min(real(1.,crm_rknd),mx(icrm,i,j,kb),mn(icrm,i,j,k)) ! Anti-diffusive flux
              enddo
            enddo
          enddo
        enddo
      endif ! nonos

      ! Update precipitation mass fraction and liquid-ice static
      ! energy using precipitation fluxes computed in this column.
      !$acc parallel loop gang vector collapse(4) async(asyncid)
      do j=1,ny
        do i=1,nx
          do k=1,nzm
            do icrm = 1 , ncrms
              kc=k+1
              ! Update precipitation mass fraction.
              ! Note that fz is the total flux, including both the
              ! upwind flux and the anti-diffusive correction.
              micro_field(icrm,i,j,k,2)=micro_field(icrm,i,j,k,2)-(fz(icrm,i,j,kc)-fz(icrm,i,j,k))*irhoadz(icrm,k)
              tmp = -(fz(icrm,i,j,kc)-fz(icrm,i,j,k))*irhoadz(icrm,k)*flagstat  ! For qp budget
              !$acc atomic update
              qpfall(icrm,k)=qpfall(icrm,k)+tmp
              lat_heat = -(lfac(icrm,i,j,kc)*fz(icrm,i,j,kc)-lfac(icrm,i,j,k)*fz(icrm,i,j,k))*irhoadz(icrm,k)
              t(icrm,i,j,k)=t(icrm,i,j,k)-lat_heat
              !$acc atomic update
              tlat(icrm,k)=tlat(icrm,k)-lat_heat            ! For energy budget
              tmp = fz(icrm,i,j,k)*flagstat
              !$acc atomic update
              precflux(icrm,k) = precflux(icrm,k) - tmp   ! For statistics
              if (k == 1) then
                precsfc(icrm,i,j) = precsfc(icrm,i,j) - fz(icrm,i,j,1)*flagstat ! For statistics
                precssfc(icrm,i,j) = precssfc(icrm,i,j) - fz(icrm,i,j,1)*(1.-omega(icrm,i,j,1))*flagstat ! For statistics
                prec_xy(icrm,i,j) = prec_xy(icrm,i,j) - fz(icrm,i,j,1)*flagstat ! For 2D output
              endif
            enddo
          enddo
        enddo
      enddo

      if (iprec.lt.nprec) then
        ! Re-compute precipitation velocity using new value of qp.
        !$acc parallel loop gang vector collapse(4) async(asyncid)
        do j=1,ny
          do i=1,nx
            do k=1,nzm
              do icrm = 1 , ncrms
                !Passing variables via first index because of PGI bug with pointers
                call term_vel_qp(ncrms,icrm,i,j,k,ind,micro_field(icrm,i,j,k,2),rho(1,1),&
                                 tabs(1,1,1,1),qp_threshold,tprmin,a_pr,vrain,crain,tgrmin,a_gr,vgrau,cgrau,vsnow,csnow,tmp)
                wp(icrm,i,j,k) = rhofac(icrm,k)*tmp
                ! Decrease precipitation velocity by factor of nprec
                wp(icrm,i,j,k) = -wp(icrm,i,j,k)*rhow(icrm,k)*dtn/dz(icrm)/real(nprec,crm_rknd)
                ! Note: Don't bother checking CFL condition at each
                ! substep since it's unlikely that the CFL will
                ! increase very much between substeps when using
                ! monotonic advection schemes.
                if (k == 1) then
                  fz(icrm,i,j,nz)=0.
                  www(icrm,i,j,nz)=0.
                  lfac(icrm,i,j,nz)=0.
                endif
              enddo
            enddo
          enddo
        enddo
      endif

    enddo
    
    deallocate( mx      )
    deallocate( mn      )
    deallocate( lfac    )
    deallocate( www     )
    deallocate( fz      )
    deallocate( wp      )
    deallocate( tmp_qp  )
    deallocate( irhoadz )
    deallocate( iwmax   )
    deallocate( rhofac  )

  end subroutine precip_fall

  !----------------------------------------------------------------------
  ! called when stepout() called

  !-----------------------------------------------------------------------
  ! Supply function that computes total water in a domain:
  !

end module microphysics
