#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_si_ref_mod
  use kinds, only: r8 => real_kind, iulog
  use dimensions_mod, only: plev => nlev, plevp => nlevp
  implicit none
  private

  type, public :: ref_state_t
     real(r8) psr                ! reference surface pressure for linearization
     real(r8) hypi(plevp)        ! reference pressures at interfaces
     real(r8) hypm(plev)         ! reference pressures at midpoints
     real(r8) hypd(plev)         ! reference pressure layer thickness

     real(r8) Tref(plev)         ! reference temperature
     real(r8) RTref(plev)        ! coefficient for ln(ps) term in velocity equation
     real(r8) Pvec(plev)         ! diagonal P matrix
     real(r8) Href(plev,plev)    ! reference hydrostatic matrix (Href)
     real(r8) Cref(plev,plev)    ! reference hydrostatic matrix (Cmat)
     real(r8) Tmat(plev,plev)    ! tau matrix (Tmat)

     real(r8) Amat(plev,plev)     ! vertical structure matrix
     real(r8) Amat_inv(plev,plev) ! inverse vertical structure matrix
     real(r8) Emat(plev,plev)     ! right eigenvector matrix
     real(r8) Emat_inv(plev,plev) ! right eigenvector matrix
     real(r8) Lambda(plev)        ! solver eigenvalues
  end type ref_state_t

  public  :: prim_si_refstate_init
  public  :: set_vert_struct_mat
  public  :: prim_set_mass
contains

  ! =====================================================
  ! prim_si_refstate_init:
  !
  ! given a hybrid vertical coordinate system, initialize
  ! the reference pressures needed by the semi-implicit.
  ! =====================================================

  function prim_si_refstate_init(lprint,masterproc,hvcoord) result(refstate)
    use hybvcoord_mod, only : hvcoord_t

    logical, intent(in) :: lprint
    logical, intent(in) :: masterproc
    type (hvcoord_t)    :: hvcoord
    type (ref_state_t)  :: refstate

    ! =============================
    ! Local variables
    ! =============================

    integer k                 ! Level index


    refstate%psr    = 1000.0_r8  ! Reference surface pressure (millibars)

    ! ================================================================
    ! Set mean temperature
    ! NOTE: Making t0 an actual function of height ***DOES NOT WORK***
    ! ================================================================

    do k=1,plev
       refstate%Tref(k) = 300.0_r8
    end do

    ! ===========================================
    ! Set reference state midpoint pressures
    ! ===========================================

    do k=1,plev
       refstate%hypm(k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*refstate%psr
    end do

    ! ============================================
    ! Reference state interface pressures
    ! ============================================

    do k=1,plevp
       refstate%hypi(k) = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*refstate%psr
    end do

    ! ============================================
    ! Reference state layer thicknesses
    ! ============================================

    do k=1,plev
       refstate%hypd(k) = refstate%hypi(k+1) - refstate%hypi(k)
    end do

    if (lprint.and.masterproc) then
       write(iulog,9820)
       do k=1,plev
          write(iulog,9830) k, refstate%hypi(k)
          write(iulog,9840) refstate%hypm(k), refstate%hypd(k)
       end do
       write(iulog,9830) plevp, refstate%hypi(plevp)
    end if

    call prim_si_refmat_init(lprint,masterproc,refstate)

9820 format(1x,'reference pressures (Pa)')
9830 format(1x,i3,f15.4)
9840 format(1x,3x,15x,2f15.4)

  end function prim_si_refstate_init

  ! ==============================================
  ! prim_si_refmat_init:
  ! 
  ! initialize the reference state hydrostatic and
  ! energy conversion matrices.
  ! ==============================================

  subroutine prim_si_refmat_init(lprint,masterproc,refstate)
    use physical_constants , only :  rgas

    logical, intent(in)  :: lprint            ! print/no print of ref matrices
    logical, intent(in)  :: masterproc        ! is the master process?
    type(ref_state_t)    :: refstate          ! reference state structure

    integer k,l,kk,kkk        ! level indices
    integer ik1,ik2,nkk       ! misc. integers

    ! =================================================================
    ! Integration matrices of hydrostatic equation(href) and conversion
    ! term(a).  href computed as in ccm0 but isothermal bottom ecref
    ! calculated to conserve energy
    ! =================================================================

    do k=1,plev
       do l=1,plev
	  refstate%Href(l,k) = 0.0_r8
	  refstate%Cref(l,k) = 0.0_r8
       end do
    end do

    ! ======================================================================
    ! Mean atmosphere energy conversion term is consistent with continiuty
    ! equation.  In Cref, 1st index = column; 2nd index = row of matrix.
    ! Mean atmosphere energy conversion term is energy conserving
    ! ======================================================================

    do k=1,plev
       refstate%Cref(k,k) = 0.50_r8/refstate%hypm(k) 
       do l=1,k-1
	  refstate%Cref(l,k) = 1.0_r8/refstate%hypm(k) 
       end do
    end do

    ! ====================================================================
    ! Reference hydrostatic integration matrix consistent with conversion
    ! term for energy conservation.  In href, 1st index = column; 
    ! 2nd index = row of matrix.
    ! ====================================================================

    do k = 1,plev
       do l = k,plev
	  refstate%Href(l,k) = refstate%Cref(k,l)*refstate%hypd(l)
       end do
    end do

    ! ==================================
    ! Print statements
    ! ==================================

    if (lprint.and.masterproc) then
       nkk = plev/13
       if (mod(plev,13).ne.0) nkk = nkk + 1
       write(iulog,*)' '
       write(iulog,*)'INITCOM: Hydrostatic matrix href'
       do kk=1,nkk
	  ik1 = 1 + (kk-1)*13
	  ik2 = min0( ik1+12, plev )
	  write(iulog,9920) (k,k=ik1,ik2)
	  do kkk=1,plev
	     write(iulog,9910) kkk,(refstate%Href(kkk,k),k=ik1,ik2)
	  end do
       end do
       write(iulog,*)' '
       write(iulog,*)'INITCOM: Thermodynamic matrix ecref'
       do kk=1,nkk
	  ik1 = 1 + (kk-1)*13
	  ik2 = MIN( ik1+12, plev )
	  write(iulog,9920) (k,k=ik1,ik2)
	  do kkk=1,plev
	     write(iulog,9910) kkk,(refstate%Cref(kkk,k),k=ik1,ik2)
	  end do
       end do
    end if

    ! =======================
    ! Multiply href by r
    ! =======================

    do k=1,plev
       do kk=1,plev
	  refstate%Href(kk,k) = refstate%Href(kk,k)*Rgas
       end do
    end do


9910 format( 1x,i3,13f9.5)
9920 format(/,      13i9)

  end subroutine prim_si_refmat_init

  ! =======================================================================
  ! set_vert_struct_mat:
  !
  ! Purpose: 
  ! Set time invariant hydrostatic matrices, which depend on the reference
  ! temperature and pressure in the semi-implicit time step. Note that
  ! this subroutine is actually called twice, because the effective time
  ! step changes between step 0 and step 1.
  !
  ! based on settau. 
  ! =======================================================================

  subroutine set_vert_struct_mat(dt, refstate, hvcoord, masterproc)
    use parallel_mod, only : abortmp
    use hybvcoord_mod, only : hvcoord_t
    use linear_algebra_mod, only : reigen_solver
    use physical_constants , only : g, rgas, kappa

    !------------------------------Arguments--------------------------------
    real(r8), intent(in)         :: dt       ! time step (or dt/2 at time 0)
    type(ref_state_t), target    :: refstate
    type(hvcoord_t) , target     :: hvcoord
    logical, intent(in)          :: masterproc ! master process
    !---------------------------Local workspace-----------------------------

    real(r8) zcr(plev)             ! gravity wave equivalent depth
    real(r8) zci(plev)             ! dummy, used to print phase speeds

    real(r8) tmp(plev,plev)

    real(r8) rcond
    real(r8) z(plev)
    real(r8) det(2) 
    real(r8) work(plev)
    integer ipvt(plev)

    real(r8) Imat(plev,plev)

    real(r8) factor                ! intermediate workspace
    real(r8) zdt0u                 ! vertical diff. of ref. temp (above)
    real(r8) zshu                  ! interface "sigma" (above)
    real(r8) zr2ds                 ! 1./(2.*hypd(k))
    real(r8) zdt0d                 ! vertical diff. of ref. temp (below)
    real(r8) zshd                  ! interface "sigma" (below)
    real(r8) ztd                   ! temporary accumulator
    real(r8) zcn                   ! sq(n)

    integer k,l,kk,kkk             ! level indices
    integer n,m                    ! n-wavenumber index
    integer nneg                   ! number of unstable mean temperatures
    integer info
    integer ik1,ik2,nkk            ! misc. integers

    real(r8), pointer :: Cref(:,:)
    real(r8), pointer :: Tmat(:,:)
    real(r8), pointer :: Href(:,:)
    real(r8), pointer :: Amat(:,:)
    real(r8), pointer :: Amat_inv(:,:)
    real(r8), pointer :: Emat(:,:)
    real(r8), pointer :: Emat_inv(:,:)

    real(r8), pointer :: RTref(:)
    real(r8), pointer :: Tref(:)
    real(r8), pointer :: Pvec(:)
    real(r8), pointer :: Lambda(:)

    real(r8), pointer :: hypd(:)
    real(r8), pointer :: hypi(:)
    real(r8), pointer :: hybi(:)


    if (masterproc) then
       print *,'Initializing vertical structure matrix'
    endif
    ! =====================================================
    ! Assign pointers to refstate structure for readibility
    ! =====================================================

    Cref     => refstate%Cref
    Tmat     => refstate%Tmat
    Href     => refstate%Href
    Amat     => refstate%Amat
    Amat_inv => refstate%Amat_inv
    Emat     => refstate%Emat
    Emat_inv => refstate%Emat_inv
    Tref     => refstate%Tref
    Pvec     => refstate%Pvec
    RTref    => refstate%RTref
    Lambda   => refstate%Lambda

    hybi => hvcoord%hybi
    hypd => refstate%hypd
    hypi => refstate%hypi

    ! =========================================
    ! Calculate hydrostatic matrix tau (Tmat)
    ! concordance with subroutine settau from CAM:
    ! --------------------------------------------
    !   here     settau 
    !   ---------------
    !   Cref  == ecref
    !   Href  == href
    !   Tref  == t0
    !   RTref == bps
    !   Tmat  == tau
    !   Amat  == zb
    !
    ! =========================================

    ! ===========================================================================
    ! This formula for Tmat (assumes constant (in vertical) reference Temperature
    ! ===========================================================================

    do k=1,plev
       do l=1,k
	  Tmat(l,k) = kappa*Tref(k)*Cref(l,k)*hypd(l)
       end do
       do l=k+1,plev
	  Tmat(l,k) = 0.0_r8
       end do
    end do

    ! ===============================================================
    ! Vector for linear surface pressure term in divergence
    ! Pressure gradient and diagonal term of hydrostatic components
    ! ===============================================================
#ifdef DEBUG
    print *,"hypi=", hypi(plevp)
    print *,"hypd=", hypd(:)
#endif
    do k=1,plev
       RTref(k) = Rgas*Tref(k)
       Pvec(k)  = hypd(k)/hypi(plevp)
    end do
#ifdef DEBUG
    print *,"Pvec=",Pvec(:)
#endif

    do k=1,plev
       do l=1,plev
	  ztd = RTref(k) * Pvec(l)
	  do kkk=1,plev
	     ztd = ztd + Href(kkk,k)*Tmat(l,kkk)
	  end do
	  Amat(l,k) = ztd
	  !       Amat(k,l) = ztd
       end do
    end do

    tmp(:,:)           = TRANSPOSE(Amat(:,:))
    Amat(:,:) = tmp(:,:)

#ifdef DEBUG
    print *," Amat"
    do k=1,plev
       do m=1,plev
	  Imat(k,m)=0.0_r8
	  print *,"Amat(",m,",",k,")=",Amat(m,k)-Amat(k,m)
       end do
       print *
    end do
#endif

    ! =========================================
    ! invert the vertical structure matrix
    ! Amat_inv in Steve Thomas's notation.
    ! =========================================

    do k=1,plev
       do l=1,plev
	  Amat_inv(k,l) = Amat(k,l)
       end do
    end do

#ifdef CAM
    call abortmp('not supported in cam')
#else   

    !  call dgeco(refstate%Amat_inv(1,1),plev,plev,ipvt(1),rcond,z(1))        !LINPACK
    !  call dgedi(refstate%Amat_inv(1,1),plev,plev,ipvt(1),det(1),work(1),01) !LINPACK
    call dgetrf(plev,plev,Amat_inv(1,1),plev,ipvt,info)           !LAPACK
    call dgetri(plev,Amat_inv(1,1),plev,ipvt,work,plev,info) !LAPACK
#endif
#ifdef DEBUG
    print *," Amat^-1*Amat"
    do k=1,plev
       do m=1,plev
	  Imat(m,k)=0.0_r8
	  do l=1,plev
	     Imat(m,k) = Imat(m,k) + Amat_inv(m,l)*Amat(l,k)
	  end do
	  print *," Imat(",m,",",k,")=",Imat(m,k)
       end do
    end do
#endif

    ! =======================================================================
    ! compute the right eigenvector matrix (Emat) and its inverse (Emat_inv)
    ! from the vertical structure matrix.
    ! =======================================================================
    Emat=Amat
    info=reigen_solver(plev,Emat,Emat_Inv, zcr)
#ifdef DEBUG
    print *," Emat^T*Emat"
    do k=1,plev
       do m=1,plev
	  Imat(m,k)=0.0_r8
	  do l=1,plev
	     Imat(m,k) = Imat(m,k) + Emat_inv(l,m)*Emat(l,k)
	  end do
	  print *," Imat(",m,",",k,")=",Imat(m,k)
       end do
    end do
#endif
    
    do k=1,plev
       do l=1,plev
	  Emat_inv(k,l) = Emat(k,l)
       end do
    end do
#ifdef CAM
    call abortmp('not supported in cam')
#else   
    ! call dgeco(Emat_inv(1,1),plev,plev,ipvt(1),rcond,z(1))
    ! call dgedi(Emat_inv(1,1),plev,plev,ipvt(1),det(1),work(1),01)
    call dgetrf(plev,plev,Emat_inv(1,1),plev,ipvt,info)           !LAPACK
    call dgetri(plev,Emat_inv(1,1),plev,ipvt,work,plev,info) !LAPACK
#endif
    tmp(:,:)      = TRANSPOSE(Emat(:,:))
    Emat(:,:)     = tmp(:,:)
    tmp(:,:)      = TRANSPOSE(Amat(:,:))
    Amat(:,:)     = tmp(:,:)

    tmp(:,:)      = TRANSPOSE(Emat_inv(:,:))
    Emat_inv(:,:) = tmp(:,:)

#ifdef DEBUG
    print *," Emat*Emat_inv"
    do k=1,plev
       do m=1,plev
	  Imat(k,m)=0.0_r8
	  do l=1,plev
	     !          Imat(k,m) = Imat(k,m) + Emat(k,l)*Emat_inv(l,m)
	     Imat(k,m) = Imat(k,m) + Emat(l,k)*Emat_inv(m,l)
	  end do
	  print *," Imat(",k,",",m,")=",Imat(k,m)
       end do
    end do
#endif

    ! =======================================================================
    ! Compute and print gravity wave equivalent depths and phase speeds
    ! co
    ! =======================================================================

    do k=1,plev
       zci(k) = sign(1.0_r8,zcr(k))*sqrt(abs(zcr(k)))
       Lambda(k) = zcr(k)*dt*dt               ! solver eigenvalues
       !    Lambda(k) = zcr(k)                     ! solver eigenvalues
    end do

    if (masterproc) then
       write(iulog,910) (Tref(k),k=1,plev)
       write(iulog,920) (zci(k),k=1,plev)
       write(iulog,930) (zcr(k)/g,k=1,plev)
    end if

    ! ========================================================================
    ! Test for unstable mean temperatures (negative phase speed and eqivalent
    ! depth) for at least one gravity wave.
    ! ========================================================================

    nneg = 0
    do k=1,plev
       if (zcr(k)<=0.0_r8) nneg = nneg + 1
    end do

    if (nneg/=0) then
       write(iulog,*)'---------------------------------------------------'
       write(iulog,*)'SET_VERT_STRUCT_MAT: UNSTABLE MEAN TEMPERATURE. STOP',zcr
       call abortmp('prim_si_ref_mod')
    end if

910 format(' REFERENCE TEMPERATURES FOR SEMI-IMPLICIT SCHEME = ',  /(1x,12f9.3))
920 format(' GRAVITY WAVE PHASE SPEEDS (M/S) FOR MEAN STATE = '    /(1x,12f9.3))
930 format(' GRAVITY WAVE EQUIVALENT DEPTHS (M) FOR MEAN STATE = ' /(1x,12f9.3))

9910 format( 1x,i3,13f9.5)
9920 format(/,      13i9)

  end subroutine set_vert_struct_mat


  subroutine prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)
  use kinds, only : real_kind
  use control_mod, only : initial_total_mass
  use physical_constants, only : g
  use element_mod, only : element_t
  use time_mod, only : timelevel_t 
  use hybvcoord_mod, only : hvcoord_t 
  use hybrid_mod, only : hybrid_t
  use dimensions_mod, only : np
  use global_norms_mod, only : global_integral 

  type (element_t), intent(inout) :: elem(:)
  type (TimeLevel_t), target, intent(in) :: tl
  type (hybrid_t),intent(in)     :: hybrid
  type (hvcoord_t), intent(in)   :: hvcoord
  integer,intent(in)             :: nets,nete
  
  ! local 
  real (kind=real_kind)  :: tmp(np,np,nets:nete)
  real (kind=real_kind)  :: scale,mass0
  integer :: n0,nm1,np1,ie

  if (initial_total_mass == 0) return;
  
  n0=tl%n0
  nm1=tl%nm1
  np1=tl%np1
  
  scale=1/g                                  ! assume code is using Pa
  if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
  ! after scaling, Energy is in J/m**2,  Mass kg/m**2
  
  do ie=nets,nete
     tmp(:,:,ie)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  mass0 = global_integral(elem, tmp(:,:,nets:nete),hybrid,np,nets,nete)
  mass0 = mass0*scale;  
  
  do ie=nets,nete
     elem(ie)%state%ps_v(:,:,n0)=elem(ie)%state%ps_v(:,:,n0)*(initial_total_mass/mass0)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,n0)
     elem(ie)%state%ps_v(:,:,nm1)=elem(ie)%state%ps_v(:,:,n0)
  enddo
  if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
     write (*,'(a,e24.15)') "Initializing Total Mass (kg/m^2) = ",initial_total_mass
  endif
  end subroutine prim_set_mass


end module prim_si_ref_mod
