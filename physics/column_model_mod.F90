#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
#define EMANUELSTATE
!#define FORCINGSTAT
!#define USE_MAXLAT
#define _DBG_ !DBG

!============================================
!
! For help or information:
! 
! Amik St-Cyr
! 
! e-mail: amik@ucar.edu
!
!============================================

module column_model_mod
#ifdef _PRIM
  use element_mod,     only : element_t
  use hybvcoord_mod,   only : hvcoord_t 
  use hybrid_mod,      only : hybrid_t
  use kinds,           only : real_kind, int_kind
  use time_mod,        only : TimeLevel_t
  use physics_mod,     only : elem_physics_t, Specific_Humidity, Saturation_Specific_Humidity, getsurfpress, Temp2PotTemp
  use dimensions_mod,  only : nlev, nlevp, np, qsize, nelemd
  use control_mod,     only : integration, columnpackage, test_case,  &
                              accumfreq, statefreq, &
                              TRACERADV_TOTAL_DIVERGENCE, &
                              tracer_advection_formulation
  use convect43c,      only : convect
  use held_suarez_mod, only : hs_forcing
  use forcing_mod,     only : Apply_Forcing,EXP_EULER,INTERP_BDF2
  use reduction_mod,   only : parallelmax,parallelmin

  use aquaplanet,      only : aquaplanet_forcing, isrf_forc
  use multicloud_mod,  only : verticalprojectionvector,verticalprojection,verticalprojection1D,&
       projvanddq,lambdaswitch,InitShearDamping,ApplyShearDamping,VerticalIntegral1D

  use parallel_mod, only : abortmp, global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use physics_types_mod, only : pelem
  use physical_constants, only : Cp, DD_PI, p0, Lc, Cp
  use gravity_wave_drag_mod, only : gravity_wave_drag_forcing
  use global_norms_mod,      only : global_integral
  use eigenmodes_mod,        only : phi1,phi2
  use column_types_mod,      only : ColumnModelMultiCloud_t, ColumnModelEmanuel_t,HeldSuarezForcing_t,&
       AquaplanetForcing_t,ColumnModelMulticloud_t,ColumnDataEmanuel_t, ColumnModel_t
  
  implicit none

  private
#ifdef USE_MAXLAT
  real(kind=real_kind), parameter :: emanuel_max_lat=1.40
#endif

  real (kind=real_kind), public :: mc_Q0R1
  real (kind=real_kind), public :: mc_TstarMinTeb
  real (kind=real_kind), public :: mc_TebMinTem
  real (kind=real_kind), public :: mc_a1
  real (kind=real_kind), public :: mc_a2
  real (kind=real_kind), public :: mc_a0
  real (kind=real_kind), public :: mc_a0prime
  real (kind=real_kind), public :: mc_gamma2
  real (kind=real_kind), public :: mc_mu
  real (kind=real_kind), public :: mc_alpha_c
  real (kind=real_kind), public :: mc_alpha_s
  real (kind=real_kind), public :: mc_tau_c
  real (kind=real_kind), public :: mc_tau_s
  real (kind=real_kind), public :: mc_tau_conv
  real (kind=real_kind), public :: mc_gamma2prime
  integer,               public :: mc_closure 
  real (kind=real_kind), public :: mc_relaxation
  integer,               public :: mc_shutconvec
  character(len=80)    , public :: mc_mask
  real (kind=real_kind), public :: mc_mask_k
  real (kind=real_kind), public :: mc_mask_pos
  character(len=80)    , public :: mc_surf
  real (kind=real_kind), public :: mc_csr

  public :: InitColumnModel
  public :: ApplyColumnModel


! used for physics_state
  real(kind=real_kind), allocatable :: precip_min(:), precip_max(:), precip_sum(:)

  real(kind=real_kind), allocatable :: tp_min(:), tp_max(:), tp_sum(:)
  real(kind=real_kind), allocatable :: qp_min(:), qp_max(:), qp_sum(:)
  real(kind=real_kind), allocatable :: wd_min(:), wd_max(:), wd_sum(:)
  real(kind=real_kind), allocatable :: cbmf_min(:), cbmf_max(:), cbmf_sum(:)
  real(kind=real_kind), allocatable :: pa_min(:), pa_max(:), pa_sum(:)
#ifdef FORCINGSTAT
  real(kind=real_kind), allocatable :: ft_min(:), ft_max(:), ft_sum(:)
  real(kind=real_kind), allocatable :: fu_min(:), fu_max(:), fu_sum(:)
  real(kind=real_kind), allocatable :: fv_min(:), fv_max(:), fv_sum(:)
  real(kind=real_kind), allocatable :: fq_min(:), fq_max(:), fq_sum(:)
#endif
contains

  subroutine InitColumnModel(elem, elem_physics, cm,hvcoord,hybrid,tl,nets,nete,runtype)
    type(element_t), intent(inout) :: elem(:)
    type(elem_physics_t), intent(inout) :: elem_physics(:)
    type (ColumnModel_t) :: cm
    type (hvcoord_t), intent(in), target     :: hvcoord
    type (TimeLevel_t), intent(in), target   :: tl
    type (hybrid_t), intent(in),target       :: hybrid
    integer, intent(in), target              :: nets,nete,runtype

    integer :: i,j,k,ie
    if(runtype.ne.1) then
    ! for all column type
       do ie=nets,nete
          do k=1,nlev
             do j=1,np	
                do i=1,np	
                   elem(ie)%state%Q(i,j,k,1,1:3)      = 0_real_kind
                   elem(ie)%state%Qdp(i,j,k,1,1:3)      = 0_real_kind
                   elem(ie)%derived%FQ(i,j,k,1,1:3)     = 0_real_kind
                   elem(ie)%derived%FM(i,j,1:2,k,1:3) = 0_real_kind
                   elem(ie)%derived%FT(i,j,k,1:3)     = 0_real_kind
                enddo
             enddo
          enddo
       enddo
    endif	
    cm%hvcoord = hvcoord
    cm%hybrid => hybrid
    cm%tl=> tl
    cm%nets = nets
    cm%nete = nete

    ! Special features are going here

    if(columnpackage == "emanuel")then
       call InitColumnModelEmanuel(cm%cm_em,nets,nete)
    endif

    !if(columnpackage == "multicloud")then
    !   call InitColumnModelEmanuel(cm%cm_em,nets,nete) ! ASC gets 30-90 N and 30 to 90 S
    !endif


    if(columnpackage == "none" .AND. test_case(1:12) == "held_suarez0")then
       call InitHeldSuarezForcing(cm%cm_hs)
    elseif(test_case(1:10) == "aquaplanet")then
       call InitAquaplanetForcing(cm%cm_aq)
    endif

    ! The MC needs the sounding
    
    if(columnpackage == "multicloud")then
       call InitColumnModelMulticloud(elem,elem_physics, cm%hybrid,cm%hvcoord,cm%cm_mc,nets,nete)!ASC gets 30S-30N
       call InitShearDamping(elem,cm%hybrid,nets,nete)
    endif

  end subroutine InitColumnModel

  subroutine ApplyColumnModel(elem, elem_physics, hybrid, hvcoord, cm,dt)
    use hybvcoord_mod, only : hvcoord_t

    type (element_t), intent(inout) :: elem(:)
    type(elem_physics_t), intent(inout) :: elem_physics(:)
    type (ColumnModel_t),intent(inout) :: cm
    real (kind=real_kind),intent(in)   :: dt
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in)       :: hybrid

    integer :: i,j,k,ie, nm1, forcing=0
    real (kind=real_kind)              :: tau_damp


    nm1 = cm%tl%nm1

    ! forcings are accumulated then applied

    do ie=cm%nets,cm%nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                elem(ie)%derived%FQ(i,j,k,:,nm1) = 0_real_kind
                elem(ie)%derived%FM(i,j,1:2,k,nm1)     = 0_real_kind
                elem(ie)%derived%FT(i,j,k,nm1)         = 0_real_kind
             enddo
          enddo
       enddo
    enddo

    if(columnpackage == "emanuel")then
      forcing=1
      call ApplyColumnModelEmanuel(elem, hybrid,cm%cm_em,dt,cm%hvcoord,cm%tl,cm%nets,cm%nete)
    endif
    if(test_case(1:12) == "held_suarez0")then
      forcing=1
      call ApplyHeldSuarezForcing(elem, cm%cm_hs,dt, cm%hvcoord,cm%tl,cm%nets,cm%nete)
    endif
    if(test_case(1:10) == "aquaplanet")then
       forcing=1
       if(columnpackage == "multicloud")then
          call ApplyAquaplanetForcing(elem, elem_physics, hybrid,cm%cm_aq,dt, cm%hvcoord,cm%tl,cm%nets,cm%nete,cm%cm_mc)          
          if(cm%cm_mc%relaxation > 0.0D0)then
             tau_damp = 1.0D0/cm%cm_mc%relaxation
             call ApplyShearDamping(elem, elem_physics, tau_damp,hybrid,cm%tl,cm%nets,cm%nete)
          end if
       else
          call ApplyAquaplanetForcing(elem, elem_physics, hybrid,cm%cm_aq,dt, cm%hvcoord,cm%tl,cm%nets,cm%nete)
       end if
    endif

    if (forcing==1) then
       do ie=cm%nets,cm%nete
          call Apply_Forcing(EXP_EULER,elem(ie), hvcoord, cm%tl,dt)
       end do
    endif

  end subroutine ApplyColumnModel





  ! The rest of the file is private


  subroutine InitAquaplanetForcing(cm)

    type (AquaplanetForcing_t)  :: cm

    if(.NOT.cm%INIT)then
       cm%INIT = .TRUE.
    endif

  end subroutine InitAquaplanetForcing
   
  subroutine InitHeldSuarezForcing(cm)

    type (HeldSuarezForcing_t) :: cm

    if(.NOT.cm%INIT)then
       cm%INIT = .TRUE.
    endif

  end subroutine InitHeldSuarezForcing

  subroutine InitColumnModelEmanuel(cm,nets,nete)

    type (ColumnModelEmanuel_t) :: cm
    integer,target              :: nets,nete

    integer :: i,j,k,ie

    if(.NOT.cm%INIT)then
       allocate(pelem(nets:nete))
       do ie=nets,nete
         do j=1,np
           do i=1,np
             pelem(ie)%state%CBMF(i,j)=0_real_kind
             pelem(ie)%surfc%precip(i,j)=0_real_kind
             pelem(ie)%surfc%wd(    i,j)=0_real_kind
             pelem(ie)%surfc%tprime(i,j)=0_real_kind
             pelem(ie)%surfc%qprime(i,j)=0_real_kind
             pelem(ie)%accum%precip(i,j)=0_real_kind
           enddo
         enddo
       enddo
       cm%INIT = .TRUE.
    endif

  end subroutine InitColumnModelEmanuel

  subroutine ConvertStateToColumnEmanuel(cm,elemin, tl, hvcoord)
    type (ColumnModelEmanuel_t) :: cm
    type (element_t)            :: elemin
    type (hvcoord_t)     :: hvcoord
    type (TimeLevel_t)  :: tl
    
    ! local
    
    real (kind=real_kind) :: ps(np,np),r,T,v1,v2
    real (kind=real_kind) :: hyam_ps0,hyai_ps0

    integer :: i,j,k

    integer :: nm1

    nm1 = tl%nm1

       
    do j=1,np
       do i=1,np
          ps(i,j) = EXP(elemin%state%lnps(i,j,nm1))
       enddo
    enddo
    if(integration == "explicit")then
       do k=1,nlev
!         hyam_ps0=hvcoord%hyam(k  )*hvcoord%ps0
!         hyai_ps0=hvcoord%hyai(k+1)*hvcoord%ps0
          do j=1,np
             do i=1,np
                T                          = elemin%state%T(i,j,k,nm1)
                cm%col(i,j)%T(nlevp-k)     = T
                cm%col(i,j)%P(nlevp-k)     = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
                cm%col(i,j)%PH(nlevp-k+1)  = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps(i,j)

                cm%col(i,j)%Q(nlevp-k)     = elemin%state%Q(i,j,k,1,nm1)

                cm%col(i,j)%QS(nlevp-k)    = Saturation_Specific_Humidity(cm%col(i,j)%P(nlevp-k),T)
                cm%col(i,j)%TRA(nlevp-k,1) = 0_real_kind 

                ! Project covariant vel. onto sphere

!               v1                         = elemin%state%v(i,j,1,k,nm1)
!               v2                         = elemin%state%v(i,j,2,k,nm1)
!               cm%col(i,j)%U(nlevp-k)     = v1*elemin%Dinv(1,1,i,j) + v2*elemin%Dinv(2,1,i,j)
!               cm%col(i,j)%V(nlevp-k)     = v1*elemin%Dinv(1,2,i,j) + v2*elemin%Dinv(2,2,i,j)                

                cm%col(i,j)%U(nlevp-k)     = elemin%state%v(i,j,1,k,nm1)
                cm%col(i,j)%V(nlevp-k)     = elemin%state%v(i,j,2,k,nm1)
             enddo
          enddo
          do j=1,np
             do i=1,np
                cm%col(i,j)%PH(1)  = hvcoord%hyai(nlevp)*hvcoord%ps0 + hvcoord%hybi(nlevp)*ps(i,j)
             enddo
          enddo
       enddo
    else
       do k=1,nlev
          hyam_ps0=hvcoord%hyam(k  )*hvcoord%ps0
          hyai_ps0=hvcoord%hyai(k+1)*hvcoord%ps0
          do j=1,np
             do i=1,np
#ifdef USE_MAXLAT
	     if(abs(elemin%sphereP(i,j)%lat).lt. emanuel_max_lat) then
#endif
                T                          = elemin%state%T(i,j,k,nm1)
                cm%col(i,j)%T(nlevp-k)     = T
                cm%col(i,j)%P(nlevp-k)     = hyam_ps0 + hvcoord%hybm(k  )*ps(i,j)
                cm%col(i,j)%PH(nlevp-k)    = hyai_ps0 + hvcoord%hybi(k+1)*ps(i,j)
                cm%col(i,j)%Q(nlevp-k)     = Specific_Humidity(elemin%state%Q(i,j,k,1,nm1))
                cm%col(i,j)%QS(nlevp-k)    = Saturation_Specific_Humidity(cm%col(i,j)%P(nlevp-k),T)
                cm%col(i,j)%TRA(nlevp-k,1) = 0_real_kind 

                ! Project cotravariant vel. onto sphere

                v1                         = elemin%state%v(i,j,1,k,nm1)
                v2                         = elemin%state%v(i,j,2,k,nm1)
                cm%col(i,j)%U(nlevp-k)     = v1*elemin%D(1,1,i,j) + v2*elemin%D(1,2,i,j)            
                cm%col(i,j)%V(nlevp-k)     = v1*elemin%D(2,1,i,j) + v2*elemin%D(2,2,i,j)                
#ifdef USE_MAXLAT
		endif
#endif
             enddo
          enddo
       enddo
    endif

  end subroutine ConvertStateToColumnEmanuel

  subroutine ConvertColumnToStateEmanuel(cm,elemin, tl)

    type (ColumnModelEmanuel_t),intent(in), target :: cm
    type (element_t),intent(inout)         :: elemin
    type (TimeLevel_t)  :: tl
    ! local
    type (ColumnDataEmanuel_t),pointer :: CD

    real (kind=real_kind)                  :: u1,v1,fu,fv,RCD
    integer                                :: k,nm1

    integer                     :: i,j

    nm1 = tl%nm1

    if(integration == "explicit" )then
       if(isrf_forc.ne.1) then
       do k=1,nlev
          do j=1,np
             do i=1,np
#ifdef USE_MAXLAT
    	     if(abs(elemin%sphereP(i,j)%lat).lt. emanuel_max_lat) then
#endif
                cd => cm%col(i,j)            ! Project back onto cubed sphere
!               u1 = cd%U( nlevp-k)          ! dry adiabatic adjustment in convect43c
!               v1 = cd%V( nlevp-k)          ! requires update altered U, V, T, Q, QS

!               elemin%state%v(i,j,1,k,nm1)  = u1*elemin%D(1,1,i,j) + v1*elemin%D(2,1,i,j)
!               elemin%state%v(i,j,2,k,nm1)  = u1*elemin%D(1,2,i,j) + v1*elemin%D(2,2,i,j)

                elemin%state%v(i,j,1,k,nm1)  = cd%U( nlevp-k)
                elemin%state%v(i,j,2,k,nm1)  = cd%V( nlevp-k)

                elemin%state%T(i,j,k,nm1)    = cd%T(nlevp-k)
                elemin%state%Q(i,j,k,1,nm1)    = cd%Q(nlevp-k)
#ifdef USE_MAXLAT
             endif
#endif
             enddo
          enddo
       enddo
       endif

       do k=1,nlev
          do j=1,np
             do i=1,np
#ifdef USE_MAXLAT
    	     if(abs(elemin%sphereP(i,j)%lat).lt. emanuel_max_lat) then
#endif
                cd => cm%col(i,j)            ! Project back onto cubed sphere

!               fu = cd%FU(nlevp-k)
!               fv = cd%FV(nlevp-k)
!               RCD = (1_real_kind+elemin%state%Q(i,j,k,nm1))**2_real_kind

!               elemin%derived%FM(i,j,1,k,nm1) = &
!                 elemin%derived%FM(i,j,1,k,nm1) + fu*elemin%D(1,1,i,j) + fv*elemin%D(2,1,i,j)
!               elemin%derived%FM(i,j,2,k,nm1) = &
!                 elemin%derived%FM(i,j,2,k,nm1) + fu*elemin%D(1,2,i,j) + fv*elemin%D(2,2,i,j)

                elemin%derived%FM(i,j,1,k,nm1) = cd%FU(nlevp-k)
                elemin%derived%FM(i,j,2,k,nm1) = cd%FV(nlevp-k)

                elemin%derived%FT(i,j,k,nm1)   = &
                  elemin%derived%FT(i,j,k,nm1)   +     cd%FT(nlevp-k)
                elemin%derived%FQ(i,j,k,nm1,1) = &
                  elemin%derived%FQ(i,j,k,nm1,1) + cd%FQ(nlevp-k)

#ifdef USE_MAXLAT
             endif
#endif
             enddo
          enddo
       enddo
    else
       if(isrf_forc.ne.1) then
       do k=1,nlev
          do j=1,np
             do i=1,np
#ifdef USE_MAXLAT
    	     if(abs(elemin%sphereP(i,j)%lat).lt. emanuel_max_lat) then
#endif
                cd => cm%col(i,j)            ! Project back onto cubed sphere
                u1 = cd%U( nlevp-k)          ! dry adiabatic adjustment in convect43c
                v1 = cd%V( nlevp-k)          ! requires update altered U, V, T, Q, QS

                elemin%state%v(i,j,1,k,nm1)  = u1*elemin%Dinv(1,1,i,j) + v1*elemin%Dinv(1,2,i,j)
                elemin%state%v(i,j,2,k,nm1)  = u1*elemin%Dinv(2,1,i,j) + v1*elemin%Dinv(2,2,i,j)
                elemin%state%T(i,j,k,nm1)    = cd%T(nlevp-k)
                elemin%state%Q(i,j,k,1,nm1)      = cd%Q(nlevp-k)/(1_real_kind - cd%Q(nlevp-k))
#ifdef USE_MAXLAT
             endif
#endif
             enddo
          enddo
       enddo
       endif

       do k=1,nlev
          do j=1,np
             do i=1,np
#ifdef USE_MAXLAT
    	     if(abs(elemin%sphereP(i,j)%lat).lt. emanuel_max_lat) then
#endif
                cd => cm%col(i,j)            ! Project back onto cubed sphere

                fu = cd%FU(nlevp-k)
                fv = cd%FV(nlevp-k)
                RCD = (1_real_kind+elemin%state%Q(i,j,k,1,nm1))**2_real_kind

                elemin%derived%FM(i,j,1,k,nm1) = &
                elemin%derived%FM(i,j,1,k,nm1) + fu*elemin%Dinv(1,1,i,j) + fv*elemin%Dinv(1,2,i,j)
                elemin%derived%FM(i,j,2,k,nm1) = &
                elemin%derived%FM(i,j,2,k,nm1) + fu*elemin%Dinv(2,1,i,j) + fv*elemin%Dinv(2,2,i,j)
                elemin%derived%FT(i,j,k,nm1)   = &
                elemin%derived%FT(i,j,k,nm1)   +     cd%FT(nlevp-k)
                elemin%derived%FQ(i,j,k,nm1,1) = &
                elemin%derived%FQ(i,j,k,nm1,1) + RCD*cd%FQ(nlevp-k)
#ifdef USE_MAXLAT
             endif
#endif
             enddo

          enddo

       enddo
    endif

  end subroutine ConvertColumnToStateEmanuel

  !==================================
  !Initializes  Qt_0 Qt_1 Qt_2
  ! also inits the convex sum
  !==================================
  !
  ! pour t1 et t2 enlever les soundings (aussi dans multicloud_mod)
  !
  subroutine InitColumnModelMulticloud(elem, elem_physics, hybrid,hvcoord,cm,nets,nete)

    type(element_t), intent(inout)  :: elem(:)
    type(elem_physics_t), intent(inout)  :: elem_physics(:)
    type (hybrid_t), intent(in)     :: hybrid
    type (hvcoord_t), intent(inout) :: hvcoord
    type (ColumnModelMultiCloud_t)  :: cm
    integer,target                  :: nets,nete

    integer :: i,j,k,ie,n0
    real (kind=real_kind)    ::  lat, lon
    real (kind=real_kind)    ::  pot(np,np,nlev)
    real (kind=real_kind)    ::  psi1trunc(nlev)
    real (kind=real_kind)    ::  psi2trunc(nlev)
    real (kind=real_kind)    ::  OneMinusPonPs(np,np,nlev)
    real (kind=real_kind)    ::  p(np,np,nlev)
    real (kind=real_kind)    ::  ps(np,np)

    ! If intial sounding is not uniform this takes care of it
    real (kind=real_kind)    ::  Qt0glob(np,np,nets:nete)
    real (kind=real_kind)    ::  Tmask(np,np,nets:nete)
    real (kind=real_kind)    ::  Qt1glob(np,np,nets:nete)
    real (kind=real_kind)    ::  Qt2glob(np,np,nets:nete)
    real (kind=real_kind)    ::  max_ps,Pbar,psi1avg,psi2avg,lsw,LcOnCp
    real (kind=real_kind)    ::  psi1avgtrunc,psi2avgtrunc
    real (kind=real_kind)    ::  phi1norm,phi2norm

    real (kind=real_kind)    ::  A(2,2),up1(2),up2(2),Adetphi,Adetpsi

    logical                  :: Debug=.FALSE.

    n0=1
    max_ps   = 0.0D0    

    ! =======================================================
    ! Known or "namelist" part of the initialization 
    ! =======================================================

    cm%alpha = 0.1D0


    
    ! exp((1 - 1./cos(x).^2).^5)
    
    ! ============================================
    ! Init the mask (keep the scheme in the tropic)
    ! exp((1 - 1./cos(x).^2).^5)
    ! ============================================

    if(mc_mask == "none")then
       do ie=nets,nete          
          do j=1,np
             do i=1,np
                lat  =  elem(ie)%sphereP(i,j)%lat
                elem_physics(ie)%mask(i,j) = 1.0D0
                Tmask(i,j,ie) = 60.0D0*(1.0D0 - cos(lat))
                elem_physics(ie)%invmask(i,j) = 1.0D0 - elem_physics(ie)%mask(i,j)
             enddo
          enddo
       enddo
       if(hybrid%par%masterproc)then
          print *,"PLOTTING 1 deg MASK FUNCTION from -pi/2 to pi/2:"
          do ie=0,180         
             lat = 3.14159D0*(real(ie,kind=real_kind)/180.0D0 - 0.5D0) 
             print *,"deg = ",180.0D0*lat/3.14159D0," H = ", 1.0D0
          enddo
       end if     
    end if

    if(mc_mask == "cos2")then
       do ie=nets,nete          
          do j=1,np
             do i=1,np
                lat  =  elem(ie)%sphereP(i,j)%lat
                elem_physics(ie)%mask(i,j) = cos(lat)*cos(lat)
                Tmask(i,j,ie) = 60.0D0*(1.0D0 - cos(lat))
                elem_physics(ie)%invmask(i,j) = 1.0D0 - elem_physics(ie)%mask(i,j)
             enddo
          enddo
       enddo
       if(hybrid%par%masterproc)then
          print *,"PLOTTING 1 deg MASK FUNCTION from -pi/2 to pi/2:"
          do ie=0,180         
             lat = 3.14159D0*(real(ie,kind=real_kind)/180.0D0 - 0.5D0) 
             print *,"deg = ",180.0D0*lat/3.14159D0," H = ", cos(lat)*cos(lat)
          end do
       end if     
    end if

    if(mc_mask == "heaviside")then
       do ie=nets,nete          
          do j=1,np
             do i=1,np
                lat  =  elem(ie)%sphereP(i,j)%lat
                elem_physics(ie)%mask(i,j) = Heaviside(lat,mc_mask_k,mc_mask_pos)
                Tmask(i,j,ie) = 60.0D0*(1.0D0 - cos(lat))
                elem_physics(ie)%invmask(i,j) = 1.0D0 - elem_physics(ie)%mask(i,j)
             enddo
          enddo
       enddo
       if(hybrid%par%masterproc)then
          print *,"PLOTTING 1 deg MASK FUNCTION from -pi/2 to pi/2:"
          do ie=0,180         
             lat = 3.14159D0*(real(ie,kind=real_kind)/180.0D0 - 0.5D0) 
             print *,"deg = ",180.0D0*lat/3.14159D0," H = ", Heaviside(lat,mc_mask_k,mc_mask_pos)
          end do
       end if     
    end if

    ! ================================================
    ! Constants present in forcings
    
    cm%closure       = mc_closure
    cm%Q0R1          = mc_Q0R1
       
    cm%TstarMinTeb   = mc_TstarMinTeb 
    
    if(mc_surf == "uniform")then
       do ie=nets,nete          
          do j=1,np
             do i=1,np
                lat  =  elem(ie)%sphereP(i,j)%lat
                lon  =  elem(ie)%sphereP(i,j)%lon
                elem_physics(ie)%delthetasurf(i,j) =  cm%TstarMinTeb
             enddo
          enddo
       enddo
    end if

    if(mc_surf == "warm_pool")then
       do ie=nets,nete          
          do j=1,np
             do i=1,np
                lat  =  elem(ie)%sphereP(i,j)%lat
                lon  =  elem(ie)%sphereP(i,j)%lon
                elem_physics(ie)%delthetasurf(i,j) = Warm_pool(lat,lon,cm%TstarMinTeb)
             enddo
          enddo
       enddo
    end if

    
    ! Variable parameters t=0 init 
    cm%TebMinTem     = mc_TebMinTem
    cm%Tau_convec    = mc_tau_conv
    cm%Tau_strat     = mc_tau_s
    cm%Tau_congest   = mc_tau_c
    cm%alpha_strat   = mc_alpha_s
    cm%alpha_congest = mc_alpha_c
    cm%mu            = mc_mu
    cm%gamma2        = mc_gamma2
    cm%gamma2prime   = mc_gamma2prime
    cm%alpha2        = 0.1D0  ! Contribution of theta2 to theta_em (tem)      

    ! Deep convective parameters

    cm%a1 = mc_a1 ! Coefficient of theta_eb
    cm%a2 = mc_a2 ! Coefficient of moisture 
    cm%a0 = mc_a0 ! Coefficient of theta1
    cm%a0prime = mc_a0prime
    cm%relaxation = mc_relaxation
    cm%shutconvec = 0
    if(mc_shutconvec >0)then
       cm%shutconvec = mc_shutconvec
    end if

    ! Switch
    if(cm%closure == 1)then       
       cm%lambdastar = 0.2D0       
    else
       cm%lambdastar = 0.0D0
    end if

    cm%Tminus        = 10.0D0                 ! Threshold
    cm%Tplus         = 20.0D0                 ! Threshold
    cm%Tau_Radiation = 50.0D0*24.0D0*3600.0D0 ! Temperature
    cm%Tau_Damping   = 75.0D0*24.0D0*3600.0D0 ! Momentum
    cm%CD0           = 0.001D0                ! Turbulent drag coefficient
    cm%U0            = 2.0D0                  ! (m/s) Turbulent velocity scale
    ! was 500 m changed to 1000m to match APE 6/3/2008
    cm%h             = 500.0D0                ! Boundary layer height
    cm%Htall         = 16000.0D0              ! Rest of the atm
    cm%Tau_r         = 1.0D0/(50.0D0*24.0D0*3600.0D0)

    cm%A  =  (1.0D0 - cm%lambdastar)/(cm%Tplus - cm%Tminus)
    cm%B  =  1.0D0 - cm%A*cm%Tplus

    ! ================================================
    ! Used to truncate for QHeating
    ! ================================================

    do k=1,nlev
       cm%vertmask(k) = 1.0D0
    end do

    cm%vertmask(1:5) = 0.0D0

 
    ! =======================================================
    ! What follows depends on the background sounding at t=0 
    ! =======================================================
    ! SOUNDING IS SUPPOSED UNIFORM in x-y

    ! BASIS INTITIALIZATION ...

    ! Initialize basis phi_1, phi_2 and psi_1 psi_2
    do k=1,nlev
       cm%D%phi(k,1)   = phi1(k)
       cm%D%phi(k,2)   = -phi2(k)
       cm%D%cstmode(k) = 1.0D0 
    enddo

    ! Normalize phi_n

    call VerticalProjection1D(cm%D%phi(:,1),phi1norm,cm%D%phi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%phi(:,2),phi2norm,cm%D%phi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)

    do k=1,nlev
       cm%D%phi(k,1)   = cm%D%phi(k,1)/dsqrt(phi1norm)
       cm%D%phi(k,2)   = cm%D%phi(k,2)/dsqrt(phi2norm)
    enddo
    
    ! Verify orthogonality with respect to our inner product

    call VerticalProjection1D(cm%D%phi(:,1),phi1norm,cm%D%phi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%phi(:,2),phi2norm,cm%D%phi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)

    if(hybrid%par%masterproc)then
       print *,"--------------------------------------------------------"
       print *,"(phi_1,phi_2) = ", phi1norm, ", (phi_2,phi_1) = ", phi2norm
    endif

    ! Create psi_n
    
    call VerticalIntegral1D(cm%D%phi(:,1),cm%D%psi(:,1),cm%D%cstmode,elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalIntegral1D(cm%D%phi(:,2),cm%D%psi(:,2),cm%D%cstmode,elem(nets)%state%lnps(1,1,n0),hvcoord)

    ! Remove average

    do k=1,nlev
       cm%D%psi(k,1)   = cm%D%psi(k,1) - cm%D%psi(nlev,1)
       cm%D%psi(k,2)   = cm%D%psi(k,2) - cm%D%psi(nlev,2)
    enddo


    ! Looks like to have O(1) norms
    ! 10/28/2009

    do k=1,nlev
       cm%D%phi(k,1)   = cm%D%phi(k,1)*6.25D0
       cm%D%phi(k,2)   = cm%D%phi(k,2)*12.5D0
       cm%D%psi(k,1)   = cm%D%psi(k,1)*6.25D0
       cm%D%psi(k,2)   = cm%D%psi(k,2)*12.5D0
    enddo



    ! Verify orthogonality with respect to our inner product (should - not - be orthogonal)

    call VerticalProjection1D(cm%D%psi(:,1),phi1norm,cm%D%psi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%psi(:,2),phi2norm,cm%D%psi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)

    if(hybrid%par%masterproc)then
       print *,"--------------------------------------------------------"
       print *,"(psi_1,psi_2) = ", phi1norm, ", (psi_2,psi_1) = ", phi2norm
    endif

    call VerticalProjection1D(cm%D%phi(:,1),A(1,1),cm%D%phi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%phi(:,1),A(1,2),cm%D%phi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%phi(:,2),A(2,1),cm%D%phi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%phi(:,2),A(2,2),cm%D%phi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)
    
    Adetphi = A(1,1)*A(2,2) - A(2,1)*A(1,2)
    cm%D%Aphi(1,1) = A(2,2)/Adetphi
    cm%D%Aphi(1,2) =-A(1,2)/Adetphi
    cm%D%Aphi(2,1) =-A(2,1)/Adetphi
    cm%D%Aphi(2,2) = A(1,1)/Adetphi

    if(hybrid%par%masterproc)then
       print *," Reproj for A (phi)"
       print *,"A11 A12", A(1,1)," ",A(1,2) 
       print *,"A21 A22", A(2,1)," ",A(2,2) 
       print *," Reproj^-1 for A (phi)"
       print *,"A11 A12", cm%D%Aphi(1,1)," ",cm%D%Aphi(1,2) 
       print *,"A21 A22", cm%D%Aphi(2,1)," ",cm%D%Aphi(2,2) 
       print *," Reproj^-1 * Reproj for A (phi)"
       print *,"1 0", cm%D%Aphi(1,1)*A(1,1)+cm%D%Aphi(1,2)*A(2,1)," ",cm%D%Aphi(1,1)*A(1,2)+cm%D%Aphi(1,2)*A(2,2) 
       print *,"0 1", cm%D%Aphi(2,1)*A(1,1)+cm%D%Aphi(2,2)*A(2,1)," ",cm%D%Aphi(2,1)*A(1,2)+cm%D%Aphi(2,2)*A(2,2) 
    endif

    call VerticalProjection1D(cm%D%psi(:,1),A(1,1),cm%D%psi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%psi(:,1),A(1,2),cm%D%psi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%psi(:,2),A(2,1),cm%D%psi(:,1),elem(nets)%state%lnps(1,1,n0),hvcoord)
    call VerticalProjection1D(cm%D%psi(:,2),A(2,2),cm%D%psi(:,2),elem(nets)%state%lnps(1,1,n0),hvcoord)

    Adetpsi = A(1,1)*A(2,2) - A(2,1)*A(1,2)

    cm%D%Apsi(1,1) = A(2,2)/Adetpsi
    cm%D%Apsi(1,2) =-A(1,2)/Adetpsi
    cm%D%Apsi(2,1) =-A(2,1)/Adetpsi
    cm%D%Apsi(2,2) = A(1,1)/Adetpsi
    
    if(hybrid%par%masterproc)then
       print *," Reproj for A (psi)"
       print *,"A11 A12", A(1,1)," ",A(1,2) 
       print *,"A21 A22", A(2,1)," ",A(2,2) 
       print *," Reproj^-1 for A (psi)"
       print *,"A11 A12", cm%D%Apsi(1,1)," ",cm%D%Apsi(1,2) 
       print *,"A21 A22", cm%D%Apsi(2,1)," ",cm%D%Apsi(2,2) 
       print *," Reproj^-1 * Reproj for A (psi)"
       print *,"1 0", cm%D%Apsi(1,1)*A(1,1)+cm%D%Apsi(1,2)*A(2,1)," ",cm%D%Apsi(1,1)*A(1,2)+cm%D%Apsi(1,2)*A(2,2) 
       print *,"0 1", cm%D%Apsi(2,1)*A(1,1)+cm%D%Apsi(2,2)*A(2,1)," ",cm%D%Apsi(2,1)*A(1,2)+cm%D%Apsi(2,2)*A(2,2) 
    endif

    do ie=nets,nete
       ! useful only if v at t=0 is not zero ...
       ! V1
       call VerticalProjectionVector(elem(ie)%state%v(:,:,:,:,n0),elem_physics(ie)%uproj1(:,:,:,n0),&
            cm%D%phi(:,1),elem(ie)%state%lnps(:,:,n0),hvcoord)
       ! V2
       call VerticalProjectionVector(elem(ie)%state%v(:,:,:,:,n0),elem_physics(ie)%uproj2(:,:,:,n0),&
            cm%D%phi(:,2),elem(ie)%state%lnps(:,:,n0),hvcoord)
       ! Vbar
       call VerticalProjectionVector(elem(ie)%state%v(:,:,:,:,n0),elem_physics(ie)%ubar(:,:,:,n0),&
            cm%D%cstmode(:),elem(ie)%state%lnps(:,:,n0),hvcoord)

       do j=1,np
          do i=1,np
             up1 = elem_physics(ie)%uproj1(i,j,:,n0)
             up2 = elem_physics(ie)%uproj2(i,j,:,n0)
             elem_physics(ie)%uproj1(i,j,1,n0)=up1(1)*cm%D%Aphi(1,1)+up2(1)*cm%D%Aphi(1,2)
             elem_physics(ie)%uproj1(i,j,2,n0)=up1(2)*cm%D%Aphi(1,1)+up2(2)*cm%D%Aphi(1,2)
             elem_physics(ie)%uproj2(i,j,1,n0)=up1(1)*cm%D%Aphi(2,1)+up2(1)*cm%D%Aphi(2,2)
             elem_physics(ie)%uproj2(i,j,2,n0)=up1(2)*cm%D%Aphi(2,1)+up2(2)*cm%D%Aphi(2,2)
          enddo
       enddo

       if(Debug) print *,ie,":AVERAGE = ",elem_physics(ie)%qmc(:,:,n0), " ===> Q = ", elem(ie)%state%Q(:,:,:,1,n0)

       ! copy other time-slabs

       ! This is needed to compute potential temperature

       ps(:,:) = getsurfpress(elem(ie)%state%lnps(:,:,n0))

       do j=1,np
          do i=1,np                   
             max_ps = max(max_ps,ps(i,j))
          enddo
       enddo

       do k=1,nlev
          do j=1,np
             do i=1,np
                p(i,j,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             enddo
          enddo
       enddo

       if(hybrid%par%masterproc)then
          print *,"   k  p	          phi1	    phi2      psi1      psi2"
          print *,"--------------------------------------------------------"
          do k=1,nlev
             print 314,k,p(1,1,k),cm%D%phi(k,1),cm%D%phi(k,2),cm%D%psi(k,1),cm%D%psi(k,2)
          enddo
314       format(1x,i3,5f10.4)
       endif
    
       ! Find Q0,Q1 and Q2
       ! WE SUPPOSE t=0 and that the soundings are 
       ! loaded in the the state variables

       ! Stores potential temp at t=0 ... 
   
       pot(:,:,:) = Temp2PotTemp(elem(ie)%state%T(:,:,:,n0),p(:,:,:))
       elem_physics(ie)%pot0(:,:,:) = pot(:,:,:) ! this is to do pot -pot(t=0) later
       do k=1,nlev
          OneMinusPonPs(:,:,k) = 1.0D0 - p(:,:,k)/ps(:,:)
       enddo

       ! Q0 global

       if(Debug)then
          if(hybrid%par%masterproc)then
             print *,"Q=",elem(ie)%state%Q(:,:,:,1,n0)
             print *,"1-p/ps=",OneMinusPonPs
             print *,"lnps=",elem(ie)%state%lnps(:,:,n0)
          end if
       end if

       call ProjVandDQ(OneMinusPonPs,elem(ie)%state%Q(:,:,:,1,n0),&
            Qt0glob(:,:,ie),elem(ie)%state%lnps(:,:,n0),hvcoord)

       do k=1,nlev
          do j=1,np
             do i=1,np
                pot(i,j,k) = cm%D%psi(k,1)  
             enddo
          enddo
       enddo

       ! Q1 global
       call ProjVandDQ(pot,elem(ie)%state%Q(:,:,:,1,n0),Qt1glob(:,:,ie),elem(ie)%state%lnps(:,:,n0),hvcoord)

       do k=1,nlev
          do j=1,np
             do i=1,np
                pot(i,j,k) = cm%D%psi(k,2)  
             enddo
          enddo
       enddo
       ! Q2 global
       call ProjVandDQ(pot,elem(ie)%state%Q(:,:,:,1,n0),Qt2glob(:,:,ie),elem(ie)%state%lnps(:,:,n0),hvcoord)
!
!To change the strength of moisture background gradient Qt1, Qt2, etc. coomment/uncomment the command line
! which defines Lc0nCp = (c) Lc/Cp. (c)=1 == original GATE sounding, (c) = 2, GATE sounding is doubled.   
!
       !  LcOnCp = Lc/Cp
       !LcOnCp multiplied by isome factor to  artificially jack up the moisture stratification constants Qti
       !50%higher
       ! LcOnCp = 1.5d0*Lc/Cp 
       !Double  
         LcOnCp = 2.d0*Lc/Cp 

       do j=1,np
          do i=1,np
             Qt0glob(i,j,ie)=Qt0glob(i,j,ie)*LcOnCp
             Qt1glob(i,j,ie)=Qt1glob(i,j,ie)*LcOnCp
             Qt2glob(i,j,ie)=Qt2glob(i,j,ie)*LcOnCp
          enddo
       enddo
    enddo

    max_ps=ParallelMax(max_ps,hybrid)

    cm%D%Qt0 = global_integral(elem,Qt0glob(:,:,nets:nete),hybrid,np,nets,nete)
    cm%D%Qt1 = global_integral(elem,Qt1glob(:,:,nets:nete),hybrid,np,nets,nete)
    cm%D%Qt2 = global_integral(elem,Qt2glob(:,:,nets:nete),hybrid,np,nets,nete)

    if(hybrid%par%masterproc)then
       print *,"---------------!SHOULD BE POSITIVE!---------------------"
       print *,"Qt0 = ", cm%D%Qt0
       print *,"Qt1 = ", cm%D%Qt1
       print *,"Qt2 = ", cm%D%Qt2

       
       print *,"Qt1 => ", Qt1glob(1,1,nets)
       print *,"Qt2 => ", Qt2glob(1,1,nets)
       print *,"--------------------------------------------------------"
    endif

    do k=1,nlev
       cm%D%psitrunc(k,1) = cm%vertmask(k)*cm%D%psi(k,1)
       cm%D%psitrunc(k,2) = cm%vertmask(k)*cm%D%psi(k,2)
    end do

    cm%D%psitrunc(6:7,2) = 0.0D0;

    call VerticalProjection1D(cm%D%psi(:,1),psi1avg,cm%D%cstmode,log(p0),hvcoord)
    call VerticalProjection1D(cm%D%psi(:,2),psi2avg,cm%D%cstmode,log(p0),hvcoord)

    cm%D%psi1avg = psi1avg
    cm%D%psi2avg = psi2avg

    call VerticalProjection1D(cm%D%psitrunc(:,1),psi1avgtrunc,cm%D%cstmode,log(p0),hvcoord)
    call VerticalProjection1D(cm%D%psitrunc(:,2),psi2avgtrunc,cm%D%cstmode,log(p0),hvcoord)

    cm%D%psi1avgtrunc = psi1avgtrunc
    cm%D%psi2avgtrunc = psi2avgtrunc

    lsw     = LambdaSwitch(cm%TebMinTem,cm%Tminus,cm%Tplus,cm%A,cm%B,cm%lambdastar) ! LAMBDA BAR in notes

    if(hybrid%par%masterproc)then
       print *, "LSW   = ", lsw
       print *, "1-LSW = ", 1.0D0-lsw
       print *, "Q0R1  = ", cm%Q0R1
    end if

    cm%D%Q0c  = ((1.0D0 - cm%lambdastar)/(1.0D0 - lsw))*cm%Q0R1
   
    if(hybrid%par%masterproc)then
       print *, "Q0c  = ", cm%D%Q0c
    end if

    if(cm%closure == 1)then
       cm%D%Q0R2 =  cm%Q0R1*(cm%alpha_congest*psi1avgtrunc*(lsw - cm%lambdastar)/(1.0D0 - cm%lambdastar) - cm%alpha_strat)/ &
            (1.0D0 - cm%alpha_congest*psi2avgtrunc*(lsw - cm%lambdastar)/(1.0D0 - cm%lambdastar))
    else
       cm%D%Q0R2 =  cm%Q0R1*(cm%alpha_congest*(lsw-cm%lambdastar)/(1.0D0 - lsw)- cm%alpha_strat)
    end if

    Pbar              = cm%Q0R1*cm%D%psi1avgtrunc  + cm%D%Q0R2*cm%D%psi2avgtrunc ! THIS IS EQUIV TO HTOT AVERAGE IN THE NOTES       
    
    ! Here we use the average relaxation
    cm%D%Tau_evap_Inv = (cm%Htall/cm%h)*Pbar*(1.0D0/cm%TstarMinTeb)

    if(cm%closure == 1)then
       cm%D%m0           = Pbar*cm%Htall/(lsw*(1.0D0 - cm%mu*cm%D%Q0R2/cm%Q0R1)*cm%TebMinTem)
    else
       cm%D%m0           = Pbar*cm%Htall/((1.0D0 - cm%mu*cm%D%Q0R2/cm%Q0R1)*cm%TebMinTem)
    end if
    
    ! ===============================================
    ! scale equations
    ! ===============================================
    
    ! SCALING NOT COMPLETED IN MULTICLOUD_MOD.F90

#if 0
    ! SCALED
    cm%sAlpha = 15.0D0
    cm%sT     = 8.33D0*60.0D0*60.0D0
    cm%sBeta  = cm%sAlpha/cm%sT
#else
    ! NOT SCALED
    cm%sAlpha = 1.0D0
    cm%sT     = 1.0D0
    cm%sBeta  = cm%sAlpha/cm%sT
#endif
    
    cm%D%m0          = cm%D%m0*cm%sT
    cm%Q0R1        = cm%Q0R1/cm%sBeta
    cm%D%Q0c         = cm%D%Q0c/cm%sBeta
    cm%TstarMinTeb = cm%TstarMinTeb/cm%sAlpha
    cm%D%Tau_evap_Inv= cm%D%Tau_evap_Inv*cm%sT
    cm%Tau_Convec  = cm%Tau_Convec/cm%sT
    cm%Tau_strat   = cm%Tau_strat/cm%sT
    cm%Tau_congest = cm%Tau_congest/cm%sT

    if(hybrid%par%masterproc)then
       print *,"MULTICLOUD parameters:" 
       print *,"lws                 = ",lsw
       print *,"Q0R2                = ",cm%D%Q0R2*86400.0D0 
       print *,"Q0C                 = ",cm%D%Q0c
       print *,"m0                  = ",cm%D%m0
       print *,"Htot average (Pbar) = ",Pbar*86400.0D0 
       print *,"psi1_avg            = ",psi1avg
       print *,"psi2_avg            = ",psi2avg
       print *,"psi1_avg truncated  = ",psi1avgtrunc
       print *,"psi2_avg truncated  = ",psi2avgtrunc
       print *,"Tau evap (per hour) = ",1.0D0/(cm%D%Tau_evap_Inv*3600.0D0)
       print *,"Hs                  = ",cm%alpha_strat*cm%Q0R1*86400.0D0 
       print *,"Hc                  = ",(cm%D%Q0R2+cm%alpha_strat*cm%Q0R1)*86400.0D0 
    end if

    do ie=nets,nete      
       ! copy for other time steps
       do j=1,np
          do i=1,np                 
             elem_physics(ie)%ubar(i,j,1:2,2)   = elem_physics(ie)%ubar(i,j,1:2,1) 
             elem_physics(ie)%ubar(i,j,1:,3)    = elem_physics(ie)%ubar(i,j,1:2,1) 
             elem_physics(ie)%uproj1(i,j,1:2,2) = elem_physics(ie)%uproj1(i,j,1:2,1) 
             elem_physics(ie)%uproj1(i,j,1:2,3) = elem_physics(ie)%uproj1(i,j,1:2,1) 
             elem_physics(ie)%uproj2(i,j,1:2,2) = elem_physics(ie)%uproj2(i,j,1:2,1) 
             elem_physics(ie)%uproj2(i,j,1:2,3) = elem_physics(ie)%uproj2(i,j,1:2,1)

             lat  =  elem(ie)%sphereP(i,j)%lat

             elem_physics(ie)%qmc(i,j,1)        = 0.0D0 !FOR TESTING ADVECTION cos(lat)*cos(lat)
             elem_physics(ie)%qmc(i,j,2)        = elem_physics(ie)%qmc(i,j,1) 
             elem_physics(ie)%qmc(i,j,3)        = elem_physics(ie)%qmc(i,j,1) 

             elem_physics(ie)%Hs(i,j,1:3)       = cm%alpha_strat*cm%Q0R1
             elem_physics(ie)%Hc(i,j,1:3)       = cm%D%Q0R2+cm%alpha_strat*cm%Q0R1
             elem_physics(ie)%QHeating(i,j,:)   = 0.0D0
             elem_physics(ie)%teb(i,j,1:3)      = 0.0D0
          enddo
       enddo

       if(Debug)print *,"QMC1 = ",elem_physics(ie)%qmc(:,:,1)
       if(Debug)print *,"QMC2 = ",elem_physics(ie)%qmc(:,:,2)
       if(Debug)print *,"QMC3 = ",elem_physics(ie)%qmc(:,:,3)

    enddo

    if(hybrid%par%masterproc)then
       print *,"   k  p	          phi1	    phi2   psitrunc1   psitrunc2"
       print *,"--------------------------------------------------------"
       do k=1,nlev
          print 315,k,p(1,1,k),cm%D%phi(k,1),cm%D%phi(k,2),cm%D%psitrunc(k,1),cm%D%psitrunc(k,2)
       enddo
315    format(1x,i3,5f10.4)
    endif


    
    cm%INIT = .TRUE.
    !    endif

  end subroutine InitColumnModelMulticloud

  function Heaviside(lat,k,xpos) result(fx)
    real (kind=real_kind),intent(in) :: lat,k,xpos
    real (kind=real_kind)            :: fx
  
    fx = .5D0*(1.0D0-tanh(k*(abs(lat)-xpos)))
  
  end function Heaviside

  function Warm_pool(lat,lon,tbg) result(wp)

    real (kind=real_kind),intent(in) :: lat,lon,tbg
    real (kind=real_kind)            :: wp

    real (kind=real_kind)            :: tmin,tmax,delt

    tmin = -5.0D0
    tmax =  5.0D0
    delt = tmax-tmin

    wp = tmin
  !  if((lon <= DD_PI*0.5D0).and.(lon >= 0.0D0))then
  !     wp = wp + delt*sin(2.0D0*lon) 
  !  end if
!Enlarge the width of warm pool
      if((lon <= DD_PI).and.(lon >= 0.0D0))then
        wp = wp + delt*sin(lon)
      end if

    
    wp = Heaviside(lat,5.0D0,mc_mask_pos)*(wp+tbg)
    
  end function Warm_pool

  subroutine ApplyHeldSuarezForcing(elem, cm,dt, hvcoord, tl, nets, nete)
    type (element_t),intent(inout) :: elem(:)
    type (HeldSuarezForcing_t),intent(inout) :: cm
    real (kind=real_kind),intent(in)             :: dt
    type (hvcoord_t)     :: hvcoord
    type (TimeLevel_t)   :: tl
    integer,target       :: nets,nete
    integer                                      :: ie

    if(cm%INIT)then
       do ie=nets,nete          
          call hs_forcing(elem(ie),hvcoord,tl,dt)
       enddo
    endif

  end subroutine ApplyHeldSuarezForcing


  subroutine ApplyAquaplanetForcing(elem, elem_physics, hybrid,cm,dt,hvcoord, tl, nets, nete,mc)
    type (element_t),intent(inout) :: elem(:)
    type(elem_physics_t), intent(inout) :: elem_physics(:)
    type (hybrid_t), intent(in)       :: hybrid
    type (AquaplanetForcing_t),intent(inout) :: cm
    real (kind=real_kind),intent(in)             :: dt
    type (hvcoord_t)     :: hvcoord
    type (TimeLevel_t)   :: tl
    integer,target       :: nets,nete
    integer                                      :: ie
    type (ColumnModelMulticloud_t), intent(in),optional :: mc
    if(cm%INIT)then
       do ie=nets,nete
          call aquaplanet_forcing(dt,ie,elem(ie),elem_physics(ie), hybrid,hvcoord,nets,nete,tl,mc)
          ! AMIK was removed on 5/29/2008: put back on 12/24/2008 
          ! Removed again on 1/13/2009
          !call gravity_wave_drag_forcing(dt,elem(ie),tl)
       enddo
    endif




  end subroutine ApplyAquaplanetForcing

  subroutine ApplyColumnModelEmanuel(elem, hybrid,cm,dt,hvcoord,tl,nets, nete)
    type (element_t), intent(inout)      :: elem(:)
    type (hybrid_t), intent(in)       :: hybrid
    type (ColumnModelEmanuel_t),intent(inout),target :: cm
    real (kind=real_kind),intent(in)                 :: dt
    type (TimeLevel_t)   :: tl
    type (hvcoord_t)     :: hvcoord
    integer,target       :: nets,nete

    ! local
    type (ColumnDataEmanuel_t),pointer :: CD
    integer                            :: iflag
    integer                            :: i,j,ie,kk

    if(cm%INIT)then
       do ie=nets,nete
          
          call ConvertStateToColumnEmanuel(cm,elem(ie), tl, hvcoord)
          
          do j=1,np
             do i=1,np
#ifdef USE_MAXLAT
               if(abs(elem(ie)%sphereP(i,j)%lat).lt. emanuel_max_lat) then                
#endif
                CD => cm%col(i,j)
                
                call convect(CD%T,CD%Q,CD%QS,CD%U,CD%V,CD%TRA,CD%P,CD%PH,nlev,nlev-1,1,dt,&
                     iflag,CD%FT,CD%FQ,CD%FU,CD%FV,CD%FTRA,CD%PRECIP,CD%WD,CD%TPRIME,CD%QPRIME,&
                     pelem(ie)%state%CBMF(i,j))


                     pelem(ie)%surfc%precip(i,j)=cd%precip
                     pelem(ie)%surfc%wd(    i,j)=cd%wd
                     pelem(ie)%surfc%tprime(i,j)=cd%tprime
                     pelem(ie)%surfc%qprime(i,j)=cd%qprime
                     pelem(ie)%accum%precip(i,j)=pelem(ie)%accum%precip(i,j)+cd%precip

#ifdef USE_MAXLAT
		 endif
#endif
             enddo
          enddo

          call ConvertColumnToStateEmanuel(cm,elem(ie),tl)

#ifdef EMANUELSTATE
          if(MODULO(tl%nstep,statefreq) == 0) then
             call physics_state(hybrid,cm,ie,nets, nete)
          endif
#endif
!          if(MODULO(tl%nstep,moviefreq) == 0) then
!             call physics_movie(cm,ie)
!          endif
!          if(MODULO(tl%nstep,accumfreq) == 0) then
!             call physics_accum(cm,ie)
!          endif

       enddo

    endif

  end subroutine ApplyColumnModelEmanuel


  subroutine physics_state(hybrid,cm,ie, nets, nete)
    type (hybrid_t), intent(in)       :: hybrid
    integer, intent(in) :: ie
    integer,intent(in)       :: nets,nete
    type (ColumnModelEmanuel_t),intent(in), target :: cm
    type (ColumnDataEmanuel_t),pointer :: CD

#ifdef FORCINGSTAT
    real(kind=real_kind) :: ftmin_p, ftmax_p, ftsum_p
    real(kind=real_kind) :: fumin_p, fumax_p, fusum_p
    real(kind=real_kind) :: fvmin_p, fvmax_p, fvsum_p
    real(kind=real_kind) :: fqmin_p, fqmax_p, fqsum_p
#endif

    real(kind=real_kind) :: pmin_p, pmax_p, psum_p
    real(kind=real_kind) :: tpmin_p, tpmax_p, tpsum_p
    real(kind=real_kind) :: qpmin_p, qpmax_p, qpsum_p
    real(kind=real_kind) :: wdmin_p, wdmax_p, wdsum_p

    real(kind=real_kind) :: cbmfmin_p, cbmfmax_p, cbmfsum_p
    real(kind=real_kind) :: pamin_p, pamax_p, pasum_p
    integer :: i,j,k,iee


    if(ie==nets) then ! first column
       allocate(precip_min(nets:nete), precip_max(nets:nete), precip_sum(nets:nete) )
       allocate(tp_min(nets:nete), tp_max(nets:nete), tp_sum(nets:nete) )
       allocate(qp_min(nets:nete), qp_max(nets:nete), qp_sum(nets:nete) )
       allocate(wd_min(nets:nete), wd_max(nets:nete), wd_sum(nets:nete) )
       allocate(cbmf_min(nets:nete), cbmf_max(nets:nete), cbmf_sum(nets:nete) )
       allocate(pa_min(nets:nete), pa_max(nets:nete), pa_sum(nets:nete) )
#ifdef FORCINGSTAT
       allocate(ft_min(nets:nete), ft_max(nets:nete), ft_sum(nets:nete) )
       allocate(fu_min(nets:nete), fu_max(nets:nete), fu_sum(nets:nete) )
       allocate(fv_min(nets:nete), fv_max(nets:nete), fv_sum(nets:nete) )
       allocate(fq_min(nets:nete), fq_max(nets:nete), fq_sum(nets:nete) )
#endif
       precip_min = 1.0e30
       precip_max = -1.0e30
       precip_sum = 0.0D0

       tp_min = 1.0e30
       tp_max = -1.0e30
       tp_sum = 0.0D0

       qp_min = 1.0e30
       qp_max = -1.0e30
       qp_sum = 0.0D0

       wd_min = 1.0e30
       wd_max = -1.0e30
       wd_sum = 0.0D0

       cbmf_min = 1.0e30
       cbmf_max = -1.0e30
       cbmf_sum = 0.0D0

       pa_min = 1.0e30
       pa_max = -1.0e30
       pa_sum = 0.0D0
#ifdef FORCINGSTAT
       ft_min = 1.0e30
       ft_max = -1.0e30
       ft_sum = 0.0D0

       fu_min = 1.0e30
       fu_max = -1.0e30
       fu_sum = 0.0D0

       fv_min = 1.0e30
       fv_max = -1.0e30
       fv_sum = 0.0D0

       fq_min = 1.0e30
       fq_max = -1.0e30
       fq_sum = 0.0D0
#endif
    endif

    cbmf_min(ie) = minval(pelem(ie)%state%cbmf)
    cbmf_max(ie) = maxval(pelem(ie)%state%cbmf)
    cbmf_sum(ie) = sum(pelem(ie)%state%cbmf)

    pa_min(ie) = minval(pelem(ie)%accum%precip)
    pa_max(ie) = maxval(pelem(ie)%accum%precip)
    pa_sum(ie) = sum(pelem(ie)%accum%precip)

    
    do j=1,np
       do i=1,np
          cd => cm%col(i,j)
          if(cd%precip < precip_min(ie)) then
             precip_min(ie)=cd%precip 
          end if
          if(cd%precip > precip_max(ie)) then
             precip_max(ie)=cd%precip 
          end if

          if(cd%tprime < tp_min(ie)) then
             tp_min(ie)=cd%tprime 
          end if
          if(cd%tprime > tp_max(ie)) then
             tp_max(ie)=cd%tprime 
          end if

          if(cd%qprime < qp_min(ie)) then
             qp_min(ie)=cd%qprime 
          end if
          if(cd%qprime > qp_max(ie)) then
             qp_max(ie)=cd%qprime 
          end if

          if(cd%wd < wd_min(ie)) then
             wd_min(ie)=cd%wd 
          end if
          if(cd%wd > wd_max(ie)) then
             wd_max(ie)=cd%wd 
          end if
#ifdef FORCINGSTAT
          do k=1,nlev
          if(cd%ft(k) < ft_min(ie)) then
             ft_min(ie)=cd%ft(k)
          end if
          if(cd%ft(k) > ft_max(ie)) then
             ft_max(ie)=cd%ft(k) 
          end if

          if(cd%fu(k) < fu_min(ie)) then
             fu_min(ie)=cd%fu(k) 
          end if
          if(cd%fu(k) > fu_max(ie)) then
             fu_max(ie)=cd%fu(k) 
          end if

          if(cd%fv(k) < fv_min(ie)) then
             fv_min(ie)=cd%fv(k) 
          end if
          if(cd%fv(k) > fv_max(ie)) then
             fv_max(ie)=cd%fv(k) 
          end if

          if(cd%fq(k) < fq_min(ie)) then
             fq_min(ie)=cd%fq(k) 
          end if
          if(cd%fq(k) > fq_max(ie)) then
             fq_max(ie)=cd%fq(k) 
          end if
          enddo
#endif
          precip_sum(ie)=precip_sum(ie)+cd%precip
          tp_sum(ie)=tp_sum(ie)+cd%tprime
          qp_sum(ie)=qp_sum(ie)+cd%qprime
          wd_sum(ie)=wd_sum(ie)+cd%wd
#ifdef FORCINGSTAT
          do k=1,nlev
          ft_sum(ie)=ft_sum(ie)+cd%ft(k)
          fu_sum(ie)=fu_sum(ie)+cd%fu(k)
          fv_sum(ie)=fv_sum(ie)+cd%fv(k)
          fq_sum(ie)=fq_sum(ie)+cd%fq(k)
          enddo
#endif
       end do
    end do

    if(ie==nete) then ! last column
       pmin_p = ParallelMin(precip_min,hybrid)
       pmax_p = ParallelMax(precip_max,hybrid)

       tpmin_p = ParallelMin(tp_min,hybrid)
       tpmax_p = ParallelMax(tp_max,hybrid)

       qpmin_p = ParallelMin(qp_min,hybrid)
       qpmax_p = ParallelMax(qp_max,hybrid)

       wdmin_p = ParallelMin(wd_min,hybrid)
       wdmax_p = ParallelMax(wd_max,hybrid)

       cbmfmin_p = ParallelMin(cbmf_min,hybrid)
       cbmfmax_p = ParallelMax(cbmf_max,hybrid)

       pamin_p = ParallelMin(pa_min,hybrid)
       pamax_p = ParallelMax(pa_max,hybrid)

       do iee = nets, nete
         global_shared_buf(iee,1) = precip_sum(iee)
         global_shared_buf(iee,2) = tp_sum(iee)
         global_shared_buf(iee,3) = qp_sum(iee)
         global_shared_buf(iee,4) = wd_sum(iee)
         global_shared_buf(iee,5) = cbmf_sum(iee)
         global_shared_buf(iee,6) = pa_sum(iee)
       enddo
       call wrap_repro_sum(nvars=6, comm=hybrid%par%comm)
       psum_p = global_shared_sum(1)
       tpsum_p = global_shared_sum(2)
       qpsum_p = global_shared_sum(3)
       wdsum_p = global_shared_sum(4)
       cbmfsum_p = global_shared_sum(5)
       pasum_p = global_shared_sum(6)

#ifdef FORCINGSTAT
       ftmin_p = ParallelMin(ft_min,hybrid)
       ftmax_p = ParallelMax(ft_max,hybrid)

       fumin_p = ParallelMin(fu_min,hybrid)
       fumax_p = ParallelMax(fu_max,hybrid)

       fvmin_p = ParallelMin(fv_min,hybrid)
       fvmax_p = ParallelMax(fv_max,hybrid)

       fqmin_p = ParallelMin(fq_min,hybrid)
       fqmax_p = ParallelMax(fq_max,hybrid)

       do iee = nets, nete
         global_shared_buf(iee,1) = ft_sum(iee)
         global_shared_buf(iee,2) = fu_sum(iee)
         global_shared_buf(iee,3) = fv_sum(iee)
         global_shared_buf(iee,4) = fq_sum(iee)
       enddo
       call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
       ftsum_p = global_shared_sum(1)
       fusum_p = global_shared_sum(2)
       fvsum_p = global_shared_sum(3)
       fqsum_p = global_shared_sum(4)
#endif

       if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
#ifdef FORCINGSTAT
          write (*,100) "Uef    = "     ,fumin_p,fumax_p,fusum_p
          write (*,100) "Vef    = "     ,fvmin_p,fvmax_p,fvsum_p
          write (*,100) "Tef    = "     ,ftmin_p,ftmax_p,ftsum_p
          write (*,100) "Qef    = "     ,fqmin_p,fqmax_p,fqsum_p
          write (*,100) "TefW   = "     ,ftmin_p*Cp,ftmax_p*Cp,ftsum_p*Cp
          write (*,100) "QefW   = "     ,fqmin_p*2.53e6,fqmax_p*2.53e6,fqsum_p*2.53e6
          write (*,*) 
#endif
          write (*,100) "precip = "     ,pmin_p,pmax_p,psum_p
          write (*,100) "accumulated = ",pamin_p,pamax_p,pasum_p
          write (*,100) "tprime = "     ,tpmin_p,tpmax_p,tpsum_p
          write (*,100) "qprime = "     ,qpmin_p,qpmax_p,qpsum_p
          write (*,100) "wd     = "     ,wdmin_p,wdmax_p,wdsum_p
          write (*,100) "cbmf     = "   ,cbmfmin_p,cbmfmax_p,cbmfsum_p
          write (*,*)
       endif
       deallocate(precip_min,precip_max, precip_sum)
       deallocate(tp_min,tp_max, tp_sum)
       deallocate(qp_min,qp_max, qp_sum)
       deallocate(wd_min,wd_max, wd_sum)
       deallocate(cbmf_min,cbmf_max, cbmf_sum)
       deallocate(pa_min,pa_max, pa_sum)
#ifdef FORCINGSTAT
       deallocate(ft_min,ft_max, ft_sum)
       deallocate(fu_min,fu_max, fu_sum)
       deallocate(fv_min,fv_max, fv_sum)
       deallocate(fq_min,fq_max, fq_sum)
#endif
    endif

100 format (A9,3(E24.15))
  end subroutine physics_state

#endif

end module column_model_mod
