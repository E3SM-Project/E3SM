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
  use physics_mod,     only : Specific_Humidity, Saturation_Specific_Humidity
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

  use parallel_mod, only : abortmp, global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use physics_types_mod, only : pelem
  use physical_constants, only : Cp, DD_PI
  use gravity_wave_drag_mod, only : gravity_wave_drag_forcing

  implicit none

  private
#ifdef USE_MAXLAT
  real(kind=real_kind), parameter :: emanuel_max_lat=1.40
#endif

  type, public :: ColumnDataEmanuel_t
     sequence
     ! in
     real (kind=real_kind)             :: T(nlev)
     real (kind=real_kind)             :: Q(nlev)
     real (kind=real_kind)             :: QS(nlev)
     real (kind=real_kind)             :: U(nlev)
     real (kind=real_kind)             :: V(nlev)
     real (kind=real_kind)             :: TRA(nlev,1)
     real (kind=real_kind)             :: P(nlev)
     real (kind=real_kind)             :: PH(nlevp)
     real (kind=real_kind)             :: R(nlev)

     ! out
     real (kind=real_kind)             :: FT(nlev)
     real (kind=real_kind)             :: FQ(nlev)
     real (kind=real_kind)             :: FU(nlev)
     real (kind=real_kind)             :: FV(nlev)
     real (kind=real_kind)             :: FTRA(nlev,1)
     real (kind=real_kind)             :: PRECIP   
     real (kind=real_kind)             :: WD
     real (kind=real_kind)             :: TPRIME
     real (kind=real_kind)             :: QPRIME
  end type ColumnDataEmanuel_t
  
  type, public :: ColumnModelEmanuel_t
     type (ColumnDataEmanuel_t) :: col(np,np)
#if 0
     type (hvcoord_t),pointer   :: hvcoord
     type (TimeLevel_t),pointer :: tl
     integer,pointer            :: nets,nete
#endif
     logical                    :: INIT = .FALSE.
  end type ColumnModelEmanuel_t
  
  type, public :: HeldSuarezForcing_t
#if 0
     type (hvcoord_t),pointer   :: hvcoord
     type (TimeLevel_t),pointer :: tl
     integer,pointer            :: nets,nete
#endif
     logical                    :: INIT = .FALSE.
  end type HeldSuarezForcing_t

  type, public :: AquaplanetForcing_t
#if 0
     type (hvcoord_t),pointer   :: hvcoord
     type (TimeLevel_t),pointer :: tl
     integer,pointer            :: nets,nete
#endif
     logical                    :: INIT = .FALSE.
  end type AquaplanetForcing_t


#if 0
  ! Polymorphism example (suggestion)
  
  ! CRCP
  type, public :: ColumnDataCRCP_t
  end type ColumnDataCRCP_t
  type, public :: ColumnModelCRCP_t
  end type ColumnModelCRCP_t
  
  ! CAM
  type, public :: ColumnDataCAM_t
  end type ColumnDataCAM_t
  type, public :: ColumnModelCAM_t
  end type ColumnModelCAM_t
#endif
  
  
  type, public :: ColumnModel_t
     type (ColumnModelEmanuel_t)    :: cm_em
     type (AquaplanetForcing_t)     :: cm_aq
     type (HeldSuarezForcing_t)     :: cm_hs
     type (hvcoord_t)            :: hvcoord
     type (TimeLevel_t),pointer  :: tl
     integer                     :: nets,nete
  end type ColumnModel_t

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

  subroutine InitColumnModel(elem, cm,hvcoord,tl,nets,nete,runtype)
    type(element_t), intent(inout) :: elem(:)
    type (ColumnModel_t) :: cm
    type (hvcoord_t), intent(in), target     :: hvcoord
    type (TimeLevel_t), intent(in), target   :: tl
    integer, intent(in), target              :: nets,nete,runtype

    integer :: i,j,k,ie
    if(runtype.ne.1) then
    ! for all column type
       do ie=nets,nete
          do k=1,nlev
             do j=1,np	
                do i=1,np	
                   elem(ie)%state%Q(i,j,k,1,1:3)      = 0_real_kind
                   elem(ie)%derived%FQ(i,j,k,1,1:3)     = 0_real_kind
                   elem(ie)%derived%FM(i,j,1:2,k,1:3) = 0_real_kind
                   elem(ie)%derived%FT(i,j,k,1:3)     = 0_real_kind
                enddo
             enddo
          enddo
       enddo
    endif	
    cm%hvcoord = hvcoord
    cm%tl=> tl
    cm%nets = nets
    cm%nete = nete

    ! Special features are going here

    if(columnpackage == "emanuel")then
       call InitColumnModelEmanuel(cm%cm_em,nets,nete)
    endif
    if(columnpackage == "none" .AND. test_case(1:12) == "held_suarez0")then
       call InitHeldSuarezForcing(cm%cm_hs)
    elseif(test_case(1:10) == "aquaplanet")then
       call InitAquaplanetForcing(cm%cm_aq)
    endif

  end subroutine InitColumnModel



  subroutine ApplyColumnModel(elem, hybrid,cm,dt)
    type (element_t), intent(inout) :: elem(:)
    type (ColumnModel_t),intent(inout) :: cm
    real (kind=real_kind),intent(in)   :: dt
    type (hybrid_t), intent(in)       :: hybrid

    integer :: i,j,k,ie, nm1, forcing=0

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
      call ApplyAquaplanetForcing(elem, hybrid,cm%cm_aq,dt, cm%hvcoord,cm%tl,cm%nets,cm%nete)
    endif

    if (forcing==1) then
       do ie=cm%nets,cm%nete
          call Apply_Forcing(EXP_EULER,elem(ie),cm%hvcoord,cm%tl,dt)
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
	     if(abs(elemin%spherev(i,j)%lat).lt. emanuel_max_lat) then
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
    	     if(abs(elemin%spherev(i,j)%lat).lt. emanuel_max_lat) then
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
    	     if(abs(elemin%spherev(i,j)%lat).lt. emanuel_max_lat) then
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
    	     if(abs(elemin%spherev(i,j)%lat).lt. emanuel_max_lat) then
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
    	     if(abs(elemin%spherev(i,j)%lat).lt. emanuel_max_lat) then
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


  subroutine ApplyAquaplanetForcing(elem, hybrid,cm,dt,hvcoord, tl, nets, nete)
    type (element_t),intent(inout) :: elem(:)
    type (hybrid_t), intent(in)       :: hybrid
    type (AquaplanetForcing_t),intent(inout) :: cm
    real (kind=real_kind),intent(in)             :: dt
    type (hvcoord_t)     :: hvcoord
    type (TimeLevel_t)   :: tl
    integer,target       :: nets,nete
    integer                                      :: ie

    if(cm%INIT)then
       do ie=nets,nete
          call aquaplanet_forcing(dt,ie,elem(ie),hybrid,hvcoord,nets,nete,tl)
          call gravity_wave_drag_forcing(dt,elem(ie),tl)
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
               if(abs(elem(ie)%spherev(i,j)%lat).lt. emanuel_max_lat) then                
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
