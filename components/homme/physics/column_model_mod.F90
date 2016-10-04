#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
#define EMANUELSTATE
!#define FORCINGSTAT
!#define USE_MAXLAT
#define _DBG_ !DBG

module column_model_mod
#ifdef _PRIM
  use element_mod,     only : element_t
  use hybvcoord_mod,   only : hvcoord_t 
  use hybrid_mod,      only : hybrid_t
  use kinds,           only : real_kind, int_kind
  use time_mod,        only : TimeLevel_t
  use physics_mod,     only : elem_physics_t, Specific_Humidity, Saturation_Specific_Humidity, getsurfpress, Temp2PotTemp
  use dimensions_mod,  only : nlev, nlevp, np, qsize, nelemd
  use control_mod,     only : integration, columnpackage, test_case,  physics, &
                              statefreq, &
                              TRACERADV_TOTAL_DIVERGENCE, &
                              tracer_advection_formulation
  use held_suarez_mod, only : hs_forcing
  use forcing_mod,     only : Apply_Forcing,EXP_EULER,INTERP_BDF2
  use reduction_mod,   only : parallelmax,parallelmin

  use parallel_mod, only : abortmp, global_shared_buf, global_shared_sum
  use global_norms_mod, only: wrap_repro_sum
  use physics_types_mod, only : pelem
  use physical_constants, only : Cp, DD_PI, p0, Cp
  use global_norms_mod,      only : global_integral
  use column_types_mod,      only : HeldSuarezForcing_t, ColumnModel_t
  
  implicit none

  private
#ifdef USE_MAXLAT
  real(kind=real_kind), parameter :: emanuel_max_lat=1.40
#endif

  public :: InitColumnModel
  public :: ApplyColumnModel

contains

  subroutine InitColumnModel(elem, cm,hvcoord,hybrid,tl,nets,nete,runtype)

    use Manager

    type(element_t), intent(inout) :: elem(:)
    type(elem_physics_t), pointer :: elem_physics(:)
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
                   !elem(ie)%state%Q(i,j,k,1,1:3)      = 0_real_kind
                   !elem(ie)%state%Qdp(i,j,k,1,1:3)      = 0_real_kind
                   elem(ie)%state%Q(i,j,k,1)      = 0_real_kind
                   elem(ie)%state%Qdp(i,j,k,1,1:2)      = 0_real_kind

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

    ! The MC needs the sounding
    if(columnpackage == "none" .AND. test_case(1:12) == "held_suarez0")then
       call InitHeldSuarezForcing(cm%cm_hs)
    endif


  end subroutine InitColumnModel

  subroutine ApplyColumnModel(elem,  hybrid, hvcoord, cm,dt)
    use hybvcoord_mod, only : hvcoord_t
    use Manager

    type (element_t), intent(inout) :: elem(:)
    type(elem_physics_t), pointer :: elem_physics(:)
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
    if(test_case(1:12) == "held_suarez0")then
      forcing=1
      call ApplyHeldSuarezForcing(elem, cm%cm_hs,dt, cm%hvcoord,cm%tl,cm%nets,cm%nete)
    endif
    if (forcing==1) then
       do ie=cm%nets,cm%nete
          call Apply_Forcing(EXP_EULER,elem(ie), hvcoord, cm%tl,dt)
       end do
    endif

  end subroutine ApplyColumnModel


  ! The rest of the file is private

  subroutine InitHeldSuarezForcing(cm)

    type (HeldSuarezForcing_t) :: cm

    if(.NOT.cm%INIT)then
       cm%INIT = .TRUE.
    endif

  end subroutine InitHeldSuarezForcing

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
#endif

end module column_model_mod
