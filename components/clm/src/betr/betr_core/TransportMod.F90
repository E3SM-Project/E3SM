module TransportMod
  !
  ! !DESCRIPTION:
  !
  ! subroutines to do 1d vertical multiphase transport in soil/water
  ! History: created by Jinyun Tang, Jun 2011

#include "shr_assert.h"
  ! !USES:
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use tracer_varcon       , only : bndcond_as_conc, bndcond_as_flux
  use clm_varctl          , only : iulog
  use abortutils          , only : endrun
  use shr_kind_mod        , only : r8 => shr_kind_r8
  implicit none
  private
  public :: DiffusTransp              !do tracer transport through diffusion, for both lake and soil
  public :: calc_interface_conductance
  public :: init_transportmod
  public :: get_cntheta
  public :: calc_col_CFL              !claculate CFL critieria
  public :: semi_lagrange_adv_backward
  interface DiffusTransp
     module procedure DiffusTransp_gw
     module procedure DiffusTransp_solid
  end interface DiffusTransp

  type, private :: Extra_type
     real(r8), pointer :: zi(:)               !interfaces
     real(r8), pointer :: us(:)               !flow velocity at the interfaces
     integer           :: nlen                !total number of interfaces
   contains
     procedure, public :: InitAllocate
     procedure, public :: DDeallocate
     procedure, public :: AAssign
  end type Extra_type

  type(Extra_type), private :: Extra_inst

  !default configuration parameters
  real(r8), private :: cntheta

contains
  !-------------------------------------------------------------------------------
  subroutine InitAllocate(this, lbj, ubj)
    !
    ! !DESCRIPTION:
    ! allocate memory for arrays of the specified data type

    ! !ARGUMENTS:
    class(Extra_type)   :: this
    integer, intent(in) :: lbj, ubj
    character(len=32) :: subname ='InitAllocate'


    allocate(this%zi(lbj:ubj))
    allocate(this%us(lbj:ubj))

  end subroutine InitAllocate
  !-------------------------------------------------------------------------------
  
  subroutine DDeallocate(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memories
    !
    ! !ARGUMENTS:
    class(Extra_type) :: this
    character(len=32) :: subname ='DDeallocate'
    deallocate(this%zi)
    deallocate(this%us)

  end subroutine DDeallocate

  !-------------------------------------------------------------------------------
  
  subroutine AAssign(this, zi_t,us_t)
    !
    ! !DESCRIPTION:
    ! Assgin values for member variables for the specified data type
    !
    ! !ARGUMENTS:
    class(Extra_type) :: this
    real(r8), dimension(:), intent(in) :: zi_t
    real(r8), dimension(:), intent(in) :: us_t

    ! !LOCAL VARIABLES:
    integer :: n1, n2
    character(len=32) :: subname ='AAssign'

    n1 = size(zi_t)
    n2 = size(us_t)
    SHR_ASSERT_ALL((n1              == n2),        errMsg(__FILE__,__LINE__))
    this%zi(1:n1) = zi_t
    this%us(1:n2) = us_t
    this%nlen = n1
  end subroutine AAssign
  !-------------------------------------------------------------------------------
  function get_cntheta()result(ans)
    !
    ! !DESCRIPTION:
    ! return the theta factor
    !
    implicit none

    ! !LOCAL VARIABLES:
    real(r8) :: ans
    character(len=32) :: subname ='get_cntheta'

    ans = cntheta
    return
  end function get_cntheta
  !-------------------------------------------------------------------------------
  subroutine init_transportmod(lcntheta)
    !
    ! !DESCRIPTION:
    ! initialize transportmod
    !
    implicit none
    ! !ARGUMENTS:
    real(r8), optional, intent(in) :: lcntheta
    character(len=32) :: subname ='init_transportmod'
    if(present(lcntheta))then
       cntheta = lcntheta
    else
       ! use implicit solver by default
       cntheta = 1._r8
    endif

  end subroutine init_transportmod
  !-------------------------------------------------------------------------------
  
  subroutine calc_interface_conductance(bounds, lbj, ubj, jtop, numfl, filter, bulkdiffus, dz, hmconductance)
    !
    ! !DESCRIPTION:
    ! calcualte conductances at the interfaces using input layered diffusivity and
    ! thickness
    !
    ! !USES:
    !
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod,  only : bounds_type
    implicit none
    ! !ARGUMENTS:
    type(bounds_type),  intent(in) :: bounds                              !bounds
    integer,     intent(in)        :: lbj, ubj                            ! lbinning and ubing level indices
    integer,     intent(in)        :: jtop(bounds%begc: )                 ! index of upper boundary, which could be variable
    integer,     intent(in)        :: numfl                               ! length of the filter
    integer,     intent(in)        :: filter(:)                           ! the actual filter
    real(r8),    intent(in)        :: bulkdiffus(bounds%begc: ,lbj: )     !weighted bulk diffusivity for dual-phase diffusion
    real(r8),    intent(in)        :: dz(bounds%begc: , lbj: )
    real(r8), intent(inout)        :: hmconductance(bounds%begc: , lbj: ) !weighted bulk conductance

    ! !LOCAL VARIABLES:
    integer :: n, c, fc

    SHR_ASSERT_ALL((ubound(jtop)              == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(dz)                == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bulkdiffus)        == (/bounds%endc, ubj/)),     errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(hmconductance)     == (/bounds%endc, ubj-1/)),   errMsg(__FILE__,__LINE__))

    do n=lbj, ubj-1
       do fc = 1, numfl
          c = filter(fc)
          if(n>=jtop(c))then
             hmconductance(c,n) = 2._r8/(dz(c,n)/bulkdiffus(c,n)+dz(c,n+1)/bulkdiffus(c,n+1))
          endif
       enddo
    enddo


  end subroutine calc_interface_conductance
  !-------------------------------------------------------------------------------
  subroutine DiffusTransp_gw_tridiag(bounds, lbj, ubj, jtop, numfl, filter, ntrcs, trcin_mobile, &
       Rfactor, hmconductance, dtime, dz, source, trc_concflx_air,condc_toplay, topbc_type,&
       bot_concflx, update_col, source_only, rt, at,bt,ct, botbc_type, condc_botlay)
    !
    ! !DESCRIPTION:
    ! Assemble the tridiagonal matrix for the multiphase diffusive transport
    !
    ! !USES
    !
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod,  only : bounds_type

    implicit none
    ! !ARGUMENTS:
    type(bounds_type),  intent(in) :: bounds                           ! bounds
    integer,  intent(in) :: lbj, ubj                                   ! lbinning and ubing level indices
    integer,  intent(in) :: jtop(bounds%begc: )                        ! index of upper boundary, which could be variable
    integer,  intent(in) :: numfl                                      ! length of the filter
    integer,  intent(in) :: filter(:)                                  ! the actual filter
    integer,  intent(in) :: ntrcs
    real(r8), intent(in) :: Rfactor(bounds%begc: , lbj: )              ! conversion parameter from the given tracer phase to bulk mobile phase
    real(r8), intent(in) :: hmconductance(bounds%begc: , lbj: )        ! weighted bulk tracer conductances
    real(r8), intent(in) :: dz(bounds%begc: , lbj: )                   ! node thickness
    real(r8), intent(in) :: dtime(bounds%begc: )                       ! time step
    real(r8), intent(in) :: condc_toplay(bounds%begc: )                ! top layer conductance
    integer,  intent(in) :: topbc_type                                 ! type of top boundary condtion: 1, concentration, 2 flux
    real(r8), intent(in) :: bot_concflx    (bounds%begc: , 1: , 1: )   ! flux or concentration at the bottom boundary
    real(r8), intent(in) :: trc_concflx_air(bounds%begc: , 1: , 1: )   ! atmospheric tracer concentration (topbc_type=1) or flux (topbc_type=2)
    real(r8), intent(in) :: trcin_mobile   (bounds%begc: , lbj: , 1: ) ! incoming mobile tracer concentration
    real(r8), intent(in) :: source         (bounds%begc: , lbj: , 1: ) ! chemical sources [mol/m3]

    logical,  intent(in) :: source_only                                ! if .true. only update the source array rt, used for explicit solver
    logical,  intent(in) :: update_col(bounds%begc: )                  ! logical switch indicating if the column is for active update

    real(r8), intent(out):: rt(bounds%begc: ,lbj: , 1: )               ! tridiagonal matrix element r
    real(r8), optional,intent(inout):: at(bounds%begc: , lbj: )        ! tridiagonal matrix element a
    real(r8), optional,intent(inout):: bt(bounds%begc: , lbj: )        ! tridiagonal matrix element b
    real(r8), optional,intent(inout):: ct(bounds%begc: , lbj: )        ! tridiagonal matrix element c
    integer,  optional,intent(in)   :: botbc_type                      ! type of bottom boundary condition
    real(r8), optional,intent(in)   :: condc_botlay(bounds%begc: )     !conductance at bottom layer

    ! !LOCAL VARIABLES:
    integer :: j, fc, c, k      !indices
    integer :: botbc_ltype      !temp. variable
    real(r8) ::Fl, Fr
    character(len=255) :: subname='DiffusTransp_gw'

    SHR_ASSERT_ALL((ubound(jtop)            == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(Rfactor)         == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(hmconductance)   == (/bounds%endc, ubj-1/)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(dz)              == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(dtime)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(condc_toplay)    == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(update_col)      == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))

    SHR_ASSERT_ALL((ubound(source)           == (/bounds%endc, ubj, ntrcs/)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(rt)              == (/bounds%endc, ubj, ntrcs/)), errMsg(__FILE__,__LINE__))

    SHR_ASSERT_ALL((ubound(trcin_mobile , 1)  == (/bounds%endc/))  , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(trcin_mobile , 2)  == (/ubj/))          , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(trcin_mobile   , 3)  == (/3/))            , errMsg(__FILE__,__LINE__))

    SHR_ASSERT_ALL((ubound(bot_concflx,1) == (/bounds%endc/))  , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(bot_concflx,2) == (/2/))            , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(bot_concflx,  3) == (/ntrcs/))        , errMsg(__FILE__,__LINE__))

    SHR_ASSERT_ALL((ubound(trc_concflx_air, 1) == (/bounds%endc/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(trc_concflx_air, 2) == (/2/))           , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(trc_concflx_air,   3) == (/3/))           , errMsg(__FILE__,__LINE__))


    if(.not. source_only) then
       SHR_ASSERT_ALL((ubound(at)            == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
       SHR_ASSERT_ALL((ubound(bt)            == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
       SHR_ASSERT_ALL((ubound(ct)            == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
    endif
    !unless specified explicitly, the bottom boundary condition is given as flux
    if(present(botbc_type))then
       botbc_ltype = botbc_type
       if(botbc_type==bndcond_as_conc)then
          SHR_ASSERT_ALL((ubound(condc_botlay)    == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
       endif
    else
       botbc_ltype = bndcond_as_flux
    endif

    do fc = 1, numfl
       !form the diffusion matrix
       c = filter(fc)
       if(update_col(c))then

          do j = jtop(c), ubj
             do k = 1, ntrcs
                if(j == jtop(c))then
                   !by default the top node is always assumed as snow surface,
                   !though when it snow free, the conductance is defined with respect to the soil
                   Fr = -hmconductance(c,j)*(trcin_mobile(c,j+1, k)/rfactor(c,j+1)-trcin_mobile(c,j, k)/rfactor(c,j))
                   if(topbc_type == bndcond_as_conc)then
                      !top boundary condition given as concentration
                      Fl = -condc_toplay(c)*(trcin_mobile(c,j, k)/rfactor(c,j)-trc_concflx_air(c, 1, k))
                   elseif(topbc_type == bndcond_as_flux)then
                      !top boundary condition given as flux, this only happens when the flux is given at soil surface
                      Fl = trc_concflx_air(c,1, k)
                   endif
                elseif(j == ubj)then
                   Fl = -hmconductance(c,j-1)*(trcin_mobile(c,j,k)/rfactor(c,j)-trcin_mobile(c,j-1,k)/rfactor(c,j-1))
                   if(botbc_ltype==bndcond_as_conc)then
                      Fr = - condc_botlay(c)*(bot_concflx(c,1,k)-trcin_mobile(c,j,k)/rfactor(c,j))
                   else
                      Fr = bot_concflx(c,1,k)
                   endif
                else
                   Fl = -hmconductance(c,j-1)*(trcin_mobile(c,j,k)/rfactor(c,j)-trcin_mobile(c,j-1,k)/rfactor(c,j-1))
                   Fr = -hmconductance(c,j)*(trcin_mobile(c,j+1,k)/rfactor(c,j+1)-trcin_mobile(c,j,k)/rfactor(c,j))
                endif
                rt(c,j,k) = Fl-Fr + source(c,j,k)*dz(c,j)
                if(j==jtop(c) .and. topbc_type == bndcond_as_conc)then
                   rt(c,j,k) = rt(c,j,k)+cntheta*condc_toplay(c)*(trc_concflx_air(c, 2,k)-trc_concflx_air(c, 1,k))
                endif
                if(j == ubj .and. botbc_ltype==bndcond_as_conc)then
                   rt(c,j,k) = rt(c,j,k) + cntheta*condc_botlay(c)*(bot_concflx(c,2,k) - bot_concflx(c,1,k))
                endif
             enddo
          enddo
       endif
    enddo

    if(source_only)return
    do fc = 1, numfl
       !form the diffusion matrix
       c = filter(fc)
       if(update_col(c))then
          do j = jtop(c), ubj
             if(j == jtop(c))then
                if(topbc_type == bndcond_as_conc)then !top boundary condition given as concentration
                   bt(c,j)=dz(c,j)/dtime(c)+cntheta*(hmconductance(c,j) &
                        +condc_toplay(c))/Rfactor(c,j)
                elseif(topbc_type == bndcond_as_flux)then   !top boundary condition given as flux
                   bt(c,j)=dz(c,j)/dtime(c)+cntheta*hmconductance(c,j)/Rfactor(c,j)
                endif
                ct(c,j)=-cntheta*hmconductance(c,j)/Rfactor(c,j+1)
             elseif(j==ubj)then
                at(c,j)=-cntheta*hmconductance(c,j-1)/rfactor(c,j-1)
                if(botbc_ltype == bndcond_as_conc)then
                   bt(c,j)=dz(c,j)/dtime(c)+cntheta*(hmconductance(c,j-1)+condc_botlay(c))/Rfactor(c,j)
                else
                   bt(c,j)=dz(c,j)/dtime(c)+cntheta*hmconductance(c,j-1)/Rfactor(c,j)
                endif
             else
                at(c,j)=-cntheta*hmconductance(c,j-1)/rfactor(c,j-1)
                ct(c,j)=-cntheta*hmconductance(c,j)/rfactor(c,j+1)
                bt(c,j)=dz(c,j)/dtime(c)+cntheta*(hmconductance(c,j-1)+hmconductance(c,j))/Rfactor(c,j)
             endif
          enddo
       endif
    enddo

  end subroutine DiffusTransp_gw_tridiag
!-------------------------------------------------------------------------------
   subroutine DiffusTransp_gw(bounds, lbj, ubj, jtop, numfl, filter, ntrcs, trcin_mobile, &
       Rfactor, hmconductance, dtime, dz, source, trc_concflx_air,condc_toplay, topbc_type,&
       bot_flux, update_col, dtracer, botbc_type, condc_botlay)
   !
   ! !DESCRIPTION:
   ! solve the dual phase transport problem.
   ! the solver returns the tracer change due to diffusive transport
   !
   ! !USES:
   use shr_kind_mod  , only : r8 => shr_kind_r8
   use decompMod     , only : bounds_type
   use TridiagonalMod, only : Tridiagonal

   implicit none
   ! !ARGUMENTS:
   type(bounds_type) , intent(in)    :: bounds                                   !bounds
   integer           , intent(in)    :: lbj, ubj                                 ! lbinning and ubing level indices
   integer           , intent(in)    :: jtop(bounds%begc: )                      ! index of upper boundary, which could be variable
   integer           , intent(in)    :: numfl                                    ! length of the filter
   integer           , intent(in)    :: filter(:)                                ! the actual filter
   integer           , intent(in)    :: ntrcs
   real(r8)          , intent(in)    :: Rfactor(bounds%begc: , lbj:  )           !conversion parameter from the given tracer phase to bulk mobile phase
   real(r8)          , intent(in)    :: hmconductance(bounds%begc: , lbj: )      !weighted bulk tracer conductances
   real(r8)          , intent(in)    :: dz(bounds%begc: , lbj: )                 !node thickness
   real(r8)          , intent(in)    :: dtime(bounds%begc: )                     !time step
   real(r8)          , intent(in)    :: condc_toplay(bounds%begc: )              !top layer conductance
   integer           , intent(in)    :: topbc_type                               !type of top boundary condtion: 1, concentration, 2 flux
   integer , optional, intent(in)    :: botbc_type
   real(r8), optional, intent(in)    :: condc_botlay(bounds%begc: )
   logical           , intent(in)    :: update_col(bounds%begc: )                !logical switch indicating if the column is for active update
   real(r8)          , intent(in)    :: trcin_mobile(bounds%begc: , lbj: ,1: )   ! incoming mobile tracer concentration
   real(r8)          , intent(in)    :: source(bounds%begc: , lbj: , 1: )        !chemical sources [mol/m3]
   real(r8)          , intent(in)    :: bot_flux(bounds%begc: , 1: , 1: )        !flux at the bottom boundary
   real(r8)          , intent(in)    :: trc_concflx_air(bounds%begc: , 1: , 1: ) !atmospheric tracer concentration (topbc_type=1) or flux (topbc_type=2)
   real(r8)          , intent(inout) :: dtracer(bounds%begc: , lbj: , 1: )       !change of tracer concentration during the time step

   ! !LOCAL VARIABLES:
   real(r8) :: rt(bounds%begc:bounds%endc, lbj:ubj, 1:ntrcs) !tridiagonal matrix element r
   real(r8) :: at(bounds%begc:bounds%endc, lbj:ubj)          !tridiagonal matrix element a
   real(r8) :: bt(bounds%begc:bounds%endc, lbj:ubj)          !tridiagonal matrix element b
   real(r8) :: ct(bounds%begc:bounds%endc, lbj:ubj)          !tridiagonal matrix element c
   real(r8) :: dtracer1(bounds%begc:bounds%endc, lbj:ubj)
   character(len=255) :: subname = 'DiffusTransp_gw'
   integer :: kk, fc, c

   SHR_ASSERT_ALL((ubound(jtop)              == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(dtime)             == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(update_col)        == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(condc_toplay)      == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(dz)                == (/bounds%endc, ubj/))  , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(Rfactor   )        == (/bounds%endc, ubj/))  , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(hmconductance)     == (/bounds%endc, ubj-1/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(source)    == (/bounds%endc, ubj,  ntrcs/)), errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(dtracer)       == (/bounds%endc, ubj, ntrcs/)) , errMsg(__FILE__,__LINE__))

   SHR_ASSERT_ALL((ubound(trcin_mobile , 1)  == (/bounds%endc/))  , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(trcin_mobile , 2)  == (/ubj/))          , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((size(trcin_mobile   , 3)  == (/ntrcs/))        , errMsg(__FILE__,__LINE__))

   SHR_ASSERT_ALL((ubound(bot_flux     , 1)   == (/bounds%endc/)) , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(bot_flux     , 2)   == (/2/))           , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((size(bot_flux       , 3)   == (/ntrcs/))       , errMsg(__FILE__,__LINE__))

   SHR_ASSERT_ALL((ubound(trc_concflx_air, 1) == (/bounds%endc/)) , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((ubound(trc_concflx_air, 2) == (/2/))           , errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((size(trc_concflx_air,   3) == (/ntrcs/))       , errMsg(__FILE__,__LINE__))

   !assemble the tridiagonal maxtrix
   if(present(botbc_type))then
     SHR_ASSERT_ALL((ubound(condc_botlay)    == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))

     call DiffusTransp_gw_tridiag(bounds, lbj, ubj, jtop, numfl, filter, ntrcs, trcin_mobile, &
        Rfactor, hmconductance, dtime, dz, source, trc_concflx_air,&
        condc_toplay, topbc_type, bot_flux, update_col, source_only=.false.,&
        rt=rt, at=at,bt=bt,ct=ct, botbc_type=botbc_type, condc_botlay=condc_botlay)
   else
     call DiffusTransp_gw_tridiag(bounds, lbj, ubj, jtop, numfl, filter, ntrcs, trcin_mobile, &
        Rfactor, hmconductance, dtime, dz, source, trc_concflx_air,&
        condc_toplay, topbc_type, bot_flux, update_col, source_only=.false.,&
        rt=rt, at=at,bt=bt,ct=ct)
   endif

   !calculate the change to tracer
   call Tridiagonal (bounds, lbj, ubj, jtop, numfl, filter, ntrcs, at, bt, ct, rt, dtracer, update_col)

   end subroutine DiffusTransp_gw

   !-------------------------------------------------------------------------------
   subroutine Diffustransp_solid_tridiag(bounds, lbj, ubj, lbn, numfl, filter, ntrcs, trcin,&
        hmconductance,  dtime_col, dz, source, update_col, at,bt,ct, rt)
     !
     ! !DESCRIPTION:
     !
     ! Do solid phase transport with tracer source
     !
     ! !USES:
     use shr_kind_mod, only: r8 => shr_kind_r8
     use decompMod,  only : bounds_type

     implicit none
     ! !ARGUMENTS:
     type(bounds_type),  intent(in) :: bounds                              !bounds
     integer  , intent(in)          :: lbj, ubj                            ! lbinning and ubing level indices
     integer  , intent(in)          :: lbn(bounds%begc: )                  !indices of top boundary
     integer  , intent(in)          :: numfl                               !filter dimension
     integer  , intent(in)          :: filter(:)                           !filter
     integer  , intent(in)          :: ntrcs
     real(r8) , intent(in)          :: trcin(bounds%begc: , lbj: ,1: )     !tracer concentration [mol/m3]
     real(r8) , intent(in)          :: hmconductance(bounds%begc: , lbj: ) !weighted conductance
     real(r8) , intent(in)          :: dtime_col(bounds%begc: )            !model time step
     real(r8) , intent(in)          :: dz(bounds%begc: , lbj: )            !layer thickness
     real(r8) , intent(in)          :: source(bounds%begc: , lbj: ,1: )    !chemical sources [mol/m3]
     logical  , intent(in)          :: update_col(bounds%begc: )           !logical switch indicating if the column is for active update
     real(r8) , intent(out)         :: at(bounds%begc: , lbj: )            !returning tridiagonal a matrix
     real(r8) , intent(out)         :: bt(bounds%begc: , lbj: )            !returning tridiagonal b matrix
     real(r8) , intent(out)         :: ct(bounds%begc: , lbj: )            !returning tridiagonal c matrix
     real(r8) , intent(out)         :: rt(bounds%begc: , lbj: ,1: )        !returning tridiagonal r matrix

     !LOCAL VARIABLES:
     real(r8) :: bot
     integer :: j, k, fc, c
     real(r8) :: Fl, Fr
     real(r8) :: dtime
     character(len=255) :: subname='DiffusTransp_solid_tridiag'

     SHR_ASSERT_ALL((ubound(lbn)           == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(hmconductance) == (/bounds%endc, ubj-1/)), errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dz)            == (/bounds%endc, ubj/))  , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dtime_col)     == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(update_col)    == (/bounds%endc/))       , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(at)            == (/bounds%endc, ubj/))  , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(bt)            == (/bounds%endc, ubj/))  , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(ct)            == (/bounds%endc, ubj/))  , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(rt)            == (/bounds%endc, ubj, ntrcs/))  , errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(source)        == (/bounds%endc, ubj, ntrcs/)), errMsg(__FILE__,__LINE__))

     SHR_ASSERT_ALL(((/ubound(trcin,1),ubound(trcin,2),size(trcin,3)/)           == (/bounds%endc, ubj, ntrcs/)), errMsg(__FILE__,__LINE__))



     !zero flux is imposed both at the top and bottom boundaries
     !set zero outgoing flux
     bot = 0._r8
     do fc = 1, numfl
        c = filter(fc)
        if(update_col(c))then
           dtime=dtime_col(c)
           do j = lbn(c), ubj
              do k = 1, ntrcs
                 if(j==lbn(c))then
                    Fr=-hmconductance(c,j)*(trcin(c,j+1,k)-trcin(c,j,k))
                    Fl=0._r8    !zero flux at top boundary for solid phase
                 elseif(j==ubj)then
                    !assume zero flux for diffusion
                    Fl=-hmconductance(c,j-1)*(trcin(c,j,k)-trcin(c,j-1,k))
                    Fr=bot
                 else
                    Fl=-hmconductance(c,j-1)*(trcin(c,j,k)-trcin(c,j-1,k))
                    Fr=-hmconductance(c,j)*(trcin(c,j+1,k)-trcin(c,j,k))
                 endif
                 rt(c,j,k) = Fl-Fr + source(c,j,k)*dz(c,j)
              enddo
           enddo

           do j = lbn(c), ubj
              if(j==lbn(c))then
                 !top boundary condition given as flux
                 at(c,j)=0._r8
                 bt(c,j)=dz(c,j)/dtime+cntheta*hmconductance(c,j)
                 ct(c,j)=-cntheta*hmconductance(c,j)
              elseif(j==ubj)then
                 at(c,j)=-cntheta*hmconductance(c,j-1)
                 bt(c,j)=dz(c,j)/dtime+cntheta*hmconductance(c,j-1)
              else
                 at(c,j)=-cntheta*hmconductance(c,j-1)
                 ct(c,j)=-cntheta*hmconductance(c,j)
                 bt(c,j)=dz(c,j)/dtime-at(c,j)-ct(c,j)
              endif
           enddo
        endif
     enddo

   end subroutine DiffusTransp_solid_tridiag
   !-------------------------------------------------------------------------------
   
   subroutine DiffusTransp_solid(bounds, lbj, ubj, lbn, numfl, filter, ntrcs, trcin,&
        hmconductance,  dtime_col, dz, source, update_col, dtracer)
     !
     ! !DESCRIPTION:
     ! Do diffusive solid phase tracer transport
     !
     ! !USES:
     use shr_kind_mod, only: r8 => shr_kind_r8
     use decompMod,  only : bounds_type
     use TridiagonalMod, only : Tridiagonal

     implicit none
     ! !ARGUMENTS:
     type(bounds_type),  intent(in) :: bounds                       ! bounds
     integer  , intent(in)   :: lbj, ubj                            ! lbinning and ubing level indices
     integer  , intent(in)   :: lbn(bounds%begc: )                  ! indices of top boundary
     integer  , intent(in)   :: numfl                               ! filter dimension
     integer  , intent(in)   :: filter(:)                           ! filter
     integer  , intent(in)   :: ntrcs
     real(r8) , intent(in)   :: hmconductance(bounds%begc: , lbj: ) ! weighted conductance
     real(r8) , intent(in)   :: dtime_col(bounds%begc: )            ! model time step
     real(r8) , intent(in)   :: dz(bounds%begc: , lbj: )            ! layer thickness
     real(r8) , intent(in)   :: trcin (bounds%begc: , lbj: , 1: )   ! tracer concentration [mol/m3]
     real(r8) , intent(in)   :: source(bounds%begc: , lbj: , 1: )   ! chemical sources [mol/m3/s]
     logical  , intent(in)   :: update_col(bounds%begc: )           ! logical switch indicating if the column is for active update
     real(r8), intent(inout) :: dtracer(bounds%begc: , lbj: ,1: )   ! update to the tracer

     ! !LOCAL VARIABLES:
     real(r8) :: at(bounds%begc:bounds%endc, lbj:ubj)             !returning tridiagonal a matrix
     real(r8) :: bt(bounds%begc:bounds%endc, lbj:ubj)             !returning tridiagonal b matrix
     real(r8) :: ct(bounds%begc:bounds%endc, lbj:ubj)             !returning tridiagonal c matrix
     real(r8) :: rt(bounds%begc:bounds%endc, lbj:ubj, 1:ntrcs)             !returning tridiagonal r matrix
     character(len=255) :: subname = 'DiffusTransp_solid'

     SHR_ASSERT_ALL((ubound(lbn)           == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(hmconductance) == (/bounds%endc, ubj-1/)), errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dz)            == (/bounds%endc, ubj/)),   errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(update_col)    == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dtime_col)     == (/bounds%endc/)),        errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dtracer)       == (/bounds%endc, ubj, ntrcs/)),   errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(source)        == (/bounds%endc, ubj, ntrcs/)),   errMsg(__FILE__,__LINE__))

     SHR_ASSERT_ALL(((/ubound(trcin,1),ubound(trcin,2),size(trcin,3)/)   == (/bounds%endc, ubj,ntrcs/)),   errMsg(__FILE__,__LINE__))


     !assemble the tridiagonal matrix
     call Diffustransp_solid_tridiag(bounds, lbj, ubj, lbn, numfl, filter, ntrcs, trcin,&
          hmconductance,  dtime_col, dz, source, update_col, at,bt,ct, rt)

     !calculate the change to tracer
     call Tridiagonal (bounds, lbj, ubj, lbn, numfl, filter, ntrcs, at, bt, ct, rt, dtracer, update_col)

   end subroutine DiffusTransp_solid
   !-------------------------------------------------------------------------------
   function calc_col_CFL(lbj, ubj, us, dx, dtime) result(cfl)
     !
     ! DESCRIPTION:
     ! calculate the CFL number for the given grid and velocity field
     ! this subroutine is now not actively used, but can be used
     ! when a Eulerian advection scheme is adopted.
     implicit none
     ! !ARGUMENTS:
     integer,  intent(in) :: lbj, ubj    !left and right bounds
     real(r8), intent(in) :: us(lbj: )   !velocity vector,    [m/s]
     real(r8), intent(in) :: dx(lbj: )   !node length,        [m]
     real(r8), intent(in) :: dtime       !imposed time step,  [s]

     ! !LOCAL VARIABLES:
     real(r8) :: cfl
     integer  :: len, j
     character(len=32) :: subname ='calc_col_CFL'

     SHR_ASSERT_ALL((ubound(us)         == (/ubj+1/)),        errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dx)         == (/ubj/)),   errMsg(__FILE__,__LINE__))

     cfl = 0._r8
     !the column cfl number is defined as the maximum over the whole domain
     do j = lbj, ubj
        if(us(j)>0._r8)then
           if(us(j)<0._r8)then
              cfl= max(dtime/dx(j)*max(abs(us(j)), abs(us(j+1))), cfl)
           else
              cfl=max(abs(dtime*us(j)/dx(j)), cfl)
           endif
        else
           if(us(j+1)>0._r8)then
              cfl=max(abs(dtime*us(j+1)/dx(j)),cfl)
           else
              cfl= max(dtime/dx(j)*max(abs(us(j)), abs(us(j+1))), cfl)
           endif
        endif
     enddo
   end function calc_col_CFL

   !-------------------------------------------------------------------------------
   subroutine semi_lagrange_adv_backward(bounds, lbj, ubj, lbn, numfl, filter, ntrcs, dtime, dz, &
        zi, us, inflx_top, inflx_bot, update_col, halfdt_col, trcin, trcou, leaching_mass)
     !
     ! DESCRIPTION:
     ! do semi-lagrangian advection for equation
     ! pu/pt+c*pu/px=0
     ! for a certain tracer group
     !
     ! !USES:
     use shr_kind_mod,     only : r8 => shr_kind_r8
     use decompMod,        only : bounds_type
     use MathfuncMod,      only : cumsum, cumdif, safe_div, dot_sum, asc_sort_vec
     use InterpolationMod, only : pchip_polycc, pchip_interp

     implicit none
     ! !ARGUMENTS:
     type(bounds_type) , intent(in)  :: bounds                             !bounds
     integer           , intent(in)  :: lbj, ubj                          ! lbinning and ubing level indices
     integer           , intent(in)  :: lbn(bounds%begc: )                !label of the top/left boundary
     integer           , intent(in)  :: numfl
     integer           , intent(in)  :: ntrcs
     integer           , intent(in)  :: filter(:)
     real(r8)          , intent(in)  :: dtime(bounds%begc: )
     real(r8)          , intent(in)  :: zi(bounds%begc: , lbj-1: )
     real(r8)          , intent(in)  :: dz(bounds%begc: , lbj: )
     real(r8)          , intent(in)  :: inflx_top(bounds%begc: , 1: )     ! incoming tracer flow at top boundary [mol/m2/s]
     real(r8)          , intent(in)  :: inflx_bot(bounds%begc: , 1: )     !incoming tracer flow at bottom boundary
     logical           , intent(in)  :: update_col(bounds%begc: )         !indicator of active clumns
     real(r8)          , intent(in)  :: us(bounds%begc: , lbj-1: )        !convective flux defined at the boundary, positive downwards, [m/s]
     logical           , intent(out) :: halfdt_col(bounds%begc:bounds%endc)
     real(r8)          , intent(in)  :: trcin(bounds%begc: , lbj: , 1: )  !input tracer concentration
     real(r8)          , intent(out) :: trcou(bounds%begc: , lbj: , 1: )
     real(r8), optional, intent(out) :: leaching_mass(bounds%begc: , 1: ) !leaching tracer mass

     ! !LOCAL VARIABLES:
     integer, parameter :: pn = 2                !first order lagrangian interpolation to avoid overshooting
     integer  :: j, fc, c, k
     integer  :: ntr                             !indices for tracer
     integer  :: length, lengthp2
     real(r8) :: mass_curve(0:ubj-lbj+5 , ntrcs) !total number of nodes + two ghost cells at each boundary
     real(r8) :: cmass_curve(0:ubj-lbj+5, ntrcs)
     real(r8) :: mass_new(1:ubj-lbj+1   , ntrcs)
     real(r8) :: cmass_new(0:ubj-lbj+1  , ntrcs)
     real(r8) :: zold(0:ubj-lbj+1)
     real(r8) :: di(0:ubj-lbj+5)
     real(r8) :: zghostl(1:2)                    !ghost grid left interface at the left boundary
     real(r8) :: zghostr(1:2)                    !ghost grid left interface at the right boundary
     real(r8) :: ughostl(1:2)                    !flow velocity at the ghost grid leff interface at the left boundary
     real(r8) :: ughostr(1:2)                    !flow velocity at the ghost grid leff interface at the right boundary
     real(r8) :: z0
     real(r8) :: zf
     real(r8) :: utmp
     real(r8) :: dinfl_mass
     character(len=32) :: subname='semi_lagrange_adv_backward'

     SHR_ASSERT_ALL((ubound(lbn)        == (/bounds%endc/)),         errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dtime)      == (/bounds%endc/)),         errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(dz)         == (/bounds%endc, ubj/)),    errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(update_col) == (/bounds%endc/)),         errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(us)         == (/bounds%endc, ubj/)),    errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(zi)         == (/bounds%endc, ubj/)),    errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(inflx_top)  == (/bounds%endc, ntrcs/)),  errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(inflx_bot)  == (/bounds%endc, ntrcs/)),  errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(leaching_mass)  == (/bounds%endc,ntrcs/)), errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((ubound(trcou)      == (/bounds%endc, ubj,ntrcs/)),    errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL(((/ubound(trcin,1),ubound(trcin,2),size(trcin,3)/)      == (/bounds%endc, ubj,ntrcs/)),    errMsg(__FILE__,__LINE__))

     call Extra_inst%InitAllocate(1,ubj-lbj+6)
     halfdt_col(:) = .false.
     do fc = 1, numfl

        c = filter(fc)
        if(.not. update_col(c))cycle
        !do backward advection for all boundaries, including leftmost (lbn(c)-1) and rightmost (ubj)
        length = ubj - lbn(c) + 1   !total number of grid cells
        lengthp2 = length + 4       ! add 2 ghost cells both at the left and right boundaries

        !define ghost boundary
        !NOTE: because of the setup, the left boundary and right boundary should have non-zero flow
        utmp = us(c,lbn(c)-1)
        zghostl(1) = -abs(utmp)*dtime(c)*2._r8 + zi(c,lbn(c)-1)-2.e-20_r8
        zghostl(2) = -abs(utmp)*dtime(c)       + zi(c,lbn(c)-1)-1.e-20_r8
        ughostl(1) = us(c,lbn(c)-1)
        ughostl(2) = us(c,lbn(c)-1)

        zghostr(1) = zi(c,ubj) + abs(us(c,ubj)) * dtime(c)        + 1.e-14_r8
        zghostr(2) = zi(c,ubj) + abs(us(c,ubj)) * dtime(c) * 2._r8+ 2.e-14_r8

        ughostr(1) = us(c,ubj)
        ughostr(2) = us(c,ubj)

        call backward_advection((/zghostl, zi(c, lbn(c)-1:ubj),zghostr/), (/ughostl, us(c, lbn(c)-1:ubj), ughostr/),  dtime(c), zold(0:length))

        if(.not. is_ascending_vec(zold(0:length)))then
           halfdt_col(c) = .true.
           cycle
        endif

        !create the cumulative mass curve
        do ntr = 1, ntrcs
           !left boundary ghost grids
           j = 0
           mass_curve(j, ntr) = 0._r8
           j = 1
           mass_curve(j, ntr) = inflx_top(c, ntr)*dtime(c)
           j = 2
           mass_curve(j, ntr) = inflx_top(c, ntr)*dtime(c)

           !regular grids
           do k = lbn(c), ubj
              j = k - lbn(c) + 3
              mass_curve(j, ntr) = trcin(c,k, ntr)*dz(c,k)
           enddo

           !right ghost grids
           if(inflx_bot(c,ntr)==0._r8)then
              j = ubj - lbn(c) + 4
              mass_curve(j, ntr) = trcin(c,ubj, ntr)*(zghostr(1)-zi(c,ubj))

              j = ubj - lbn(c) + 5
              mass_curve(j, ntr) = trcin(c,ubj, ntr)*(zghostr(2)-zghostr(1))
           else
              j = ubj - lbn(c) + 4
              mass_curve(j, ntr) = inflx_bot(c, ntr) * dtime(c)

              j = ubj - lbn(c) + 5
              mass_curve(j, ntr) = inflx_bot(c, ntr) * dtime(c)
           endif
        enddo
        !compute cumulative mass curve
        call cumsum(mass_curve(0:lengthp2,1:ntr), cmass_curve(0:lengthp2, 1:ntr),idim=1)

        !do mass interpolation
        do ntr = 1, ntrcs
           call pchip_polycc((/zghostl,zi(c,lbn(c)-1:ubj),zghostr/), cmass_curve(0:lengthp2, ntr), di(0:lengthp2))

           call pchip_interp((/zghostl,zi(c,lbn(c)-1:ubj),zghostr/), cmass_curve(0:lengthp2, ntr), di(0:lengthp2),&
                zold(0:length), cmass_new(0:length, ntr))

           !ensure mass is increasing monotonically
           call asc_sort_vec(cmass_new(0:length,ntr))

           !ensure no negative leaching
           call cmass_mono_smoother(cmass_new(0:length, ntr),cmass_curve(ubj-lbn(c)+3, ntr))

           !diagnose the leaching flux
           if(present(leaching_mass))then
              leaching_mass(c, ntr) = cmass_curve(ubj-lbn(c)+3, ntr)-cmass_new(length, ntr) !add the numerical error to leaching
           endif

           !obtain the grid concentration
           call cumdif(cmass_new(0:length, ntr), mass_new(0:length, ntr))
           do k = lbn(c), ubj
              j = k - lbn(c) + 1
              !correct for small negative values
              if(mass_new(j, ntr)<0._r8)then
                 write(iulog,*)j,mass_new(j, ntr),cmass_new(j, ntr),cmass_new(j-1, ntr)
                 call endrun('negative tracer '//errMsg(__FILE__, __LINE__))
                 if(present(leaching_mass))then
                    leaching_mass(c, ntr) = leaching_mass(c, ntr)+mass_new(j, ntr) !add the numerical error to leaching
                 endif
                 mass_new(j, ntr)=mass_curve_correct(mass_new(j, ntr))
              endif
              trcou(c,k, ntr)=mass_new(j, ntr)/dz(c,k)
           enddo
        enddo

     enddo
     call Extra_inst%DDeallocate()
   end subroutine semi_lagrange_adv_backward
   !-------------------------------------------------------------------------------
   subroutine cmass_mono_smoother(cmass,mass_thc)
     !
     ! !DESCRIPTION:
     ! assuming cmass is sorted as ascending vector, make sure no mass is greater than mass_thc
     !
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:), intent(inout) :: cmass
     real(r8), intent(in)                  :: mass_thc
     ! !LOCAL VARIABLES:
     integer :: n , j
     character(len=32) :: subname = 'cmass_mono_smoother'

     n = size(cmass)
     do j = n, 1
        if(cmass(j)>=mass_thc)then
           cmass(j) = mass_thc
        else
           exit
        endif
     enddo
   end subroutine cmass_mono_smoother
   !-------------------------------------------------------------------------------
   
   function is_ascending_vec(zcor)result(ans)
     !
     ! DESCRIPTION:
     ! check if it is an ascending array

     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:), intent(in) :: zcor

     ! !LOCAL VARIABLES:
     logical :: ans
     integer :: j, n
     character(len=32) :: subname= 'is_ascending_vec'

     n = size(zcor)
     ans = .true.
     do j = 2 , n
        if(zcor(j)<zcor(j-1))then
           ans=.false.
           exit
        endif
     enddo
   end function is_ascending_vec

   !-------------------------------------------------------------------------------
   function mass_curve_correct(mass_curve)result(ans)
     !
     ! !DESCRIPTION:
     ! Correct random truncation error induced negative values

     implicit none
     ! !ARGUMENTS:
     real(r8), intent(in) :: mass_curve

     ! !LOCAL VARIABLES:
     real(r8) :: ans
     real(r8) :: eps = 1.e-15_r8
     character(len=32) :: subname = 'mass_curve_correct'

     if(abs(mass_curve)<eps)then
        ans = 0._r8
     else
        ans = mass_curve
     endif
     return
   end function mass_curve_correct
   !-------------------------------------------------------------------------------
   subroutine backward_advection(zi, us,  dtime, zold)
     !
     ! !DESCRIPTION:
     ! do backward trajectory track of the boundaries
     !
     ! !USES:
     use ODEMod, only : ode_rk4, ode_rk2
     implicit none
     ! !ARGUMENTS:
     real(r8), dimension(:), intent(in) :: zi     !boundary interfaces,  [0 : n+4], 2:n+2 are to be advected
     real(r8), dimension(:), intent(in) :: us     !interface velocities, [0 : n+4] including ghost cells
     real(r8),               intent(in) :: dtime  !time stepping
     real(r8), dimension(:), intent(out):: zold   !the starting point of the interfaces, [0:n]

     ! !LOCAL VARIABLES:
     integer :: neq  ! number of equations
     real(r8):: time
     character(len=32) :: subname = 'backward_advection'

     SHR_ASSERT_ALL((size(zi)        == size(us)),            errMsg(__FILE__,__LINE__))
     SHR_ASSERT_ALL((size(zi)        == size(zold)+4),        errMsg(__FILE__,__LINE__))

     neq = size(zold)
     call Extra_inst%AAssign(zi,us)
     time =0._r8
     call ode_rk2(trajectory, zi(3:neq+2), neq, time, dtime, zold)

   end subroutine backward_advection
   
   !-------------------------------------------------------------------------------
   
   subroutine trajectory(y0, dt, ti, neq, dxdt)
     !
     ! !DESCRIPTION:
     ! update the trajectory

     ! !USES:
     use InterpolationMod, only : Lagrange_interp
     implicit none
     ! !ARGUMENTS:
     real(r8), intent(in)  :: y0(neq)
     real(r8), intent(in)  :: dt
     real(r8), intent(in)  :: ti
     integer,  intent(in)  :: neq
     real(r8), intent(out) :: dxdt(neq)

     ! !LOCAL VARIABLES:
     integer            :: j
     integer, parameter :: pn = 1
     real(r8)           :: ui(neq)
     character(len=32)  :: subname ='trajectory'

     call Lagrange_interp(pn, Extra_inst%zi(1:Extra_inst%nlen), Extra_inst%us(1:Extra_inst%nlen), y0, ui)
     do j = 1, neq
        dxdt(j) = -ui(j)
     enddo

   end subroutine trajectory


end module TransportMod
