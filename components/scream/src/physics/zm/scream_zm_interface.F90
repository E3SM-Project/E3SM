#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_zm_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR, c_float
  use physics_utils, only: r8 => rtype, rtype8, itype, btype
  use zm_conv,       only: zm_convr, zm_conv_evap, zm_convi, convtran, momtran,zmconv_readnl
 
  implicit none
#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#else
# define c_real c_float
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type #endif
#endif



  public :: zm_init_f90
  public :: zm_main_f90
  public :: zm_finalize_f90

  real   :: test
  
  integer(kind=c_int) :: pcols  = 32
  integer(kind=c_int) :: pver   = 73
  integer(kind=c_int) :: pverp  = 75
  
 ! real(kind=c_real) :: pcols = 32
 ! real(kind=c_real) :: pver = 75
  real(kind=c_real) :: cpair  !=    1004.64000000000
  real(kind=c_real) :: gravit !=    9.80616000000000
  real(kind=c_real) :: latvap !=    2501000.00000000
  real(kind=c_real) :: plevp = 0.0 

!  real(kind=c_real) :: pverp = 0

contains

  !====================================================================!

  subroutine zm_init_f90 (limcnv_in, no_deep_pbl_in) bind(c)
    integer, intent(in)                :: limcnv_in
    logical, intent(in), optional      :: no_deep_pbl_in    
    call zm_convi(limcnv_in, no_deep_pbl_in)
    call zmconv_readnl()

  end subroutine zm_init_f90
  !====================================================================!
subroutine zm_main_f90(lchnk   ,ncol    , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd,     & 
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,landfrac,hu_nm1  , &
                    cnv_nm1 ,tm1     ,qm1     ,t_star  ,q_star  , &
                    dcape   ,q       ,tend_s  ,tend_q  ,cld     , &
                    snow    ,ntprprd ,ntsnprd ,flxprec ,flxsnow, &
                    ztodt   , pguall , pgdall , icwu   , ncnst, fracis) bind(c)
   integer :: lchnk
   integer lengath
   integer, intent(in) :: ncol
   real(kind=c_real), intent(inout) :: t(pcols,pver) ! State array  kg/kg
   real(kind=c_real) :: qh(pcols,pver) 
   real(kind=c_real), intent(out) :: prec(pcols)
   real(kind=c_real), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(kind=c_real), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer maxg(pcols)                        ! wg gathered values of maxi.
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer mx(pcols) 
!   
!   real(r8) landfracg(pcols)            ! wg grid slice of landfrac  
   real(kind=c_real) delt                     ! length of model time-step in seconds.
!   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)
!   
   real(kind=c_real), intent(in) :: pblh(pcols)
   real(kind=c_real), intent(in) :: geos(pcols)
   real(kind=c_real), intent(in) :: zi(pcols,pver+1)
   real(kind=c_real), intent(in) :: zm(pcols,pver)
   real(kind=c_real), intent(in) :: pap(pcols,pver)     
   real(kind=c_real), intent(in) :: paph(pcols,pver+1)
   real(kind=c_real), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(kind=c_real), intent(in) :: tpert(pcols)
   real(kind=c_real), intent(in) :: tm1(pcols,pver)       ! grid slice of temperature at mid-layer.
   real(kind=c_real), intent(in) :: qm1(pcols,pver)       ! grid slice of specific humidity.
   real(kind=c_real), intent(in) :: landfrac(pcols) ! RBN Landfrac
   

   real(kind=c_real), intent(in) :: t_star(pcols,pver) ! intermediate T between n and n-1 time step
   real(kind=c_real), intent(in) :: q_star(pcols,pver) ! intermediate q between n and n-1 time step
   
   real(kind=c_real), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(kind=c_real), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(kind=c_real), intent(out) :: mcon(pcols,pverp)
   real(kind=c_real), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(kind=c_real), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(kind=c_real), intent(out) :: cme(pcols,pver)
   real(kind=c_real), intent(out) :: cape(pcols)        ! w  convective available potential energy.
   real(kind=c_real), intent(out) :: zdu(pcols,pver)
   real(kind=c_real), intent(out) :: rprd(pcols,pver)     ! rain production rate
   real(kind=c_real), intent(out) :: mu(pcols,pver)
   real(kind=c_real), intent(out) :: md(pcols,pver)
   real(kind=c_real), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(kind=c_real), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(kind=c_real), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(kind=c_real), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(kind=c_real), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   real(kind=c_real), intent(out) :: ql(pcols,pver)
   real(kind=c_real), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(kind=c_real), intent(out) :: dcape(pcols)           ! output dynamical CAPE
!   
   real(kind=c_real), intent(inout) :: hu_nm1 (pcols,pver)
   real(kind=c_real), intent(inout) :: cnv_nm1 (pcols,pver)
!
!
   real(kind=c_real) q(pcols,pver)              
!   
!!   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
!!Used for convtran exclusively 
   real(kind=c_real), intent(in) :: fracis(pcols,pver,ncnst)  
   real(kind=c_real) :: fake_dpdry(pcols,pver)       
   real(kind=c_real) :: fake_dqdt(pcols,pver,ncnst)  ! Tracer tendency array
!
   integer :: il1g
   integer :: nstep             
!!
!
   real(kind=c_real) :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(kind=c_real) :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(kind=c_real),intent(inout), dimension(pcols,pver) :: tend_s    ! heating rate (J/kg/s)
   real(kind=c_real),intent(inout), dimension(pcols,pver) :: tend_q    ! heating rate (J/kg/s)
   real(kind=c_real), intent(inout), dimension(pcols,pver) :: cld
   real(kind=c_real), intent(inout), dimension(pcols)   :: snow         ! snow from ZM convection 
   real(kind=c_real) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(kind=c_real) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(kind=c_real), intent(out), dimension(pcols,pver) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(kind=c_real), intent(out), dimension(pcols,pver) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(kind=c_real), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(kind=c_real) :: pguall(pcols, pver, 2)
   real(kind=c_real) :: pgdall(pcols, pver, 2)
   real(kind=c_real) :: icwu(pcols,pver, 2)
   real(kind=c_real) :: icwd(pcols,pver, 2)
   real(kind=c_real) :: seten(pcols, pver)
   integer, intent(in) :: ncnst 
!
   logical :: domomtran(ncnst)
   logical :: doconvtran(ncnst)
   fake_dpdry(:,:) = 0._r8
   il1g = 1
   
   call zm_convr(lchnk   ,ncol    , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,landfrac,hu_nm1  , &
                    cnv_nm1 ,tm1     ,qm1     ,t_star  ,q_star, dcape)
   
   call zm_conv_evap(ncol    ,lchnk, &
                     t       ,pap     ,dpp     ,q     , &
                     tend_s,     tend_s_snwprd ,tend_s_snwevmlt , &
                     tend_q   ,rprd,       cld ,ztodt            , &
                     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow )



  call momtran(lchnk, ncol, &
                   domomtran,q       ,ncnst   ,mu      ,md    , &
                   du      ,eu      ,ed      ,dp      ,dsubcld , &
                   jt      ,mx      ,ideep   ,il1g    ,lengath    , &
                   nstep   ,fake_dqdt    ,pguall     ,pgdall, icwu, icwd, ztodt, seten    )

  call convtran(lchnk   , &
                   doconvtran,q       ,ncnst   ,mu      ,md      , &
                   du      ,eu      ,ed      ,dp      ,dsubcld , &
                   jt      ,mx      ,ideep   ,il1g    , lengath    , &
                   1   ,fracis  ,fake_dqdt,  fake_dpdry   )

   end subroutine zm_main_f90
  !====================================================================!
  subroutine zm_finalize_f90 () bind(c)


  end subroutine zm_finalize_f90
  !====================================================================!

end module scream_zm_interface_mod
