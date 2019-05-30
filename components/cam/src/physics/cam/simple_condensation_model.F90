!-----------------------------------------------------------------------------------------
! This version of the Rasch-Kristjansson macro-physics scheme includes
!  - the in-cloud condensation term;
!  - the clear-sky liquid advection term;
!  - the cloud expansion term.
!
! Multiple versions of the second and third terms are included for testing.
!
! History: 
!  - first version by Hui Wan (PNNL, 2014-2015)
!  - revised for the SciDAC Convergence project, Hui Wan (PNNL, 2017)
!-----------------------------------------------------------------------------------------
module simple_condensation_model

  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use ppgrid,       only: pver, pcols

  implicit none
  private
  public :: simple_RKZ_tend
  public :: simple_RKZ_init

contains

  !------------------------------------------------------------------------
  ! Initialization of the RKZ simple condensation model.
  ! Currently this contains just registration of output variables.
  !------------------------------------------------------------------------
  subroutine simple_RKZ_init()

    use cam_history,    only: addfld

    implicit none

    ! There are two steps of code changes for adding an variable (physical quantity) to the model output file.
    !  1. add a "call addfld" following examples in this subroutine to register the variable during model initialization.
    !  2. add a "call outfld" in subroutine simple_RKZ_tend following examples there in to send the values at each time
    !     step to the corresponding data structure for further processing (e.g., time averaging if applicable) and 
    !     for output.
    ! After these changes are made, re-compile the model, and include the newly added variable name to 
    ! the namelist variable fincl1 (or fincl2, fincl3, ...). This latter step is usually done by revising
    ! your run script.

    ! The arguements of subroutine addfld are:
    !  -  variable name (character, max length = 24)
    !  -  dimension names (character, for variables that have the shape of (pcols,pver), use "(/'lev'/)")
    !  -  time averaging flag (character, 'I' for instantaneous, 'A' for average)
    !  -  unit (character)
    !  -  long name (character)

    call addfld ('RKZ_ql',  (/'lev'/), 'I','kg/kg','grid-box mean ql used in the simple RKZ scheme')
    call addfld ('RKZ_qv',  (/'lev'/), 'I','kg/kg','grid-box mean qv used in the simple RKZ scheme')
    call addfld ('RKZ_Tbf', (/'lev'/), 'I','K','grid-box mean temperature after the simple RKZ scheme')
    call addfld ('RKZ_Taf', (/'lev'/), 'I','K','grid-box mean temperature before the simple RKZ scheme')

    call addfld ('RKZ_qsat',    (/'lev'/), 'I', 'kg/kg',   'saturation specific humidity used in the simple RKZ scheme')
    call addfld ('RKZ_dqsatdT', (/'lev'/), 'I', 'kg/kg/K', 'derivative of saturation specific humidity wrt temperature used in the simple RKZ scheme')

    call addfld ('RKZ_RH',   (/'lev'/), 'I','1','grid-box mean relative humidity')
    call addfld ('RKZ_f',    (/'lev'/), 'I','1','cloud fraction')
    call addfld ('RKZ_fac',  (/'lev'/), 'I','1','cloud fraction after macro- and microphysics calculation in the previous time step')
    call addfld ('RKZ_dfacdRH',  (/'lev'/), 'I','1','derivative of cloud fraction wrt relative humidity after macro- and microphysics calculation in the previous time step')

    call addfld ('RKZ_dfdRH',(/'lev'/), 'I','1','derivative of cloud fraction wrt relative humidity')
    call addfld ('RKZ_dlnfdRH',(/'lev'/), 'I','1','derivative of logrithmic cloud fraction wrt relative humidity')

    call addfld ('RKZ_dfdt', (/'lev'/), 'I','1','estimated cloud fraction tendency')

    call addfld ('RKZ_ql_incld',(/'lev'/), 'I','kg/kg','in-cloud ql used for calculating term C')

    call addfld ('RKZ_term_A',(/'lev'/), 'I','kg/kg/s','grid-box mean condensation rate, term A')
    call addfld ('RKZ_term_B',(/'lev'/), 'I','kg/kg/s','grid-box mean condensation rate, term B')
    call addfld ('RKZ_term_C',(/'lev'/), 'I','kg/kg/s','grid-box mean condensation rate, term C')

    call addfld ('RKZ_Al',(/'lev'/), 'I','kg/kg/s', 'grid-box mean ql tendency caused by other processes')
    call addfld ('RKZ_Av',(/'lev'/), 'I','kg/kg/s', 'grid-box mean qv tendency caused by other processes')
    call addfld ('RKZ_AT',(/'lev'/), 'I','K/s',     'grid-box mean temperature tendency caused by other processes')

    call addfld ('RKZ_zqme', (/'lev'/), 'I', 'kg/kg/s', 'condensation rate before limiters in the simple RKZ scheme')
    call addfld ('RKZ_qme',  (/'lev'/), 'I', 'kg/kg/s', 'condensation rate after limiters in the simple RKZ scheme')

    call addfld ('RKZ_qme_lm4_qv',  (/'lev'/), 'I', 'kg/kg/s', 'condensation rate difference before and after limiter 4 for qv in the simple RKZ scheme')
    call addfld ('RKZ_qme_lm4_ql',  (/'lev'/), 'I', 'kg/kg/s', 'condensation rate difference before and after limiter 4 for ql in the simple RKZ scheme')
    call addfld ('RKZ_qme_lm5_ps',  (/'lev'/), 'I', 'kg/kg/s', 'condensation rate difference before and after limiter 5 with positive sign in the simple RKZ scheme')
    call addfld ('RKZ_qme_lm5_ng',  (/'lev'/), 'I', 'kg/kg/s', 'condensation rate difference before and after limiter 5 with negative sign in the simple RKZ scheme')

    call addfld ('RKZ_qvneg_lm4',  (/'lev'/), 'I', 'kg/kg', 'negative qv clipped by limiter 4 in the simple RKZ scheme')
    call addfld ('RKZ_qlneg_lm4',  (/'lev'/), 'I', 'kg/kg', 'negative ql clipped by limiter 4 in the simple RKZ scheme')

    call addfld ('RKZ_stend_lm4_qv',  (/'lev'/), 'I', 'J/kg', 'dry static energy tendency associated with limiter 4 for clipping negative qv in the simple RKZ scheme')
    call addfld ('RKZ_stend_lm4_ql',  (/'lev'/), 'I', 'J/kg', 'dry static energy tendency associated with limiter 4 for clipping negative ql in the simple RKZ scheme')

    call addfld ('RKZ_lmt4_flg',  (/'lev'/), 'I', '1', 'flag indicating cells in which condition is met to trigger limiter 4')
    call addfld ('RKZ_lmt5_flg',  (/'lev'/), 'I', '1', 'flag indicating cells in which condition is met to trigger limiter 5')
    call addfld ('RKZ_lmt45_flg', (/'lev'/), 'I', '1', 'flag indicating cells in which condition is met to trigger both limiter 4 and limiter 5')

    call addfld ('RKZ_gam',  (/'lev'/), 'I', '1', 'coefficient for (dqsat/dT)*(lv/cp) or alpha*beta used in simple condensation model')



  end subroutine simple_RKZ_init

  !------------------------------------------------------------------------
  ! Calculate condensation rate and the resulting tendencies of the model 
  ! state variables
  !------------------------------------------------------------------------
  subroutine simple_RKZ_tend(state, ptend, tcwat, qcwat, lcwat, ast, qmeold, astwat, dfacdRH,&
                             dtime, ixcldliq, &
                             rkz_cldfrc_opt, &
                             rkz_term_A_opt, &
                             rkz_term_B_opt, &
                             rkz_term_C_opt, &
                             rkz_term_C_ql_opt, & 
                             rkz_term_C_fmin, & 
                             rkz_zsmall_opt, &
                             rkz_lmt5_opt, &
                             l_rkz_lmt_2, &
                             l_rkz_lmt_3, &
                             l_rkz_lmt_4, &
                             l_rkz_lmt_5  &
                             )

  use shr_kind_mod, only: r8=>shr_kind_r8
  use constituents, only: pcnst
  use physics_types,only: physics_state, physics_ptend, physics_ptend_init
  use time_manager,  only: get_nstep
  use wv_saturation, only: qsat_water
  use physconst,     only: latvap, cpair
  use simple_cloud_fraction, only: smpl_frc
  use cam_history,   only: outfld

  implicit none
  !
  ! arguments
  !
  type(physics_state), intent(in), target    :: state       ! State variables
  type(physics_ptend), intent(out)           :: ptend       ! Package tendencies

  real(r8), intent(inout) :: tcwat(:,:)       ! temperature after macro- and microphysics calculation in the previous time step
  real(r8), intent(inout) :: qcwat(:,:)       ! qv          after macro- and microphysics calculation in the previous time step
  real(r8), intent(inout) :: lcwat(:,:)       ! ql          after macro- and microphysics calculation in the previous time step
  real(r8), intent(inout) :: astwat(:,:)      ! f           after macro- and microphysics calculation in the previous time step 
  real(r8), intent(inout) :: dfacdRH(:,:)     ! df/dRH      after macro- and microphysics calculation in the previous time step 
                                              !             where f is the cloud fraction and RH the relative humidity

  real(r8), intent(inout) ::   ast(:,:)       ! cloud fraction
  real(r8), intent(inout) ::   qmeold(:,:)    ! total condensation rate calculation in previous time step 

  real(r8), intent(in) :: dtime               ! Set model physics timestep
  integer,  intent(in) :: ixcldliq            ! constituent index 

  integer,  intent(in) :: rkz_cldfrc_opt       ! cloud fraction scheme 
  integer,  intent(in) :: rkz_term_A_opt
  integer,  intent(in) :: rkz_term_B_opt
  integer,  intent(in) :: rkz_term_C_opt
  integer,  intent(in) :: rkz_term_C_ql_opt
  real(r8), intent(in) :: rkz_term_C_fmin

  integer,  intent(in) :: rkz_zsmall_opt       
  integer,  intent(in) :: rkz_lmt5_opt

  logical,  intent(in) :: l_rkz_lmt_2
  logical,  intent(in) :: l_rkz_lmt_3
  logical,  intent(in) :: l_rkz_lmt_4
  logical,  intent(in) :: l_rkz_lmt_5

  ! tmp work arrays

  real(r8) :: rdtime                 ! 1/dtime

  integer  :: lchnk                  ! index of (grid) chunk handled by this call of the subroutine
  integer  :: ncol                   ! number of active columns in this chunk for which calculations will be done
  integer  :: nstep                  ! time step index
  integer  :: i,k                    ! loop indices for column and vertical layer
  logical  :: lq(pcnst)              ! logical array used when calling subroutine physics_ptend_init indicating 
                                     ! which tracers are affected by this parameterization

  real(r8) :: ast_old(pcols,pver)    ! cloud fraction of previous time step before condensation

  real(r8) :: qsat(pcols,pver)       ! saturation specific humidity
  real(r8) :: esl(pcols,pver)        ! saturation vapor pressure (output from subroutine qsat_water, not used)
  real(r8) :: dqsatdT(pcols,pver)    ! dqsat/dT
  real(r8) :: gam(pcols,pver)        ! L/cpair * dqsat/dT

  real(r8) :: rhu00                  ! threshold grid-box-mean RH used in the diagnostic cldfrc scheme
  real(r8) :: rhgbm(pcols,pver)      ! grid box mean relative humidity
  real(r8) :: dastdRH(pcols,pver)    ! df/dRH where f is the cloud fraction and RH the relative humidity
  real(r8) :: dlnastdRH(pcols,pver)  ! dlnf/dRH where lnf is the logrithm of the cloud fraction and RH the relative humidity

  real(r8) :: qtend(pcols,pver)      ! Moisture tendency caused by other processes
  real(r8) :: ltend(pcols,pver)      ! liquid condensate tendency caused by other processes
  real(r8) :: ttend(pcols,pver)      ! Temperature tendency caused by other processes

  real(r8) :: qme   (pcols,pver)     ! total condensation rate

  real(r8) :: qmebf   (pcols,pver)   ! total condensation rate before appling limiter 
  real(r8) :: qmedf   (pcols,pver)   ! difference of total condensation rate before and after limiter application
  real(r8) :: tmp     (pcols,pver)   ! temporary work array 

  real(r8) :: term_A (pcols,pver)
  real(r8) :: term_B (pcols,pver)
  real(r8) :: term_C (pcols,pver)
  real(r8) :: rdenom (pcols,pver)    ! 1/denominator in the final expression of total grid-box mean condensation rate 
                                     ! when option 2 is used for term C. Note that in this case, 1/denominator will be 
                                     ! multipled to all three terms (A, B, and C).

  real(r8) :: ql_incld(pcols,pver)   ! in-cloud liquid concentration
  real(r8) :: dfdt    (pcols,pver)   ! df/dt where f is the cloud fraction
  real(r8) :: zforcing (pcols,pver)
  real(r8) :: zc3      (pcols,pver)

  real(r8) :: zqvnew(pcols,pver)     ! qv at new time step if total condenation rate is not limited. Might be negative.
  real(r8) :: zqlnew(pcols,pver)     ! ql at new time step if total condenation rate is not limited. Might be negative.
  real(r8) :: zlim (pcols,pver)
  real(r8) :: zsmall                 ! bottom boundary for ql and qv used in limiter 4 and 5

  real(r8) :: zqvneg(pcols,pver)     ! Save the negative qv for analysis.
  real(r8) :: zqlneg(pcols,pver)     ! Save the negative ql for analysis.
  real(r8) :: zstend(pcols,pver)     ! Save the dry static energy change before and after limiter 

  real(r8) :: flag_lmt4         (pcols,pver)  ! 1 = condition is met for triggering limiter 4
  real(r8) :: flag_lmt5         (pcols,pver)  ! 1 = condition is met for triggering limiter 5
  real(r8) :: flag_lmt4_and_lmt5(pcols,pver)  ! 1 = condition is met for triggering both limiter 4 and limiter 5

  logical :: lcondition4_qv(pcols,pver)  ! condition is met for triggering limiter 4 to prevent negative qv
  logical :: lcondition4_ql(pcols,pver)  ! condition is met for triggering limiter 4 to prevent negative ql
  logical :: lcondition5   (pcols,pver)  ! condition is met for triggering limiter 5
  real(r8),parameter :: pi = 3.141592653589793
  real(r8),parameter :: fmax= 1._r8 ! upper limit for the cloud fraction
  real(r8) :: dastdb(pcols,pver)
  real(r8) :: rhlim

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  rdtime = 1._r8/dtime
  ncol  = state%ncol
  nstep = get_nstep()
  lchnk = state%lchnk  ! needed by "call outfld" for model output

  call outfld('RKZ_fac',     astwat,    pcols, lchnk)
  call outfld('RKZ_dfacdRH', dfacdRH,   pcols, lchnk)

  ! Calculate saturation specific humidity (qsat) and its derivative wrt temperature (dqsatdT)

  do k=1,pver
  do i=1,ncol
     call qsat_water( state%t(i,k), state%pmid(i,k), esl(i,k), qsat(i,k), gam(i,k), dqsatdT(i,k) )
  end do
  end do
  call outfld('RKZ_qsat',    qsat,    pcols, lchnk)
  call outfld('RKZ_dqsatdT', dqsatdT, pcols, lchnk)
  call outfld('RKZ_gam',     gam,     pcols, lchnk)

  ! Cloud fraction (ast): save old values, diagnose new values

  if (nstep > 1) ast_old(:ncol,:pver) = ast(:ncol,:pver)

  if ( rkz_term_C_opt.eq.4 .and. ((rkz_cldfrc_opt.ne.1).and.(rkz_cldfrc_opt.ne.3)) ) then
     write(iulog,*) "rkz_term_C_opt = 4 can be used only when rkz_cldfrc_opt = 1 or 3."
     write(iulog,*) "Your choice was: rkz_term_C_opt = ",rkz_term_C_opt, ", rkz_cldfrc_opt = ", rkz_cldfrc_opt, ".Abort."
     call endrun
  end if

  call  smpl_frc( state%q(:,:,1), state%q(:,:,ixcldliq), qsat,      &! all in
                  ast, rhu00, rhgbm, dastdRH, dlnastdRH,            &! inout, out, out
                  rkz_cldfrc_opt, 0.5_r8, 0.5_r8, pcols, pver, ncol )! all in
 
 !!!add bounded condition for dln(f)/dt (case1)!!!
 ! where ( ast(:ncol,:pver) .le. rkz_term_C_fmin )
 !   dlnastdRH(:ncol,:pver) = 300._r8
 ! end where
 !!!add bounded condition for dln(f)/dt (case2)!!!
  where ( ast(:ncol,:pver) .le. rkz_term_C_fmin )
    dlnastdRH(:ncol,:pver) = 100._r8
  end where

 !Output the grid-box-mean relative humidity (rhgbm)
  call outfld('RKZ_RH', rhgbm, pcols, lchnk)

  call outfld('RKZ_qv',    state%q(:,:,1),        pcols, lchnk)
  call outfld('RKZ_ql',    state%q(:,:,ixcldliq), pcols, lchnk)
  call outfld('RKZ_Taf',   state%t(:,:), pcols, lchnk)
  call outfld('RKZ_Tbf',   tcwat,        pcols, lchnk)

  call outfld('RKZ_f',     ast,      pcols, lchnk)
  call outfld('RKZ_dfdRH', dastdRH,  pcols, lchnk)
  call outfld('RKZ_dlnfdRH', dlnastdRH,  pcols, lchnk)

  !===================================================================
  ! Condensation/evaporation rate
  !===================================================================
  if (nstep > 1) then

     !--------------------------------------------------------------------------------------------------
     ! Diagose "other forcing" on qv, ql, and temperature. I.e., tendencies caused by other processes.

     qtend(:ncol,:pver) = ( state%q(:ncol,:pver,1)        - qcwat(:ncol,:pver) )*rdtime
     ttend(:ncol,:pver) = ( state%t(:ncol,:pver)          - tcwat(:ncol,:pver) )*rdtime
     ltend(:ncol,:pver) = ( state%q(:ncol,:pver,ixcldliq) - lcwat(:ncol,:pver) )*rdtime

     call outfld('RKZ_Av', qtend, pcols, lchnk)
     call outfld('RKZ_Al', ltend, pcols, lchnk)
     call outfld('RKZ_AT', ttend, pcols, lchnk)
     
     !------------------------------------------------------------------------------------
     ! Term A: condensation in cloudy portion of a grid box, weighted by cloud fraction

     select case (rkz_term_A_opt) 
     case(0) ! omit this term
        term_A(:ncol,:pver) = 0._r8

     case(1)
        term_A(:ncol,:pver) = ast(:ncol,:pver)     &
                           *( qtend(:ncol,:pver) - dqsatdT(:ncol,:pver)*ttend(:ncol,:pver) ) &
                           /( 1._r8 + gam(:ncol,:pver) )
     case(21)
     !use the RK98 formula,the first part on the right side of equation.           
     !note that gam=alpha*beta
        term_A(:ncol,:pver) = ast(:ncol,:pver)     &
                           *( qtend(:ncol,:pver) - dqsatdT(:ncol,:pver)*ast(:ncol,:pver)*ttend(:ncol,:pver) ) &
                           /( 1._r8 + ast(:ncol,:pver)*gam(:ncol,:pver))

     case default
         write(iulog,*) "Unrecognized value of rkz_term_A_opt:",rkz_term_A_opt,". Abort."
         call endrun
     end select

     !--------------------------------------------------------------------------------------
     ! Term B: condensation in cloud-free portion of a grid box in response to cloud liquid tencendy
     ! caused by other processes.

     select case (rkz_term_B_opt) 
     case(0) ! omit this term
        term_B(:ncol,:pver) = 0._r8

     case(1) 
     ! Following Zhang et al. (2003), assume the in-cloud A_l equals the grid-box mean A_l.
     ! Hence, - (\overline{A_l} - f \hat{A_l}) = - (1-f)\overline{A_l}.

        term_B(:ncol,:pver) = - (1._r8 - ast(:ncol,:pver))*ltend(:ncol,:pver)

     case(2) 
     ! For testing only: term B = - 0.5*\overline{A_l}

        term_B(:ncol,:pver) = - 0.5_r8*ltend(:ncol,:pver)

     case(3)
     ! For testing only: Following Zhang et al. (2003), assume the in-cloud A_l
     ! equals the grid-box mean A_l.! Hence, - (\overline{A_l} - f \hat{A_l}) =
     ! - (1-f)\overline{A_l}.However, if overline{A_l} <0, term_B=0.0

        term_B(:ncol,:pver) = - (1._r8 - ast(:ncol,:pver))*ltend(:ncol,:pver)
  
        where( ltend(:ncol,:pver) .lt. 0._r8)
          term_B(:ncol,:pver) = 0._r8
        end where

     case(20) 
     !!Origionally, there are only two terms on the right side of RK98 formula, 
     !!thus, term B should be always zero at any time
        term_B(:ncol,:pver) = 0._r8

     case default
         write(iulog,*) "Unrecognized value of rkz_term_B_opt:",rkz_term_B_opt,". Abort."
         call endrun
     end select

     !---------------------------------------------------------------------------------------------
     ! Term C: condensation in cloud-free portion of a grid box, related to cloud expansion/erosion 

     ! Calculate in-cloud liquid concentration:
     ! 1 or 11: ql_incld=constant      , only for testing
     ! 3 or 13: ql_incld=qv            , only for testing
     ! 4 or 14: ql_incld=ql            , only for testing
     ! 7 or 17: ql_incld=ql/max(f,fmin)

     SELECT CASE (rkz_term_C_ql_opt)
     CASE (1,11)
       ql_incld(:ncol,:pver) = 1e-5_r8

     CASE (3,13)
       ql_incld(:ncol,:pver) = state%q(:ncol,:pver,1)

     CASE (4,14)
       ql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)

     CASE (5,15)
       ql_incld(:ncol,:pver) = 0.5_r8*state%q(:ncol,:pver,ixcldliq)

     CASE (6,16)
       ql_incld(:ncol,:pver) = 2.0_r8*state%q(:ncol,:pver,ixcldliq)

     CASE (7,17)
       ql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)/max(ast(:ncol,:pver),rkz_term_C_fmin)

     CASE (8,18)
       ql_incld(:ncol,:pver) = lcwat(:ncol,:pver)/max(ast(:ncol,:pver),rkz_term_C_fmin)

     CASE (9,19)
       ql_incld(:ncol,:pver) = lcwat(:ncol,:pver)/max(astwat(:ncol,:pver),rkz_term_C_fmin)

     CASE DEFAULT
       write(iulog,*) "Unrecognized value of rkz_term_C_ql_opt:",rkz_term_C_ql_opt,". Abort."
       call endrun
     END SELECT

     ! term C = ql_incld * df/dt (or, term C = ql * dln(f)/dt if rkz_term_C_opt = 4)

     select case (rkz_term_C_opt)
     case(0)  ! omit this term
        term_C(:ncol,:pver) = 0._r8

     case(1) ! Use simple finite-difference to approximate df/dt
        term_C(:ncol,:pver) = ql_incld(:ncol,:pver)* rdtime*( ast(:ncol,:pver) - ast_old(:ncol,:pver) )

     case(2) ! Use chain rule: df/dt = df/dRH * dRH/dt = df/dRH * (dRH/dqv + dRH/dT)

        ! zforcing is the "forcing" term, i.e., grid box cooling and/or moistening caused by 
        ! processes other than condensation. It appears on the nominator of the expression 
        ! for the grid-box-mean condensation rate. Using the notation of Zhang et al. (2003),
        !
        !     zforcing = C_alpha * A_v - C_beta * A_T
        !
        ! Using the definition of C_alpha and C_beta, we have
        !
        !     zforcing = 1/qsat * ( A_v - RH * d(qsat)/dT * A_T)
        !
        ! Note that A_v = qtend in this subr., and A_T = ttend. 
        
        zforcing(:ncol,:pver) = ( qtend(:ncol,:pver)  &
                                 -ttend(:ncol,:pver)*rhgbm(:ncol,:pver)*dqsatdT(:ncol,:pver) ) &
                               /qsat(:ncol,:pver)

        ! zc3 is the term C_gamma in Zhang et al. (2003). It appears as part of the denominator
        ! of the final expression for total grid-box-mean condensation rate (i.e., qme in this subr.)
        !
        ! C_gamma = 1/qsat + Lv/Cp * (qv/qsat^2) * d(qsat)/dT
        !         = 1/qsat * [ 1 + (Lv/Cp) * RH * d(qsat)/dT ]
        !         = 1/qsat * [ 1 + RH * gam ]                   ! gam in this subr. = Lv/Cp * d(qsat)/dT

        zc3(:ncol,:pver) = ( 1._r8 + rhgbm(:ncol,:pver)*gam(:ncol,:pver) )/qsat(:ncol,:pver)

        ! Now calculate the denominator of the grid-box mean condensation rate.
        ! rdenom = 1/denominator.

        rdenom(:ncol,:pver) = 1._r8/( 1._r8 + ql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)*zc3(:ncol,:pver) )

        ! Calculate term C, then use the same denominator to re-scale all three terms (A, B, and C). 

        term_C(:ncol,:pver) = ql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)*zforcing(:ncol,:pver)

        term_C(:ncol,:pver) = term_C(:ncol,:pver)*rdenom(:ncol,:pver) 
        term_A(:ncol,:pver) = term_A(:ncol,:pver)*rdenom(:ncol,:pver) 
        term_B(:ncol,:pver) = term_B(:ncol,:pver)*rdenom(:ncol,:pver) 

     case(21) !!term C for RK98 scheme with the finite difference method the same as case (1)
        ! use the finite difference to estimate df/dt, which is (fnew-fold)/(tn+1-tn)
        ! term_C = ql_incld*(df/dt)/(1+f*alpha*beta)
         term_C(:ncol,:pver) = ql_incld(:ncol,:pver)*rdtime*( ast(:ncol,:pver) - ast_old(:ncol,:pver) ) &
                              /( 1._r8 + ast(:ncol,:pver)*gam(:ncol,:pver) )

     case(22) !!term C for RK98 scheme with analytical method the same as case (2) 
        ! df/dt = df/dRH * dRH/dt = df/dRH * (dRH/dqv + dRH/dT) 
        ! The zforcing and zc3 are exactly the same as those in case (2) Copy and use them here
        ! Refer to case(2) part for detailed information for these two variables 
        zforcing(:ncol,:pver) = ( qtend(:ncol,:pver)  &
                                 -ttend(:ncol,:pver)*rhgbm(:ncol,:pver)*dqsatdT(:ncol,:pver)) &
                               /qsat(:ncol,:pver)

        zc3(:ncol,:pver) = ( 1._r8 + rhgbm(:ncol,:pver)*gam(:ncol,:pver) )/qsat(:ncol,:pver)

        ! Now calculate the denominator of the grid-box mean condensation rate.
        ! rdenom = 1/denominator. Note that this term is different to that in case(2)
        rdenom(:ncol,:pver) = 1._r8/( 1._r8 + ql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)*zc3(:ncol,:pver) &
                                              /(1._r8 + ast(:ncol,:pver)*gam(:ncol,:pver)) )
        ! Calculate term C, then use the same denominator to re-scale all
        ! three terms (A, B, and C). 

        term_C(:ncol,:pver) = ql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)*zforcing(:ncol,:pver) &
                             /(1._r8 + ast(:ncol,:pver)*gam(:ncol,:pver))

        term_C(:ncol,:pver) = term_C(:ncol,:pver)*rdenom(:ncol,:pver)
        term_A(:ncol,:pver) = term_A(:ncol,:pver)*rdenom(:ncol,:pver)
        term_B(:ncol,:pver) = term_B(:ncol,:pver)*rdenom(:ncol,:pver)

     case default
         write(iulog,*) "Unrecognized value of rkz_term_C_opt:",rkz_term_C_opt,". Abort."
         call endrun
     end select

     !------------------------------------------------------------------
     ! Sum up all three contributors to the grid-box mean condensation.
     !------------------------------------------------------------------
     qme(:ncol,:pver) = term_A(:ncol,:pver) + term_B(:ncol,:pver) + term_C(:ncol,:pver)

     !------------------------------------------------------------------
     ! Send diagnostics to output
     !------------------------------------------------------------------
     call outfld('RKZ_zqme',   qme,    pcols, lchnk)
     call outfld('RKZ_term_A', term_A, pcols, lchnk)
     call outfld('RKZ_term_B', term_B, pcols, lchnk)
     call outfld('RKZ_term_C', term_C, pcols, lchnk)
     call outfld('RKZ_ql_incld', ql_incld, pcols, lchnk)

     select case (rkz_term_C_opt)
     case(0)
         dfdt(:ncol,:pver) = 0._r8 

     case(1)
         dfdt(:ncol,:pver) = rdtime*( ast(:ncol,:pver) - ast_old(:ncol,:pver) ) 
  
     case(2)
         dfdt(:ncol,:pver) = dastdRH(:ncol,:pver) *( zforcing(:ncol,:pver) - zc3(:ncol,:pver)*qme(:ncol,:pver) )

     case(21)
         dfdt(:ncol,:pver) = rdtime*( ast(:ncol,:pver) - ast_old(:ncol,:pver) )

     case(22)
         dfdt(:ncol,:pver) = dastdRH(:ncol,:pver) *( zforcing(:ncol,:pver) - zc3(:ncol,:pver)*qme(:ncol,:pver) )

     case default
         write(iulog,*) "Unrecognized value of rkz_term_C_opt:",rkz_term_C_opt,". Abort."
         call endrun
     end select

     call outfld('RKZ_dfdt', dfdt, pcols, lchnk)


     !------------------------------------------------------------------
     !set up choices of q and ql bottom boundaries to apply limiters
     !------------------------------------------------------------------

     select case (rkz_zsmall_opt)
     case(0) 
        zsmall = 0._r8

     case(1)
        zsmall = 1e-12_r8

     case default
         write(iulog,*) "Unrecognized value of rkz_zsmall_opt:",rkz_zsmall_opt,". Abort."
         call endrun
     end select

     !-----------
     ! limiters
     !-----------
     if (l_rkz_lmt_2) then

        ! Check the sign of qme:
        ! If other forcing leads to strong cooling or moistening,
        ! then condensation might happen, evaporate should not happen

        where ( rhgbm(:ncol,:pver) > 1._r8 ) 
          qme(:ncol,:pver) = max( qme(:ncol,:pver), 0._r8 )

        ! If other forcing dries/warms the atmos a lot, 
        ! evaporation might happen, condensation should not happen. 

        elsewhere ( rhgbm(:ncol,:pver) < rhu00 ) 
          qme(:ncol,:pver) = min( qme(:ncol,:pver), 0._r8 )

        end where
     end if

     !-----------------------------
     ! Check the magnitude of qme:
     if (l_rkz_lmt_3) then

        ! What qme would be needed to bring the whole grid cell to a new saturation equilibrium?
        zlim(:ncol,:pver) = ( state%q(:ncol,:pver,1)-qsat(:ncol,:pver) ) &
                           /( 1._r8 + (latvap/cpair)*dqsatdT(:ncol,:pver) )*rdtime

        ! Condensation should not lead to sub-saturation
        where ( qme(:ncol,:pver) > 0._r8 ) 
          qme(:ncol,:pver) = min( qme(:ncol,:pver), zlim(:ncol,:pver) )

        ! Evaporation should not lead to supersaturation
        elsewhere ( qme(:ncol,:pver) < 0._r8 ) 
          qme(:ncol,:pver) = max( qme(:ncol,:pver), zlim(:ncol,:pver) )
        end where

     end if

     !---------------------------------------------------------------------------
     ! Limit condensation/evaporation rate to avoid negative water concentrations
     !---------------------------------------------------------------------------
     if (l_rkz_lmt_4) then

        ! Avoid negative qv

        qmebf(:ncol,:pver)  = qme(:ncol,:pver)

        zqvneg(:ncol,:pver) = 0._r8
        zstend(:ncol,:pver) = 0._r8
        zqvnew(:ncol,:pver) = state%q(:ncol,:pver,1) - dtime*qme(:ncol,:pver)
        lcondition4_qv(:ncol,:pver) = .false.
        where( zqvnew(:ncol,:pver).lt.zsmall )
           qme(:ncol,:pver) = ( state%q(:ncol,:pver,1) - zsmall )*rdtime
           zqvneg(:ncol,:pver) = zqvnew(:ncol,:pver)
           zstend(:ncol,:pver) = qme(:ncol,:pver)*latvap - qmebf(:ncol,:pver)*latvap
           lcondition4_qv(:ncol,:pver) = .true.
        end where
        qmedf(:ncol,:pver)  = qme(:ncol,:pver) - qmebf(:ncol,:pver)

        call outfld('RKZ_qme_lm4_qv',    qmedf,    pcols, lchnk)
        call outfld('RKZ_qvneg_lm4',    zqvneg,    pcols, lchnk)
        call outfld('RKZ_stend_lm4_qv', zstend,    pcols, lchnk)

        ! Avoid negative ql (note that qme could be negative, which would mean evaporation)

        qmebf(:ncol,:pver)  = qme(:ncol,:pver)

        zqlnew(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq) + dtime*qme(:ncol,:pver)
        zqlneg(:ncol,:pver) = 0._r8
        zstend(:ncol,:pver) = 0._r8
        lcondition4_ql(:ncol,:pver) = .false.
        where( zqlnew(:ncol,:pver).lt.zsmall )
           qme(:ncol,:pver) = ( zsmall - state%q(:ncol,:pver,ixcldliq) )*rdtime
           zqlneg(:ncol,:pver) = zqlnew(:ncol,:pver)
           zstend(:ncol,:pver) = qme(:ncol,:pver)*latvap - qmebf(:ncol,:pver)*latvap
           lcondition4_ql(:ncol,:pver) = .true.
        end where
        qmedf(:ncol,:pver)  = qme(:ncol,:pver) - qmebf(:ncol,:pver)

        call outfld('RKZ_qme_lm4_ql',    qmedf,    pcols, lchnk)
        call outfld('RKZ_qlneg_lm4',    zqlneg,    pcols, lchnk)
        call outfld('RKZ_stend_lm4_ql', zstend,    pcols, lchnk)
       
        ! identify cells that would be "corrected" by limiter 5 

        lcondition5(:ncol,:pver) = (ast(:ncol,:pver) == 0._r8).and.(dfdt(:ncol,:pver) == 0._r8).and.(ltend(:ncol,:pver) /= 0._r8)
        flag_lmt5         (:ncol,:pver) = 0._r8
        flag_lmt4         (:ncol,:pver) = 0._r8
        flag_lmt4_and_lmt5(:ncol,:pver) = 0._r8

        where (lcondition5(:ncol,:pver))
          flag_lmt5(:ncol,:pver) = 1._r8
        end where
        call outfld('RKZ_lmt5_flg', flag_lmt5, pcols, lchnk)

        where ( lcondition4_qv(:ncol,:pver).or.lcondition4_ql(:ncol,:pver)  )
          flag_lmt4(:ncol,:pver) = 1._r8
        end where
        call outfld('RKZ_lmt4_flg', flag_lmt4, pcols, lchnk)

        where ( lcondition5(:ncol,:pver) .and. (lcondition4_qv(:ncol,:pver).or.lcondition4_ql(:ncol,:pver))  )
          flag_lmt4_and_lmt5(:ncol,:pver) = 1._r8
        end where
        call outfld('RKZ_lmt45_flg', flag_lmt4_and_lmt5, pcols, lchnk)

     end if

     !---------------------------------------------------------------------------
     ! Limit condensation/evaporation rate to avoid negative water
     ! concentrations. Also, let qme=0 when f=0 and df/dt=0.
     ! (say, without any changes in cloud fraction and its tendency, the
     ! condensation does not play any roles, we think qme should be zero, thus,
     ! we add this limiter to avoid nonzero qme with f=0 and df/dt=0) . 
     !---------------------------------------------------------------------------
     if (l_rkz_lmt_5) then

        ! Avoid nonzero qme if f=0 and df/dt=0 

        qmebf(:ncol,:pver)  = qme(:ncol,:pver)

        select case (rkz_lmt5_opt)
        case(0)
           where ((ast(:ncol,:pver) == 0._r8).and.(dfdt(:ncol,:pver) == 0._r8))
             qme(:ncol,:pver) = 0._r8
           end where
        case(1)
           where ((ast(:ncol,:pver) == 0._r8).and.(dfdt(:ncol,:pver) == 0._r8).and.(ltend(:ncol,:pver) < 0._r8))
             qme(:ncol,:pver) = 0._r8
           end where
        case(2)
           where ((ast(:ncol,:pver) == 0._r8).and.(dfdt(:ncol,:pver) == 0._r8).and.(ltend(:ncol,:pver) > 0._r8))
             qme(:ncol,:pver) = 0._r8
           end where
        case default
           write(iulog,*) "Unrecognized value of rkz_lmt5_opt:",rkz_lmt5_opt,". Abort."
           call endrun
        end select

        qmedf(:ncol,:pver)  = qme(:ncol,:pver) - qmebf(:ncol,:pver)

        tmp(:ncol,:pver) = 0._r8
        where (qmedf(:ncol,:pver) < 0._r8)
          tmp(:ncol,:pver) = qmedf(:ncol,:pver)
        end where
        call outfld('RKZ_qme_lm5_ng',   tmp,    pcols, lchnk)

        tmp(:ncol,:pver) = 0._r8
        where (qmedf(:ncol,:pver) > 0._r8)
          tmp(:ncol,:pver) = qmedf(:ncol,:pver)
        end where
        call outfld('RKZ_qme_lm5_ps',   tmp,    pcols, lchnk)

     end if

     ! Send limited condensation rate to output
     call outfld('RKZ_qme', qme, pcols, lchnk)

  else

    qme(:ncol,:pver) = 0._r8

  end if !nstep > 1

  !===================================================================
  ! ptend to be returned to the calling subroutine
  !===================================================================
   lq(:) = .FALSE.
   lq(1) = .TRUE.
   lq(ixcldliq) = .TRUE.
   call physics_ptend_init(ptend,state%psetcols, "RKZ macro simplified", ls=.true., lq=lq)

   ptend%q(:ncol,:pver,1)        = -qme(:ncol,:pver)
   ptend%q(:ncol,:pver,ixcldliq) =  qme(:ncol,:pver)
   ptend%s(:ncol,:pver)          =  qme(:ncol,:pver) *latvap
   qmeold(:ncol,:pver)           =  qme(:ncol,:pver) ! save qme for next step 

  end subroutine simple_RKZ_tend

end module simple_condensation_model
