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

    call addfld ('RKZ_qsat',(/'lev'/), 'I','kg/kg','saturation specific humidity used in the simple RKZ scheme')

    call addfld ('RKZ_RH',   (/'lev'/), 'I','1','grid-box mean relative humidity')
    call addfld ('RKZ_f',    (/'lev'/), 'I','1','cloud fraction')
    call addfld ('RKZ_dfdRH',(/'lev'/), 'I','1','derivative of cloud fraction wrt relative humidity')

    call addfld ('RKZ_dfdt', (/'lev'/), 'I','1','estimated cloud fraction tendency')

    call addfld ('RKZ_ql_incld',(/'lev'/), 'I','kg/kg','in-cloud ql used for calculating term C')

    call addfld ('RKZ_term_A',(/'lev'/), 'I','kg/kg/s','grid-box mean condensation rate, term A')
    call addfld ('RKZ_term_B',(/'lev'/), 'I','kg/kg/s','grid-box mean condensation rate, term B')
    call addfld ('RKZ_term_C',(/'lev'/), 'I','kg/kg/s','grid-box mean condensation rate, term C')

    call addfld ('RKZ_Al',(/'lev'/), 'I','kg/kg/s','grid-box mean ql tendency caused by other processes')
    call addfld ('RKZ_Av',(/'lev'/), 'I','kg/kg/s','grid-box mean qv tendency caused by other processes')
    call addfld ('RKZ_AT',(/'lev'/), 'I','kg/kg/s','grid-box mean temperature tendency caused by other processes')

  end subroutine simple_RKZ_init

  !------------------------------------------------------------------------
  ! Calculate condensation rate and the resulting tendencies of the model 
  ! state variables
  !------------------------------------------------------------------------
  subroutine simple_RKZ_tend(state, ptend, tcwat, qcwat, lcwat, ast,  &
                             dtime, ixcldliq, &
                             rkz_cldfrc_opt, &
                             rkz_term_A_opt, &
                             rkz_term_B_opt, &
                             rkz_term_C_opt, &
                             rkz_term_C_ql_opt, & 
                             rkz_term_C_fmin,   & 
                             l_rkz_lmt_2, &
                             l_rkz_lmt_3, &
                             l_rkz_lmt_4  &
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
  real(r8), intent(inout) :: lcwat(:,:)       ! ql          macro- and microphysics calculation in the previous time step

  real(r8), intent(inout) ::   ast(:,:)       ! cloud fraction

  real(r8), intent(in) :: dtime               ! Set model physics timestep
  integer,  intent(in) :: ixcldliq            ! constituent index 

  integer,  intent(in) :: rkz_cldfrc_opt       ! cloud fraction scheme 
  integer,  intent(in) :: rkz_term_A_opt
  integer,  intent(in) :: rkz_term_B_opt
  integer,  intent(in) :: rkz_term_C_opt
  integer,  intent(in) :: rkz_term_C_ql_opt
  real(r8), intent(in) :: rkz_term_C_fmin

  logical,  intent(in) :: l_rkz_lmt_2
  logical,  intent(in) :: l_rkz_lmt_3
  logical,  intent(in) :: l_rkz_lmt_4

  ! tmp work arrays

  real(r8) :: rdtime                 ! 1/dtime

  integer  :: lchnk                  ! index of (grid) chunk handled by this call of the subroutine
  integer  :: ncol                   ! number of active columns in this chunk for which calculations will be done
  integer  :: nstep                  ! time step index
  integer  :: i,k                    ! loop indices for column and vertical layer
  logical  :: lq(pcnst)              ! logical array used when calling subroutine physics_ptend_init indicating 
                                     ! which tracers are affected by this parameterization

  real(r8) ::     ast_old(pcols,pver)! cloud fraction of previous time step

  real(r8) ::        qsat(pcols,pver)! saturation specific humidity
  real(r8) ::         esl(pcols,pver)! saturation vapor pressure (output from subroutine qsat_water, not used)
  real(r8) ::     dqsatdT(pcols,pver)! dqsat/dT
  real(r8) ::         gam(pcols,pver)! L/cpair * dqsat/dT

  real(r8) :: rhu00                  ! threshold grid-box-mean RH used in the diagnostic cldfrc scheme
  real(r8) :: rhgbm(pcols,pver)      ! grid box mean relative humidity
  real(r8) ::     dastdRH(pcols,pver)! df/dRH where f is the cloud fraction and RH the relative humidity

  real(r8) :: qtend(pcols,pver)      ! Moisture tendency caused by other processes
  real(r8) :: ltend(pcols,pver)      ! liquid condensate tendency caused by other processes
  real(r8) :: ttend(pcols,pver)      ! Temperature tendency caused by other processes

  real(r8) :: qme   (pcols,pver)     ! total condensation rate

  real(r8) :: term_A (pcols,pver)
  real(r8) :: term_B (pcols,pver)
  real(r8) :: term_C (pcols,pver)
  real(r8) :: rdenom (pcols,pver)    ! 1/denominator in the final expression of total grid-box mean condensation rate 
                                     ! when option 2 is used for term C. Note that in this case, 1/denominator will be 
                                     ! multipled to all three terms (A, B, and C).

  real(r8) ::  ql_incld(pcols,pver)  ! in-cloud liquid concentration
  real(r8) ::  dfdt    (pcols,pver)  ! df/dt where f is the cloud fraction
  real(r8) :: zforcing (pcols,pver)
  real(r8) :: zc3      (pcols,pver)

  real(r8) :: zqvnew(pcols,pver)     ! qv at new time step if total condenation rate is not limited. Might be negative.
  real(r8) :: zqlnew(pcols,pver)     ! ql at new time step if total condenation rate is not limited. Might be negative.
  real(r8) :: zlim (pcols,pver)
  real(r8),parameter :: zsmall = 1e-12_r8

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  rdtime = 1._r8/dtime
  ncol  = state%ncol
  nstep = get_nstep()
  lchnk = state%lchnk  ! needed by "call outfld" for model output

  ! Calculate saturation specific humidity (qsat) and its derivative wrt temperature (dqsatdT)

  do k=1,pver
  do i=1,ncol
     call qsat_water( state%t(i,k), state%pmid(i,k), esl(i,k), qsat(i,k), gam(i,k), dqsatdT(i,k) )
  end do
  end do
  call outfld('RKZ_qsat', qsat, pcols, lchnk)

  ! Calculate the grid-box-mean relative humidity (rhgbm)

  rhgbm(:ncol,:pver) = state%q(:ncol,:pver,1)/qsat(:ncol,:pver)
  call outfld('RKZ_RH', rhgbm, pcols, lchnk)

  ! Cloud fraction (ast): save old values, diagnose new values

  if (nstep > 1) ast_old(:ncol,:pver) = ast(:ncol,:pver)

  call  smpl_frc( state%q(:,:,1), state%q(:,:,ixcldliq), qsat,      &! all in
                  ast, rhu00, dastdRH,                              &! inout, out, out
                  rkz_cldfrc_opt, 0.5_r8, 0.5_r8, pcols, pver, ncol )! all in

  call outfld('RKZ_qv',    state%q(:,:,1),        pcols, lchnk)
  call outfld('RKZ_ql',    state%q(:,:,ixcldliq), pcols, lchnk)

  call outfld('RKZ_f',     ast,      pcols, lchnk)
  call outfld('RKZ_dfdRH', dastdRH,  pcols, lchnk)

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
        term_A = 0._r8

     case(1)
        term_A(:ncol,:pver) = ast(:ncol,:pver)     &
                           *( qtend(:ncol,:pver) - dqsatdT(:ncol,:pver)*ttend(:ncol,:pver) ) &
                           /( 1._r8 + gam(:ncol,:pver) )
     case default
         write(iulog,*) "Unrecognized value of rkz_term_A_opt:",rkz_term_A_opt,". Abort."
         call endrun
     end select

     !--------------------------------------------------------------------------------------
     ! Term B: condensation in cloud-free portion of a grid box in response to cloud liquid tencendy
     ! caused by other processes.

     select case (rkz_term_B_opt) 
     case(0) ! omit this term
        term_B = 0._r8

     case(1) 
     ! Following Zhang et al. (2003), assume the in-cloud A_l equals the grid-box mean A_l.
     ! Hence, - (\overline{A_l} - f \hat{A_l}) = - (1-f)\overline{A_l}.

        term_B(:ncol,:pver) = - (1._r8 - ast(:ncol,:pver))*ltend(:ncol,:pver)

     case(2) 
     ! For testing only: term B = - 0.5*\overline{A_l}

        term_B(:ncol,:pver) = - 0.5_r8*ltend(:ncol,:pver)

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

     CASE (7,17)
       ql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)/max(ast(:ncol,:pver),rkz_term_C_fmin)

     CASE DEFAULT
       write(iulog,*) "Unrecognized value of rkz_term_C_ql_opt:",rkz_term_C_ql_opt,". Abort."
       call endrun
     END SELECT

     ! term C = ql_incld * df/dt

     select case (rkz_term_C_opt)
     case(0)  ! omit this term
        term_C = 0._r8

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

        rdenom(:ncol,:pver) = 1._r8/( 1.+ ql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)*zc3(:ncol,:pver) )

        ! Calculate term C, then use the same denominator to re-scale all three terms (A, B, and C). 

        term_C(:ncol,:pver) = ql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)*zforcing(:ncol,:pver)

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
     end select
     call outfld('RKZ_dfdt', dfdt, pcols, lchnk)

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

        zqvnew(:ncol,:pver) = state%q(:ncol,:pver,1) - dtime*qme(:ncol,:pver)
        where( zqvnew(:ncol,:pver).lt.zsmall )
           qme(:ncol,:pver) = ( state%q(:ncol,:pver,1) - zsmall )*rdtime
        end where

        ! Avoid negative ql (note that qme could be negative, which would mean evaporation)

        zqlnew(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq) + dtime*qme(:ncol,:pver)
        where( zqlnew(:ncol,:pver).lt.zsmall )
           qme(:ncol,:pver) = ( zsmall - state%q(:ncol,:pver,ixcldliq) )*rdtime
        end where

     end if

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

  end subroutine simple_RKZ_tend

end module simple_condensation_model
