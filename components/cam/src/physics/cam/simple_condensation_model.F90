!-----------------------------------------------------------------------------------------
! This version of the Rasch-Kristjansson macro-physics scheme
!  - includes the in-cloud condensation term;
!  - omits the clear-sky liquid advection term;
!  - use the semi-analytical method of Zhang et al. (2003) for the cloud expansion term.
!
! Multiple versions of the third term was included for testing.
!
! History: 
!  - first version by Hui Wan (PNNL, 2014-2015)
!  - revised for the SciDAC Convergence project, Hui Wan (PNNL, 2017)
!-----------------------------------------------------------------------------------------
module simple_condensation_model

  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none
  private
  public :: simple_RKZ_tend

contains

  subroutine simple_RKZ_tend(state, ptend, tcwat, qcwat, lcwat, ast,  &
                             dtime, ixcldliq, &
                             rkz_cldfrc_opt, &
                             rkz_term_A_opt, &
                             rkz_term_B_opt, &
                             rkz_term_C_opt, &
                             rkz_term_C_ql_opt, & 
                             l_rkz_lmt_2, &
                             l_rkz_lmt_3, &
                             l_rkz_lmt_4  &
                             )

  use shr_kind_mod, only: r8=>shr_kind_r8
  use ppgrid,       only: pver, pcols
  use constituents, only: pcnst
  use physics_types,only: physics_state, physics_ptend, physics_ptend_init
  use time_manager,  only: get_nstep
  use wv_saturation, only: qsat_water
  use physconst,     only: latvap, cpair
  use simple_cloud_fraction, only: smpl_frc
 !use cam_history,   only: outfld

  implicit none
  !
  ! arguments
  !
  type(physics_state), intent(in), target    :: state       ! State variables
  type(physics_ptend), intent(out)   :: ptend       ! Package tendencies

  real(r8), intent(inout) :: tcwat(:,:)
  real(r8), intent(inout) :: qcwat(:,:)
  real(r8), intent(inout) :: lcwat(:,:)
  real(r8), intent(inout) ::   ast(:,:)

  real(r8), intent(in) :: dtime               ! Set model physics timestep
  integer,  intent(in) :: ixcldliq            ! constituent index 

  integer,  intent(in) :: rkz_cldfrc_opt       ! cloud fraction scheme 
  integer,  intent(in) :: rkz_term_A_opt
  integer,  intent(in) :: rkz_term_B_opt
  integer,  intent(in) :: rkz_term_C_opt
  integer,  intent(in) :: rkz_term_C_ql_opt

  logical,  intent(in) :: l_rkz_lmt_2
  logical,  intent(in) :: l_rkz_lmt_3
  logical,  intent(in) :: l_rkz_lmt_4

  real(r8) :: rdtime
  integer  :: ncol, itim, ifld, nstep, lchnk
  integer  :: i,k
  logical  :: lq(pcnst)


  real(r8) ::     ast_old(pcols,pver)
  real(r8) ::         esl(pcols,pver)        ! not used
  real(r8) ::        qsat(pcols,pver)
  real(r8) ::       dqsdt(pcols,pver)        ! dqsat/dT
  real(r8) ::         gam(pcols,pver)        ! L/cpair * dqsat/dT
  real(r8) ::     dastdRH(pcols,pver)

  real(r8) :: qtend(pcols,pver)        ! Moisture tendencies
  real(r8) :: ltend(pcols,pver)        ! liquid condensate tendencies
  real(r8) :: ttend(pcols,pver)        ! Temperature tendencies
  real(r8) :: zforcing(pcols,pver)
  real(r8) :: zc3(pcols,pver)
  real(r8) :: zqlf(pcols,pver)

  real(r8) :: qme   (pcols,pver)
  real(r8) :: zqvnew(pcols,pver)
  real(r8) :: zqlnew(pcols,pver)
  real(r8) :: ztmp  (pcols,pver)
  real(r8) :: zql_incld(pcols,pver)
  real(r8),parameter :: zsmall = 1e-12_r8

  real(r8) :: rhgbm(pcols,pver)
  real(r8) :: zlim (pcols,pver)

  real(r8) :: rhu00  ! threshold grid-box-mean RH used in the diagnostic cldfrc scheme

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  rdtime = 1._r8/dtime
  ncol  = state%ncol
  nstep = get_nstep()

 !lchnk = state%lchnk  ! needed for "call outfld"

  ! Calculate qsat and its derivative wrt temperature

  do k=1,pver
  do i=1,ncol
     call qsat_water( state%t(i,k), state%pmid(i,k), esl(i,k), qsat(i,k), gam(i,k), dqsdt(i,k) )
  end do
  end do

  ! Calculate the grid-box-mean relative humidity

  rhgbm(:ncol,:pver) = state%q(:ncol,:pver,1)/qsat(:ncol,:pver)

  ! Cloud fraction: save old values, diagnose new values

  if (nstep > 1) ast_old(:ncol,:pver) = ast(:ncol,:pver)

  call  smpl_frc( state%q(:,:,1), state%q(:,:,ixcldliq), qsat,     &! in
                  ast, rhu00, dastdRH,                             &! inout, out
                  rkz_cldfrc_opt, 0.5_r8, 0.5_r8, pcols, pver, ncol )! in

  !===================================================================
  ! Condensation/evaporation rate
  !===================================================================
  ! qme is the array to store the grid-box mean condensation/evaporation rate.
  ! Here it is initialized with zeros.
  ! Later we will accumulate the contribution from various terms.

  qme(:ncol,:pver) = 0._r8

  if (nstep > 1) then
     !------------------------------------------------------------------
     ! Diagose "other forcing" on qv and temperature

     qtend(:ncol,:pver) = ( state%q(:ncol,:pver,1) - qcwat(:ncol,:pver) )*rdtime
     ttend(:ncol,:pver) = ( state%t(:ncol,:pver)   - tcwat(:ncol,:pver) )*rdtime

     !------------------------------------------------------------------
     ! Term A: qme in cloudy portion of a grid box, weighted by cloud fraction

     select case (rkz_term_A_opt) 
     case(0)
        continue  ! omit this term
     case(1)
        qme(:ncol,:pver) = qme(:ncol,:pver) + ast(:ncol,:pver)     &
                           *( qtend(:ncol,:pver) - dqsdt(:ncol,:pver)*ttend(:ncol,:pver) ) &
                           /( 1._r8 + gam(:ncol,:pver) )
     case default
         write(iulog,*) "Unrecognized value of rkz_term_A_opt:",rkz_term_A_opt,". Abort."
         call endrun
     end select

     !--------------------------------------------------------------------------------------
     ! Term B: qme in cloud-free portion of a grid box in response to cloud liquid tencendy
     ! caused by other processes.

     select case (rkz_term_B_opt) 
     case(0)
        continue  ! omit this term

     case(1)
     ! Following Zhang et al. (2003), assume the in-cloud A_l equals the grid-box mean A_l.
     ! Hence - (\overline{A_l} - f \hat{A_l}) = - (1-f)\overline{A_l}.
     ! First derive A_l:

        ltend(:ncol,:pver) = ( state%q(:ncol,:pver,ixcldliq) - lcwat(:ncol,:pver) )*rdtime

       ztmp(:ncol,:pver) = - (1._r8 - ast(:ncol,:pver))*ltend(:ncol,:pver)
        qme(:ncol,:pver) = qme(:ncol,:pver) + ztmp(:ncol,:pver)

     case default
         write(iulog,*) "Unrecognized value of rkz_term_B_opt:",rkz_term_B_opt,". Abort."
         call endrun
     end select

     !------------------------------------------------------------------------------------
     ! Term C: qme in cloud-free portion of a grid box, related to cloud expansion/erosion 

     ztmp(:ncol,:pver) = 0._r8

     ! Do the calculation only if this term is turned on and cloud fraction is time-dependent
     select case (rkz_term_C_opt) 
     case(0)
        continue  ! omit this term

     case(2)
     if (rkz_cldfrc_opt.ne.0) then  

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
        
        zforcing(:ncol,:pver) = 0._r8 
        where ( ast(:ncol,:pver) > 0._r8 )
              zforcing(:ncol,:pver) = ( qtend(:ncol,:pver)  &
                                       -ttend(:ncol,:pver)*rhgbm(:ncol,:pver)*dqsdt(:ncol,:pver) ) &
                                      /qsat(:ncol,:pver)
        end where

        ! zc3 is the term C_gamma in Zhang et al. 2003). It appears as part of the denominator
        ! of the final expression for total grid-box-mean condensation rate (i.e., qme in this subr.)
        !
        ! C_gamma = 1/qsat + Lv/Cp * (qv/qsat^2) * d(qsat)/dT
        !         = 1/qsat * [ 1 + (Lv/Cp) * RH * d(qsat)/dT ]
        !         = 1/qsat * [ 1 + RH * gam ]                   !(with gam in this subr. = Lv/Cp * d(qsat)/dT

        zc3(:ncol,:pver) = ( 1._r8 + rhgbm(:ncol,:pver)*gam(:ncol,:pver) )/qsat(:ncol,:pver)

        ! The next block of code tests different expressions for the (ql/f * df/dRH) term. 
        ! The value of that term is stored in the array zqlf.
        ! 
        ! 1,11 (ql_incloud=constant, df/dRH=constant)       
        ! 2,12 (ql_incloud=constant, df/dRH=calculated)       
        ! 3,14 (ql_incloud=qv,       df/dRH=calculated)      
        ! 4,14 (ql_incloud=ql,       df/dRH=calculated)      
        ! 5,15 (ql_incloud=ql/f,     df/dRH=calculated)             
        ! 6,16 (ql_incloud=0 if f<0.01%)   
        ! 7,17 (ql_incloud=ql/max(f,0.01%) 

        SELECT CASE (rkz_term_C_ql_opt)
        CASE (1,11)
          zqlf(:ncol,:pver) = 1e-5_r8

        CASE (2,12)
          zqlf(:ncol,:pver) = 1e-5_r8*dastdRH(:ncol,:pver)

        CASE (3,13)
          zqlf(:ncol,:pver) = state%q(:ncol,:pver,1)*dastdRH(:ncol,:pver)

        CASE (4,14)
          zqlf(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)*dastdRH(:ncol,:pver)

        CASE (5,15)

          zql_incld(:,:) = 0._r8
          where ( state%q(:ncol,:pver,ixcldliq) > zsmall .and. ast(:ncol,:pver) > 1e-8_r8)
           zql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)/ast(:ncol,:pver)
          end where
          zqlf(:ncol,:pver) = zql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)

        CASE (6,16)

          zql_incld(:,:) = 0._r8
          where ( state%q(:ncol,:pver,ixcldliq) > zsmall .and. ast(:ncol,:pver) > 1e-4_r8)
           zql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)/ast(:ncol,:pver)
          end where
          zqlf(:ncol,:pver) = zql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)

        CASE (7,17)

          zql_incld(:,:) = 0._r8
          where ( state%q(:ncol,:pver,ixcldliq) > zsmall .and. ast(:ncol,:pver) > 1e-8_r8)
          !zql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)/max(ast(:ncol,:pver), 1e-4)
           zql_incld(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)/max(ast(:ncol,:pver), 1e-3)
          end where
          zqlf(:ncol,:pver) = zql_incld(:ncol,:pver)*dastdRH(:ncol,:pver)

        END SELECT

        ! Now calculate the denominator of the grid-box mean condensation rate

        if (rkz_term_C_ql_opt < 10 ) then
           ztmp(:ncol,:pver) = 1._r8
        else
           ztmp(:ncol,:pver) = 1._r8/( 1.+ zqlf(:ncol,:pver)*zc3(:ncol,:pver) )
        end if

        ! Add the contribution of cloud growth to the grid-box mean condensation rate

        qme(:ncol,:pver) = ( qme(:ncol,:pver) + zqlf(:ncol,:pver)*zforcing(:ncol,:pver) ) &
                          *ztmp(:ncol,:pver) 

     end if

     case default
         write(iulog,*) "Unrecognized value of rkz_term_C_opt:",rkz_term_C_opt,". Abort."
         call endrun
     end select

     !-----------
     ! limiting
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
                           /( 1._r8 + (latvap/cpair)*dqsdt(:ncol,:pver) )*rdtime

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

        ! Avoid negative ql

        zqlnew(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq) + dtime*qme(:ncol,:pver)
        where( zqlnew(:ncol,:pver).lt.zsmall )
           qme(:ncol,:pver) = ( zsmall - state%q(:ncol,:pver,ixcldliq) )*rdtime
        end where

     end if

  end if !nstep > 1

  !===================================================================
  ! ptend
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
