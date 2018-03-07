module reed_jablonowski_condensation_model

  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use ppgrid,       only: pver, pcols

  implicit none
  private
  public :: reed_jablonowski_sat_adj_tend
  public :: reed_jablonowski_init

contains

  !------------------------------------------------------------------------
  ! Initialization of the reed_jablonowski simple condensation model.
  ! Currently this contains just registration of output variables.
  !------------------------------------------------------------------------

  subroutine reed_jablonowski_init()

    use cam_history,    only: addfld

    implicit none
!!!!add history fileds!!!!!!!!!!!!!!!!!!!!!
    call addfld ('RKZ_dqdt', (/'lev'/), 'I', 'kg/kg/s', 'condensation rate in the Reed-Jablonowski scheme')
    call addfld ('RKZ_dsdt', (/'lev'/), 'I', 'J/kg/s', 'dry static energy tendency in the Reed-Jablonowski scheme')

  end subroutine reed_jablonowski_init


subroutine reed_jablonowski_sat_adj_tend( state, ptend, dtime )
!----------------------------------------------------------------------- 
! 
! Purpose: Simple parameterization of large-scale condensation
!
! Based on the simple physics suite of K. A. Reed and C. Jablonowski
! Technical modifications by H. Wan (PNNL, 2014-06).
!
! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone 
!            simulations of intermediate complexity: A test case for AGCMs, 
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
!-----------------------------------------------------------------------
  ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
  ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM


   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols
   use constituents, only: pcnst
   use physics_types,only: physics_state, physics_ptend, physics_ptend_init
   use cam_history,    only: addfld
   use cam_history,   only: outfld

   implicit none
!
! arguments
!
   real(r8), intent(in) :: dtime        ! Set model physics timestep

   type(physics_state), intent(in),  target   :: state   ! State variables
   type(physics_ptend), intent(out), target   :: ptend   ! Package tendencies

   real(r8) :: precl(pcols)         ! Precipitation rate (m_water / s)

!
!---------------------------Local workspace-----------------------------
!

! Integers for loops

   integer  i,k                         ! Longitude, level indices
   integer  ncol

   logical  :: lq(pcnst)

   real(r8),pointer :: t(:,:)
   real(r8),pointer :: q(:,:)
   real(r8),pointer :: pmid(:,:)
   real(r8),pointer :: pdel(:,:)

   
! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air 
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)
   real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
   real(r8) pi                          ! pi

! Simple Physics Specific Constants 


   real(r8) SST_tc                      ! Sea Surface Temperature for tropical cyclone test
   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation
   real(r8) T00                         ! Horizontal mean T at surface for moist baro test
   real(r8) u0                          ! Zonal wind constant for moist baro test
   real(r8) latw                        ! halfwidth for  for baro test
   real(r8) eta0                        ! Center of jets (hybrid) for baro test
   real(r8) etav                        ! Auxiliary variable for baro test
   real(r8) q0                          ! Maximum specific humidity for baro test

! Physics Tendency Arrays
  real(r8) dsdt(pcols,pver)             ! Dry static energy tendency 
  real(r8) dqdt(pcols,pver)             ! Specific humidity tendency


! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

  integer lchnk                  ! index of (grid) chunk handled by this call of the subroutine


!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
   pi     = 4._r8*atan(1._r8)            ! pi

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
      C        = 0.0011_r8      ! From Smith and Vogl 2008
      SST_tc   = 302.15_r8      ! Constant Value for SST for tropical cyclone test
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Smith and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Smith and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Smith and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      T00      = 288.0_r8         ! Horizontal mean T at surface for moist baro test
      u0       = 35.0_r8          ! Zonal wind constant for moist baro test
      latw     = 2.0_r8*pi/9.0_r8 ! Halfwidth for  for baro test
      eta0     = 0.252_r8         ! Center of jets (hybrid) for baro test
      etav     = (1._r8-eta0)*0.5_r8*pi ! Auxiliary variable for baro test
      q0       = 0.021            ! Maximum specific humidity for baro test


!===============================================================================
!===============================================================================
    ncol  = state%ncol
    lchnk = state%lchnk  ! needed by "call outfld" for model output

! input

       t => state%t
       q => state%q(:,:,1)
    pmid => state%pmid
    pdel => state%pdel

! local/output

    lq(:) = .FALSE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, "Reed-Jablonowski saturation adjustment", ls=.true., lq=lq)

    precl(:ncol) = 0._r8                  ! initialize precipitation rate with zero
    dsdt(:ncol,:pver)  = 0._r8            ! initialize temperature tendency with zero
    dqdt(:ncol,:pver)  = 0._r8            ! initialize specific humidity tendency with zero
 
!===============================================================================
!
! Large-Scale Condensation and Precipitation Rate
!
!===============================================================================
!
! Calculate Tendencies
!
      do k=1,pver
         do i=1,ncol
            qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/t(i,k))-1._r8/T0))  ! saturation specific humidity
            if (q(i,k) > qsat) then                                                 ! saturated?
               tmp  = 1._r8/dtime*(q(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(i,k)**2)))
               dsdt(i,k) = dsdt(i,k)+latvap*tmp
               dqdt(i,k) = dqdt(i,k)-tmp
               precl(i)  = precl(i) + tmp*pdel(i,k)/(gravit*rhow)                    ! precipitation rate, computed via a vertical integral
                                                                                    ! corrected in version 1.3
            end if
         end do
      end do

      ptend%q(:ncol,:pver,1) = dqdt(:ncol,:pver)
      ptend%s(:ncol,:pver)   = dsdt(:ncol,:pver)

   call outfld('RKZ_dqdt', dqdt, pcols, lchnk)
   call outfld('RKZ_dsdt', dsdt, pcols, lchnk)
   
!===============================================================================
! Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the history file
!       after the moist physics process.  This ensures that the moisture fields
!       are somewhat in equilibrium.
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM

!===============================================================================
!
! Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC 
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. Add your own routine here.
!===============================================================================
  !  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

   return
end subroutine reed_jablonowski_sat_adj_tend 

end module reed_jablonowski_condensation_model
