
module MOSART_RES_type

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: MOSART_RES_type
!
! !DESCRIPTION:
! Module containing data structure for reservoir dynamics
!
!! !REVISION HISTORY:
! Hongyi Li: Created 03/2015
!USES:
  use RunoffMod     , only : Tctl, TUnit, rtmCTL
  use RtmSpmd       , only : masterproc
  use RtmVar        , only : iulog
  use RtmIO         
  use rof_cpl_indices, only : nt_rtm
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_flush, shr_sys_abort
  use netcdf
  use pio
  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit 

! !PUBLIC TYPES:
  implicit none
  private

! control information for reservoirs
  type Treservoir_control
     integer :: NDam             ! number of dams
     integer :: localNumDam      ! number of dams on decomposition

     integer :: month            ! month of the simulation
     integer :: year             ! year of the simulation
     integer :: RESFlag          ! Flag for invoking reservoir processes or not

     character(len=350) :: paraFile         ! the path of the parameter files
     character(len=350) :: paraPath         ! the path of the parameter files

  end type Treservoir_control

! reservoirs information: geometry etc.
    type Tpara_reservoir

        ! reservoirs or lakes
        real(r8), pointer :: height(:)       ! reservoir height, [m]
        real(r8), pointer :: length(:)       ! reservoir length, [m]
        real(r8), pointer :: area(:)         ! reservoir surface area, [m2]
        real(r8), pointer :: storage(:,:)    ! reservoir storage, [m3] for water, [kg] for sediment
        real(r8), pointer :: Tres(:)         ! Change of water residence time due to a reservoir on main channel, [yrs]
                                             ! here the reservoirs both regulate flow and trap sediment
        real(r8), pointer :: Tres_t(:)       ! Change of water residence time due to a reservoir on sub-network channel, [yrs]
        real(r8), pointer :: Tres_r(:)       ! Change of water residence time due to a reservoir on main channel, [yrs]
                                             ! here the reservoirs only trap sediment 
        real(r8), pointer :: Eff_trapping(:) ! reservoir trapping efficiency on main channel, [-]
        real(r8), pointer :: Eff_trapping_t(:) ! reservoir trapping efficiency on sub-network channel, [-]
        real(r8), pointer :: Eff_trapping_r(:) ! reservoir trapping efficiency on main channel, [-]
    end type Tpara_reservoir

! reservoirs status and fluxes
    type TstatusFlux_reservoir
     ! reservoirs, lakes on the main channel
     !! states
     real(r8), pointer :: wres(:,:)    ! MOSART reservoir storage (m3 for water, kg for mud and sand sediment)
     real(r8), pointer :: dwres(:,:)   ! MOSART reservoir storage change (m3 for water, kg for mud and sand sediment)
     !! exchange fluxes
     real(r8), pointer :: eres_in(:,:) ! MOSART reservoir inflow  (m3/s for water, kg/s for mud and sand sediment) 
     real(r8), pointer :: eres_out(:,:)! MOSART reservoir outflow  (m3/s for water, kg/s for mud and sand sediment) 

     ! reservoirs, lakes on the subnetwork channel
     !! states
     real(r8), pointer :: wres_t(:,:)    ! MOSART reservoir storage (m3 for water, kg for mud and sand sediment)
     real(r8), pointer :: dwres_t(:,:)   ! MOSART reservoir storage change (m3 for water, kg for mud and sand sediment)
     !! exchange fluxes
     real(r8), pointer :: eres_in_t(:,:) ! MOSART reservoir inflow  (m3/s for water, kg/s for mud and sand sediment) 
     real(r8), pointer :: eres_out_t(:,:)! MOSART reservoir outflow  (m3/s for water, kg/s for mud and sand sediment) 

    end type TstatusFlux_reservoir

    type (Treservoir_control)   ,   public :: Tres_ctl
    type (Tpara_reservoir),   public :: Tres_para
    type (TstatusFlux_reservoir),   public :: Tres

  public :: MOSART_reservoir_sed_init

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------

  subroutine MOSART_reservoir_sed_init
! !DESCRIPTION:
! initialize MOSART-reservoir variables
! 
! !USES:
! !ARGUMENTS:
!
! !REVISION HISTORY:
! Author: Hongyi Li
!

     ! !DESCRIPTION: initilization of reservoir_sed module
     implicit none

    integer :: ier                  ! error code
    integer :: begr, endr, iunit, nn, n, cnt, nr, nt
    integer  :: damID
    character(len=*),parameter :: subname = '(MOSART_reservoir_sed_init)'
    character(len=*),parameter :: FORMI = '(2A,2i10)'
    character(len=*),parameter :: FORMR = '(2A,2g15.7)'
    type(file_desc_t):: ncid       ! netcdf file
    type(var_desc_t) :: vardesc    ! netCDF variable description
    type(io_desc_t)    :: iodesc_dbl ! pio io desc

    begr = rtmCTL%begr
    endr = rtmCTL%endr  

    if(endr >= begr) then
        allocate(Tres_para%Tres(begr:endr))
        Tres_para%Tres = 0._r8
        allocate(Tres_para%Eff_trapping(begr:endr))
        Tres_para%Eff_trapping = 0._r8
        do iunit=begr,endr
            damID = WRMUnit%INVicell(iunit) 
            if (.not.(damID > ctlSubwWRM%LocalNumDam .OR. damID <= 0 .or. WRMUnit%MeanMthFlow(damID,13) <= 0.01_r8)) then
                Tres_para%Tres(iunit) = CRTres(WRMUnit%StorCap(damID), WRMUnit%MeanMthFlow(damID,13))
                Tres_para%Eff_trapping(iunit) = CREff_trapping(Tres_para%Tres(iunit))
                !write(iulog,*) ' Reservoir Trapping ', iunit, WRMUnit%StorCap(damID), WRMUnit%MeanMthFlow(damID,13), Tres_para%Tres(iunit), Tres_para%Eff_trapping(iunit)
            else
                Tres_para%Tres(iunit) = 0._r8
                Tres_para%Eff_trapping(iunit) = 0._r8
            end if
        end do    

        allocate (Tres%wres(begr:endr,nt_rtm))
        Tres%wres = 0._r8

        allocate (Tres%dwres(begr:endr,nt_rtm))
        Tres%dwres = 0._r8

        allocate (Tres%eres_in(begr:endr,nt_rtm))
        Tres%eres_in = 0._r8

        allocate (Tres%eres_out(begr:endr,nt_rtm))
        Tres%eres_out = 0._r8    

        if(0>1) then ! comment out temporily, please do not delete. We may need to use this block in future refinement
            Tres_para%Eff_trapping_t = 0._r8
            do iunit=begr,endr
			    if(Tres_para%Tres_t(iunit)>0._r8) then
                    Tres_para%Eff_trapping_t(iunit) = CREff_trapping(Tres_para%Tres_t(iunit))
				else
                    Tres_para%Eff_trapping_t(iunit) = 0._r8
				end if
			end do
        
            allocate(Tres_para%Eff_trapping_r(begr:endr))
            Tres_para%Eff_trapping_r = 0._r8
            do iunit=begr,endr
			    if(Tres_para%Tres_r(iunit)>0._r8) then
                   Tres_para%Eff_trapping_r(iunit) = CREff_trapping(Tres_para%Tres_r(iunit))
				else
                    Tres_para%Eff_trapping_r(iunit) = 0._r8
				end if
			end do
        end if

        allocate (Tres%wres_t(begr:endr,nt_rtm))
        Tres%wres_t = 0._r8

        allocate (Tres%dwres_t(begr:endr,nt_rtm))
        Tres%dwres_t = 0._r8

        allocate (Tres%eres_in_t(begr:endr,nt_rtm))
        Tres%eres_in_t = 0._r8

        allocate (Tres%eres_out_t(begr:endr,nt_rtm))
        Tres%eres_out_t = 0._r8    

    end if  

  end subroutine MOSART_reservoir_sed_init  

    function CRTres(V_, Q_) result(Tres_)
      ! !DESCRIPTION: calculate large reservoir trapping efficiency based on Eqn (1) and Fig. 2 in Vorosmarty et al. (2003)
      !! Vorosmarty et al, Anthropogenic sediment retention: Major global impact from registered river impoundments, Glob. Planet. Change, 39, 169-190
      implicit none
      real(r8), intent(in) :: V_       ! reservoir maximum reported storage capacity, [m3]; 
      real(r8), intent(in) :: Q_       ! Mean annual inflow, [m3/s]; 
      real(r8)             :: Tres_        ! approximated residence time of regulated portion of basin [yrs]

      real(r8) :: V_operational  ! proportion of maximum storage at which reservoirs are assumed to operate routinely [m3]
      real(r8) :: Eeff  ! Utilization factor [-]

      Eeff = 0.67_r8
      V_operational = Eeff * V_
      Tres_ = V_operational / Q_ 
      Tres_ = Tres_ / (365*24*3600)  !seconds --> yrs

      return
    end function CRTres

    function CREff_trapping(Tres_) result(Eff_trapping_)
      ! !DESCRIPTION: calculate large reservoir trapping efficiency based on Eqn (1) and Fig. 2 in Vorosmarty et al. (2003)
      !! Vorosmarty et al, Anthropogenic sediment retention: Major global impact from registered river impoundments, Glob. Planet. Change, 39, 169-190
      implicit none
      real(r8), intent(in) :: Tres_        ! approximated residence time of regulated portion of basin [yrs]
      real(r8)             :: Eff_trapping_        ! trapping efficiency [-]


      Eff_trapping_ = 1 - 0.05/sqrt(Tres_)

      if(Eff_trapping_ < 0._r8) then
          Eff_trapping_ = 0._r8
      end if

      return
    end function CREff_trapping


!
!
!EOP
!-----------------------------------------------------------------------


end module MOSART_RES_type