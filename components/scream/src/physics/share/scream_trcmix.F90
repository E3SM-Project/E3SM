
module scream_trcmix

  !------------------------------------------------------------------------------------------------
  ! Purpose:
  ! Provide default distributions of CH4, N2O, CFC11 and CFC12 to the radiation routines.
  ! **NOTE** CO2 is assumed by the radiation to a be constant value.  This value is
  !          currently supplied directly by the chem_surfvals module.
  !
  ! Much of this code was copied from components/eam/src/physics/cam/ghg_data.F90
  !------------------------------------------------------------------------------------------------

  ! use physconst,      only: mwdry, mwch4, mwn2o, mwf11, mwf12, mwco2
  ! use chem_surfvals,  only: chem_surfvals_get, chem_surfvals_co2_rad

  use physics_utils, only: rtype

  implicit none
  save

  !================================================================================================
contains
  !================================================================================================

  function chem_surfvals_get(name, rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr)
    ! Copied from components/eam/src/physics/chem_surfvals.F90

    character(len=*), intent(in) :: name
    real(rtype),      intent(in) :: rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr

    real(rtype) :: chem_surfvals_get

    select case (name)
    case ('CO2VMR')
       chem_surfvals_get = co2vmr
    case ('CO2MMR')
       chem_surfvals_get = rmwco2 * co2vmr
    case ('N2OVMR')
       chem_surfvals_get = n2ovmr
    case ('CH4VMR')
       chem_surfvals_get = ch4vmr
    case ('F11VMR')
       chem_surfvals_get = f11vmr
    case ('F12VMR')
       chem_surfvals_get = f12vmr
    case ('O2MMR')
       chem_surfvals_get = o2mmr
    case default

    end select

  end function chem_surfvals_get

  function chem_surfvals_co2_rad(rmwco2, co2vmr, co2vmr_rad)
    ! Copied from components/eam/src/physics/chem_surfvals.F90

    ! Return the value of CO2 (as mmr) that is radiatively active.

    ! This method is used by ghg_data to set the prescribed value of CO2 in
    ! the physics buffer.  If the user has set the co2vmr_rad namelist
    ! variable then that value will override either the value set by the
    ! co2vmr namelist variable, or the values time interpolated from a
    ! dataset.

    ! This method is also used by cam_history to write the radiatively active
    ! CO2 to the history file.  The optional argument allows returning the
    ! value as vmr.

    real(rtype), intent(in) :: rmwco2, co2vmr, co2vmr_rad

    ! Return value
    real(rtype) :: chem_surfvals_co2_rad

    !-----------------------------------------------------------------------

    if (co2vmr_rad > 0._rtype) then
       chem_surfvals_co2_rad = rmwco2 * co2vmr_rad
    else
       chem_surfvals_co2_rad = rmwco2 * co2vmr
    end if

  end function chem_surfvals_co2_rad

  subroutine trcmix( &
       name, ncol, pcols, pver, clat, pmid, q, &
       mwdry, mwco2, mwn2o, mwch4, mwf11, mwf12, &
       o2mmr, co2vmr_rad, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr)

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
    ! CFC12
    !
    ! Method:
    ! Distributions assume constant mixing ratio in the troposphere
    ! and a decrease of mixing ratio in the stratosphere. Tropopause
    ! defined by ptrop. The scale height of the particular trace gas
    ! depends on latitude. This assumption produces a more realistic
    ! stratospheric distribution of the various trace gases.
    !
    ! Author: J. Kiehl
    !
    !-----------------------------------------------------------------------

    ! Arguments
    character(len=*), intent(in)  :: name              ! constituent name
    integer,          intent(in)  :: ncol              ! number of columns
    integer,          intent(in)  :: pcols
    integer,          intent(in)  :: pver
    real(rtype),      intent(in)  :: clat(pcols)       ! latitude in radians for columns
    real(rtype),      intent(in)  :: pmid(pcols,pver)  ! model pressures
    real(rtype),      intent(out) :: q(pcols,pver)     ! constituent mass mixing ratio

    real(rtype), intent(in) :: mwdry, mwco2, mwn2o, mwch4, mwf11, mwf12

    real(rtype), intent(in) :: o2mmr, co2vmr_rad, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr

    integer i                ! longitude loop index
    integer k                ! level index

    real(rtype) coslat(pcols)   ! cosine of latitude
    real(rtype) dlat            ! latitude in degrees
    real(rtype) ptrop           ! pressure level of tropopause
    real(rtype) pratio          ! pressure divided by ptrop
    real(rtype) trop_mmr        ! tropospheric mass mixing ratio
    real(rtype) scale           ! pressure scale height

    real(rtype) :: rmwn2o
    real(rtype) :: rmwch4
    real(rtype) :: rmwf11
    real(rtype) :: rmwf12
    real(rtype) :: rmwco2

    rmwn2o = mwn2o/mwdry
    rmwch4 = mwch4/mwdry
    rmwf11 = mwf11/mwdry
    rmwf12 = mwf12/mwdry
    rmwco2 = mwco2/mwdry

    !-----------------------------------------------------------------------

    do i = 1, ncol
       coslat(i) = cos(clat(i))
    end do

    if (name == 'O2') then

       q = chem_surfvals_get('O2MMR', rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr)

    else if (name == 'CO2') then

       q = chem_surfvals_co2_rad(rmwco2, co2vmr, co2vmr_rad)

    else if (name == 'CH4') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwch4 * chem_surfvals_get('CH4VMR', rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr)

       do k = 1,pver
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_rtype * clat(i))
             if(dlat.le.45.0_rtype) then
                scale = 0.2353_rtype
             else
                scale = 0.2353_rtype + 0.0225489_rtype * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_rtype - 150.0e2_rtype*coslat(i)**2.0_rtype

             ! determine output mass mixing ratios
             if (pmid(i,k) >= ptrop) then
                q(i,k) = trop_mmr
             else
                pratio = pmid(i,k)/ptrop
                q(i,k) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else if (name == 'N2O') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwn2o * chem_surfvals_get('N2OVMR', rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr)

       do k = 1,pver
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_rtype * clat(i))
             if(dlat.le.45.0_rtype) then
                scale = 0.3478_rtype + 0.00116_rtype * dlat
             else
                scale = 0.4000_rtype + 0.013333_rtype * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_rtype - 150.0e2_rtype*coslat(i)**2.0_rtype

             ! determine output mass mixing ratios
             if (pmid(i,k) >= ptrop) then
                q(i,k) = trop_mmr
             else
                pratio = pmid(i,k)/ptrop
                q(i,k) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else if (name == 'CFC11') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwf11 * chem_surfvals_get('F11VMR', rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr)

       do k = 1,pver
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_rtype * clat(i))
             if(dlat.le.45.0_rtype) then
                scale = 0.7273_rtype + 0.00606_rtype * dlat
             else
                scale = 1.00_rtype + 0.013333_rtype * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_rtype - 150.0e2_rtype*coslat(i)**2.0_rtype

             ! determine output mass mixing ratios
             if (pmid(i,k) >= ptrop) then
                q(i,k) = trop_mmr
             else
                pratio = pmid(i,k)/ptrop
                q(i,k) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else if (name == 'CFC12') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwf12 * chem_surfvals_get('F12VMR', rmwco2, co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, o2mmr)

       do k = 1,pver
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_rtype * clat(i))
             if(dlat.le.45.0_rtype) then
                scale = 0.4000_rtype + 0.00222_rtype * dlat
             else
                scale = 0.50_rtype + 0.024444_rtype * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_rtype - 150.0e2_rtype*coslat(i)**2.0_rtype

             ! determine output mass mixing ratios
             if (pmid(i,k) >= ptrop) then
                q(i,k) = trop_mmr
             else
                pratio = pmid(i,k)/ptrop
                q(i,k) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    end if

  end subroutine trcmix

end module scream_trcmix
