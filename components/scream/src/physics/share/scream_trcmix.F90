
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

  ! Copied from physconst
  real(rtype), public, parameter :: mwdry       = 28.966_rtype ! molecular weight dry air ~ kg/kmole
  real(rtype), public, parameter :: mwco2       =  44._rtype   ! molecular weight co2
  real(rtype), public, parameter :: mwn2o       =  44._rtype   ! molecular weight n2o
  real(rtype), public, parameter :: mwch4       =  16._rtype   ! molecular weight ch4
  real(rtype), public, parameter :: mwf11       = 136._rtype   ! molecular weight cfc11
  real(rtype), public, parameter :: mwf12       = 120._rtype   ! molecular weight cfc12
  real(rtype), public, parameter :: mwo3        =  48._rtype   ! molecular weight O3
  real(rtype), public, parameter :: mwso2       =  64._rtype
  real(rtype), public, parameter :: mwso4       =  96._rtype
  real(rtype), public, parameter :: mwh2o2      =  34._rtype
  real(rtype), public, parameter :: mwdms       =  62._rtype
  real(rtype), public, parameter :: mwnh4       =  18._rtype

  real(rtype), public, parameter :: rmwn2o = mwn2o/mwdry ! ratio of molecular weight n2o   to dry air
  real(rtype), public, parameter :: rmwch4 = mwch4/mwdry ! ratio of molecular weight ch4   to dry air
  real(rtype), public, parameter :: rmwf11 = mwf11/mwdry ! ratio of molecular weight cfc11 to dry air
  real(rtype), public, parameter :: rmwf12 = mwf12/mwdry ! ratio of molecular weight cfc12 to dry air
  real(rtype), public, parameter :: rmwco2 = mwco2/mwdry ! ratio of molecular weights of co2 to dry air


  ! Default values for namelist variables -- now set by build-namelist
  real(rtype) :: o2mmr = .23143_rtype               ! o2 mass mixing ratio
  real(rtype) :: co2vmr_rad = -1.0_rtype            ! co2 vmr override for radiation
  real(rtype) :: co2vmr = -1.0_rtype                ! co2   volume mixing ratio
  real(rtype) :: n2ovmr = -1.0_rtype                ! n2o   volume mixing ratio
  real(rtype) :: ch4vmr = -1.0_rtype                ! ch4   volume mixing ratio
  real(rtype) :: f11vmr = -1.0_rtype                ! cfc11 volume mixing ratio
  real(rtype) :: f12vmr = -1.0_rtype                ! cfc12 volume mixing ratio

  !================================================================================================
contains
  !================================================================================================

  function chem_surfvals_get(name)
    ! Copied from components/eam/src/physics/chem_surfvals.F90

    character(len=*), intent(in) :: name

    real(rtype) :: rmwco2
    real(rtype) :: chem_surfvals_get

    rmwco2 = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
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

  function chem_surfvals_co2_rad()
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

    ! Return value
    real(rtype) :: chem_surfvals_co2_rad

    ! Local variables
    real(rtype) :: convert_vmr      ! convert vmr to desired output
    !-----------------------------------------------------------------------

    ! by default convert vmr to mmr
    convert_vmr = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air

    if (co2vmr_rad > 0._rtype) then
       chem_surfvals_co2_rad = convert_vmr * co2vmr_rad
    else
       chem_surfvals_co2_rad = convert_vmr * co2vmr
    end if

  end function chem_surfvals_co2_rad

  subroutine trcmix(name, ncol, pcols, pver, clat, pmid, q)
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

    integer i                ! longitude loop index
    integer k                ! level index

    real(rtype) coslat(pcols)   ! cosine of latitude
    real(rtype) dlat            ! latitude in degrees
    real(rtype) ptrop           ! pressure level of tropopause
    real(rtype) pratio          ! pressure divided by ptrop
    real(rtype) trop_mmr        ! tropospheric mass mixing ratio
    real(rtype) scale           ! pressure scale height

    !-----------------------------------------------------------------------

    do i = 1, ncol
       coslat(i) = cos(clat(i))
    end do

    if (name == 'O2') then

       q = chem_surfvals_get('O2MMR')

    else if (name == 'CO2') then

       q = chem_surfvals_co2_rad()

    else if (name == 'CH4') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwch4 * chem_surfvals_get('CH4VMR')

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
       trop_mmr = rmwn2o * chem_surfvals_get('N2OVMR')

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
       trop_mmr = rmwf11 * chem_surfvals_get('F11VMR')

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
       trop_mmr = rmwf12 * chem_surfvals_get('F12VMR')

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
