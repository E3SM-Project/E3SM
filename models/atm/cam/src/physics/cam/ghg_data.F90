
module ghg_data

!------------------------------------------------------------------------------------------------
! Purpose:
! Provide default distributions of CH4, N2O, CFC11 and CFC12 to the radiation routines.
! **NOTE** CO2 is assumed by the radiation to a be constant value.  This value is
!          currently supplied directly by the chem_surfvals module.
!
! Revision history:
! 2004-08-29  B. Eaton        Create CAM interface to trcmix.
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pcols, pver, begchunk, endchunk
use physics_types,  only: physics_state
use physconst,      only: mwdry, mwch4, mwn2o, mwf11, mwf12, mwco2
use chem_surfvals,  only: chem_surfvals_get, chem_surfvals_co2_rad
use abortutils,     only: endrun
use error_messages, only: handle_err


implicit none
private
save

! Public interfaces
public ::&
   ghg_data_register, &! register ghg's with pbuf2d
   ghg_data_timestep_init    ! place data model of ghg's in pbuf2d

! Private variables

real(r8) :: rmwn2o ! = mwn2o/mwdry ! ratio of molecular weight n2o   to dry air
real(r8) :: rmwch4 ! = mwch4/mwdry ! ratio of molecular weight ch4   to dry air
real(r8) :: rmwf11 ! = mwf11/mwdry ! ratio of molecular weight cfc11 to dry air
real(r8) :: rmwf12 ! = mwf12/mwdry ! ratio of molecular weight cfc12 to dry air
real(r8) :: rmwco2 ! = mwco2/mwdry ! ratio of molecular weights of co2 to dry air

integer, parameter :: ncnst = 6                        ! number of constituents
character(len=8), dimension(ncnst), parameter :: &
   cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12', 'CO2  ', 'O2   '/) ! constituent names
integer  :: pbuf_idx(ncnst)

!================================================================================================
contains
!================================================================================================

subroutine ghg_data_register()
!-------------------------------------------------------------------------------
! register ghg's with pbuf2d
!-------------------------------------------------------------------------------
  use physics_buffer, only : pbuf_add_field, dtype_r8

  integer iconst

 
  do iconst = 1,ncnst
     call pbuf_add_field(cnst_names(iconst),'physpkg',dtype_r8,(/pcols,pver/),pbuf_idx(iconst))
  enddo

end subroutine ghg_data_register

subroutine ghg_data_timestep_init(pbuf2d, state)
!-------------------------------------------------------------------------------
! place data model of ghg's in pbuf2d at each timestep
!-------------------------------------------------------------------------------
  use ppgrid,              only: begchunk, endchunk, pcols, pver
  use physics_types,       only: physics_state
  use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

  
  type(physics_state), intent(in), dimension(begchunk:endchunk) :: state
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
 
  type(physics_buffer_desc), pointer :: pbuf_chnk(:)
  real(r8), pointer :: tmpptr(:,:)

  integer iconst
  integer lchnk

  rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
  rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
  rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
  rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
  rmwco2 = mwco2/mwdry      ! ratio of molecular weights of co2 to dry air

   do iconst = 1,ncnst
!$OMP PARALLEL DO PRIVATE (LCHNK,tmpptr,pbuf_chnk)
     do lchnk = begchunk, endchunk
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call pbuf_get_field(pbuf_chnk, pbuf_idx(iconst), tmpptr) 
       call trcmix(cnst_names(iconst), state(lchnk)%ncol, &
                   state(lchnk)%lat, state(lchnk)%pmid, &
                   tmpptr)
     enddo
  enddo

end subroutine ghg_data_timestep_init


!================================================================================================

subroutine trcmix(name, ncol, clat, pmid, q)
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
   real(r8),         intent(in)  :: clat(pcols)       ! latitude in radians for columns
   real(r8),         intent(in)  :: pmid(pcols,pver)  ! model pressures
   real(r8),         intent(out) :: q(pcols,pver)     ! constituent mass mixing ratio

   integer i                ! longitude loop index
   integer k                ! level index

   real(r8) coslat(pcols)   ! cosine of latitude
   real(r8) dlat            ! latitude in degrees
   real(r8) ptrop           ! pressure level of tropopause
   real(r8) pratio          ! pressure divided by ptrop
   real(r8) trop_mmr        ! tropospheric mass mixing ratio
   real(r8) scale           ! pressure scale height
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
            dlat = abs(57.2958_r8 * clat(i))
            if(dlat.le.45.0_r8) then
               scale = 0.2353_r8
            else
               scale = 0.2353_r8 + 0.0225489_r8 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

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
            dlat = abs(57.2958_r8 * clat(i))
            if(dlat.le.45.0_r8) then
               scale = 0.3478_r8 + 0.00116_r8 * dlat
            else
               scale = 0.4000_r8 + 0.013333_r8 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

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
            dlat = abs(57.2958_r8 * clat(i))
            if(dlat.le.45.0_r8) then
               scale = 0.7273_r8 + 0.00606_r8 * dlat
            else
               scale = 1.00_r8 + 0.013333_r8 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

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
            dlat = abs(57.2958_r8 * clat(i))
            if(dlat.le.45.0_r8) then
               scale = 0.4000_r8 + 0.00222_r8 * dlat
            else
               scale = 0.50_r8 + 0.024444_r8 * (dlat - 45)
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

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

end module ghg_data
