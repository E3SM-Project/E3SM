module co2_diagnostics

!------------------------------------------------------------------------------------------------
 
! Purpose:
! Write out mass of CO2 in each tracer (total, fossil fuel, land, and ocean)
! I wanted this to go in co2_cyle.F90, but it created a circular dependency
! with camsrfexch.F90, so this became its own module.
!
! Author: Bryce Harrop
! 
!                                              
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8, cxx =>SHR_KIND_CXX, cl =>SHR_KIND_CL
use camsrfexch,     only: cam_in_t
use co2_cycle,      only: c_i, co2_transport

implicit none
private
save

public co2_gmean_check_wflux         ! printout co2 global means
public co2_gmean_check2_wflux        ! higher level co2 checker

! Number of CO2 tracers
integer, parameter :: ncnst = 4                      ! number of constituents implemented

character(len=7), dimension(ncnst), parameter :: & ! constituent names
     c_names = (/'CO2_OCN', 'CO2_FFF', 'CO2_LND', 'CO2    '/)

integer :: co2_ocn_glo_ind ! global index of 'CO2_OCN'
integer :: co2_fff_glo_ind ! global index of 'CO2_FFF'
integer :: co2_lnd_glo_ind ! global index of 'CO2_LND'
integer :: co2_glo_ind     ! global index of 'CO2'

contains

   subroutine co2_gmean_check_wflux(title, state, cam_in)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Prints diagnostics to log file.
!
! Authors: B. Eaton (based on gavglook) & Bryce Harrop
!
!-----------------------------------------------------------------------
      use ppgrid,         only: pver, pcols, begchunk, endchunk
      use physconst,      only: gravit
      use phys_grid,      only: get_ncols_p
      use physics_types,  only: physics_state, physics_ptend
      use constituents,   only: pcnst, cnst_name
      use cam_logfile,    only: iulog
      use phys_gmean,     only: gmean
      use spmd_utils,     only: masterproc
!
! Arguments
!
      character(len=*),    intent(in) :: title    ! location of this call
      type(physics_state), intent(in) :: state(begchunk:endchunk)
      type(cam_in_t),      intent(in) :: cam_in(begchunk:endchunk)
!
! Local workspace
!
      character(len=*), parameter :: sub_name='co2_gmean_check: '

      integer :: c, i, k, m
      integer :: ierr
      integer :: ncols

      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
      real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
      real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
      real(r8), pointer :: sfc_flux(:,:,:) ! constituent surface flux
      real(r8) :: sfc_flux_mean(pcnst)     ! global mean constituent surface flux

!
!-----------------------------------------------------------------------
!
      allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_dry'

      allocate(sfc_flux(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate sfc_flux'

      do m = 1, pcnst
         do c = begchunk, endchunk
            ncols = get_ncols_p(c)
            do i = 1, ncols

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
                                    state(c)%pdel(i,k) * state(c)%q(i,k,m)
               end do
               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

               mass_dry(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
                                    state(c)%pdeldry(i,k) * state(c)%q(i,k,m)
               end do
               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

               sfc_flux(i,c,m) = 0.0_r8
               sfc_flux(i,c,m) = sfc_flux(i,c,m) + &
                                 cam_in(c)%cflx(i,m)

            end do
         end do
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)
      call gmean(sfc_flux, sfc_flux_mean, pcnst)

      ! report to log file
      if (masterproc) then

         do m = 1, ncnst
               write (6,66) trim(title)//' m=',c_i(m), &
                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc ', &
                  mass_dry_mean(c_i(m)), mass_wet_mean(c_i(m)), &
                  sfc_flux_mean(c_i(m))
66             format (a32,i2,a36,1p,3e25.13)
         end do

      endif

      deallocate(mass_wet)
      deallocate(mass_dry)
      deallocate(sfc_flux)

   end subroutine co2_gmean_check_wflux

   subroutine co2_gmean_check2_wflux(title, state, ptend, cam_in, pbuf)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Prints diagnostics to log file.
!
! Authors: B. Eaton (based on gavglook) & Bryce Harrop
!
!-----------------------------------------------------------------------
      use ppgrid,         only: pver, pcols, begchunk, endchunk
      use physconst,      only: gravit
      use phys_grid,      only: get_ncols_p
      use physics_types,  only: physics_state, physics_ptend
      use constituents,   only: pcnst, cnst_name
      use cam_logfile,    only: iulog
      use phys_gmean,     only: gmean
      use spmd_utils,     only: masterproc
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
!
! Arguments
!
      character(len=*),    intent(in) :: title    ! location of this call
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(physics_ptend), intent(in) :: ptend
      type(cam_in_t),      intent(in) :: cam_in
      !integer, dimension(:), intent(in) :: cnst_ind_arr
!
! Local workspace
!
      character(len=*), parameter :: sub_name='co2_gmean_check2: '

      integer :: c, i, k, m
      integer :: ierr
      integer :: ncol
      integer :: ifld, idx
      real(r8)          :: ac_CO2_tot(pcols)
      real(r8), pointer :: ac_CO2(:,:)

      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
      real(r8) :: mass_wet_mean(pcnst)     ! global mean constituent masses assuming moist mmr
      real(r8) :: mass_dry_mean(pcnst)     ! global mean constituent masses assuming dry mmr
      real(r8), pointer :: sfc_flux(:,:,:) ! constituent surface flux
      real(r8) :: sfc_flux_mean(pcnst)     ! global mean constituent surface flux
      real(r8), pointer :: air_flux(:,:,:) ! aircraft emission flux
      real(r8) :: air_flux_mean(pcnst)     ! global mean aircraft emission flux

!
!-----------------------------------------------------------------------
!
      allocate(mass_wet(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(mass_dry(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

      allocate(sfc_flux(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate sfc_flux'

      allocate(air_flux(pcols,begchunk:endchunk,pcnst), stat=ierr)
      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate air_flux'

      ifld = pbuf_get_index('ac_CO2')   
      call pbuf_get_field(pbuf, ifld, ac_CO2)

      ! Set CO2 global indices
      do idx = 1, ncnst
         select case (trim(c_names(idx)))
         case ('CO2_OCN')
            co2_ocn_glo_ind = c_i(idx)
         case ('CO2_FFF')
            co2_fff_glo_ind = c_i(idx)
         case ('CO2_LND')
            co2_lnd_glo_ind = c_i(idx)
         case ('CO2')
            co2_glo_ind     = c_i(idx)
         end select
      end do

      ncol = state%ncol
      c    = state%lchnk
      !ncnst = size(cnst_ind_arr)
      do m = 1, pcnst
!         do c = begchunk, endchunk
!            ncols = get_ncols_p(c)
            do i = 1, ncol

               ! Compute column masses assuming both dry and wet mixing ratios

               mass_wet(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
                                    state%pdel(i,k)*state%q(i,k,m)
               end do
               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

               mass_dry(i,c,m) = 0.0_r8
               do k = 1, pver
                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
                                    state%pdeldry(i,k)*state%q(i,k,m)
               end do
               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

               sfc_flux(i,c,m) = 0.0_r8
               sfc_flux(i,c,m) = sfc_flux(i,c,m) + &
                                 cam_in%cflx(i,m)
               
               ac_CO2_tot(i) = 0.0_r8
               do k = 1, pver
                  ac_CO2_tot(i) = ac_CO2_tot(i) + ac_CO2(i,k)
               end do

            end do
!         end do
            if (m .eq. co2_glo_ind .or. m .eq. co2_fff_glo_ind) then
               air_flux(:ncol,c,m) = ac_CO2_tot(:ncol)
            else
               air_flux(:ncol,c,m) = 0.0_r8
            end if
      end do

      ! compute global mean mass
      call gmean(mass_wet, mass_wet_mean, pcnst)
      call gmean(mass_dry, mass_dry_mean, pcnst)
      call gmean(sfc_flux, sfc_flux_mean, pcnst)
      call gmean(air_flux, air_flux_mean, pcnst)

      ! report to log file
      if (masterproc) then

         do m = 1, ncnst
               write (6,66) trim(title)//' m=',c_i(m), &
                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc, air ', &
                  mass_dry_mean(c_i(m)), mass_wet_mean(c_i(m)), &
                  sfc_flux_mean(c_i(m)), air_flux_mean(c_i(m))
66             format (a32,i2,a36,1p,4e25.13)
         end do

      endif

      deallocate(mass_wet)
      deallocate(mass_dry)
      deallocate(sfc_flux)
      deallocate(air_flux)

   end subroutine co2_gmean_check2_wflux

end module co2_diagnostics
