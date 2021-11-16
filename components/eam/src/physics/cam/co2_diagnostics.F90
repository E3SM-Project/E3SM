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
use ppgrid,         only: pver, pcols, begchunk, endchunk
use physics_types,  only: physics_state, physics_tend, physics_ptend, physics_ptend_init
use constituents,   only: pcnst, cnst_name
use cam_logfile,    only: iulog
use spmd_utils,     only: masterproc

implicit none
private
save

public co2_diags_register            ! setup co2 pbuf fields
public check_co2_change_pr2
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

   subroutine co2_diags_register()
      !-------------------------------------------------
      ! Purpose: register co2 fields into pbuf
      !-------------------------------------------------
      use physics_buffer, only: pbuf_add_field, dtype_r8

      integer :: m,idx

      if (co2_transport()) then
         do m = 1,4
!            call outfld(trim(cnst_name(c_i(m)))//'_BOT', state%q(1,pver,c_i(m)), pcols, lchnk)
!            call outfld(sflxnam(c_i(m)), cam_in%cflx(:,c_i(m)), pcols, lchnk)
!            call pbuf_add_field(rflx_pname(i), 'physpkg', dtype_r8, (/pcols/), idx)
!
            call pbuf_add_field(trim(cnst_name(c_i(m)))//'_cur', 'physpkg', dtype_r8, (/pcols/), idx)
         enddo
      endif

   end subroutine co2_diags_register

   subroutine check_co2_change_pr2(state, tend, pbuf2d, cam_in, wet_or_dry)
      use physconst,      only: gravit
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, pbuf_set_field 
      use phys_control,   only: ieflx_opt
      use phys_gmean,     only: gmean

!      integer , intent(in) :: nstep        ! current timestep number
      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
      type(physics_tend ), intent(in   ), dimension(begchunk:endchunk) :: tend
      type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
      type(physics_buffer_desc),          pointer :: pbuf2d(:,:)
      character(len=3),    intent(in   ) :: wet_or_dry    ! is co2 wet or dry at this point
 
      integer :: ncol                      ! number of active columns
      integer :: lchnk                     ! chunk index

      integer :: i, k, m

      character(len=*), parameter :: title='CO2 end of phys run2:)'
      character(len=*), parameter :: sub_name='co2_gmean_check: '


      real(r8) :: mass_wet_mean(ncnst)
      real(r8) :: mass_dry_mean(ncnst)
      real(r8) :: sfc_flux_mean(ncnst)
      real(r8) :: air_flux_mean(ncnst)

      real(r8) :: mass_wet(pcols, begchunk:endchunk, ncnst)
      real(r8) :: mass_dry(pcols, begchunk:endchunk, ncnst)
      real(r8) :: sfc_flux(pcols, begchunk:endchunk, ncnst)
      real(r8) :: air_flux(pcols, begchunk:endchunk, ncnst)
      

      if ( .not. co2_transport() ) return

      do lchnk = begchunk, endchunk
!         qflx(:ncol,lchnk) = cam_in(lchnk)%cflx(:ncol,1)
         mass_wet(:ncol, lchnk, :) = 0._r8
         mass_dry(:ncol, lchnk, :) = 0._r8
         sfc_flux(:ncol, lchnk, :) = 0._r8
         air_flux(:ncol, lchnk, :) = 0._r8
         do i = 1, ncol
            do k = 1, pver
               do m = 1, ncnst
                  mass_wet(i, lchnk, m) = mass_wet(i, lchnk, m) + &
                       state(lchnk)%pdel(i, k)*state(lchnk)%q(i, k, c_i(m))
                  mass_dry(i, lchnk, m) = mass_dry(i, lchnk, m) + &
                       state(lchnk)%pdel(i, k)*state(lchnk)%q(i, k, c_i(m))
               end do ! m = 1, ncnst
!               air_flux(i, lchnk) = air_flux(i, lchnk) + ac_CO2(i, k)
            end do ! k = 1, pver
            do m = 1, ncnst
               mass_wet(i, lchnk, m) = mass_wet(i, lchnk, m) / gravit
               mass_dry(i, lchnk, m) = mass_dry(i, lchnk, m) / gravit
               sfc_flux(i, lchnk, m) = sfc_flux(i, lchnk, m) + &
                    cam_in(lchnk)%cflx(i, c_i(m))
               air_flux_mean(m) = 0._r8  ! DONT KEEP THIS!!!
            end do ! m = 1, ncnst
         end do ! i = 1, ncol
      end do ! lchnk = begchunk, endchunk

      call gmean(mass_wet, mass_wet_mean, ncnst)
      call gmean(mass_dry, mass_dry_mean, ncnst)
      call gmean(sfc_flux, sfc_flux_mean, ncnst)
!      call gmean(air_flux, air_flux_mean, ncnst)

      if (begchunk .le. endchunk) then
         if (masterproc) then
            do m = 1, ncnst
               write (6,'(a32,i2,a36,1p,4e25.15)') trim(title)//' m=',c_i(m), &
                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc, air ', &
                  mass_dry_mean(m), mass_wet_mean(m), &
                  sfc_flux_mean(m), air_flux_mean(m)
            end do
         end if
      end if
      
   end subroutine check_co2_change_pr2

!   subroutine check_co2_change(wet_or_dry, state, ptend, cam_in, pbuf)
!-----------------------------------------------------------------------
!
! Purpose:
! Computes global mean mass, max and min mmr, of constituents on the 
! physics decomposition. Checks for conservation from previous timestep
!
! Author: Bryce Harrop
!
!-----------------------------------------------------------------------
!      use physconst,      only: gravit
!      use phys_grid,      only: get_ncols_p
!      use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
!      use phys_gmean,     only: gmean
!
! Arguments
!
!      character(len=3),    intent(in) :: wet_or_dry    ! is co2 wet or dry at this point
!      type(physics_state), intent(in) :: state
!      type(physics_buffer_desc), pointer :: pbuf(:)
!      type(physics_ptend), intent(in) :: ptend
!      type(cam_in_t),      intent(in) :: cam_in
      !integer, dimension(:), intent(in) :: cnst_ind_arr
!
! Local workspace
!
!      character(len=*), parameter :: sub_name='check_co2_change: '

!      integer :: c, i, k, m
!      integer :: ierr, ncol, idx
!      integer :: index_co2_ocn, index_co2_fff, index_co2_lnd, index_co2
!      real(r8)          :: ac_CO2_tot(pcols)
!      real(r8), pointer :: ac_CO2(:,:)

!      real(r8), pointer :: mass_wet(:,:,:) ! constituent masses assuming moist mmr
!      real(r8), pointer :: mass_dry(:,:,:) ! constituent masses assuming dry mmr
!      real(r8) :: mass_wet_mean(ncnst)     ! global mean constituent masses assuming moist mmr
!      real(r8) :: mass_dry_mean(ncnst)     ! global mean constituent masses assuming dry mmr
!      real(r8), pointer :: sfc_flux(:,:,:) ! constituent surface flux
!      real(r8) :: sfc_flux_mean(ncnst)     ! global mean constituent surface flux
!      real(r8), pointer :: air_flux(:,:)   ! aircraft emission flux
!      real(r8) :: air_flux_mean            ! global mean aircraft emission flux

!      real(r8), pointer :: mass_cur(:,:,:) ! previous timestep mass
!      real(r8) :: mass_mean_cur(ncnst)

!      real(r8), pointer :: co2_cur(:,:,:)  ! previous co2 mixing ratios
!
!-----------------------------------------------------------------------
!
!      if ( .not. co2_transport() ) return

!      allocate(mass_wet(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

!      allocate(mass_dry(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate mass_wet'

!      allocate(sfc_flux(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate sfc_flux'

!      allocate(air_flux(pcols,begchunk:endchunk), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate air_flux'

!      allocate(mass_cur(pcols,begchunk:endchunk,ncnst), stat=ierr)
!      if (ierr /= 0) write(iulog,*) sub_name // 'FAIL to allocate air_flux'

!      index_co2_ocn = pbuf_get_index(trim(cnst_name(c_i(1)))//'_cur')
!      index_co2_fff = pbuf_get_index(trim(cnst_name(c_i(2)))//'_cur')
!      index_co2_lnd = pbuf_get_index(trim(cnst_name(c_i(3)))//'_cur')
!      index_co2     = pbuf_get_index(trim(cnst_name(c_i(4)))//'_cur')
!      index_acco2   = pbuf_get_index('ac_CO2')  

!      call pbuf_get_field(pbuf, index_co2_ocn, co2_ocn_cur)
!      call pbuf_get_field(pbuf, index_co2_fff, co2_fff_cur)
!      call pbuf_get_field(pbuf, index_co2_lnd, co2_lnd_cur)
!      call pbuf_get_field(pbuf, index_co2    , co2_cur)
!      call pbuf_get_field(pbuf, index_acco2  , ac_CO2)

      ! Set CO2 global indices
!      do idx = 1, ncnst
!         select case (trim(c_names(idx)))
!         case ('CO2_OCN')
!            co2_ocn_glo_ind = c_i(idx)
!         case ('CO2_FFF')
!            co2_fff_glo_ind = c_i(idx)
!         case ('CO2_LND')
!            co2_lnd_glo_ind = c_i(idx)
!         case ('CO2')
!            co2_glo_ind     = c_i(idx)
!         end select
!      end do

!      ncol = state%ncol
!      c    = state%lchnk
!      !ncnst = size(cnst_ind_arr)
!      do m = 1, ncnst
!!         do c = begchunk, endchunk
!!            ncols = get_ncols_p(c)
!            do i = 1, ncol

!               ! Compute column masses assuming both dry and wet mixing ratios

!               mass_wet(i,c,m) = 0.0_r8
!               do k = 1, pver
!                  mass_wet(i,c,m) = mass_wet(i,c,m) + &
!                                    state%pdel(i,k)*state%q(i,k,c_i(m))
!               end do
!               mass_wet(i,c,m) = mass_wet(i,c,m)/gravit

!               mass_dry(i,c,m) = 0.0_r8
!               do k = 1, pver
!                  mass_dry(i,c,m) = mass_dry(i,c,m) + &
!                                    state%pdeldry(i,k)*state%q(i,k,c_i(m))
!               end do
!               mass_dry(i,c,m) = mass_dry(i,c,m)/gravit

!               sfc_flux(i,c,m) = 0.0_r8
!               sfc_flux(i,c,m) = sfc_flux(i,c,m) + &
!                                 cam_in%cflx(i,c_i(m))
               
!!               ac_CO2_tot(i) = 0.0_r8
!!               do k = 1, pver
!!                  ac_CO2_tot(i) = ac_CO2_tot(i) + ac_CO2(i,k)
!!               end do

!            end do
!!         end do
!!            if (m .eq. co2_glo_ind .or. m .eq. co2_fff_glo_ind) then
!!               air_flux(:ncol,c,m) = ac_CO2_tot(:ncol)
!!            else
!!               air_flux(:ncol,c,m) = 0.0_r8
!!            end if
!      end do

!      do i = 1, ncol
!         ac_CO2_tot(i) = 0.0_r8
!         do k = 1, pver
!            ac_CO2_tot(i) = ac_CO2_tot(i) + ac_CO2(i,k)
!         end do
!         air_flux(:ncol,c) = ac_CO2_tot(:ncol)
!      end do

!      do i = 1, ncol
!         mass_cur(i,c,:) = 0._r8
!         if (wet_or_dry == 'wet') then
!            mass_cur(i,c,1) = mass_wet(i,c,1) + &
!                                 state%pdel(i,k)*co2_ocn_cur(i,k)
!            mass_cur(i,c,2) = mass_wet(i,c,2) + &
!                                 state%pdel(i,k)*co2_fff_cur(i,k)
!            mass_cur(i,c,3) = mass_wet(i,c,3) + &
!                                 state%pdel(i,k)*co2_lnd_cur(i,k)
!            mass_cur(i,c,4) = mass_wet(i,c,4) + &
!                                 state%pdel(i,k)*co2_cur(i,k)
!         end if
!      end do

!      ! compute global mean mass
!      call gmean(mass_wet, mass_wet_mean, ncnst)
!      call gmean(mass_dry, mass_dry_mean, ncnst)
!      call gmean(sfc_flux, sfc_flux_mean, ncnst)
!      call gmean(air_flux, air_flux_mean)

!      call gmean(mass_cur, mass_mean_cur, ncnst)

!      ! report to log file
!      if (masterproc) then

!         do m = 1, ncnst
!               write (6,66) trim(title)//' m=',c_i(m), &
!                  'name='//trim(cnst_name(c_i(m)))//' gavg dry, wet, sfc, air ', &
!                  mass_dry_mean(c_i(m)), mass_wet_mean(c_i(m)), &
!                  sfc_flux_mean(c_i(m)), air_flux_mean(c_i(m))
!66             format (a32,i2,a36,1p,4e25.15)
!!old version66             format (a32,i2,a36,1p,4e25.13)
!         end do

!      endif

!      deallocate(air_flux)
!      deallocate(sfc_flux)
!      deallocate(mass_dry)
!      deallocate(mass_wet)

!   end subroutine check_co2_change

!   subroutine check_co2_change_nope(state, pbuf2d)
      !-----------------------------------------------
      ! Purpose: check timestep level co2 conservation
      !
      ! Author: Bryce Harrop
      !-----------------------------------------------
!      use physics_types,  only: physics_state
!      use ppgrid,         only: begchunk, endchunk, pver, pcols
!      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
!                                pbuf_get_chunk, pbuf_get_index
!      use physconst,      only: rga
!      use co2_cycle,      only: c_i, co2_transport

!      implicit none
      
!      type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
!      type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

!      !local vars
!      type(physics_buffer_desc), pointer     :: pbuf_chnk(:)

!      real(r8), pointer :: ac_CO2(:,:)

!      real(r8) :: ftem      (pcols,pver) ! tmp space
!      real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
!      real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor

!      real(r8) :: ftem(pcols,pver)         ! tmp space
!      real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
!      real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor
!      real(r8) :: CIDiff(pcols)            ! Difference in vertically integrated static energy

!      integer  :: index_co2_ocn, index_co2_fff, index_co2_lnd, index_co2
!      integer  :: ichnk, ncol, i, k, m, lchnk

!      if ( .not. co2_transport() ) return

!      do ichnk = begchunk, endchunk
!         ncol       =  state(ichnk)%ncol
!         pbuf_chunk => pbuf_get_chunk(pbuf2d, ichnk)

!         index_co2_ocn = pbuf_get_index(trim(cnst_name(c_i(1)))//'_cur')
!         index_co2_fff = pbuf_get_index(trim(cnst_name(c_i(2)))//'_cur')
!         index_co2_lnd = pbuf_get_index(trim(cnst_name(c_i(3)))//'_cur')
!         index_co2     = pbuf_get_index(trim(cnst_name(c_i(4)))//'_cur')
!         index_acco2   = pbuf_get_index('ac_CO2')   

!         call pbuf_get_field(pbuf_chnk, index_co2_ocn, co2_ocn_cur)
!         call pbuf_get_field(pbuf_chnk, index_co2_fff, co2_fff_cur)
!         call pbuf_get_field(pbuf_chnk, index_co2_lnd, co2_lnd_cur)
!         call pbuf_get_field(pbuf_chnk, index_co2    , co2_cur)
!         call pbuf_get_field(pbuf_chnk, index_acco2  , ac_CO2)

!         do i = 1, ncol
!            flux_ocn(i) = cam_in(ichnk)%cflx(i, c_i(1))
!            flux_fff(i) = cam_in(ichnk)%cflx(i, c_i(2)) + 
!            flux_lnd(i) = cam_in(ichnk)%cflx(i, c_i(3))
!            flux(i)     = flux_ocn(i) + flux_fff(i) + flux_lnd(i)
!         end do ! i = 1, ncol

!      end do ! ichnk = begchunk, endchunk


!      static_ener_ac_idx = pbuf_get_index('static_ener_ac')
!      call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
!      water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
!      call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

      !Integrate column static energy
!      ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
!    do k=2,pver
!       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
!    end do
!    static_ener_ac_2d(:ncol) = ftem(:ncol,1)

!    static_ener_ac_idx = pbuf_get_index('static_ener_ac')
!    call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
!    water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
!    call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

    ! Integrate and compute the difference
    ! CIDiff = difference of column integrated values
!    if( nstep == 0 ) then
!       CIDiff(:ncol) = 0.0_r8
!       call outfld('DTENDTH', CIDiff, pcols, lchnk )
!       call outfld('DTENDTQ', CIDiff, pcols, lchnk )
!    else
!       ! MSE first
!       ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
!       do k=2,pver
!          ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
!       end do
!       CIDiff(:ncol) = (ftem(:ncol,1) - static_ener_ac_2d(:ncol))*rtdt!

!       call outfld('DTENDTH', CIDiff, pcols, lchnk )
!       ! Water vapor second
!       ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
!       do k=2,pver
!          ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
!       end do
!       CIDiff(:ncol) = (ftem(:ncol,1) - water_vap_ac_2d(:ncol))*rtdt

!       call outfld('DTENDTQ', CIDiff, pcols, lchnk )
!    end if
      

!   end subroutine check_co2_change_nope

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

      deallocate(sfc_flux)
      deallocate(mass_dry)
      deallocate(mass_wet)

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
66             format (a32,i2,a36,1p,4e25.15)
!old version66             format (a32,i2,a36,1p,4e25.13)
         end do

      endif

      deallocate(air_flux)
      deallocate(sfc_flux)
      deallocate(mass_dry)
      deallocate(mass_wet)

   end subroutine co2_gmean_check2_wflux

end module co2_diagnostics
