! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

MODULE MOD_COSP_RADAR
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  USE MOD_COSP_UTILS
  use radar_simulator_types
  use array_lib
  use atmos_lib
  use format_input
  IMPLICIT NONE
  
  INTERFACE
    subroutine radar_simulator(freq,k2,do_ray,use_gas_abs,use_mie_table,mt, &
        nhclass,hp,nprof,ngate,nsizes,D,hgt_matrix,hm_matrix,re_matrix,p_matrix,t_matrix, &
        rh_matrix,Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe, &
        g_to_vol_in,g_to_vol_out)
  
        use m_mrgrnk 
        use array_lib
        use math_lib
        use optics_lib
        use radar_simulator_types
        implicit none
        ! ----- INPUTS -----  
        type(mie), intent(in) :: mt
        type(class_param) :: hp
        real*8, intent(in) :: freq,k2
        integer, intent(in) ::  do_ray,use_gas_abs,use_mie_table, &
            nhclass,nprof,ngate,nsizes
        real*8, dimension(nsizes), intent(in) :: D
        real*8, dimension(nprof,ngate), intent(in) :: hgt_matrix, p_matrix, &
            t_matrix,rh_matrix
        real*8, dimension(nhclass,nprof,ngate), intent(in) :: hm_matrix
        real*8, dimension(nhclass,nprof,ngate), intent(inout) :: re_matrix
        ! ----- OUTPUTS -----
        real*8, dimension(nprof,ngate), intent(out) :: Ze_non,Ze_ray, &
            g_atten_to_vol,dBZe,h_atten_to_vol    
        ! ----- OPTIONAL -----
        real*8, optional, dimension(ngate,nprof) :: &
            g_to_vol_in,g_to_vol_out
     end subroutine radar_simulator
  END INTERFACE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_RADAR ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_RADAR(gbx,sgx,sghydro,z)
  IMPLICIT NONE

  ! Arguments
  type(cosp_gridbox),intent(inout) :: gbx  ! Gridbox info
  type(cosp_subgrid),intent(in) :: sgx  ! Subgrid info
  type(cosp_sghydro),intent(in) :: sghydro  ! Subgrid info for hydrometeors
  type(cosp_sgradar),intent(inout) :: z ! Output from simulator, subgrid

  ! Local variables 
  real*8 :: &
  freq, &           ! radar frequency (GHz)
  k2                ! |K|^2, -1=use frequency dependent default
  
  real*8, dimension(:,:), allocatable :: &
  g_to_vol ! integrated atten due to gases, r>v (dB)
  
  real*8, dimension(:,:), allocatable :: &
  Ze_non, &         ! radar reflectivity withOUT attenuation (dBZ)
  Ze_ray, &         ! Rayleigh reflectivity (dBZ)
  h_atten_to_vol, &     ! attenuation by hydromets, radar to vol (dB)
  g_atten_to_vol, &     ! gaseous atteunation, radar to vol (dB)
  dBZe, &           ! effective radar reflectivity factor (dBZ)
  hgt_matrix, &         ! height of hydrometeors (km)
  t_matrix, &                   !temperature (k)
  p_matrix, &                   !pressure (hPa)
  rh_matrix                     !relative humidity (%)
  
  real*8, dimension(:,:,:), allocatable :: &
  hm_matrix, &          ! hydrometeor mixing ratio (g kg^-1)
  re_matrix

  integer, parameter :: one = 1
  logical :: hgt_reversed
  integer :: pr,i

! ----- main program settings ------

  freq = gbx%radar_freq
  k2 = gbx%k2
 
  !
  ! note:  intitialization section that was here has been relocated to SUBROUTINE CONSTRUCT_COSP_GRIDBOX by roj, Feb 2008
  !
  mt_ttl=gbx%mt_ttl  ! these variables really should be moved into the mt structure rather than kept as global arrays.
  mt_tti=gbx%mt_tti

  ! Inputs to Quickbeam
  allocate(hgt_matrix(gbx%Npoints,gbx%Nlevels),p_matrix(gbx%Npoints,gbx%Nlevels), &
           t_matrix(gbx%Npoints,gbx%Nlevels),rh_matrix(gbx%Npoints,gbx%Nlevels))
  allocate(hm_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels)) 
  allocate(re_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))

  ! Outputs from Quickbeam
  allocate(Ze_non(gbx%Npoints,gbx%Nlevels))
  allocate(Ze_ray(gbx%Npoints,gbx%Nlevels))
  allocate(h_atten_to_vol(gbx%Npoints,gbx%Nlevels))
  allocate(g_atten_to_vol(gbx%Npoints,gbx%Nlevels))
  allocate(dBZe(gbx%Npoints,gbx%Nlevels))
  
  ! Optional argument. It is computed and returned in the first call to
  ! radar_simulator, and passed as input in the rest
  allocate(g_to_vol(gbx%Nlevels,gbx%Npoints))
  
  p_matrix   = gbx%p/100.0     ! From Pa to hPa
  hgt_matrix = gbx%zlev/1000.0 ! From m to km
  t_matrix   = gbx%T-273.15    ! From K to C
  rh_matrix  = gbx%q
  re_matrix  = 0.0
  
  ! Quickbeam assumes the first row is closest to the radar
  call order_data(hgt_matrix,hm_matrix,p_matrix,t_matrix, &
      rh_matrix,gbx%surface_radar,hgt_reversed)
  
  ! ----- loop over subcolumns -----
  do pr=1,sgx%Ncolumns
      !  atmospheric profiles are the same within the same gridbox
      !  only hydrometeor profiles will be different
      if (hgt_reversed) then  
         do i=1,gbx%Nhydro  
            hm_matrix(i,:,:) = sghydro%mr_hydro(:,pr,gbx%Nlevels:1:-1,i)*1000.0 ! Units from kg/kg to g/kg
            if (gbx%use_reff) then
              re_matrix(i,:,:) = sghydro%Reff(:,pr,gbx%Nlevels:1:-1,i)*1.e6     ! Units from m to micron
            endif
         enddo  
      else  
         do i=1,gbx%Nhydro
            hm_matrix(i,:,:) = sghydro%mr_hydro(:,pr,:,i)*1000.0 ! Units from kg/kg to g/kg
            if (gbx%use_reff) then
              re_matrix(i,:,:) = sghydro%Reff(:,pr,:,i)*1.e6       ! Units from m to micron
            endif
         enddo
      endif  

      !   ----- call radar simulator -----
      if (pr == 1) then ! Compute gaseous attenuation for all profiles
         call radar_simulator(freq,k2,gbx%do_ray,gbx%use_gas_abs,gbx%use_mie_tables,gbx%mt, &    !  v0.2: mt changed to gbx%mt, roj
           gbx%Nhydro,gbx%hp,gbx%Npoints,gbx%Nlevels,gbx%nsizes,gbx%D, &                         !  v0.2: hp->gbx%hp, D->gbx%d, nsizes->gbx%nsizes, roj
           hgt_matrix,hm_matrix,re_matrix,p_matrix,t_matrix,rh_matrix, &
           Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe,g_to_vol_out=g_to_vol)
      else ! Use gaseous atteunuation for pr = 1
         call radar_simulator(freq,k2,gbx%do_ray,gbx%use_gas_abs,gbx%use_mie_tables,gbx%mt, &
           gbx%Nhydro,gbx%hp,gbx%Npoints,gbx%Nlevels,gbx%nsizes,gbx%D, &
           hgt_matrix,hm_matrix,re_matrix,p_matrix,t_matrix,rh_matrix, &
           Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe,g_to_vol_in=g_to_vol)
      endif
      ! ----- BEGIN output section -----
      ! spaceborne radar : from TOA to SURFACE
      if (gbx%surface_radar == 1) then
        z%Ze_tot(:,pr,:)=dBZe(:,:)
      else if (gbx%surface_radar == 0) then ! Spaceborne
        z%Ze_tot(:,pr,:)=dBZe(:,gbx%Nlevels:1:-1)
      endif

  enddo !pr
  
  ! Change undefined value to one defined in COSP
  where (z%Ze_tot == -999.0) z%Ze_tot = R_UNDEF
  
  deallocate(hgt_matrix,p_matrix,t_matrix,rh_matrix)
  deallocate(hm_matrix,re_matrix, &
      Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe)
  deallocate(g_to_vol)
 
  ! deallocate(mt_ttl,mt_tti)   !v0.2: roj feb 2008 can not be done here,
                                !these variables now part of gbx structure and dealocated later

END SUBROUTINE COSP_RADAR

END MODULE MOD_COSP_RADAR
