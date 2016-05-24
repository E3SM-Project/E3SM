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
    subroutine radar_simulator(hp,nprof,ngate,undef, &
        hgt_matrix,hm_matrix,re_matrix,Np_matrix, &
        p_matrix,t_matrix,rh_matrix, &
        Ze_non,Ze_ray,g_to_vol,a_to_vol,dBZe, &
        g_to_vol_in,g_to_vol_out)

        use m_mrgrnk
        use array_lib
        use math_lib
        use optics_lib
        use radar_simulator_types
        implicit none

        ! ----- INPUTS -----  
        type(class_param) :: hp

        integer, intent(in) :: nprof,ngate

        real undef
        real*8, dimension(nprof,ngate), intent(in) :: hgt_matrix, p_matrix, &
            t_matrix,rh_matrix
        real*8, dimension(hp%nhclass,nprof,ngate), intent(in) :: hm_matrix
        real*8, dimension(hp%nhclass,nprof,ngate), intent(inout) :: re_matrix
        real*8, dimension(hp%nhclass,nprof,ngate), intent(inout) :: Np_matrix

        ! ----- OUTPUTS -----
        real*8, dimension(nprof,ngate), intent(out) :: Ze_non,Ze_ray, &
            g_to_vol,dBZe,a_to_vol
        ! ----- OPTIONAL -----
        real*8, optional, dimension(nprof,ngate) :: &
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
  integer :: & 
  nsizes            ! num of discrete drop sizes

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
  re_matrix, &          ! effective radius (microns).   Optional. 0 ==> use Np_matrix or defaults
  Np_matrix         ! total number concentration (kg^-1).   Optional 0==> use defaults 

  integer, parameter :: one = 1
  ! logical :: hgt_reversed
  logical :: hgt_descending
  integer :: pr,i,j,k,unt,ngate

! ----- main program settings ------

  ! Inputs to Quickbeam
  allocate(hgt_matrix(gbx%Npoints,gbx%Nlevels),p_matrix(gbx%Npoints,gbx%Nlevels), &
           t_matrix(gbx%Npoints,gbx%Nlevels),rh_matrix(gbx%Npoints,gbx%Nlevels))
  allocate(hm_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))
  allocate(re_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))
  allocate(Np_matrix(gbx%Nhydro,gbx%Npoints,gbx%Nlevels))

  ! Outputs from Quickbeam
  allocate(Ze_non(gbx%Npoints,gbx%Nlevels))
  allocate(Ze_ray(gbx%Npoints,gbx%Nlevels))
  allocate(h_atten_to_vol(gbx%Npoints,gbx%Nlevels))
  allocate(g_atten_to_vol(gbx%Npoints,gbx%Nlevels))
  allocate(dBZe(gbx%Npoints,gbx%Nlevels))

  ! Optional argument. It is computed and returned in the first call to
  ! radar_simulator, and passed as input in the rest
  allocate(g_to_vol(gbx%Npoints,gbx%Nlevels))

  ! Even if there is no unit conversion, they are needed for type conversion
  p_matrix   = gbx%p/100.0     ! From Pa to hPa
  hgt_matrix = gbx%zlev/1000.0 ! From m to km
  t_matrix   = gbx%T
  rh_matrix  = gbx%q
  re_matrix  = 0.0


  ! set flag denoting position of radar relative to hgt_matrix orientation
	  ngate = size(hgt_matrix,2)

	  hgt_descending = hgt_matrix(1,1) > hgt_matrix(1,ngate)

	  if ( &
	     (gbx%surface_radar == 1 .and. hgt_descending) .or.  &
	     (gbx%surface_radar == 0 .and. (.not. hgt_descending)) &
	     ) &
	  then
	    gbx%hp%radar_at_layer_one = .false.
	  else
	    gbx%hp%radar_at_layer_one = .true.
	  endif

  ! ----- loop over subcolumns -----
  do pr=1,sgx%Ncolumns

      !  NOTE:
      !  atmospheric profiles are the same within the same gridbox
      !  only hydrometeor profiles will be different for each subgridbox

         do i=1,gbx%Nhydro
            hm_matrix(i,:,:) = sghydro%mr_hydro(:,pr,:,i)*1000.0 ! Units from kg/kg to g/kg
            if (gbx%use_reff) then
              re_matrix(i,:,:) = sghydro%Reff(:,pr,:,i)*1.e6       ! Units from m to micron
              Np_matrix(i,:,:) = sghydro%Np(:,pr,:,i)              ! Units [#/kg]
            endif
         enddo

      !   ----- call radar simulator -----
      if (pr == 1) then ! Compute gaseous attenuation for all profiles
         call radar_simulator(gbx%hp,gbx%Npoints,gbx%Nlevels,R_UNDEF, &
           hgt_matrix,hm_matrix,re_matrix,Np_matrix, &
           p_matrix,t_matrix,rh_matrix, &
           Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe,g_to_vol_out=g_to_vol)
      else ! Use gaseous atteunuation for pr = 1
         call radar_simulator(gbx%hp,gbx%Npoints,gbx%Nlevels,R_UNDEF, &
           hgt_matrix,hm_matrix,re_matrix,Np_matrix, &
           p_matrix,t_matrix,rh_matrix, &
           Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe,g_to_vol_in=g_to_vol)
      endif

      ! store caluculated dBZe values for later output/processing
      z%Ze_tot(:,pr,:)=dBZe(:,:)
  enddo !pr

  deallocate(hgt_matrix,p_matrix,t_matrix,rh_matrix)
  deallocate(hm_matrix,re_matrix, &
      Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe)
  deallocate(g_to_vol)
END SUBROUTINE COSP_RADAR

END MODULE MOD_COSP_RADAR
