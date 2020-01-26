! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! 11/2005: John Haynes - Created
! 09/2006  placed into subroutine form (Roger Marchand,JMH)
! 08/2007  added equivalent volume spheres, Z and N scalling most distrubtion types (Roger Marchand)
! 01/2008  'Do while' to determine if hydrometeor(s) present in volume
!           changed for vectorization purposes (A. Bodas-Salcedo)
!
! 07/2010  V3.0 ... Modified to load or save scale factors to disk as a Look-Up Table (LUT)
!  ... All hydrometeor and radar simulator properties now included in hp structure
!  ... hp structure should be initialized by call to radar_simulator_init prior 
!  ... to calling this subroutine.  
!     Also ... Support of Morrison 2-moment style microphyscis (Np_matrix) added 
!  ... Changes implement by Roj Marchand following work by Laura Fowler
!
!   10/2011  Modified ngate loop to go in either direction depending on flag 
!     hp%radar_at_layer_one.  This affects the direction in which attenuation is summed.
!
!     Also removed called to AVINT for gas and hydrometeor attenuation and replaced with simple
!     summation. (Roger Marchand)
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module quickbeam
  USE COSP_KINDS,           ONLY: wp
  USE MOD_COSP_CONFIG,      ONLY: DBZE_BINS,DBZE_MIN,DBZE_MAX,CFAD_ZE_MIN,CFAD_ZE_WIDTH, &
                                  R_UNDEF,cloudsat_histRef,use_vgrid,vgrid_zl,vgrid_zu
  USE MOD_COSP_STATS,       ONLY: COSP_LIDAR_ONLY_CLOUD,hist1D,COSP_CHANGE_VERTICAL_GRID
  implicit none

  integer,parameter :: &
       maxhclass     = 20,  & ! Qucikbeam maximum number of hydrometeor classes.
       nRe_types     = 550, & ! Quickbeam maximum number or Re size bins allowed in N and Z_scaled look up table.
       nd            = 85,  & ! Qucikbeam number of discrete particles used in construction DSDs.
       mt_ntt        = 39,  & ! Quickbeam number of temperatures in mie LUT.
       Re_BIN_LENGTH = 10,  & ! Quickbeam minimum Re interval in scale LUTs  
       Re_MAX_BIN    = 250    ! Quickbeam maximum Re interval in scale LUTs
  real(wp),parameter :: &
       dmin          = 0.1, & ! Quickbeam minimum size of discrete particle
       dmax          = 10000. ! Quickbeam maximum size of discrete particle
  
  !djs logical :: radar_at_layer_one   ! If true radar is assume to be at the edge 
                                  ! of the first layer, if the first layer is the
                                  ! surface than a ground-based radar.   If the
                                  ! first layer is the top-of-atmosphere, then
                                  ! a space borne radar.

  ! ##############################################################################################
  type radar_cfg
     ! Radar properties
     real(wp) :: freq,k2
     integer  :: nhclass               ! Number of hydrometeor classes in use
     integer  :: use_gas_abs, do_ray
     logical  :: radar_at_layer_one    ! If true radar is assume to be at the edge 
                                       ! of the first layer, if the first layer is the
                                       ! surface than a ground-based radar.   If the
                                       ! first layer is the top-of-atmosphere, then
                                       ! a space borne radar.
     
     ! Variables used to store Z scale factors
     character(len=240)                             :: scale_LUT_file_name
     logical                                        :: load_scale_LUTs, update_scale_LUTs
     logical, dimension(maxhclass,nRe_types)        :: N_scale_flag
     logical, dimension(maxhclass,mt_ntt,nRe_types) :: Z_scale_flag,Z_scale_added_flag
     real(wp),dimension(maxhclass,mt_ntt,nRe_types) :: Ze_scaled,Zr_scaled,kr_scaled
     real(wp),dimension(maxhclass,nd,nRe_types)     :: fc, rho_eff
     real(wp),dimension(Re_MAX_BIN)                 :: base_list,step_list

  end type radar_cfg

contains
  ! ######################################################################################
  ! SUBROUTINE quickbeam_subcolumn
  ! ######################################################################################
  !subroutine quickbeam_subcolumn(rcfg,nprof,ngate,hgt_matrix,z_vol,kr_vol,g_vol,&
  !                               a_to_vol,g_to_vol,dBZe,Ze_non,Ze_ray)
  subroutine quickbeam_subcolumn(rcfg,nprof,ngate,hgt_matrix,z_vol,kr_vol,g_vol,dBZe)

    ! INPUTS
    type(radar_cfg),intent(inout) :: &
         rcfg             ! Derived type for radar simulator setup
    integer,intent(in) :: &
         nprof,         & ! Number of hydrometeor profiles
         ngate            ! Number of vertical layers
    real(wp),intent(in),dimension(nprof,ngate) :: &
         hgt_matrix,    & ! Height of hydrometeors (km)
         z_vol,         & ! Effective reflectivity factor (mm^6/m^3)
         kr_vol,        & ! Attenuation coefficient hydro (dB/km)
         g_vol            ! Attenuation coefficient gases (dB/km)
    
    ! OUTPUTS
    real(wp), intent(out),dimension(nprof,ngate) :: &
!         Ze_non,        & ! Radar reflectivity without attenuation (dBZ)
!         Ze_ray,        & ! Rayleigh reflectivity (dBZ)
!         g_to_vol,      & ! Gaseous atteunation, radar to vol (dB)
!         a_to_vol,      & ! Hydromets attenuation, radar to vol (dB)
         dBZe             ! Effective radar reflectivity factor (dBZ)

    ! LOCAL VARIABLES
    integer :: k,pr,start_gate,end_gate,d_gate
    real(wp),dimension(nprof,ngate) :: &
         Ze_non,        & ! Radar reflectivity without attenuation (dBZ)
         Ze_ray,        & ! Rayleigh reflectivity (dBZ)
         g_to_vol,      & ! Gaseous atteunation, radar to vol (dB)
         a_to_vol,      & ! Hydromets attenuation, radar to vol (dB) 
         z_ray            ! Reflectivity factor, Rayleigh only (mm^6/m^3)

    ! Load scaling matricies from disk -- but only the first time this subroutine is called
    if(rcfg%load_scale_LUTs) then
       call load_scale_LUTs(rcfg)
       rcfg%load_scale_LUTs=.false.
       rcfg%Z_scale_added_flag = .false. ! will be set true if scaling Look Up Tables are modified during run
    endif

    ! Initialization
    g_to_vol = 0._wp
    a_to_vol = 0._wp

    ! Loop over each range gate (ngate) ... starting with layer closest to the radar !
    if(rcfg%radar_at_layer_one) then
       start_gate = 1
       end_gate   = ngate
       d_gate     = 1
    else
       start_gate = ngate
       end_gate   = 1
       d_gate     = -1
    endif
    do k=start_gate,end_gate,d_gate
       ! Loop over each profile (nprof)
       do pr=1,nprof
          ! Attenuation due to hydrometeors between radar and volume
          
          ! NOTE old scheme integrates attenuation only for the layers ABOVE
          ! the current layer ... i.e. 1 to k-1 rather than 1 to k ...
          ! which may be a problem.   ROJ
          ! in the new scheme I assign half the attenuation to the current layer
          if(d_gate==1) then
             ! dheight calcuations assumes hgt_matrix points are the cell mid-points.
             if (k>2) then
                ! add to previous value to half of above layer + half of current layer
                a_to_vol(pr,k)=  a_to_vol(pr,k-1) + &
                     (kr_vol(pr,k-1)+kr_vol(pr,k))*(hgt_matrix(pr,k-1)-hgt_matrix(pr,k))
             else
                a_to_vol(pr,k)=  kr_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k+1))
             endif
          else   ! d_gate==-1
             if(k<ngate) then
                ! Add to previous value half of above layer + half of current layer
                a_to_vol(pr,k) = a_to_vol(pr,k+1) + &
                     (kr_vol(pr,k+1)+kr_vol(pr,k))*(hgt_matrix(pr,k+1)-hgt_matrix(pr,k))
             else
                a_to_vol(pr,k)= kr_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k-1))
             endif
          endif
          
          ! Attenuation due to gaseous absorption between radar and volume
          if ((rcfg%use_gas_abs == 1) .or. (rcfg%use_gas_abs == 2 .and. pr .eq. 1)) then
             if (d_gate==1) then
                if (k>1) then
                   ! Add to previous value to half of above layer + half of current layer
                   g_to_vol(pr,k) =  g_to_vol(pr,k-1) + &
                        0.5*(g_vol(pr,k-1)+g_vol(pr,k))*(hgt_matrix(pr,k-1)-hgt_matrix(pr,k))
                else
                   g_to_vol(pr,k)=  0.5_wp*g_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k+1))
                endif
             else   ! d_gate==-1
                if (k<ngate) then
                   ! Add to previous value to half of above layer + half of current layer
                   g_to_vol(pr,k) = g_to_vol(pr,k+1) + &
                        0.5_wp*(g_vol(pr,k+1)+g_vol(pr,k))*(hgt_matrix(pr,k+1)-hgt_matrix(pr,k))
                else
                   g_to_vol(pr,k)= 0.5_wp*g_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k-1))
                endif
             endif
          elseif(rcfg%use_gas_abs == 2) then
             ! Using value calculated for the first column
             g_to_vol(pr,k) = g_to_vol(1,k)
          elseif (rcfg%use_gas_abs == 0) then
             g_to_vol(pr,k) = 0._wp
          endif
       enddo   ! End loop over pr (profile)
    enddo ! End loop of k (range gate)
    
    ! Compute Rayleigh reflectivity, and full, attenuated reflectivity
    if(rcfg%do_ray == 1) then
       where(z_ray(1:nprof,1:ngate) > 0._wp)
          Ze_ray(1:nprof,1:ngate) = 10._wp*log10(z_ray(1:nprof,1:ngate))
       elsewhere
          Ze_Ray(1:nprof,1:ngate) = 0._wp
       endwhere
!djs       Ze_ray(1:nprof,1:ngate) = merge(10._wp*log10(z_ray(1:nprof,1:ngate)), 1._wp*R_UNDEF, z_ray(1:nprof,1:ngate) > 0._wp)
    else 
      Ze_ray(1:nprof,1:ngate) = R_UNDEF
    end if

    where(z_vol(1:nprof,1:ngate) > 0._wp) 
      Ze_non(1:nprof,1:ngate) = 10._wp*log10(z_vol(1:nprof,1:ngate))
      dBZe(1:nprof,1:ngate) = Ze_non(1:nprof,1:ngate)-a_to_vol(1:nprof,1:ngate)-g_to_vol(1:nprof,1:ngate)
    elsewhere
      dBZe(1:nprof,1:ngate) = R_UNDEF
      Ze_non(1:nprof,1:ngate) = R_UNDEF
    end where 

    ! Save any updates made 
    if (rcfg%update_scale_LUTs) call save_scale_LUTs(rcfg)
 
  end subroutine quickbeam_subcolumn
  ! ######################################################################################
  ! SUBROUTINE quickbeam_column
  ! ######################################################################################
  subroutine quickbeam_column(npoints,ncolumns,nlevels,llm,Ze_tot,zlev,zlev_half,cfad_ze)
    ! Inputs
    integer,intent(in) :: &
         npoints,    & ! Number of horizontal grid points
         ncolumns,   & ! Number of subcolumns
         nlevels,    & ! Number of vertical layers in OLD grid
         llm           ! NUmber of vertical layers in NEW grid
    real(wp),intent(in),dimension(npoints,ncolumns,Nlevels) :: &
         Ze_tot        ! 
    real(wp),intent(in),dimension(npoints,Nlevels) :: &
         zlev          ! Model full levels
    real(wp),intent(in),dimension(npoints,Nlevels+1) :: &
         zlev_half     ! Model half levels
         
    ! Outputs
    real(wp),intent(inout),dimension(npoints,DBZE_BINS,llm) :: &
         cfad_ze    !

    ! Local variables
    integer :: i,j
    real(wp),dimension(npoints,ncolumns,llm) :: ze_totFlip
    
    if (use_vgrid) then
       ! Regrid in the vertical
       call cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,zlev(:,nlevels:1:-1),&
            zlev_half(:,nlevels:1:-1),Ze_tot(:,:,nlevels:1:-1),llm,vgrid_zl(llm:1:-1),&
            vgrid_zu(llm:1:-1),Ze_totFlip(:,:,llm:1:-1),log_units=.true.)

       ! Effective reflectivity histogram
       do i=1,Npoints
          do j=1,llm
             cfad_ze(i,:,j) = hist1D(Ncolumns,Ze_totFlip(i,:,j),DBZE_BINS,cloudsat_histRef)
          enddo
       enddo
       where(cfad_ze .ne. R_UNDEF) cfad_ze = cfad_ze/Ncolumns

    else
       ! Effective reflectivity histogram
       do i=1,Npoints
          do j=1,llm
             cfad_ze(i,:,j) = hist1D(Ncolumns,Ze_tot(i,:,j),DBZE_BINS,cloudsat_histRef)
          enddo
       enddo
       where(cfad_ze .ne. R_UNDEF) cfad_ze = cfad_ze/Ncolumns
    endif   

  end subroutine quickbeam_column
  ! ##############################################################################################
  ! ##############################################################################################

  
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine load_scale_LUTs(rcfg)
    
    type(radar_cfg), intent(inout) :: rcfg
    logical                        :: LUT_file_exists
    integer                        :: i,j,k,ind
    
    ! Load scale LUT from file 
    inquire(file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
         exist=LUT_file_exists)
    
    if(.not.LUT_file_exists) then  
       write(*,*) '*************************************************'
       write(*,*) 'Warning: Could NOT FIND radar LUT file: ', &
            trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'        
       write(*,*) 'Will calculated LUT values as needed'
       write(*,*) '*************************************************'
       return
    else
       OPEN(unit=12,file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
            form='unformatted', &
            err= 89, &
            access='DIRECT',&
            recl=28)
       write(*,*) 'Loading radar LUT file: ', &
            trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
       
       do i=1,maxhclass
          do j=1,mt_ntt
             do k=1,nRe_types
                ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
                read(12,rec=ind) rcfg%Z_scale_flag(i,j,k), &
                     rcfg%Ze_scaled(i,j,k), &
                     rcfg%Zr_scaled(i,j,k), &
                     rcfg%kr_scaled(i,j,k)
             enddo
          enddo
       enddo
       close(unit=12)
       return 
    endif
    
89  write(*,*) 'Error: Found but could NOT READ radar LUT file: ', &
         trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
  end subroutine load_scale_LUTs
  
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine save_scale_LUTs(rcfg)
    type(radar_cfg), intent(inout) :: rcfg
    logical                        :: LUT_file_exists
    integer                        :: i,j,k,ind
    
    inquire(file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
         exist=LUT_file_exists)
    
    OPEN(unit=12,file=trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
         form='unformatted',err= 99,access='DIRECT',recl=28)
    
    write(*,*) 'Creating or Updating radar LUT file: ', &
         trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
    do i=1,maxhclass
       do j=1,mt_ntt
          do k=1,nRe_types
             ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
             if(.not.LUT_file_exists .or. rcfg%Z_scale_added_flag(i,j,k)) then
                rcfg%Z_scale_added_flag(i,j,k)=.false.
                write(12,rec=ind) rcfg%Z_scale_flag(i,j,k), &
                     rcfg%Ze_scaled(i,j,k), &
                     rcfg%Zr_scaled(i,j,k), &
                     rcfg%kr_scaled(i,j,k)
             endif
          enddo
       enddo
    enddo
    close(unit=12)
    return 
    
99  write(*,*) 'Error: Unable to create/update radar LUT file: ', &
         trim(rcfg%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    return  
    
  end subroutine save_scale_LUTs
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine quickbeam_init()

    
  end subroutine quickBeam_init
  ! ##############################################################################################
  ! ##############################################################################################


end module quickbeam


