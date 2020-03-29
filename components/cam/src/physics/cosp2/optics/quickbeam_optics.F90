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
! History:
! May 2015:  Dustin Swales - Initial version
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
module mod_quickbeam_optics
  USE COSP_KINDS,          ONLY: wp,dp
  USE array_lib,           ONLY: infind
  USE math_lib,            ONLY: path_integral,avint,gamma
  USE optics_lib,          ONLY: m_wat,m_ice,MieInt
  USE cosp_math_constants, ONLY: pi
  USE cosp_phys_constants, ONLY: rhoice  
  use quickbeam,           ONLY: radar_cfg,dmin,dmax,Re_BIN_LENGTH,  &
                                 Re_MAX_BIN,nRe_types,nd,maxhclass,save_scale_LUTs
  use mod_cosp_config,     ONLY: N_HYDRO
  use mod_cosp_error,      ONLY: errorMessage
  implicit none

  ! Derived type for particle size distribution
  TYPE size_distribution
     real(wp),dimension(maxhclass) :: p1,p2,p3,dmin,dmax,apm,bpm,rho
     integer, dimension(maxhclass) :: dtype,phase
  END TYPE size_distribution
  
  ! Parameters
  integer,parameter ::   & !
       cnt_liq       = 19, & ! Liquid temperature count
       cnt_ice       = 20    ! Lce temperature count
  
  ! Initialization variables
  real(wp),dimension(cnt_ice) :: mt_tti 
  real(wp),dimension(cnt_liq) :: mt_ttl
  real(wp),dimension(nd)      :: D
  
contains
  ! ######################################################################################
  ! SUBROUTINE quickbeam_optics_init
  ! ######################################################################################
  subroutine quickbeam_optics_init()
    integer :: j
    
    mt_tti = (/ ((j-1)*5-90 + 273.15, j = 1, cnt_ice) /) 
    mt_ttl = (/ ((j-1)*5-60 + 273.15, j = 1, cnt_liq) /)
    D(1) = dmin
    do j=2,nd
       D(j) = D(j-1)*exp((log(dmax)-log(dmin))/(nd-1))
    enddo
  end subroutine quickbeam_optics_init
  
  ! ######################################################################################
  ! SUBROUTINE QUICKBEAM_OPTICS
  ! ######################################################################################
  subroutine quickbeam_optics(sd, rcfg, nprof, ngate, undef, hm_matrix, re_matrix,       &
       Np_matrix, p_matrix, t_matrix, sh_matrix,z_vol,kr_vol)
    
    ! INPUTS
    type(size_distribution),intent(inout) :: &
         sd               !
    type(radar_cfg),intent(inout) :: &
         rcfg             !
    integer,intent(in) :: &
         nprof,         & ! Number of hydrometeor profiles
         ngate            ! Number of vertical layers
    real(wp),intent(in) :: &
         undef            ! Missing data value
    real(wp),intent(in),dimension(nprof,ngate) :: &
         p_matrix,      & ! Pressure profile (hPa)
         t_matrix,      & ! Temperature profile (K)
         sh_matrix        ! Specific humidity profile (%) -- only needed if gaseous aborption calculated.
    real(wp),intent(in),dimension(nprof,ngate,rcfg%nhclass) :: &
         re_matrix,     & ! Table of hydrometeor effective radii.       0 ==> use defaults. (units=microns)
         hm_matrix        ! Table of hydrometeor mixing ratios (g/kg)
    real(wp),intent(inout),dimension(nprof,ngate,rcfg%nhclass) :: &
         Np_matrix        ! Table of hydrometeor number concentration.  0 ==> use defaults. (units = 1/kg)

    ! OUTPUTS
    real(wp),intent(out), dimension(nprof, ngate) :: &
         z_vol,         & ! Effective reflectivity factor (mm^6/m^3)
         kr_vol           ! Attenuation coefficient hydro (dB/km)

    ! INTERNAL VARIABLES   
    integer :: &
         phase, ns,tp,j,k,pr,itt,iRe_type,n 
    logical :: &
         hydro
    real(wp) :: &
         t_kelvin,Re_internal
    real(wp) :: &
         rho_a,kr,ze,zr,scale_factor,Re,Np,base,step 

    real(wp),dimension(:),allocatable :: &
         Deq,     & ! Discrete drop sizes (um)
         Ni,      & ! Discrete concentrations (cm^-3 um^-1)
         rhoi,    & ! Discrete densities (kg m^-3)
         xxa,     & !
         Di         ! Discrete drop sizes (um)
    
    real(wp), dimension(nprof, ngate) :: &
         z_ray      ! Reflectivity factor, Rayleigh only (mm^6/m^3)
    
    ! PARAMETERS    
    logical, parameter ::       & !
         DO_LUT_TEST = .false., & !
         DO_NP_TEST  = .false.    !
    real(wp), parameter :: &
         one_third   = 1._wp/3._wp    !

    ! Initialization
    z_vol    = 0._wp
    z_ray    = 0._wp
    kr_vol   = 0._wp

    do k=1,ngate       ! Loop over each profile (nprof)
       do pr=1,nprof

          ! Determine if hydrometeor(s) present in volume
          hydro = .false.
          do j=1,rcfg%nhclass
             if ((hm_matrix(pr,k,j) > 1E-12) .and. (sd%dtype(j) > 0)) then
                hydro = .true.
                exit
             endif
          enddo

          t_kelvin = t_matrix(pr,k)
          ! If there is hydrometeor in the volume
          if (hydro) then
             rho_a = (p_matrix(pr,k))/(287._wp*(t_kelvin))
             
             ! Loop over hydrometeor type
             do tp=1,rcfg%nhclass
                Re_internal = re_matrix(pr,k,tp)

                if (hm_matrix(pr,k,tp) <= 1E-12) cycle
                
                ! Index into temperature dimension of scaling tables
                !   These tables have regular steps -- exploit this and abandon infind
                phase = sd%phase(tp)
                if (phase==0) then
                   itt = infind(mt_ttl,t_kelvin)
                else
                   itt = infind(mt_tti,t_kelvin) 
                endif
                
                ! Compute effective radius from number concentration and distribution parameters
                if (Re_internal .eq. 0) then
                   call calc_Re(hm_matrix(pr,k,tp),Np_matrix(pr,k,tp),rho_a, &
                        sd%dtype(tp),sd%apm(tp),sd%bpm(tp),sd%rho(tp),sd%p1(tp),sd%p2(tp),sd%p3(tp),Re)
                   Re_internal=Re
                   !re_matrix(pr,k,tp)=Re
                else
                   if (Np_matrix(pr,k,tp) > 0) then
                      call errorMessage('WARNING(optics/quickbeam_optics.f90): '//&
                         'Re and Np set for the same volume & hydrometeor type.  Np is being ignored.')
                   endif
                   Re = Re_internal
                   !Re = re_matrix(pr,k,tp)
                endif
                
                ! Index into particle size dimension of scaling tables 
                iRe_type=1
                if(Re.gt.0) then
                   ! Determine index in to scale LUT
                   ! Distance between Re points (defined by "base" and "step") for
                   ! each interval of size Re_BIN_LENGTH
                   ! Integer asignment, avoids calling floor intrinsic
                   n=Re/Re_BIN_LENGTH
                   if (n>=Re_MAX_BIN) n=Re_MAX_BIN-1
                   step = rcfg%step_list(n+1)
                   base = rcfg%base_list(n+1)
                   iRe_type=Re/step
                   if (iRe_type.lt.1) iRe_type=1
                   Re=step*(iRe_type+0.5_wp)    ! set value of Re to closest value allowed in LUT.
                   iRe_type=iRe_type+base-int(n*Re_BIN_LENGTH/step)
                   
                   ! Make sure iRe_type is within bounds
                   if (iRe_type.ge.nRe_types) then
                      !write(*,*) 'Warning: size of Re exceed value permitted ', &
                      !            'in Look-Up Table (LUT).  Will calculate. '
                      ! No scaling allowed
                      iRe_type=nRe_types
                      rcfg%Z_scale_flag(tp,itt,iRe_type)=.false.
                   else
                      ! Set value in re_matrix to closest values in LUT
                      if (.not. DO_LUT_TEST) re_internal=Re
                      !if (.not. DO_LUT_TEST) re_matrix(pr,k,tp)=Re
                   endif
                endif
                
                ! Use Ze_scaled, Zr_scaled, and kr_scaled ... if know them
                ! if not we will calculate Ze, Zr, and Kr from the distribution parameters
!                if( rcfg%Z_scale_flag(tp,itt,iRe_type) .and. .not. DO_LUT_TEST)  then
!                   ! can use z scaling
!                   scale_factor=rho_a*hm_matrix(pr,k,tp)
!                   zr = rcfg%Zr_scaled(tp,itt,iRe_type) * scale_factor
!                   ze = rcfg%Ze_scaled(tp,itt,iRe_type) * scale_factor
!                   kr = rcfg%kr_scaled(tp,itt,iRe_type) * scale_factor
!                else
                if( (.not. rcfg%Z_scale_flag(tp,itt,iRe_type)) .or. DO_LUT_TEST)  then
                   ! Create a discrete distribution of hydrometeors within volume
                   select case(sd%dtype(tp))
                   case(4)
                      ns = 1
                      allocate(Di(ns),Ni(ns),rhoi(ns),xxa(ns),Deq(ns))
                      Di = sd%p1(tp)
                      Ni = 0._wp
                   case default
                      ns = nd   ! constant defined in simulator/quickbeam.f90
                      allocate(Di(ns),Ni(ns),rhoi(ns),xxa(ns),Deq(ns))
                      Di = D
                      Ni = 0._wp
                   end select
                   call dsd(hm_matrix(pr,k,tp),re_internal,Np_matrix(pr,k,tp), &
                        Di,Ni,ns,sd%dtype(tp),rho_a,t_kelvin, &
                        sd%dmin(tp),sd%dmax(tp),sd%apm(tp),sd%bpm(tp), &
                        sd%rho(tp),sd%p1(tp),sd%p2(tp),sd%p3(tp))
                   
                   ! Calculate particle density
                   if (phase == 1) then
                      if (sd%rho(tp) < 0) then
                         ! Use equivalent volume spheres.
                         rcfg%rho_eff(tp,1:ns,iRe_type) = rhoice ! solid ice == equivalent volume approach
                         Deq = ( ( 6/pi*sd%apm(tp)/rhoice) ** one_third ) * ( (Di*1E-6) ** (sd%bpm(tp)/3._wp) )  * 1E6
                         ! alternative is to comment out above two lines and use the following block
                         ! MG Mie approach - adjust density of sphere with D = D_characteristic to match particle density
                         !
                         ! rcfg%rho_eff(tp,1:ns,iRe_type) = (6/pi)*sd%apm(tp)*(Di*1E-6)**(sd%bpm(tp)-3)   !MG Mie approach
                         
                         ! as the particle size gets small it is possible that the mass to size relationship of 
                         ! (given by power law in hclass.data) can produce impossible results 
                         ! where the mass is larger than a solid sphere of ice.  
                         ! This loop ensures that no ice particle can have more mass/density larger than an ice sphere.
                         ! do i=1,ns
                         ! if(rcfg%rho_eff(tp,i,iRe_type) > 917 ) then
                         ! rcfg%rho_eff(tp,i,iRe_type) = 917
                         ! endif
                         ! enddo
                      else
                         ! Equivalent volume sphere (solid ice rhoice=917 kg/m^3).
                         rcfg%rho_eff(tp,1:ns,iRe_type) = rhoice
                         Deq=Di * ((sd%rho(tp)/rhoice)**one_third)
                         ! alternative ... coment out above two lines and use the following for MG-Mie
                         ! rcfg%rho_eff(tp,1:ns,iRe_type) = sd%rho(tp)   !MG Mie approach
                      endif
                   else
                      ! I assume here that water phase droplets are spheres.
                      ! sd%rho should be ~ 1000  or sd%apm=524 .and. sd%bpm=3
                      Deq = Di
                   endif

                   ! Calculate effective reflectivity factor of volume
                   ! xxa are unused (Mie scattering and extinction efficiencies)
                   xxa(1:ns) = -9.9_wp
                   rhoi = rcfg%rho_eff(tp,1:ns,iRe_type)
                   call zeff(rcfg%freq,Deq,Ni,ns,rcfg%k2,t_kelvin,phase,rcfg%do_ray, &
                        ze,zr,kr,xxa,xxa,rhoi)

                   ! Test compares total number concentration with sum of discrete samples 
                   ! The second test, below, compares ab initio and "scaled" computations 
                   !    of reflectivity
                   !  These should get broken out as a unit test that gets called on 
                   !    data. That routine could write to std out. 
                   
                   ! Test code ... compare Np value input to routine with sum of DSD
                   ! NOTE: if .not. DO_LUT_TEST, then you are checking the LUT approximation 
                   ! not just the DSD representation given by Ni
                   if(Np_matrix(pr,k,tp)>0 .and. DO_NP_TEST ) then
                      Np = path_integral(Ni,Di,1,ns-1)/rho_a*1.E6_wp
                      ! Note: Representation is not great or small Re < 2 
                      if( (Np_matrix(pr,k,tp)-Np)/Np_matrix(pr,k,tp)>0.1 ) then
                         call errorMessage('ERROR(optics/quickbeam_optics.f90): Error: Np input does not match sum(N)')
                      endif
                   endif

                   ! Clean up space
                   deallocate(Di,Ni,rhoi,xxa,Deq)

                   ! LUT test code
                   ! This segment of code compares full calculation to scaling result
                   if ( rcfg%Z_scale_flag(tp,itt,iRe_type) .and. DO_LUT_TEST )  then
                      scale_factor=rho_a*hm_matrix(pr,k,tp)
                      ! if more than 2 dBZe difference print error message/parameters.
                      if ( abs(10*log10(ze) - 10*log10(rcfg%Ze_scaled(tp,itt,iRe_type) * &
                           scale_factor)) > 2 ) then
                         call errorMessage('ERROR(optics/quickbeam_optics.f90): ERROR: Roj Error?')
                      endif
                   endif
                else
                   ! Use z scaling
                   scale_factor=rho_a*hm_matrix(pr,k,tp)
                   zr = rcfg%Zr_scaled(tp,itt,iRe_type) * scale_factor
                   ze = rcfg%Ze_scaled(tp,itt,iRe_type) * scale_factor
                   kr = rcfg%kr_scaled(tp,itt,iRe_type) * scale_factor
                endif  ! end z_scaling
                
                kr_vol(pr,k) = kr_vol(pr,k) + kr
                z_vol(pr,k)  = z_vol(pr,k)  + ze
                z_ray(pr,k)  = z_ray(pr,k)  + zr
                
                ! Construct Ze_scaled, Zr_scaled, and kr_scaled ... if we can
                if ( .not. rcfg%Z_scale_flag(tp,itt,iRe_type) ) then
                   if (iRe_type>1) then
                      scale_factor=rho_a*hm_matrix(pr,k,tp)
                      rcfg%Ze_scaled(tp,itt,iRe_type) = ze/ scale_factor
                      rcfg%Zr_scaled(tp,itt,iRe_type) = zr/ scale_factor
                      rcfg%kr_scaled(tp,itt,iRe_type) = kr/ scale_factor
                      rcfg%Z_scale_flag(tp,itt,iRe_type) = .true.
                      rcfg%Z_scale_added_flag(tp,itt,iRe_type)=.true.
                   endif
                endif
             enddo   ! end loop of tp (hydrometeor type)
          endif
       enddo
    enddo

    where(kr_vol(:,:) <= EPSILON(kr_vol)) 
       ! Volume is hydrometeor-free	
       !z_vol(:,:)  = undef
       z_ray(:,:)  = undef
    end where

    ! Save any updates made 
    if (rcfg%update_scale_LUTs) call save_scale_LUTs(rcfg)

  end subroutine quickbeam_optics
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine calc_Re(Q,Np,rho_a,dtype,apm,bpm,rho_c,p1,p2,p3,Re)  
    ! ##############################################################################################
    ! Purpose:
    !   Calculates Effective Radius (1/2 distribution 3rd moment / 2nd moment). 
    !
    !   For some distribution types, the total number concentration (per kg), Np
    !   may be optionally specified.   Should be set to zero, otherwise.
    !
    !   Roj Marchand July 2010
    !
    ! Inputs:
    !
    !   [Q]        hydrometeor mixing ratio (g/kg)  ! not needed for some distribution types
    !   [Np]       Optional Total number concentration (per kg).  0 = use defaults (p1, p2, p3)
    !   [rho_a]    ambient air density (kg m^-3)   
    !
    !   Distribution parameters as per quickbeam documentation.
    !   [dtype]    distribution type
    !   [apm]      a parameter for mass (kg m^[-bpm])
    !   [bmp]      b params for mass 
    !   [p1],[p2],[p3]  distribution parameters
    !
    ! Outputs:
    !   [Re]       Effective radius, 1/2 the 3rd moment/2nd moment (um)
    !
    ! Created:
    !   July 2010  Roj Marchand
    ! Modified:
    !   12/18/14  Dustin Swales: Define type REALs as double precision (dustin.swales@noaa.gov)
    !
    ! ##############################################################################################
    ! ##############################################################################################
    
    ! Inputs
    real(wp), intent(in)    :: Q,Np,rho_a,rho_c,p1,p2,p3
    integer,  intent(in)    :: dtype
    real(wp), intent(inout) :: apm,bpm  
    
    ! Outputs
    real(wp), intent(out) :: Re
    
    ! Internal
    integer  :: local_dtype
    real(wp) :: local_p3,local_Np,tmp1,tmp2
    real(wp) :: N0,D0,vu,dm,ld,rg,log_sigma_g        ! gamma, exponential variables
    
    
    ! If density is constant, set equivalent values for apm and bpm
    if ((rho_c > 0) .and. (apm < 0)) then
       apm = (pi/6)*rho_c
       bpm = 3._wp
    endif
    
    ! Exponential is same as modified gamma with vu =1
    ! if Np is specified then we will just treat as modified gamma
    if(dtype .eq. 2 .and. Np .gt. 0) then
       local_dtype = 1
       local_p3    = 1
    else
       local_dtype = dtype
       local_p3    = p3
    endif
    select case(local_dtype)
       
       ! ---------------------------------------------------------!
       ! Modified gamma                                           !
       ! Np = total number concentration (1/kg) = Nt / rho_a      !
       ! D0 = characteristic diameter (um)                        !
       ! dm = mean diameter (um) - first moment over zeroth moment!
       ! vu = distribution width parameter                        !
       ! ---------------------------------------------------------!
    case(1)  
       
       if( abs(local_p3+2) < 1E-8) then
          if(Np>1E-30) then
             ! Morrison scheme with Martin 1994 shape parameter (NOTE: vu = pc +1)
             ! fixed Roj. Dec. 2010 -- after comment by S. Mcfarlane
             vu = (1/(0.2714_wp + 0.00057145_wp*Np*rho_a*1E-6))**2 ! units of Nt = Np*rhoa = #/cm^3
          else
             call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
                'Must specify a value for Np in each volume with Morrison/Martin Scheme.')
             return
          endif
       elseif (abs(local_p3+1) > 1E-8) then
          ! vu is fixed in hp structure  
          vu = local_p3 
       else
          ! vu isn't specified
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
             'Must specify a value for vu for Modified Gamma distribution')
          return
       endif
       
       if( Np.eq.0 .and. p2+1 > 1E-8) then     ! use default value for MEAN diameter as first default  
          dm = p2             ! by definition, should have units of microns
          D0 = gamma(vu)/gamma(vu+1)*dm
       else   ! use value of Np
          if(Np.eq.0) then
             if( abs(p1+1) > 1E-8 ) then  !   use default number concentration   
                local_Np = p1 ! total number concentration / pa --- units kg^-1
             else
                call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
                   'Must specify Np or default value (p1=Dm [um] or p2=Np [1/kg]) for Modified Gamma distribution')
                return
             endif
          else
             local_Np=Np;    
          endif
          D0 = 1E6 * ( Q*1E-3*gamma(vu)/(apm*local_Np*gamma(vu+bpm)) )**(1/bpm)  ! units = microns
       endif
       Re = 0.5_wp*D0*gamma(vu+3)/gamma(vu+2)
       
       ! ---------------------------------------------------------!
       ! Exponential                                              !
       ! N0 = intercept parameter (m^-4)                          !
       ! ld = slope parameter (um)                                !
       ! ---------------------------------------------------------!
    case(2)
       
       ! Np not specified (see if statement above) 
       if((abs(p1+1) > 1E-8) ) then   ! N0 has been specified, determine ld
          N0   = p1
          tmp1 = 1._wp/(1._wp+bpm)
          ld   = ((apm*gamma(1.+bpm)*N0)/(rho_a*Q*1E-3))**tmp1
          ld   = ld/1E6                     ! set units to microns^-1
       elseif (abs(p2+1) > 1E-8) then  ! lambda=ld has been specified as default
          ld = p2     ! should have units of microns^-1 
       else
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
             'Must specify Np or default value (p1=No or p2=lambda) for Exponential distribution')
          return
       endif
       Re = 1.5_wp/ld 
       
       ! ---------------------------------------------------------!
       ! Power law                                                !
       ! ahp = Ar parameter (m^-4 mm^-bhp)                        !
       ! bhp = br parameter                                       !
       ! dmin_mm = lower bound (mm)                               !
       ! dmax_mm = upper bound (mm)                               !
       ! ---------------------------------------------------------!
    case(3)
       
       Re=0._wp  ! Not supporting LUT approach for power-law ...
       if(Np>0) then
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
             'Variable Np not supported for Power Law distribution')
          return
       endif
       
       ! ---------------------------------------------------------!
       ! Monodisperse                                             !
       ! D0 = particle diameter (um) == Re                        !
       ! ---------------------------------------------------------!
    case(4)
       
       Re = p1
       if(Np>0) then
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
             'Variable Np not supported for Monodispersed distribution')
          return
       endif
       
       ! ---------------------------------------------------------!
       ! Lognormal                                                !
       ! N0 = total number concentration (m^-3)                   !
       ! np = fixed number concentration (kg^-1)                  !
       ! rg = mean radius (um)                                    !
       ! log_sigma_g = ln(geometric standard deviation)           !
       ! ---------------------------------------------------------!
    case(5)
       
       if( abs(local_p3+1) > 1E-8 ) then
          !set natural log width
          log_sigma_g = local_p3 
       else
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
             'Must specify a value for sigma_g when using a Log-Normal distribution')
          return
       endif
       
       ! get rg ... 
       if( Np.eq.0 .and. (abs(p2+1) > 1E-8) ) then ! use default value of rg
          rg = p2     
       else
          if(Np>0) then
             local_Np=Np;  
          elseif(abs(p2+1) < 1E-8) then
             local_Np=p1
          else
             call errorMessage('ERROR(optics/quickbeam_optics.f90:Calc_Re): '//&
                'Must specify Np or default value (p2=Rg or p1=Np) for Log-Normal distribution')
          endif
          log_sigma_g = p3
          tmp1        = (Q*1E-3)/(2._wp**bpm*apm*local_Np)
          tmp2        = exp(0.5_wp*bpm*bpm*(log_sigma_g))*exp(0.5_wp*bpm*bpm*(log_sigma_g))    
          rg          = ((tmp1/tmp2)**(1._wp/bpm))*1E6
       endif
       Re = rg*exp(2.5_wp*(log_sigma_g*log_sigma_g))    
    end select
  end subroutine calc_Re
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine dsd(Q,Re,Np,D,N,nsizes,dtype,rho_a,tk,dmin,dmax,apm,bpm,rho_c,p1,p2,p3)
    ! ##############################################################################################
    ! Purpose:
    !   Create a discrete drop size distribution
    !
    !   Starting with Quickbeam V3, this routine now allows input of 
    !   both effective radius (Re) and total number concentration (Nt)
    !   Roj Marchand July 2010
    !
    !   The version in Quickbeam v.104 was modified to allow Re but not Nt 
    !   This is a significantly modified form for the version      
    !
    !   Originally Part of QuickBeam v1.03 by John Haynes
    !   http://reef.atmos.colostate.edu/haynes/radarsim
    !
    ! Inputs:
    !
    !   [Q]        hydrometeor mixing ratio (g/kg)
    !   [Re]       Optional Effective Radius (microns).  0 = use defaults (p1, p2, p3)
    !
    !   [D]        array of discrete drop sizes (um) where we desire to know the number concentraiton n(D).
    !   [nsizes]   number of elements of [D]
    !
    !   [dtype]    distribution type
    !   [rho_a]    ambient air density (kg m^-3)
    !   [tk]       temperature (K)
    !   [dmin]     minimum size cutoff (um)
    !   [dmax]     maximum size cutoff (um)
    !   [rho_c]    alternate constant density (kg m^-3)
    !   [p1],[p2],[p3]  distribution parameters
    !
    ! Input/Output:
    !   [apm]      a parameter for mass (kg m^[-bpm])
    !   [bmp]      b params for mass
    !
    ! Outputs:
    !   [N]        discrete concentrations (cm^-3 um^-1)
    !              or, for monodisperse, a constant (1/cm^3)
    !
    ! Requires:
    !   function infind
    !
    ! Created:
    !   11/28/05  John Haynes (haynes@atmos.colostate.edu)
    ! Modified:
    !   01/31/06  Port from IDL to Fortran 90
    !   07/07/06  Rewritten for variable DSD's
    !   10/02/06  Rewritten using scaling factors (Roger Marchand and JMH), Re added V1.04
    !   July 2020 "N Scale factors" (variable fc) removed (Roj Marchand).
    !   12/18/14  Define type REALs as double precision (dustin.swales@noaa.gov)
    ! ##############################################################################################
    
    ! Inputs
    integer, intent(in)                   :: &
         nsizes,& ! Number of elements of [D]
         dtype    !  distribution type
    real(wp),intent(in),dimension(nsizes) :: &
         D        ! Array of discrete drop sizes (um) where we desire to know the number concentraiton n(D).
    real(wp),intent(in) :: &
         Q,     & ! Hydrometeor mixing ratio (g/kg)
         Np,    & !
         rho_a, & ! Ambient air density (kg m^-3)
         tk,    & ! Temperature (K)
         dmin,  & ! Minimum size cutoff (um)
         dmax,  & ! Maximum size cutoff (um)
         rho_c, & ! Alternate constant density (kg m^-3)
         p1,    & ! Distribution parameter 1
         p2,    & ! Distribution parameter 2
         p3       ! Distribution parameter 3
    real(wp),intent(inout) :: &
         apm,   & ! a parameter for mass (kg m^[-bpm])
         bpm,   & ! b params for mass
         Re       ! Optional Effective Radius (microns)
    
    ! Outputs
    real(wp),intent(out),dimension(nsizes) :: &
         N        ! Discrete concentrations (cm^-3 um^-1)
                  ! or, for monodisperse, a constant (1/cm^3)
    
    ! Internal Variables
    real(wp),dimension(nsizes) :: &
         fc
    real(wp)                   :: &
         N0,D0,vu,local_np,dm,ld, & ! gamma, exponential variables
         dmin_mm,dmax_mm,ahp,bhp, & ! power law variables
         rg,log_sigma_g,          & ! lognormal variables
         rho_e,                   & ! particle density (kg m^-3)
         tmp1,tmp2,tc
    integer :: &
         k,lidx,uidx
    
    ! Convert temperature from Kelvin to Celsius
    tc = tk - 273.15_wp
    
    ! If density is constant, store equivalent values for apm and bpm
    if ((rho_c > 0) .and. (apm < 0)) then
       apm = (pi/6)*rho_c
       bpm = 3._wp
    endif
    
    ! Will preferentially use Re input over Np.
    ! if only Np given then calculate Re
    ! if neigher than use other defaults (p1,p2,p3) following quickbeam documentation
    if(Re==0 .and. Np>0) then
       call calc_Re(Q,Np,rho_a,dtype,apm,bpm,rho_c,p1,p2,p3,Re)
    endif
    select case(dtype)
       
       ! ---------------------------------------------------------!
       ! Modified gamma                                           !
       ! np = total number concentration                          !
       ! D0 = characteristic diameter (um)                        !
       ! dm = mean diameter (um) - first moment over zeroth moment!
       ! vu = distribution width parameter                        !
       ! ---------------------------------------------------------!
    case(1)  
       
       if( abs(p3+2) < 1E-8) then
          if( Np>1E-30) then
             ! Morrison scheme with Martin 1994 shape parameter (NOTE: vu = pc +1)
             ! fixed Roj. Dec. 2010 -- after comment by S. Mcfarlane
             vu = (1/(0.2714_wp + 0.00057145_wp*Np*rho_a*1E-6))**2._wp ! units of Nt = Np*rhoa = #/cm^3
          else
             call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:dsd): '//&
                'Must specify a value for Np in each volume with Morrison/Martin Scheme.')
             return
          endif
       elseif (abs(p3+1) > 1E-8) then
          ! vu is fixed in hp structure  
          vu = p3 
       else
          ! vu isn't specified
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:dsd): '//&
             'Must specify a value for vu for Modified Gamma distribution')
          return
       endif
       
       if(Re>0) then
          D0 = 2._wp*Re*gamma(vu+2)/gamma(vu+3)
          fc = (((D*1E-6)**(vu-1)*exp(-1*D/D0)) / &
               (apm*((D0*1E-6)**(vu+bpm))*gamma(vu+bpm))) * 1E-12
          N  = fc*rho_a*(Q*1E-3)
       elseif( p2+1 > 1E-8) then     ! use default value for MEAN diameter
          dm = p2
          D0 = gamma(vu)/gamma(vu+1)*dm
          fc = (((D*1E-6)**(vu-1)*exp(-1*D/D0)) / &
               (apm*((D0*1E-6)**(vu+bpm))*gamma(vu+bpm))) * 1E-12
          N  = fc*rho_a*(Q*1E-3)
       elseif(abs(p3+1) > 1E-8)  then! use default number concentration
          local_np = p1 ! total number concentration / pa check
          tmp1     = (Q*1E-3)**(1./bpm)
          fc       = (D*1E-6 / (gamma(vu)/(apm*local_np*gamma(vu+bpm)))**(1._wp/bpm))**vu
          N        = ((rho_a*local_np*fc*(D*1E-6)**(-1._wp))/(gamma(vu)*tmp1**vu) * &
               exp(-1._wp*fc**(1._wp/vu)/tmp1)) * 1E-12
       else
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:dsd): '//&
             'No default value for Dm or Np provided!')
          return
       endif
       
       ! ---------------------------------------------------------!
       ! Exponential                                              !
       ! N0 = intercept parameter (m^-4)                          !
       ! ld = slope parameter (um)                                !
       ! ---------------------------------------------------------!
    case(2)
       
       if(Re>0) then
          ld = 1.5_wp/Re   ! units 1/um
          fc = (ld*1E6)**(1.+bpm)/(apm*gamma(1+bpm))*exp(-1._wp*(ld*1E6)*(D*1E-6))*1E-12
          N  = fc*rho_a*(Q*1E-3)
       elseif (abs(p1+1) > 1E-8) then
          ! Use N0 default value
          N0   = p1
          tmp1 = 1._wp/(1._wp+bpm)
          fc   = ((apm*gamma(1.+bpm)*N0)**tmp1)*(D*1E-6)
          N    = (N0*exp(-1._wp*fc*(1._wp/(rho_a*Q*1E-3))**tmp1)) * 1E-12
       elseif (abs(p2+1) > 1E-8) then
          ! Use default value for lambda 
          ld = p2
          fc = (ld*1E6)**(1._wp+bpm)/(apm*gamma(1+bpm))*exp(-1._wp*(ld*1E6)*(D*1E-6))*1E-12
          N  = fc*rho_a*(Q*1E-3)
       else
          ! ld "parameterized" from temperature (carry over from original Quickbeam).
          ld = 1220._wp*10._wp**(-0.0245_wp*tc)*1E-6
          N0 = ((ld*1E6)**(1._wp+bpm)*Q*1E-3*rho_a)/(apm*gamma(1+bpm))
          N  = (N0*exp(-ld*D)) * 1E-12
       endif
       
       ! ---------------------------------------------------------!
       ! Power law                                                !
       ! ahp = Ar parameter (m^-4 mm^-bhp)                        !
       ! bhp = br parameter                                       !
       ! dmin_mm = lower bound (mm)                               !
       ! dmax_mm = upper bound (mm)                               !
       ! ---------------------------------------------------------!
    case(3)
       
       if(Re>0) then
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:dsd): '//&
             'Variable Re not supported for Power-Law distribution')
          return
       elseif(Np>0) then
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:dsd): '//&
             'Variable Np not supported for Power-Law distribution')
          return
       endif
       
       ! br parameter
       if (abs(p1+2) < 1E-8) then
          ! if p1=-2, bhp is parameterized according to Ryan (2000),
          ! applicatable to cirrus clouds
          if (tc < -30) then
             bhp = -1.75_wp+0.09_wp*((tc+273._wp)-243.16_wp)
          elseif ((tc >= -30) .and. (tc < -9)) then
             bhp = -3.25_wp-0.06_wp*((tc+273._wp)-265.66_wp)
          else
             bhp = -2.15_wp
          endif
       elseif (abs(p1+3) < 1E-8) then      
          ! if p1=-3, bhp is parameterized according to Ryan (2000),
          ! applicable to frontal clouds
          if (tc < -35) then
             bhp = -1.75_wp+0.09_wp*((tc+273._wp)-243.16_wp)
          elseif ((tc >= -35) .and. (tc < -17.5)) then
             bhp = -2.65_wp+0.09_wp*((tc+273._wp)-255.66_wp)
          elseif ((tc >= -17.5) .and. (tc < -9)) then
             bhp = -3.25_wp-0.06_wp*((tc+273._wp)-265.66_wp)
          else
             bhp = -2.15_wp
          endif
       else
          ! Otherwise the specified value is used
          bhp = p1
       endif
       
       ! Ar parameter
       dmin_mm = dmin*1E-3
       dmax_mm = dmax*1E-3
       
       ! Commented lines are original method with constant density
       ! rc = 500.       ! (kg/m^3)
       ! tmp1 = 6*rho_a*(bhp+4)
       ! tmp2 = pi*rc*(dmax_mm**(bhp+4))*(1-(dmin_mm/dmax_mm)**(bhp+4))
       ! ahp = (Q*1E-3)*1E12*tmp1/tmp2
       
       ! New method is more consistent with the rest of the distributions
       ! and allows density to vary with particle size
       tmp1 = rho_a*(Q*1E-3)*(bhp+bpm+1)
       tmp2 = apm*(dmax_mm**bhp*dmax**(bpm+1)-dmin_mm**bhp*dmin**(bpm+1))
       ahp  = tmp1/tmp2 * 1E24
       ! ahp = tmp1/tmp2 
       lidx = infind(D,dmin)
       uidx = infind(D,dmax)    
       do k=lidx,uidx
          N(k) = (ahp*(D(k)*1E-3)**bhp) * 1E-12    
       enddo
       
       ! ---------------------------------------------------------!
       ! Monodisperse                                             !
       ! D0 = particle diameter (um)                              !
       ! ---------------------------------------------------------!
    case(4)
       
       if (Re>0) then
          D0 = Re
       else
          D0 = p1
       endif
       
       rho_e = (6._wp/pi)*apm*(D0*1E-6)**(bpm-3)
       fc(1) = (6._wp/(pi*D0*D0*D0*rho_e))*1E12
       N(1)  = fc(1)*rho_a*(Q*1E-3)
       
       ! ---------------------------------------------------------!
       ! Lognormal                                                !
       ! N0 = total number concentration (m^-3)                   !
       ! np = fixed number concentration (kg^-1)                  !
       ! rg = mean radius (um)                                    !
       ! og_sigma_g = ln(geometric standard deviation)            !
       ! ---------------------------------------------------------!
    case(5)
       if (abs(p1+1) < 1E-8 .or. Re>0 ) then
          ! rg, log_sigma_g are given
          log_sigma_g = p3
          tmp2 = (bpm*log_sigma_g)*(bpm*log_sigma_g)
          if(Re.le.0) then 
             rg = p2
          else
             !rg = Re*exp(-2.5*(log_sigma_g*log_sigma_g))
             rg =Re*exp(-2.5_wp*(log_sigma_g**2))
             
          endif
          
          fc = 0.5_wp*((1._wp/((2._wp*rg*1E-6)**(bpm)*apm*(2._wp*pi)**(0.5_wp) * &
               log_sigma_g*D*0.5_wp*1E-6))*exp(-0.5_wp*((log(0.5_wp*D/rg)/log_sigma_g)**2._wp+tmp2)))*1E-12
          N = fc*rho_a*(Q*1E-3)
          
       elseif (abs(p2+1) < 1E-8 .or. Np>0) then
          ! Np, log_sigma_g are given    
          if(Np>0) then
             local_Np = Np
          else
             local_Np = p1
          endif
          
          log_sigma_g = p3
          N0   = local_np*rho_a
          tmp1 = (rho_a*(Q*1E-3))/(2._wp**bpm*apm*N0)
          tmp2 = exp(0.5_wp*bpm*bpm*(log_sigma_g))*exp(0.5_wp*bpm*bpm*(log_sigma_g))
          rg   = ((tmp1/tmp2)**(1/bpm))*1E6
          
          N = 0.5_wp*(N0 / ((2._wp*pi)**(0.5_wp)*log_sigma_g*D*0.5_wp*1E-6) * &
               exp((-0.5_wp*(log(0.5_wp*D/rg)/log_sigma_g)**2._wp)))*1E-12      
       else
          call errorMessage('FATAL ERROR(optics/quickbeam_optics.f90:dsd): '//&
             'Must specify a value for sigma_g')
          return
       endif
    end select
  end subroutine dsd
  ! ##############################################################################################
  ! ##############################################################################################
  subroutine zeff(freq,D,N,nsizes,k2,tt,ice,xr,z_eff,z_ray,kr,qe,qs,rho_e)
    ! ##############################################################################################
    ! Purpose:
    !   Simulates radar return of a volume given DSD of spheres
    !   Part of QuickBeam v1.03 by John Haynes
    !   http://reef.atmos.colostate.edu/haynes/radarsim
    !
    ! Inputs:
    !   [freq]      radar frequency (GHz)
    !   [D]         discrete drop sizes (um)
    !   [N]         discrete concentrations (cm^-3 um^-1)
    !   [nsizes]    number of discrete drop sizes
    !   [k2]        |K|^2, -1=use frequency dependent default 
    !   [tt]        hydrometeor temperature (K)
    !   [ice]       indicates volume consists of ice
    !   [xr]        perform Rayleigh calculations?
    !   [qe]        if using a mie table, these contain ext/sca ...
    !   [qs]        ... efficiencies; otherwise set to -1
    !   [rho_e]     medium effective density (kg m^-3) (-1 = pure)
    !
    ! Outputs:
    !   [z_eff]     unattenuated effective reflectivity factor (mm^6/m^3)
    !   [z_ray]     reflectivity factor, Rayleigh only (mm^6/m^3)
    !   [kr]        attenuation coefficient (db km^-1)
    !
    ! Created:
    !   11/28/05  John Haynes (haynes@atmos.colostate.edu)
    ! Modified:
    !   12/18/14  Dustin Swales: Define type REALs as double precision (dustin.swales@noaa.gov)
    ! ##############################################################################################
    ! Inputs
    integer,  intent(in) :: &
         ice,        & ! Indicates volume consists of ice
         xr,         & ! Perform Rayleigh calculations?
         nsizes        ! Number of discrete drop sizes
    real(wp), intent(in),dimension(nsizes) :: &
         D,     & ! Discrete drop sizes (um)
         N,     & ! Discrete concentrations (cm^-3 um^-1)
         rho_e, & ! Medium effective density (kg m^-3) (-1 = pure)
         qe,    & ! Extinction efficiency, when using Mie tables
         qs       ! Scatering efficiency, when using Mie tables
    real(wp),intent(in) :: &
         freq,  & ! Radar frequency (GHz)
         tt       ! Hydrometeor temperature (K)
    real(wp), intent(inout) :: &
         k2            ! |K|^2, -1=use frequency dependent default 
    
    ! Outputs
    real(wp), intent(out) :: &
         z_eff,      & ! Unattenuated effective reflectivity factor (mm^6/m^3)
         z_ray,      & ! Reflectivity factor, Rayleigh only (mm^6/m^3)
         kr            ! Attenuation coefficient (db km^-1)
    
    ! Internal Variables
    integer                     :: correct_for_rho ! Correct for density flag
    real(wp), dimension(nsizes) :: &
         D0,        &    ! D in (m)
         N0,        &    ! N in m^-3 m^-1
         sizep,     &    ! Size parameter
         qext,      &    ! Extinction efficiency
         qbsca,     &    ! Backscatter efficiency
         f,         &    ! Ice fraction
         xtemp           !
    real(wp) :: &
         wl, cr,eta_sum,eta_mie,const,z0_eff,z0_ray,k_sum,n_r,n_i,dqv(1),dqsc,dg,dph(1)
    complex(wp)         :: &
         m,                  &    ! Complex index of refraction of bulk form
         Xs1(1), Xs2(1)           !
    integer          :: &
         i, err                   !
    integer, parameter :: &
         one=1                    !
    real(wp),parameter :: &
         conv_d  = 1e-6,     &    ! Conversion factor for drop sizes (to m)
         conv_n  = 1e12,     &    ! Conversion factor for drop concentrations (to m^-3)
         conv_f  = 0.299792458    ! Conversion for radar frequency (to m)
    complex(wp),dimension(nsizes) ::&
         m0             ! Complex index of refraction
    
    ! Initialize
    z0_ray = 0._wp
    
    ! Conversions
    D0 = d*conv_d
    N0 = n*conv_n
    wl = conv_f/freq
    
    ! // dielectric constant |k^2| defaults
    if (k2 < 0) then
       k2 = 0.933_wp
       if (abs(94.-freq) < 3.) k2=0.75_wp
       if (abs(35.-freq) < 3.) k2=0.88_wp
       if (abs(13.8-freq) < 3.) k2=0.925_wp
    endif
    
    if (qe(1) < -9) then
       
       ! Get the refractive index of the bulk hydrometeors
       if (ice == 0) then
          call m_wat(freq,tt,n_r,n_i)
       else
          call m_ice(freq,tt,n_r,n_i)
       endif
       m = cmplx(n_r,-n_i)
       m0(1:nsizes) = m
       
       correct_for_rho = 0
       if ((ice == 1) .and. (minval(rho_e) >= 0)) correct_for_rho = 1
       
       ! Correct refractive index for ice density if needed
       if (correct_for_rho == 1) then
          f  = rho_e/rhoice
          m0 = sqrt((2+(m0*m0)+2*f*((m0*m0)-1))/(2+(m0*m0)+f*(1-(m0*m0))))
       endif
       
       ! Mie calculations
       sizep = (pi*D0)/wl
       dqv(1) = 0._wp
       do i=1,nsizes
          call mieint(sizep(i), m0(i), one, dqv, qext(i), dqsc, qbsca(i), &
               dg, xs1, xs2, dph, err)
       end do

    else
       ! Mie table used
       qext  = qe
       qbsca = qs
    endif
    
    ! eta_mie = 0.25*sum[qbsca*pi*D^2*N(D)*deltaD]
    ! <--------- eta_sum --------->
    ! z0_eff = (wl^4/!pi^5)*(1./k2)*eta_mie
    eta_sum = 0._wp
    if (size(D0) == 1) then
       eta_sum = qbsca(1)*(n(1)*1E6)*D0(1)*D0(1)
    else
       xtemp = qbsca*N0*D0*D0
       call avint(xtemp,D0,nsizes,D0(1),D0(size(D0,1)),eta_sum)
    endif
    
    eta_mie = eta_sum*0.25_wp*pi
    const   = ((wl*wl*wl*wl)/(pi*pi*pi*pi*pi))*(1._wp/k2)
    
    z0_eff  = const*eta_mie
    
    ! kr = 0.25*cr*sum[qext*pi*D^2*N(D)*deltaD]
    ! <---------- k_sum --------->  
    k_sum = 0._wp
    if (size(D0) == 1) then
       k_sum = qext(1)*(n(1)*1E6)*D0(1)*D0(1)
    else
       xtemp = qext*N0*D0*D0
       call avint(xtemp,D0,nsizes,D0(1),D0(size(D0,1)),k_sum)
    endif
    ! DS2014 START: Making this calculation in double precision results in a small 
    !               amount of very small errors in the COSP output field,dBZE94,
    !               so it will be left as is.
    !cr = 10._wp/log(10._wp)
    cr = 10./log(10.)
    ! DS2014 STOP
    kr = k_sum*0.25_wp*pi*(1000._wp*cr)
    
    ! z_ray = sum[D^6*N(D)*deltaD]
    if (xr == 1) then
       z0_ray = 0._wp
       if (size(D0) == 1) then
          z0_ray = (n(1)*1E6)*D0(1)*D0(1)*D0(1)*D0(1)*D0(1)*D0(1)
       else
          xtemp = N0*D0*D0*D0*D0*D0*D0
          call avint(xtemp,D0,nsizes,D0(1),D0(size(D0)),z0_ray)
       endif
    endif
    
    ! Convert to mm^6/m^3
    z_eff = z0_eff*1E18 !  10.*alog10(z0_eff*1E18)
    z_ray = z0_ray*1E18 !  10.*alog10(z0_ray*1E18)
    
  end subroutine zeff
  ! ##############################################################################################
  ! ##############################################################################################
  function gases(PRES_mb,T,SH,f)
    ! ##############################################################################################
    ! Purpose:
    !   Compute 2-way gaseous attenuation through a volume in microwave
    !
    ! Inputs:
    !   [PRES_mb]   pressure (mb) (hPa)
    !   [T]         temperature (K)
    !   [RH]        relative humidity (%)
    !   [f]         frequency (GHz), < 300 GHz
    !
    ! Returns:
    !   2-way gaseous attenuation (dB/km)
    !
    ! Reference:
    !   Uses method of Liebe (1985)
    !
    ! Created:
    !   12/09/05  John Haynes (haynes@atmos.colostate.edu)
    ! Modified:
    !   01/31/06  Port from IDL to Fortran 90
    !   12/19/14  Dustin Swales: Define type REALs as double precision (dustin.swales@noaa.gov)
    ! ##############################################################################################
    
    ! INPUTS
    real(wp), intent(in) :: & !
         PRES_mb,           & ! Pressure (mb) (hPa)
         T,                 & ! Temperature (K)
         SH,                & ! Specific humidity
         f                    ! Frequency (GHz), < 300 GHz
    
    ! PARAMETERS
    integer, parameter   :: & !
         nbands_o2  = 48,   & ! Number of O2 bands
         nbands_h2o = 30      ! Number of h2o bands
    ! LOCAL VARIABLES
    real(wp) :: &
         gases, th, e, p, sumo, gm0, a0, ap, term1,    &
         term2, term3, bf, be, term4, npp,e_th,one_th, &
         pth3,eth35,aux1,aux2,aux3, aux4,gm,delt,x,y,  &
         gm2,fpp_o2,fpp_h2o,s_o2,s_h2o
    integer :: i

    ! Table1 parameters  v0, a1, a2, a3, a4, a5, a6  
    real(wp),dimension(nbands_o2),parameter ::                                          &
         v0 = (/49.4523790,49.9622570,50.4742380,50.9877480,51.5033500,                 &
                52.0214090,52.5423930,53.0669060,53.5957480,54.1299999,54.6711570,      &
                55.2213650,55.7838000,56.2647770,56.3378700,56.9681000,57.6124810,      &
                58.3238740,58.4465890,59.1642040,59.5909820,60.3060570,60.4347750,      &
                61.1505580,61.8001520,62.4112120,62.4862530,62.9979740,63.5685150,      &
                64.1277640,64.6789000,65.2240670,65.7647690,66.3020880,66.8368270,      &
                67.3695950,67.9008620,68.4310010,68.9603060,69.4890210,70.0173420,      &
                118.7503410,368.4983500,424.7631200,487.2493700,715.3931500,            &
                773.8387300, 834.1453300/),                                             &
         a1 = (/0.0000001,0.0000003,0.0000009,0.0000025,0.0000061,0.0000141,            &
                0.0000310,0.0000641,0.0001247,0.0002280,0.0003918,0.0006316,0.0009535,  &
                0.0005489,0.0013440,0.0017630,0.0000213,0.0000239,0.0000146,0.0000240,  &
                0.0000211,0.0000212,0.0000246,0.0000250,0.0000230,0.0000193,0.0000152,  &
                0.0000150,0.0000109,0.0007335,0.0004635,0.0002748,0.0001530,0.0000801,  &
                0.0000395,0.0000183,0.0000080,0.0000033,0.0000013,0.0000005,0.0000002,  &
                0.0000094,0.0000679,0.0006380,0.0002350,0.0000996,0.0006710,0.0001800/),&
         a2 = (/11.8300000,10.7200000,9.6900000,8.8900000,7.7400000,6.8400000,          &
                6.0000000,5.2200000,4.4800000,3.8100000,3.1900000,2.6200000,2.1150000,  &
                0.0100000,1.6550000,1.2550000,0.9100000,0.6210000,0.0790000,0.3860000,  &
                0.2070000,0.2070000,0.3860000,0.6210000,0.9100000,1.2550000,0.0780000,  &
                1.6600000,2.1100000,2.6200000,3.1900000,3.8100000,4.4800000,5.2200000,  &
                6.0000000,6.8400000,7.7400000,8.6900000,9.6900000,10.7200000,11.8300000,&
                0.0000000,0.0200000,0.0110000,0.0110000,0.0890000,0.0790000,0.0790000/),&
         a3 = (/0.0083000,0.0085000,0.0086000,0.0087000,0.0089000,0.0092000,            &
                0.0094000,0.0097000,0.0100000,0.0102000,0.0105000,0.0107900,0.0111000,  &
                0.0164600,0.0114400,0.0118100,0.0122100,0.0126600,0.0144900,0.0131900,  &
                0.0136000,0.0138200,0.0129700,0.0124800,0.0120700,0.0117100,0.0146800,  &
                0.0113900,0.0110800,0.0107800,0.0105000,0.0102000,0.0100000,0.0097000,  &
                0.0094000,0.0092000,0.0089000,0.0087000,0.0086000,0.0085000,0.0084000,  &
                0.0159200,0.0192000,0.0191600,0.0192000,0.0181000,0.0181000,0.0181000/),&
         a4 = (/0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,            &
                0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,  &
                0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,  &
                0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,  &
                0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,  &
                0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000,  &
                0.0000000,0.6000000,0.6000000,0.6000000,0.6000000,0.6000000,0.6000000/),&
         a5 = (/0.0056000,0.0056000,0.0056000,0.0055000,0.0056000,0.0055000,            &
                0.0057000,0.0053000,0.0054000,0.0048000,0.0048000,0.0041700,0.0037500,  &
                0.0077400,0.0029700,0.0021200,0.0009400,-0.0005500,0.0059700,-0.0024400,&
                0.0034400,-0.0041300,0.0013200,-0.0003600,-0.0015900,-0.0026600,        &
                -0.0047700,-0.0033400,-0.0041700,-0.0044800,-0.0051000,-0.0051000,      &
                -0.0057000,-0.0055000,-0.0059000,-0.0056000,-0.0058000,-0.0057000,      &
                -0.0056000,-0.0056000,-0.0056000,-0.0004400,0.0000000,0.0000000,        &
                0.0000000,0.0000000,0.0000000,0.0000000/),                              & 
         a6 = (/1.7000000,1.7000000,1.7000000,1.7000000,1.8000000,1.8000000,            &
                1.8000000,1.9000000,1.8000000,2.0000000,1.9000000,2.1000000,2.1000000,  &
                0.9000000,2.3000000,2.5000000,3.7000000,-3.1000000,0.8000000,0.1000000, &
                0.5000000,0.7000000,-1.0000000,5.8000000,2.9000000,2.3000000,0.9000000, &
                2.2000000,2.0000000,2.0000000,1.8000000,1.9000000,1.8000000,1.8000000,  &
                1.7000000,1.8000000,1.7000000,1.7000000,1.7000000,1.7000000,1.7000000,  &
                0.9000000,1.0000000,1.0000000,1.0000000,1.0000000,1.0000000,1.0000000/)
    
    ! Table2 parameters  v1, b1, b2, b3
    real(wp),dimension(nbands_h2o),parameter ::                                          &
         v1 = (/22.2350800,67.8139600,119.9959400,183.3101170,321.2256440,               &
                325.1529190,336.1870000,380.1973720,390.1345080,437.3466670,439.1508120, &
                443.0182950,448.0010750,470.8889740,474.6891270,488.4911330,503.5685320, &
                504.4826920,556.9360020,620.7008070,658.0065000,752.0332270,841.0735950, &
                859.8650000,899.4070000,902.5550000,906.2055240,916.1715820,970.3150220, &
                987.9267640/),                                                           &
         b1 = (/0.1090000,0.0011000,0.0007000,2.3000000,0.0464000,1.5400000,             &
                0.0010000,11.9000000,0.0044000,0.0637000,0.9210000,0.1940000,10.6000000, &
                0.3300000,1.2800000,0.2530000,0.0374000,0.0125000,510.0000000,5.0900000, &
                0.2740000,250.0000000,0.0130000,0.1330000,0.0550000,0.0380000,0.1830000, &
                8.5600000,9.1600000,138.0000000/),                                       &
         b2 = (/2.1430000,8.7300000,8.3470000,0.6530000,6.1560000,1.5150000,             &
                9.8020000,1.0180000,7.3180000,5.0150000,3.5610000,5.0150000,1.3700000,   &
                3.5610000,2.3420000,2.8140000,6.6930000,6.6930000,0.1140000,2.1500000,   &
                7.7670000,0.3360000,8.1130000,7.9890000,7.8450000,8.3600000,5.0390000,   &
                1.3690000,1.8420000,0.1780000/),                                         &
         b3 = (/0.0278400,0.0276000,0.0270000,0.0283500,0.0214000,0.0270000,             &
                0.0265000,0.0276000,0.0190000,0.0137000,0.0164000,0.0144000,0.0238000,   &
                0.0182000,0.0198000,0.0249000,0.0115000,0.0119000,0.0300000,0.0223000,   &
                0.0300000,0.0286000,0.0141000,0.0286000,0.0286000,0.0264000,0.0234000,   &
                0.0253000,0.0240000,0.0286000/)

    ! Conversions
    th     = 300._wp/T                                             ! unitless

    ! DS2014 START: Using _wp for the exponential in the denominator results in slight errors
    !               for dBze94. 0.01 % of values differ, relative range: 1.03e-05 to  1.78e-04
    !e      = (RH*th*th*th*th*th)/(41.45_wp*10**(9.834_wp*th-10))   ! kPa
    !e = (RH*th*th*th*th*th)/(41.45_wp*10**(9.834_wp*th-10))   ! kPa
    e = SH*PRES_mb/(SH+0.622_wp)/1000._wp !kPa
    ! DS2014 END

    p      = PRES_mb/1000._wp-e                                      ! kPa
    e_th   = e*th
    one_th = 1 - th
    pth3   = p*th*th*th
    eth35  = e*th**(3.5)
    
    ! Term1
    sumo = 0._wp
    aux1 = 1.1_wp*e_th
    do i=1,nbands_o2
       aux2   = f/v0(i)
       aux3   = v0(i)-f
       aux4   = v0(i)+f
       gm     = a3(i)*(p*th**(0.8_wp-a4(i))+aux1)
       gm2    = gm*gm
       delt   = a5(i)*p*th**a6(i)
       x      = aux3*aux3+gm2
       y      = aux4*aux4+gm2
       fpp_o2 = (((1._wp/x)+(1._wp/y))*(gm*aux2) - (delt*aux2)*((aux3/(x))-(aux4/(x))))
       s_o2   = a1(i)*pth3*exp(a2(i)*one_th)
       sumo   = sumo + fpp_o2 * s_o2
    enddo
    term1 = sumo

    ! Term2
    gm0   = 5.6E-3_wp*(p+1.1_wp*e)*th**(0.8_wp)
    a0    = 3.07E-4_wp
    ap    = 1.4_wp*(1-1.2_wp*f**(1.5_wp)*1E-5)*1E-10
    term2 = (2*a0*(gm0*(1+(f/gm0)*(f/gm0))*(1+(f/60._wp)**2))**(-1) + ap*p*th**(2.5_wp))*f*p*th*th

    ! Term3
    sumo = 0._wp
    aux1 = 4.8_wp*e_th
    do i=1,nbands_h2o
       aux2    = f/v1(i)
       aux3    = v1(i)-f
       aux4    = v1(i)+f
       gm      = b3(i)*(p*th**(0.8)+aux1)
       gm2     = gm*gm
       x       = aux3*aux3+gm2
       y       = aux4*aux4+gm2
       fpp_h2o = ((1._wp/x)+(1._wp/y))*(gm*aux2) ! - (delt*aux2)*((aux3/(x))-(aux4/(x)))
       s_h2o   = b1(i)*eth35*exp(b2(i)*one_th)
       sumo    = sumo + fpp_h2o * s_h2o
    enddo
    term3 = sumo

    ! Term4
    bf    = 1.4E-6_wp
    be    = 5.41E-5_wp
    term4 = (bf*p+be*e*th*th*th)*f*e*th**(2.5_wp)

    ! Summation and result
    npp   = term1 + term2 + term3 + term4
    gases = 0.182_wp*f*npp
    
  end function gases
 subroutine hydro_class_init(lsingle,ldouble,sd)
    ! ##############################################################################################
    ! Purpose:
    !
    !   Initialize variables used by the radar simulator.
    !   Part of QuickBeam v3.0 by John Haynes and Roj Marchand
    !   
    ! Inputs:  
    !   NAME            SIZE        DESCRIPTION
    !   [lsingle]       (1)         Logical flag to use single moment
    !   [ldouble]       (1)         Logical flag to use two moment
    ! Outputs:
    !   [sd]                        Structure that define hydrometeor types
    !
    ! Local variables:
    !   [n_hydro]       (1)         Number of hydrometeor types
    !   [hclass_type]   (nhclass)   Type of distribution (see quickbeam documentation)
    !   [hclass_phase]  (nhclass)   1==ice, 0=liquid
    !   [hclass_dmin]   (nhclass)   Minimum diameter allowed is drop size distribution N(D<Dmin)=0
    !   [hclass_dmax]   (nhclass)   Maximum diameter allowed is drop size distribution N(D>Dmax)=0
    !   [hclass_apm]    (nhclass)   Density of partical apm*D^bpm or constant = rho
    !   [hclass_bpm]    (nhclass)   Density of partical apm*D^bpm or constant = rho
    !   [hclass_rho]    (nhclass)   Density of partical apm*D^bpm or constant = rho
    !   [hclass_p1]     (nhclass)   Default values of DSD parameters (see quickbeam documentation)
    !   [hclass_p2]     (nhclass)   Default values of DSD parameters (see quickbeam documentation)
    !   [hclass_p3]     (nhclass)   Default values of DSD parameters (see quickbeam documentation)    
    ! Modified:
    !   08/23/2006  placed into subroutine form (Roger Marchand)
    !   June 2010   New interface to support "radar_simulator_params" structure
    !   12/22/2014  Moved radar simulator (CLOUDSAT) configuration initialization to cloudsat_init
    ! ##############################################################################################

    ! ####################################################################################
    ! NOTES on HCLASS variables
    !
    ! TYPE - Set to
    ! 1 for modified gamma distribution,
    ! 2 for exponential distribution,
    ! 3 for power law distribution,
    ! 4 for monodisperse distribution,
    ! 5 for lognormal distribution.
	!
    ! PHASE - Set to 0 for liquid, 1 for ice.
    ! DMIN  - The minimum drop size for this class (micron), ignored for monodisperse.
    ! DMAX  - The maximum drop size for this class (micron), ignored for monodisperse.
    ! Important note: The settings for DMIN and DMAX are
    ! ignored in the current version for all distributions except for power
    ! law. Except when the power law distribution is used, particle size
    ! is fixed to vary from zero to infinity, a restriction that is expected
    ! to be lifted in future versions. A placeholder must still be specified
    ! for each.
    ! Density of particles is given by apm*D^bpm or a fixed value rho. ONLY specify ONE of these two!!
    ! APM - The alpha_m coefficient in equation (1) (kg m**-beta_m )
    ! BPM - The beta_m coefficient in equation (1), see section 4.1.
    ! RHO - Hydrometeor density (kg m-3 ).
    ! 
    ! P1, P2, P3 - are default distribution parameters that depend on the type
    ! of distribution (see quickmbeam documentation for more information)
    !
    ! Modified Gamma (must set P3 and one of P1 or P2)
    ! P1 - Set to the total particle number concentration Nt /rho_a (kg-1 ), where
    ! rho_a is the density of air in the radar volume.
    ! P2 - Set to the particle mean diameter D (micron).
    ! P3 - Set to the distribution width nu.
    !
    ! Exponetial (set one of)
    ! P1 - Set to a constant intercept parameter N0 (m-4).
    ! P2 - Set to a constant lambda (micron-1).
    !
    ! Power Law
    ! P1 - Set this to the value of a constant power law parameter br
    !
    ! Monodisperse
    ! P1 - Set to a constant diameter D0 (micron) = Re.
    !
    ! Log-normal (must set P3 and one of P1 or P2)
    ! P1 - Set to the total particle number concentration Nt /rho_a (kg-1 )
    ! P2 - Set to the geometric mean particle radius rg (micron).
    ! P3 - Set to the natural logarithm of the geometric standard deviation.
    ! ####################################################################################
    ! INPUTS
    logical,intent(in) :: &
       lsingle, & ! True -> use single moment
       ldouble    ! True -> use two moment 
                     
    ! OUTPUTS
    type(size_distribution),intent(out) ::&
         sd              !

   ! SINGLE MOMENT PARAMETERS
   integer,parameter,dimension(N_HYDRO) :: &
                    ! LSL  LSI  LSR  LSS  CVL  CVI  CVR  CVS  LSG    
       HCLASS1_TYPE  = (/5,   1,   2,   2,   5,   1,   2,   2,   2/), & ! 
       HCLASS1_PHASE = (/0,   1,   0,   1,   0,   1,   0,   1,   1/)    ! 
   real(wp),parameter,dimension(N_HYDRO) ::&
                      ! LSL   LSI    LSR    LSS    CVL   CVI    CVR    CVS    LSG    
       HCLASS1_DMIN = (/ -1.,  -1.,   -1.,   -1.,   -1.,  -1.,   -1.,   -1.,   -1.  /),  &
       HCLASS1_DMAX = (/ -1.,  -1.,   -1.,   -1.,   -1.,  -1.,   -1.,   -1.,   -1.  /),  &
       HCLASS1_APM  = (/524., 110.8, 524.,   -1.,  524., 110.8, 524.,   -1.,   -1.  /),  &
       HCLASS1_BPM  = (/  3.,   2.91,  3.,   -1.,    3.,   2.91,  3.,   -1.,   -1.  /),  &
       HCLASS1_RHO  = (/ -1.,  -1.,   -1.,  100.,   -1.,  -1.,   -1.,  100.,  400.  /),  &
       HCLASS1_P1   = (/ -1.,  -1.,    8.e6,  3.e6, -1.,  -1.,    8.e6,  3.e6,  4.e6/),  & 
       HCLASS1_P2   = (/  6.,  40.,   -1.,   -1.,    6.,  40.,   -1.,   -1.,   -1.   /), & 
       HCLASS1_P3   = (/  0.3,  2.,   -1.,   -1.,    0.3,  2.,   -1.,   -1.,   -1.   /)

    ! TWO MOMENT PARAMETERS
    integer,parameter,dimension(N_HYDRO) :: &
                      ! LSL  LSI  LSR  LSS  CVL  CVI  CVR  CVS  LSG
       HCLASS2_TYPE  = (/ 1,   1,   1,   1,   1,   1,   1,   1,   1/), &
       HCLASS2_PHASE = (/ 0,   1,   0,   1,   0,   1,   0,   1,   1/)

    real(wp),parameter,dimension(N_HYDRO) :: &
                      ! LSL    LSI      LSR     LSS   CVL    CVI   CVR     CVS    LSG
       HCLASS2_DMIN = (/ -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/), &
       HCLASS2_DMAX = (/ -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/), &        
       HCLASS2_APM  = (/524,     -1,    524,     -1,   524,    -1,  524,     -1,   -1/), &
       HCLASS2_BPM  = (/  3,     -1,      3,     -1,     3,    -1,    3,     -1,   -1/), &
       HCLASS2_RHO  = (/ -1,    500,     -1,    100,    -1,   500,   -1,    100,  900/), &
       HCLASS2_P1   = (/ -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/), &
       HCLASS2_P2   = (/ -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/), &
       HCLASS2_P3   = (/ -2,      1,      1,      1,    -2,     1,    1,      1,    1/) 
    
    if (lsingle) then    
       sd%dtype(1:N_HYDRO) = HCLASS1_TYPE(1:N_HYDRO)
       sd%phase(1:N_HYDRO) = HCLASS1_PHASE(1:N_HYDRO)
       sd%dmin(1:N_HYDRO)  = HCLASS1_DMIN(1:N_HYDRO)
       sd%dmax(1:N_HYDRO)  = HCLASS1_DMAX(1:N_HYDRO)
       sd%apm(1:N_HYDRO)   = HCLASS1_APM(1:N_HYDRO)
       sd%bpm(1:N_HYDRO)   = HCLASS1_BPM(1:N_HYDRO)
       sd%rho(1:N_HYDRO)   = HCLASS1_RHO(1:N_HYDRO)
       sd%p1(1:N_HYDRO)    = HCLASS1_P1(1:N_HYDRO)
       sd%p2(1:N_HYDRO)    = HCLASS1_P2(1:N_HYDRO)
       sd%p3(1:N_HYDRO)    = HCLASS1_P3(1:N_HYDRO)
    endif
    if (ldouble) then    
       sd%dtype(1:N_HYDRO) = HCLASS2_TYPE(1:N_HYDRO)
       sd%phase(1:N_HYDRO) = HCLASS2_PHASE(1:N_HYDRO)
       sd%dmin(1:N_HYDRO)  = HCLASS2_DMIN(1:N_HYDRO)
       sd%dmax(1:N_HYDRO)  = HCLASS2_DMAX(1:N_HYDRO)
       sd%apm(1:N_HYDRO)   = HCLASS2_APM(1:N_HYDRO)
       sd%bpm(1:N_HYDRO)   = HCLASS2_BPM(1:N_HYDRO)
       sd%rho(1:N_HYDRO)   = HCLASS2_RHO(1:N_HYDRO)
       sd%p1(1:N_HYDRO)    = HCLASS2_P1(1:N_HYDRO)
       sd%p2(1:N_HYDRO)    = HCLASS2_P2(1:N_HYDRO)
       sd%p3(1:N_HYDRO)    = HCLASS2_P3(1:N_HYDRO)
    endif    
  end subroutine hydro_class_init    
end module mod_quickbeam_optics
