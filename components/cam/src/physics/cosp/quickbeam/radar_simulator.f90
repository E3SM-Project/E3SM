  subroutine radar_simulator( &
    hp, &
    nprof,ngate, &
    undef, &
    hgt_matrix,hm_matrix,re_matrix,Np_matrix, &
    p_matrix,t_matrix,rh_matrix, &
    Ze_non,Ze_ray,a_to_vol,g_to_vol,dBZe, &
    g_to_vol_in,g_to_vol_out)
 
  use m_mrgrnk 
  use array_lib
  use math_lib
  use optics_lib
  use radar_simulator_types
  use scale_LUTs_io
  implicit none
  
! Purpose:
!
!   Simulates a vertical profile of radar reflectivity
!   Originally Part of QuickBeam v1.04 by John Haynes & Roger Marchand.
!   but has been substantially modified since that time by
!   Laura Fowler and Roger Marchand (see modifications below).
!
! Inputs:
!
!   [hp]              structure that defines hydrometeor types and other radar properties
!
!   [nprof]           number of hydrometeor profiles
!   [ngate]           number of vertical layers
!
!   [undef]           missing data value
!   (The following 5 arrays must be in order from closest to the radar
!    to farthest...)
!
!   [hgt_matrix]      height of hydrometeors (km)
!   [p_matrix]        pressure profile (hPa)
!   [t_matrix]        temperature profile (K)
!   [rh_matrix]       relative humidity profile (%) -- only needed if gaseous aborption calculated.
!
!   [hm_matrix]       table of hydrometeor mixing rations (g/kg)
!   [re_matrix]       table of hydrometeor effective radii.  0 ==> use defaults. (units=microns)    
!   [Np_matrix]	      table of hydrometeor number concentration.  0 ==> use defaults. (units = 1/kg)
!
! Outputs:
!
!   [Ze_non]          radar reflectivity without attenuation (dBZ)
!   [Ze_ray]          Rayleigh reflectivity (dBZ)
!   [h_atten_to_vol]  attenuation by hydromets, radar to vol (dB)
!   [g_atten_to_vol]  gaseous atteunation, radar to vol (dB)
!   [dBZe]            effective radar reflectivity factor (dBZ)
!
! Optional:
!   [g_to_vol_in]     integrated atten due to gases, r>v (dB).
!                     If present then is used as gaseous absorption, independently of the
!                     value in use_gas_abs
!   [g_to_vol_out]    integrated atten due to gases, r>v (dB).
!                     If present then gaseous absorption for each profile is returned here.
!
! Created:
!   11/28/2005  John Haynes (haynes@atmos.colostate.edu)
!
! Modified:
!   09/2006  placed into subroutine form (Roger Marchand,JMH)
!   08/2007  added equivalent volume spheres, Z and N scalling most distrubtion types (Roger Marchand)
!   01/2008  'Do while' to determine if hydrometeor(s) present in volume
!             changed for vectorization purposes (A. Bodas-Salcedo)
!
!   07/2010  V3.0 ... Modified to load or save scale factors to disk as a Look-Up Table (LUT)
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
!
!
! ----- INPUTS -----  
 
  logical, parameter  ::  DO_LUT_TEST = .false. 
  logical, parameter  ::  DO_NP_TEST = .false. 

  type(class_param), intent(inout) :: hp

  integer, intent(in) ::  nprof,ngate

  real undef
  real*8, dimension(nprof,ngate), intent(in) :: &
    hgt_matrix, p_matrix,t_matrix,rh_matrix

  real*8, dimension(hp%nhclass,nprof,ngate), intent(in) :: hm_matrix
  real*8, dimension(hp%nhclass,nprof,ngate), intent(inout) :: re_matrix
  real*8, dimension(hp%nhclass,nprof,ngate), intent(in)    :: Np_matrix

! ----- OUTPUTS -----
  real*8, dimension(nprof,ngate), intent(out) :: Ze_non,Ze_ray, &
       g_to_vol,dBZe,a_to_vol

! ----- OPTIONAL -----
  real*8, optional, dimension(nprof,ngate) :: &
  g_to_vol_in,g_to_vol_out ! integrated atten due to gases, r>v (dB). This allows to output and then input
                           ! the same gaseous absorption in different calls. Optional to allow compatibility
                           ! with original version. A. Bodas April 2008.
!  real*8, dimension(nprof,ngate) :: kr_matrix 

! ----- INTERNAL -----

  real, parameter :: one_third = 1.0/3.0
  real*8 :: t_kelvin
  integer :: &
  phase, & ! 0=liquid, 1=ice
  ns       ! number of discrete drop sizes

  logical :: hydro	! true=hydrometeor in vol, false=none
  real*8 :: &
  rho_a, &   ! air density (kg m^-3)
  gases      ! function: 2-way gas atten (dB/km)

  real*8, dimension(:), allocatable :: &
  Di, Deq, &   ! discrete drop sizes (um)
  Ni, &        ! discrete concentrations (cm^-3 um^-1)
  rhoi         ! discrete densities (kg m^-3)

  real*8, dimension(nprof, ngate) :: &
  z_vol, &      ! effective reflectivity factor (mm^6/m^3)
  z_ray, &                      ! reflectivity factor, Rayleigh only (mm^6/m^3)
  kr_vol, &     ! attenuation coefficient hydro (dB/km)
  g_vol         ! attenuation coefficient gases (dB/km)


  integer,parameter :: KR8 = selected_real_kind(15,300)
  real*8, parameter :: xx = -1.0_KR8
  real*8,  dimension(:), allocatable :: xxa
  real*8 :: kr, ze, zr, pi, scale_factor, tc, Re, ld, tmp1, ze2, kr2, apm, bpm
  real*8 :: half_a_atten_current,half_a_atten_above
  real*8 :: half_g_atten_current,half_g_atten_above
  integer*4 :: tp, i, j, k, pr, itt, iff

  real*8    step,base, Np
  integer*4 iRe_type,n,max_bin

  integer   start_gate,end_gate,d_gate

  logical :: g_to_vol_in_present, g_to_vol_out_present

  ! Logicals to avoid calling present within the loops
  g_to_vol_in_present  = present(g_to_vol_in)
  g_to_vol_out_present = present(g_to_vol_out)

  !
  ! load scaling matricies from disk -- but only the first time this subroutine is called
  !
  if(hp%load_scale_LUTs) then
    call load_scale_LUTs(hp)
    hp%load_scale_LUTs=.false.
    hp%Z_scale_added_flag = .false. ! will be set true if scaling Look Up Tables are modified during run
  endif

  pi = acos(-1.0)

!   ----- Initialisation -----
  g_to_vol = 0.0
  a_to_vol = 0.0
  z_vol    = 0.0
  z_ray    = 0.0
  kr_vol   = 0.0

!   // loop over each range gate (ngate) ... starting with layer closest to the radar !
  if(hp%radar_at_layer_one) then
    start_gate=1
    end_gate=ngate
    d_gate=1
  else
    start_gate=ngate
    end_gate=1
    d_gate=-1
  endif
  do k=start_gate,end_gate,d_gate
  ! // loop over each profile (nprof)
    do pr=1,nprof
      t_kelvin = t_matrix(pr,k)
!     :: determine if hydrometeor(s) present in volume
      hydro = .false.
      do j=1,hp%nhclass
        if ((hm_matrix(j,pr,k) > 1E-12) .and. (hp%dtype(j) > 0)) then
          hydro = .true.
          exit
        endif
      enddo

!     :: if there is hydrometeor in the volume
      if (hydro) then

        rho_a = (p_matrix(pr,k)*100.)/(287.0*(t_kelvin))
!       :: loop over hydrometeor type
        do tp=1,hp%nhclass
          if (hm_matrix(tp,pr,k) <= 1E-12) cycle
          phase = hp%phase(tp)
          if (phase==0) then
            itt = infind(hp%mt_ttl,t_kelvin)
          else
            itt = infind(hp%mt_tti,t_kelvin)
          endif
          if (re_matrix(tp,pr,k).eq.0) then
            call calc_Re(hm_matrix(tp,pr,k),Np_matrix(tp,pr,k),rho_a, &
              hp%dtype(tp),hp%dmin(tp),hp%dmax(tp),hp%apm(tp),hp%bpm(tp), &
              hp%rho(tp),hp%p1(tp),hp%p2(tp),hp%p3(tp),Re)
            re_matrix(tp,pr,k)=Re
          else
            if (Np_matrix(tp,pr,k)>0) then
               print *, 'Warning: Re and Np set for the same ', &
                        'volume & hydrometeor type.  Np is being ignored.'
            endif
            Re = re_matrix(tp,pr,k)
          endif

          iRe_type=1
          if(Re.gt.0) then
            ! determine index in to scale LUT
            !
            ! distance between Re points (defined by "base" and "step") for
            ! each interval of size Re_BIN_LENGTH
            ! Integer asignment, avoids calling floor intrinsic
            n=Re/Re_BIN_LENGTH
            if (n>=Re_MAX_BIN) n=Re_MAX_BIN-1
            step=hp%step_list(n+1)
            base=hp%base_list(n+1)
            iRe_type=Re/step
            if (iRe_type.lt.1) iRe_type=1

            Re=step*(iRe_type+0.5)	! set value of Re to closest value allowed in LUT.
            iRe_type=iRe_type+base-int(n*Re_BIN_LENGTH/step)

            ! make sure iRe_type is within bounds
            if (iRe_type.ge.nRe_types) then
!               write(*,*) 'Warning: size of Re exceed value permitted ', &
!                    'in Look-Up Table (LUT).  Will calculate. '
               ! no scaling allowed
               iRe_type=nRe_types
               hp%Z_scale_flag(tp,itt,iRe_type)=.false.
            else
               ! set value in re_matrix to closest values in LUT
              if (.not. DO_LUT_TEST) re_matrix(tp,pr,k)=Re
            endif
          endif
          ! use Ze_scaled, Zr_scaled, and kr_scaled ... if know them
          ! if not we will calculate Ze, Zr, and Kr from the distribution parameters
          if( (.not. hp%Z_scale_flag(tp,itt,iRe_type)) .or. DO_LUT_TEST)  then
            ! :: create a distribution of hydrometeors within volume
            select case(hp%dtype(tp))
              case(4)
                ns = 1
                allocate(Di(ns),Ni(ns),rhoi(ns),xxa(ns),Deq(ns))
                Di = hp%p1(tp)
                Ni = 0.
              case default
                ns = nd   ! constant defined in radar_simulator_types.f90
                allocate(Di(ns),Ni(ns),rhoi(ns),xxa(ns),Deq(ns))
                Di = hp%D
                Ni = 0.
            end select
            call dsd(hm_matrix(tp,pr,k),re_matrix(tp,pr,k),Np_matrix(tp,pr,k), &
                     Di,Ni,ns,hp%dtype(tp),rho_a,t_kelvin, &
                     hp%dmin(tp),hp%dmax(tp),hp%apm(tp),hp%bpm(tp), &
                     hp%rho(tp),hp%p1(tp),hp%p2(tp),hp%p3(tp))

            ! calculate particle density
            if (phase == 1) then
              if (hp%rho(tp) < 0) then
                ! Use equivalent volume spheres.
                hp%rho_eff(tp,1:ns,iRe_type) = 917  				! solid ice == equivalent volume approach
                Deq = ( ( 6/pi*hp%apm(tp)/917 ) ** (1.0/3.0) ) * ( (Di*1E-6) ** (hp%bpm(tp)/3.0) )  * 1E6
                ! alternative is to comment out above two lines and use the following block
                ! MG Mie approach - adjust density of sphere with D = D_characteristic to match particle density
                !
                ! hp%rho_eff(tp,1:ns,iRe_type) = (6/pi)*hp%apm(tp)*(Di*1E-6)**(hp%bpm(tp)-3)   !MG Mie approach

                ! as the particle size gets small it is possible that the mass to size relationship of 
                ! (given by power law in hclass.data) can produce impossible results 
                ! where the mass is larger than a solid sphere of ice.  
                ! This loop ensures that no ice particle can have more mass/density larger than an ice sphere.
                ! do i=1,ns
                ! if(hp%rho_eff(tp,i,iRe_type) > 917 ) then
                ! hp%rho_eff(tp,i,iRe_type) = 917
                ! endif
                ! enddo
              else
                ! Equivalent volume sphere (solid ice rho_ice=917 kg/m^3).
                hp%rho_eff(tp,1:ns,iRe_type) = 917
                Deq=Di * ((hp%rho(tp)/917)**(1.0/3.0))
                ! alternative ... coment out above two lines and use the following for MG-Mie
                ! hp%rho_eff(tp,1:ns,iRe_type) = hp%rho(tp)   !MG Mie approach
              endif
            else
              ! I assume here that water phase droplets are spheres.
              ! hp%rho should be ~ 1000  or hp%apm=524 .and. hp%bpm=3
              Deq = Di
            endif

            ! calculate effective reflectivity factor of volume
            xxa = -9.9
            rhoi = hp%rho_eff(tp,1:ns,iRe_type)
            call zeff(hp%freq,Deq,Ni,ns,hp%k2,t_kelvin,phase,hp%do_ray, &
                      ze,zr,kr,xxa,xxa,rhoi)

            ! test code ... compare Np value input to routine with sum of DSD
            ! NOTE: if .not. DO_LUT_TEST, then you are checking the LUT approximation 
            ! not just the DSD representation given by Ni
            if(Np_matrix(tp,pr,k)>0 .and. DO_NP_TEST ) then
              Np = path_integral(Ni,Di,1,ns-1)/rho_a*1E6
              ! Note: Representation is not great or small Re < 2 
              if( (Np_matrix(tp,pr,k)-Np)/Np_matrix(tp,pr,k)>0.1 ) then
                write(*,*) 'Error: Np input does not match sum(N)'
                write(*,*) tp,pr,k,Re,Ni(1),Ni(ns),10*log10(ze)
                write(*,*) Np_matrix(tp,pr,k),Np,(Np_matrix(tp,pr,k)-Np)/Np_matrix(tp,pr,k)
                write(*,*)
              endif
            endif

            deallocate(Di,Ni,rhoi,xxa,Deq)

            ! LUT test code
            ! This segment of code compares full calculation to scaling result
            if ( hp%Z_scale_flag(tp,itt,iRe_type) .and. DO_LUT_TEST )  then
              scale_factor=rho_a*hm_matrix(tp,pr,k)
              ! if more than 2 dBZe difference print error message/parameters.
              if ( abs(10*log10(ze) - 10*log10(hp%Ze_scaled(tp,itt,iRe_type) * &
                   scale_factor)) > 2 ) then
                write(*,*) 'Roj Error: ',tp,itt,iRe_type,hp%Z_scale_flag(tp,itt,iRe_type),n,step,base
                write(*,*) 10*log10(ze),10*log10(hp%Ze_scaled(tp,itt,iRe_type) * scale_factor)
                write(*,*) hp%Ze_scaled(tp,itt,iRe_type),scale_factor
                write(*,*) re_matrix(tp,pr,k),Re
                write(*,*)
              endif
            endif

          else ! can use z scaling
            scale_factor=rho_a*hm_matrix(tp,pr,k)
            zr = hp%Zr_scaled(tp,itt,iRe_type) * scale_factor
            ze = hp%Ze_scaled(tp,itt,iRe_type) * scale_factor
            kr = hp%kr_scaled(tp,itt,iRe_type) * scale_factor
          endif  ! end z_scaling

          kr_vol(pr,k) = kr_vol(pr,k) + kr
          z_vol(pr,k)  = z_vol(pr,k)  + ze
          z_ray(pr,k)  = z_ray(pr,k)  + zr

          ! construct Ze_scaled, Zr_scaled, and kr_scaled ... if we can
          if ( .not. hp%Z_scale_flag(tp,itt,iRe_type) ) then
            if (iRe_type>1) then
              scale_factor=rho_a*hm_matrix(tp,pr,k)
              hp%Ze_scaled(tp,itt,iRe_type) = ze/ scale_factor
              hp%Zr_scaled(tp,itt,iRe_type) = zr/ scale_factor
              hp%kr_scaled(tp,itt,iRe_type) = kr/ scale_factor
              hp%Z_scale_flag(tp,itt,iRe_type) = .true.
              hp%Z_scale_added_flag(tp,itt,iRe_type)=.true.
            endif
          endif

        enddo	! end loop of tp (hydrometeor type)

      else
!     :: volume is hydrometeor-free	
        kr_vol(pr,k) = 0
        z_vol(pr,k)  = undef
        z_ray(pr,k)  = undef
      endif

      !     :: attenuation due to hydrometeors between radar and volume
      !
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
          ! add to previous value half of above layer + half of current layer
          a_to_vol(pr,k) = a_to_vol(pr,k+1) + &
              (kr_vol(pr,k+1)+kr_vol(pr,k))*(hgt_matrix(pr,k+1)-hgt_matrix(pr,k))
        else
          a_to_vol(pr,k)= kr_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k-1))
        endif
      endif

      !     :: attenuation due to gaseous absorption between radar and volume
      if (g_to_vol_in_present) then
        g_to_vol(pr,k) = g_to_vol_in(pr,k)
      else
        if ( (hp%use_gas_abs == 1) .or. ((hp%use_gas_abs == 2) .and. (pr == 1)) ) then
          g_vol(pr,k) = gases(p_matrix(pr,k),t_kelvin,rh_matrix(pr,k),hp%freq)
          if (d_gate==1) then
            if (k>1) then
              ! add to previous value to half of above layer + half of current layer
              g_to_vol(pr,k) =  g_to_vol(pr,k-1) + &
                  0.5*(g_vol(pr,k-1)+g_vol(pr,k))*(hgt_matrix(pr,k-1)-hgt_matrix(pr,k))
            else
              g_to_vol(pr,k)=  0.5*g_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k+1))
            endif
          else   ! d_gate==-1
            if (k<ngate) then
              ! add to previous value to half of above layer + half of current layer
              g_to_vol(pr,k) = g_to_vol(pr,k+1) + &
                 0.5*(g_vol(pr,k+1)+g_vol(pr,k))*(hgt_matrix(pr,k+1)-hgt_matrix(pr,k))
            else
              g_to_vol(pr,k)= 0.5*g_vol(pr,k)*(hgt_matrix(pr,k)-hgt_matrix(pr,k-1))
            endif
          endif
        elseif(hp%use_gas_abs == 2) then
          ! using value calculated for the first column
          g_to_vol(pr,k) = g_to_vol(1,k)
        elseif (hp%use_gas_abs == 0) then
          g_to_vol(pr,k) = 0
        endif
      endif

      ! Compute Rayleigh reflectivity, and full, attenuated reflectivity
      if ((hp%do_ray == 1) .and. (z_ray(pr,k) > 0)) then
        Ze_ray(pr,k) = 10*log10(z_ray(pr,k))
      else
        Ze_ray(pr,k) = undef
      endif
      if (z_vol(pr,k) > 0) then
        Ze_non(pr,k) = 10*log10(z_vol(pr,k))
        dBZe(pr,k) = Ze_non(pr,k)-a_to_vol(pr,k)-g_to_vol(pr,k)
      else
        dBZe(pr,k) = undef
        Ze_non(pr,k) = undef
      endif

    enddo   ! end loop over pr (profile)

  enddo ! end loop of k (range gate)

  ! Output array with gaseous absorption
  if (g_to_vol_out_present) g_to_vol_out = g_to_vol

  ! save any updates made 
  if (hp%update_scale_LUTs) call save_scale_LUTs(hp)

  end subroutine radar_simulator
