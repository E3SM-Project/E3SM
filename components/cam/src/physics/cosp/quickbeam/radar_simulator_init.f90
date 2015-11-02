  subroutine radar_simulator_init(freq,k2,use_gas_abs,do_ray,undef, &
                  nhclass, &
                  hclass_type,hclass_phase, &
                      hclass_dmin,hclass_dmax, &
                          hclass_apm,hclass_bpm,hclass_rho, &
                          hclass_p1,hclass_p2,hclass_p3, &
                          load_scale_LUTs_flag,update_scale_LUTs_flag,LUT_file_name, &
                  hp &      ! output
                  )
  use radar_simulator_types
  implicit none
  
! Purpose:
!
!   Initialize variables used by the radar simulator.
!   Part of QuickBeam v3.0 by John Haynes and Roj Marchand
!   
!
! Inputs:  
!   []   from data in hydrometeor class input 
! 
!   [freq]            radar frequency (GHz)
!
!   [k2]              |K|^2, the dielectric constant, set to -1 to use the
!                     frequency dependent default
!
!   [use_gas_abs]     1=do gaseous abs calcs, 0=no gasesous absorbtion calculated,
!                     2=calculate (and use) absorption for first profile on all profiles
!
!   [undef]           mising data value
!   [nhclass]         number of hydrometeor types
!
!   For each hydrometero type:
!       hclass_type     Type of distribution (see quickbeam documentation)
!       hclass_phase            1==ice, 0=liquid
!
!       hclass_dmin         minimum diameter allowed is drop size distribution N(D<Dmin)=0
!       hclass_dmax         maximum diameter allowed is drop size distribution N(D>Dmax)=0
!
!   hclass_apm,hclass_bpm   Density of partical apm*D^bpm or constant = rho
!   hclass_rho, 
!   
!       hclass_p1,hclass_p2,    Default values of DSD parameters (see quickbeam documentation)
!   hclass_p3, 
!
!   load_scale_LUTs_flag    Flag, load scale factors Look Up Table from file at start of run
!   update_scale_LUTs_flag  Flag, save new scale factors calculated during this run to LUT
!   LUT_file_name       Name of file containing LUT
!
! Outputs:
!   [hp]            structure that define hydrometeor types
!
! Modified:
!   08/23/2006  placed into subroutine form (Roger Marchand)
!   June 2010   New interface to support "radar_simulator_params" structure
   
! ----- INPUT -----

   real, intent(in)    :: freq,k2
   integer, intent(in) :: nhclass       ! number of hydrometeor classes in use
   integer, intent(in) :: use_gas_abs,do_ray
   real :: undef
   real,dimension(nhclass),intent(in)     ::    hclass_dmin,hclass_dmax, &
                                hclass_apm,hclass_bpm,hclass_rho, &
                                hclass_p1,hclass_p2,hclass_p3 
   integer,dimension(nhclass),intent(in)  ::    hclass_type,hclass_phase
  
   logical, intent(in)       :: load_scale_LUTs_flag,update_scale_LUTs_flag
   character*240, intent(in) :: LUT_file_name

! ----- OUTPUTS -----  
  type(class_param), intent(out) :: hp

! ----- INTERNAL -----  
  integer :: i,j
  real*8  :: delt, deltp
        
    !
    ! set radar simulation properites
    !
    hp%freq=freq
    hp%k2=k2
    hp%use_gas_abs=use_gas_abs
    hp%do_ray=do_ray
    hp%nhclass=nhclass
    
    hp%load_scale_LUTs=load_scale_LUTs_flag
    hp%update_scale_LUTs=update_scale_LUTs_flag
    hp%scale_LUT_file_name=LUT_file_name
    
    ! 
        ! Store settings for hydrometeor types in hp (class_parameter) structure.   
        !
        do i = 1,nhclass
        hp%dtype(i) = hclass_type(i)
        hp%phase(i) = hclass_phase(i)
        hp%dmin(i)  = hclass_dmin(i)
        hp%dmax(i)  = hclass_dmax(i)
        hp%apm(i)   = hclass_apm(i)
        hp%bpm(i)   = hclass_bpm(i)
        hp%rho(i)   = hclass_rho(i)
        hp%p1(i)    = hclass_p1(i)
        hp%p2(i)    = hclass_p2(i)
        hp%p3(i)    = hclass_p3(i)
        enddo
   
        ! 
        ! initialize scaling array
        !
        hp%N_scale_flag = .false.
        hp%fc = undef
        hp%rho_eff = undef
    
        hp%Z_scale_flag = .false.
        hp%Ze_scaled = 0.0
        hp%Zr_scaled = 0.0
        hp%kr_scaled = 0.0
  
        !
    ! set up Re bin "structure" for z_scaling
        !
    hp%base_list(1)=0;
    do j=1,Re_MAX_BIN
        hp%step_list(j)=0.1+0.1*((j-1)**1.5);
        if(hp%step_list(j)>Re_BIN_LENGTH) then
            hp%step_list(j)=Re_BIN_LENGTH;
        endif
        if(j>1) then
            hp%base_list(j)=hp%base_list(j-1)+floor(Re_BIN_LENGTH/hp%step_list(j-1));
        endif
    enddo

    !
    ! set up Temperature bin structure used for z scaling
    !
    do i=1,cnt_ice
        hp%mt_tti(i)=(i-1)*5-90 + 273.15
    enddo
    
    do i=1,cnt_liq
        hp%mt_ttl(i)=(i-1)*5-60 + 273.15
    enddo 
  
    !
        ! define array discrete diameters used in mie calculations
        !
        delt = (log(dmax)-log(dmin))/(nd-1)
        deltp = exp(delt)

        hp%D(1) = dmin
        do i=2,nd
          hp%D(i) = hp%D(i-1)*deltp
        enddo   
 
  
  end subroutine radar_simulator_init
