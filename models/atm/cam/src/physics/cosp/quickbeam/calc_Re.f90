  subroutine calc_Re(Q,Np,rho_a, &
             dtype,dmin,dmax,apm,bpm,rho_c,p1,p2,p3, &
             Re)
  use math_lib 
  implicit none

! Purpose:
!   Calculates Effective Radius (1/2 distribution 3rd moment / 2nd moment). 
!
!   For some distribution types, the total number concentration (per kg), Np
!   may be optionally specified.   Should be set to zero, otherwise.
!
!   Roj Marchand July 2010


! Inputs:
!
!   [Q]        hydrometeor mixing ratio (g/kg)  ! not needed for some distribution types
!   [Np]       Optional Total number concentration (per kg).  0 = use defaults (p1, p2, p3)
!
!   [rho_a]    ambient air density (kg m^-3)   
!
!   Distribution parameters as per quickbeam documentation.
!   [dtype]    distribution type
!   [dmin]     minimum size cutoff (um)
!   [dmax]     maximum size cutoff (um)
!   [apm]      a parameter for mass (kg m^[-bpm])
!   [bmp]      b params for mass 
!   [p1],[p2],[p3]  distribution parameters
!
!
! Outputs:
!   [Re]       Effective radius, 1/2 the 3rd moment/2nd moment (um)
!
! Created:
!   July 2010  Roj Marchand
!

 
! ----- INPUTS -----  
  
  real*8, intent(in) :: Q,Np,rho_a
 
  integer, intent(in):: dtype
  real*8, intent(in) :: dmin,dmax,rho_c,p1,p2,p3
    
  real*8, intent(inout) :: apm,bpm  
    
! ----- OUTPUTS -----

  real*8, intent(out) :: Re
  
! ----- INTERNAL -----
 
  integer :: local_dtype
  real*8  :: local_p3,local_Np

  real*8 :: pi, &
  N0,D0,vu,dm,ld, &         ! gamma, exponential variables
  rg,log_sigma_g
  
  real*8 :: tmp1,tmp2

  pi = acos(-1.0)

  ! // if density is constant, set equivalent values for apm and bpm
  if ((rho_c > 0) .and. (apm < 0)) then
    apm = (pi/6)*rho_c
    bpm = 3.
  endif

  ! Exponential is same as modified gamma with vu =1
  ! if Np is specified then we will just treat as modified gamma
  if(dtype.eq.2 .and. Np>0) then
    local_dtype=1;
    local_p3=1;
  else
    local_dtype=dtype;
    local_p3=p3;
  endif
  
  select case(local_dtype)
  
! ---------------------------------------------------------!
! // modified gamma                                        !
! ---------------------------------------------------------!
! :: Np = total number concentration (1/kg) = Nt / rho_a
! :: D0 = characteristic diameter (um)
! :: dm = mean diameter (um) - first moment over zeroth moment
! :: vu = distribution width parameter 

  case(1)  
  
    if( abs(local_p3+2) < 1E-8) then
  
    if(Np>1E-30) then
        ! Morrison scheme with Martin 1994 shape parameter (NOTE: vu = pc +1)
        ! fixed Roj. Dec. 2010 -- after comment by S. Mcfarlane
        vu = (1/(0.2714 + 0.00057145*Np*rho_a*1E-6))**2 ! units of Nt = Np*rhoa = #/cm^3
    else
        print *, 'Error: Must specify a value for Np in each volume', &
             ' with Morrison/Martin Scheme.'
            stop    
    endif
    
    elseif (abs(local_p3+1) > 1E-8) then

      ! vu is fixed in hp structure  
      vu = local_p3 

    else

      ! vu isn't specified
      
      print *, 'Error: Must specify a value for vu for Modified Gamma distribution'
      stop    
      
    endif
    

    if( Np.eq.0 .and. p2+1 > 1E-8) then     ! use default value for MEAN diameter as first default
      
        dm = p2             ! by definition, should have units of microns
    D0 = gamma(vu)/gamma(vu+1)*dm
        
    else   ! use value of Np
        
        if(Np.eq.0) then
        
            if( abs(p1+1) > 1E-8 ) then  !   use default number concentration 
            
                local_Np = p1 ! total number concentration / pa --- units kg^-1
            else
            print *, 'Error: Must specify Np or default value ', &
                 '(p1=Dm [um] or p2=Np [1/kg]) for ', &
                 'Modified Gamma distribution'
                stop
            endif
        else
            local_Np=Np;    
    endif
    
    D0 = 1E6 * ( Q*1E-3*gamma(vu)/(apm*local_Np*gamma(vu+bpm)) )**(1/bpm)  ! units = microns

    endif  
      
    Re = 0.5*D0*gamma(vu+3)/gamma(vu+2)
    
    
! ---------------------------------------------------------!
! // exponential                                           !
! ---------------------------------------------------------!
! :: N0 = intercept parameter (m^-4)
! :: ld = slope parameter (um)

  case(2)
  
    ! Np not specified (see if statement above) 
  
    if((abs(p1+1) > 1E-8) ) then   ! N0 has been specified, determine ld
    
        N0 = p1
    tmp1 = 1./(1.+bpm)
    ld = ((apm*gamma(1.+bpm)*N0)/(rho_a*Q*1E-3))**tmp1
    ld = ld/1E6                     ! set units to microns^-1
        
    elseif (abs(p2+1) > 1E-8) then  ! lambda=ld has been specified as default

        ld = p2     ! should have units of microns^-1 
    
    else
    
    print *, 'Error: Must specify Np or default value ', &
         '(p1=No or p2=lambda) for Exponential distribution'
        stop
        
    endif

    Re = 1.5/ld ;
  
! ---------------------------------------------------------!
! // power law                                             !
! ---------------------------------------------------------!
! :: ahp = Ar parameter (m^-4 mm^-bhp)
! :: bhp = br parameter
! :: dmin_mm = lower bound (mm)
! :: dmax_mm = upper bound (mm)

  case(3)

    Re=0;  ! Not supporting LUT approach for power-law ...

    if(Np>0) then
    
        print *, 'Variable Np not supported for ', &
         'Power Law distribution'
        stop
    endif
! ---------------------------------------------------------!
! // monodisperse                                          !
! ---------------------------------------------------------!
! :: D0 = particle diameter (um) == Re

  case(4)
  
        Re = p1
    
        if(Np>0) then
        print *, 'Variable Np not supported for ', &
         'Monodispersed distribution'
        stop
    endif
    
! ---------------------------------------------------------!
! // lognormal                                             !
! ---------------------------------------------------------!
! :: N0 = total number concentration (m^-3)
! :: np = fixed number concentration (kg^-1)
! :: rg = mean radius (um)
! :: log_sigma_g = ln(geometric standard deviation)

  case(5)
  
        if( abs(local_p3+1) > 1E-8 ) then

            !set natural log width
            log_sigma_g = local_p3 
        else
            print *, 'Error: Must specify a value for sigma_g ', &
             'when using a Log-Normal distribution'
            stop
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
                print *, 'Error: Must specify Np or default value ', &
                 '(p2=Rg or p1=Np) for Log-Normal distribution'
            endif
           
            log_sigma_g = p3
            tmp1 = (Q*1E-3)/(2.**bpm*apm*local_Np)
            tmp2 = exp(0.5*bpm**2.*(log_sigma_g))**2.      
            rg = ((tmp1/tmp2)**(1/bpm))*1E6
        endif
                
        Re = rg*exp(+2.5*(log_sigma_g**2)) 
        
  end select
  
  end subroutine calc_Re
