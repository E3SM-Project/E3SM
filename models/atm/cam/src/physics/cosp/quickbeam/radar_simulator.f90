  subroutine radar_simulator(freq,k2,do_ray,use_gas_abs,use_mie_table,mt, &
    nhclass,hp,nprof,ngate,nsizes,D,hgt_matrix,hm_matrix,re_matrix,p_matrix,t_matrix, &
    rh_matrix,Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe, &
    g_to_vol_in,g_to_vol_out)

!     rh_matrix,Ze_non,Ze_ray,kr_matrix,g_atten_to_vol,dBZe)
 
  use m_mrgrnk 
  use array_lib
  use math_lib
  use optics_lib
  use radar_simulator_types
  implicit none
  
! Purpose:
!   Simulates a vertical profile of radar reflectivity
!   Part of QuickBeam v1.04 by John Haynes & Roger Marchand
!
! Inputs:
!   [freq]            radar frequency (GHz), can be anything unless
!                     use_mie_table=1, in which case one of 94,35,13.8,9.6,3
!   [k2]              |K|^2, the dielectric constant, set to -1 to use the
!                     frequency dependent default
!   [do_ray]          1=do Rayleigh calcs, 0=not
!   [use_gas_abs]     1=do gaseous abs calcs, 0=not,
!                     2=use same as first profile (undocumented)
!   [use_mie_table]   1=use Mie tables, 0=not
!   [mt]              Mie look up table
!   [nhclass]         number of hydrometeor types
!   [hp]              structure that defines hydrometeor types
!   [nprof]           number of hydrometeor profiles
!   [ngate]           number of vertical layers
!   [nsizes]          number of discrete particles in [D]
!   [D]               array of discrete particles (um)
!
!   (The following 5 arrays must be in order from closest to the radar
!    to farthest...)
!   [hgt_matrix]      height of hydrometeors (km)
!   [hm_matrix]       table of hydrometeor mixing rations (g/kg)
!   [re_matrix]       OPTIONAL table of hydrometeor effective radii (microns)
!   [p_matrix]        pressure profile (hPa)
!   [t_matrix]        temperature profile (C)
!   [rh_matrix]       relative humidity profile (%)
!
! Outputs:
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
! Modified:
!   09/2006  placed into subroutine form, scaling factors (Roger Marchand,JMH)
!   08/2007  added equivalent volume spheres, Z and N scalling most distrubtion types (Roger Marchand)
!   01/2008  'Do while' to determine if hydrometeor(s) present in volume
!             changed for vectorization purposes (A. Bodas-Salcedo)

! ----- INPUTS -----  
  type(mie), intent(in) :: mt
  type(class_param), intent(inout) :: hp
!!+jek
!!  real*8, intent(in) :: freq,k2
  real*8, intent(in) :: freq
  real*8, intent(inout) :: k2
!!-jek
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
  g_to_vol_in,g_to_vol_out ! integrated atten due to gases, r>v (dB). This allows to output and then input
                           ! the same gaseous absorption in different calls. Optional to allow compatibility
                           ! with original version. A. Bodas April 2008.
        
!  real*8, dimension(nprof,ngate) :: kr_matrix 

! ----- INTERNAL -----
  integer :: &
  phase, &               ! 0=liquid, 1=ice
  ns                     ! number of discrete drop sizes

  integer*4, dimension(ngate) :: &
  hydro                  ! 1=hydrometeor in vol, 0=none
  real*8 :: &
  rho_a, &               ! air density (kg m^-3)
  gases                  ! function: 2-way gas atten (dB/km)

  real*8, dimension(:), allocatable :: &
  Di, Deq, &             ! discrete drop sizes (um)
  Ni, Ntemp, &           ! discrete concentrations (cm^-3 um^-1)
  rhoi                   ! discrete densities (kg m^-3)
  
  real*8, dimension(ngate) :: &
  z_vol, &               ! effective reflectivity factor (mm^6/m^3)
  z_ray, &                      ! reflectivity factor, Rayleigh only (mm^6/m^3)
  kr_vol, &              ! attenuation coefficient hydro (dB/km)
  g_vol, &               ! attenuation coefficient gases (dB/km)
  a_to_vol, &            ! integrated atten due to hydometeors, r>v (dB)
  g_to_vol               ! integrated atten due to gases, r>v (dB)
   
 
  integer,parameter :: KR8 = selected_real_kind(15,300)
  real*8, parameter :: xx = -1.0_KR8
  real*8, dimension(:), allocatable :: xxa
  real*8 :: kr, ze, zr, pi, scale_factor, Re, ld, tmp1, apm,bpm
  integer*4 :: tp, i, j, k, pr, itt, iff

  real*8 bin_length,step,base,step_list(25),base_list(25)
  integer*4 iRe_type,n,max_bin
  
  logical :: g_to_vol_in_present, g_to_vol_out_present
     
  ! Logicals to avoid calling present within the loops
  g_to_vol_in_present  = present(g_to_vol_in)
  g_to_vol_out_present = present(g_to_vol_out)
  
    ! set up Re bins for z_scalling
     bin_length=50;
     max_bin=25

     step_list(1)=1
     base_list(1)=75 
     do j=2,max_bin
          step_list(j)=3*(j-1);
          if(step_list(j)>bin_length) then
               step_list(j)=bin_length;
          endif
          base_list(j)=base_list(j-1)+floor(bin_length/step_list(j-1));
     enddo


  pi = acos(-1.0)
  if (use_mie_table == 1) iff = infind(mt%freq,freq,sort=1)

     
  ! // loop over each profile (nprof)
  do pr=1,nprof

!   ----- calculations for each volume ----- 
    z_vol(:) = 0
    z_ray(:) = 0
    kr_vol(:) = 0
    hydro(:) = 0    

!   // loop over eacho range gate (ngate)
    do k=1,ngate
  
!     :: determine if hydrometeor(s) present in volume
      hydro(k) = 0
      do j=1,nhclass ! Do while changed for vectorization purposes (A. B-S)
        if ((hm_matrix(j,pr,k) > 1E-12) .and. (hp%dtype(j) > 0)) then
          hydro(k) = 1
          exit
        endif
      enddo

      if (hydro(k) == 1) then
!     :: if there is hydrometeor in the volume            

        rho_a = (p_matrix(pr,k)*100.)/(287*(t_matrix(pr,k)+273.15))

!       :: loop over hydrometeor type
        do tp=1,nhclass

          if (hm_matrix(tp,pr,k) <= 1E-12) cycle

       phase = hp%phase(tp)
       if(phase==0) then
          itt = infind(mt_ttl,t_matrix(pr,k))
       else
          itt = infind(mt_tti,t_matrix(pr,k))
      endif

       ! calculate Re if we have an exponential distribution with fixed No ... precipitation type particle
       if( hp%dtype(tp)==2 .and. abs(hp%p2(tp)+1) < 1E-8)  then

          apm=hp%apm(tp)
          bpm=hp%bpm(tp)

          if ((hp%rho(tp) > 0) .and. (apm < 0)) then
               apm = (pi/6)*hp%rho(tp)
               bpm = 3.
          endif

          tmp1 = 1./(1.+bpm)
          ld = ((apm*gamma(1.+bpm)*hp%p1(tp))/(rho_a*hm_matrix(tp,pr,k)*1E-3))**tmp1
          
          Re = 1.5E6/ld 
          
          re_matrix(tp,pr,k) = Re;

       endif
  
       if(re_matrix(tp,pr,k).eq.0) then

          iRe_type=1
          Re=0
       else
          iRe_type=1
          Re=re_matrix(tp,pr,k)
          
          n=floor(Re/bin_length)
          if(n==0) then
               if(Re<25) then
                    step=0.5
                    base=0
               else           
                    step=1
                    base=25
               endif
          else
               if(n>max_bin) then
                    n=max_bin 
               endif

               step=step_list(n)
               base=base_list(n)
          endif

          iRe_type=floor(Re/step)

          if(iRe_type.lt.1) then  
               iRe_type=1               
          endif

          Re=step*(iRe_type+0.5)
          iRe_type=iRe_type+base-floor(n*bin_length/step)

          ! make sure iRe_type is within bounds
          if(iRe_type.ge.nRe_types) then  

               ! print *, tp, re_matrix(tp,pr,k), Re, iRe_type

               ! no scaling allowed
               Re=re_matrix(tp,pr,k)

               iRe_type=nRe_types
               hp%z_flag(tp,itt,iRe_type)=.false.
               hp%scaled(tp,iRe_type)=.false.               
          endif
       endif
     
       ! use Ze_scaled, Zr_scaled, and kr_scaled ... if know them
       ! if not we will calculate Ze, Zr, and Kr from the distribution parameters
       if( .not. hp%z_flag(tp,itt,iRe_type) )  then
      
!         :: create a distribution of hydrometeors within volume   
       select case(hp%dtype(tp))
          case(4)
         ns = 1
         allocate(Di(ns),Ni(ns),rhoi(ns),xxa(ns),Deq(ns))
         if (use_mie_table == 1) allocate(mt_qext(ns),mt_qbsca(ns),Ntemp(ns))
         Di = hp%p1(tp)
         Ni = 0.
       case default
         ns = nsizes            
         allocate(Di(ns),Ni(ns),rhoi(ns),xxa(ns),Deq(ns))
         if (use_mie_table == 1) allocate(mt_qext(ns),mt_qbsca(ns),Ntemp(ns))       
         Di = D
         Ni = 0.
       end select

!         :: create a DSD (using scaling factor if applicable)
       ! hp%scaled(tp,iRe_type)=.false.   ! turn off N scaling

       call dsd(hm_matrix(tp,pr,k),Re,Di,Ni,ns,hp%dtype(tp),rho_a, &
         t_matrix(pr,k),hp%dmin(tp),hp%dmax(tp),hp%apm(tp),hp%bpm(tp), &
         hp%rho(tp),hp%p1(tp),hp%p2(tp),hp%p3(tp),hp%fc(tp,1:ns,iRe_type), &
         hp%scaled(tp,iRe_type))

!         :: calculate particle density 
          ! if ((hp%rho_eff(tp,1,iRe_type) < 0) .and. (phase == 1)) then
       if (phase == 1) then
         if (hp%rho(tp) < 0) then
                
          ! MG Mie approach - adjust density of sphere with D = D_characteristic to match particle density         
          ! hp%rho_eff(tp,1:ns,iRe_type) = (6/pi)*hp%apm(tp)*(Di*1E-6)**(hp%bpm(tp)-3)   !MG Mie approach
          
          ! as the particle size gets small it is possible that the mass to size relationship of 
          ! (given by power law in hclass.data) can produce impossible results 
          ! where the mass is larger than a solid sphere of ice.  
          ! This loop ensures that no ice particle can have more mass/density larger than an ice sphere.
          ! do i=1,ns
          ! if(hp%rho_eff(tp,i,iRe_type) > 917 ) then
          !    hp%rho_eff(tp,i,iRe_type) = 917
          !endif
          !enddo

          ! alternative is to use equivalent volume spheres.
          hp%rho_eff(tp,1:ns,iRe_type) = 917                     ! solid ice == equivalent volume approach
               Deq = ( ( 6/pi*hp%apm(tp)/917 ) ** (1.0/3.0) ) * &
                  ( (Di*1E-6) ** (hp%bpm(tp)/3.0) )  * 1E6       ! Di now really Deq in microns.
          
            else

               ! hp%rho_eff(tp,1:ns,iRe_type) = hp%rho(tp)   !MG Mie approach
               
          ! Equivalent volume sphere (solid ice rho_ice=917 kg/m^3).
               hp%rho_eff(tp,1:ns,iRe_type) = 917
               Deq=Di * ((hp%rho(tp)/917)**(1.0/3.0))  

         endif

          ! if using equivalent volume spheres
          if (use_mie_table == 1) then

               Ntemp=Ni

               ! Find N(Di) from N(Deq) which we know
               do i=1,ns
                              j=infind(Deq,Di(i))
                    Ni(i)=Ntemp(j)
               enddo
          else
               ! just use Deq and D variable input to mie code
               Di=Deq;
          endif

       endif
       rhoi = hp%rho_eff(tp,1:ns,iRe_type)
       
!         :: calculate effective reflectivity factor of volume
       if (use_mie_table == 1) then
       
         if ((hp%dtype(tp) == 4) .and. (hp%idd(tp) < 0)) then
              hp%idd(tp) = infind(mt%D,Di(1))
         endif
         
         if (phase == 0) then
         
           ! itt = infind(mt_ttl,t_matrix(pr,k))
              select case(hp%dtype(tp))
           case(4)
          mt_qext(1) = mt%qext(hp%idd(tp),itt,1,iff)
             mt_qbsca(1) = mt%qbsca(hp%idd(tp),itt,1,iff)
              case default
             mt_qext = mt%qext(:,itt,1,iff)
             mt_qbsca = mt%qbsca(:,itt,1,iff)
           end select

          call zeff(freq,Di,Ni,ns,k2,mt_ttl(itt),0,do_ray, &
             ze,zr,kr,mt_qext,mt_qbsca,spread(xx,1,ns))
         
         else

           ! itt = infind(mt_tti,t_matrix(pr,k))
           select case(hp%dtype(tp))
           case(4)
                if (hp%ifc(tp,1,iRe_type) < 0) then
                  hp%ifc(tp,1,iRe_type) = infind(mt%f,rhoi(1)/917.)
             endif             
                mt_qext(1) = &
            mt%qext(hp%idd(tp),itt+cnt_liq,hp%ifc(tp,1,iRe_type),iff)
             mt_qbsca(1) = &
            mt%qbsca(hp%idd(tp),itt+cnt_liq,hp%ifc(tp,1,iRe_type),iff)           
           case default
             do i=1,ns
               if (hp%ifc(tp,i,iRe_type) < 0) then
                    hp%ifc(tp,i,iRe_type) = infind(mt%f,rhoi(i)/917.)
               endif           
                    mt_qext(i) = mt%qext(i,itt+cnt_liq,hp%ifc(tp,i,iRe_type),iff)
            mt_qbsca(i) = mt%qbsca(i,itt+cnt_liq,hp%ifc(tp,i,iRe_type),iff)
             enddo
           end select

             call zeff(freq,Di,Ni,ns,k2,mt_tti(itt),1,do_ray, &
             ze,zr,kr,mt_qext,mt_qbsca,spread(xx,1,ns))

         endif

       else
       
         xxa = -9.9
         call zeff(freq,Di,Ni,ns,k2,t_matrix(pr,k),phase,do_ray, &
           ze,zr,kr,xxa,xxa,rhoi)

           
       endif  ! end of use mie table 

          ! xxa = -9.9
          !call zeff(freq,Di,Ni,ns,k2,t_matrix(pr,k),phase,do_ray, &
               !    ze2,zr,kr2,xxa,xxa,rhoi)

          ! if(abs(ze2-ze)/ze2 > 0.1) then
          ! if(abs(kr2-kr)/kr2 > 0.1) then
          
          ! write(*,*) pr,k,tp,ze2,ze2-ze,abs(ze2-ze)/ze2,itt+cnt_liq,iff
          ! write(*,*) pr,k,tp,ze2,kr2,kr2-kr,abs(kr2-kr)/kr2
          ! stop

          !endif

       deallocate(Di,Ni,rhoi,xxa,Deq)
       if (use_mie_table == 1) deallocate(mt_qext,mt_qbsca,Ntemp)

       else ! can use z scaling
       
          if( hp%dtype(tp)==2 .and. abs(hp%p2(tp)+1) < 1E-8 )  then
           
               ze = hp%Ze_scaled(tp,itt,iRe_type)
               zr = hp%Zr_scaled(tp,itt,iRe_type)
               kr = hp%kr_scaled(tp,itt,iRe_type)

          else
               scale_factor=rho_a*hm_matrix(tp,pr,k) 

               zr = hp%Zr_scaled(tp,itt,iRe_type) * scale_factor 
               ze = hp%Ze_scaled(tp,itt,iRe_type) * scale_factor
               kr = hp%kr_scaled(tp,itt,iRe_type) * scale_factor 
          endif

       endif  ! end z_scaling
 
       ! kr=0 

       kr_vol(k) = kr_vol(k) + kr
       z_vol(k) = z_vol(k) + ze
       z_ray(k) = z_ray(k) + zr
     
       ! construct Ze_scaled, Zr_scaled, and kr_scaled ... if we can
       if( .not. hp%z_flag(tp,itt,iRe_type) .and. 1.eq.1 ) then

          if( ( (hp%dtype(tp)==1 .or. hp%dtype(tp)==5 .or.  hp%dtype(tp)==2)  .and. abs(hp%p1(tp)+1) < 1E-8  ) .or. &
              (  hp%dtype(tp)==3 .or. hp%dtype(tp)==4 )  &
          ) then

               scale_factor=rho_a*hm_matrix(tp,pr,k) 

               hp%Ze_scaled(tp,itt,iRe_type) = ze/ scale_factor
               hp%Zr_scaled(tp,itt,iRe_type) = zr/ scale_factor
               hp%kr_scaled(tp,itt,iRe_type) = kr/ scale_factor

               hp%z_flag(tp,itt,iRe_type)=.True.

          elseif( hp%dtype(tp)==2 .and. abs(hp%p2(tp)+1) < 1E-8 ) then 
           
               hp%Ze_scaled(tp,itt,iRe_type) = ze
               hp%Zr_scaled(tp,itt,iRe_type) = zr
               hp%kr_scaled(tp,itt,iRe_type) = kr

               hp%z_flag(tp,itt,iRe_type)=.True.
          endif

       endif

        enddo  ! end loop of tp (hydrometeor type)

      else
!     :: volume is hydrometeor-free
     
        kr_vol(k) = 0
     z_vol(k) = -999
        z_ray(k) = -999
     
      endif

!     :: attenuation due to hydrometeors between radar and volume
      a_to_vol(k) = 2*path_integral(kr_vol,hgt_matrix(pr,:),1,k-1)
      
!     :: attenuation due to gaseous absorption between radar and volume
      if (g_to_vol_in_present) then
        g_to_vol(k) = g_to_vol_in(k,pr)
      else
        if ( (use_gas_abs == 1) .or. ((use_gas_abs == 2) .and. (pr == 1)) )  then
            g_vol(k) = gases(p_matrix(pr,k),t_matrix(pr,k)+273.15, &
            rh_matrix(pr,k),freq)
            g_to_vol(k) = path_integral(g_vol,hgt_matrix(pr,:),1,k-1)
        elseif (use_gas_abs == 0) then
            g_to_vol(k) = 0
        endif  
      endif
    
!      kr_matrix(pr,:)=kr_vol

!     :: store results in matrix for return to calling program
      h_atten_to_vol(pr,k)=a_to_vol(k)
      g_atten_to_vol(pr,k)=g_to_vol(k)
      if ((do_ray == 1) .and. (z_ray(k) > 0)) then
        Ze_ray(pr,k) = 10*log10(z_ray(k))
      else
        Ze_ray(pr,k) = -999
      endif
      if (z_vol(k) > 0) then
        dBZe(pr,k) = 10*log10(z_vol(k))-a_to_vol(k)-g_to_vol(k)
        Ze_non(pr,k) = 10*log10(z_vol(k))
      else
        dBZe(pr,k) = -999
        Ze_non(pr,k) = -999
      endif
      
    enddo ! end loop of k (range gate)
    ! Output array with gaseous absorption
    if (g_to_vol_out_present) g_to_vol_out(:,pr) = g_to_vol
  enddo        ! end loop over pr (profile)  

  end subroutine radar_simulator
  
