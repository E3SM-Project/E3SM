module hetfrz_classnuc

!-----------------------------------------------------------------------
!
! Purpose: Calculate heterogeneous freezing rates from classical nucleation theory
!
! Public interfaces:
!
! hetfrz_classnuc_init
! hetfrz_classnuc_calc
!
! Author: 
!   Corinna Hoose, UiO, May 2009
!   Yong Wang and Xiaohong Liu, UWyo, 12/2012, 
!   implement in CAM5 and constrain uncertain parameters using natural dust and
!   BC(soot) datasets. 
!   Yong Wang and Xiaohong Liu, UWyo, 05/2013, implement the PDF-contact angle
!   approach: Y. Wang et al., Atmos. Chem. Phys., 2014.
!   Jack Chen, NCAR, 09/2015, modify calculation of dust activation fraction.
!
!-----------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use wv_saturation, only: svp_water, svp_ice
use shr_spfn_mod,  only: erf => shr_spfn_erf

implicit none
private
save

public :: hetfrz_classnuc_init, hetfrz_classnuc_calc

real(r8) :: rair
real(r8) :: cpair
real(r8) :: rh2o
real(r8) :: rhoh2o
real(r8) :: mwh2o
real(r8) :: tmelt
real(r8) :: pi

integer  :: iulog

!===================================================================================================
contains
!===================================================================================================

subroutine hetfrz_classnuc_init( &
   rair_in, cpair_in, rh2o_in, rhoh2o_in, mwh2o_in, &
   tmelt_in, pi_in, iulog_in)

   real(r8), intent(in) :: rair_in
   real(r8), intent(in) :: cpair_in
   real(r8), intent(in) :: rh2o_in
   real(r8), intent(in) :: rhoh2o_in
   real(r8), intent(in) :: mwh2o_in
   real(r8), intent(in) :: tmelt_in
   real(r8), intent(in) :: pi_in
   integer,  intent(in) :: iulog_in

   rair   = rair_in
   cpair  = cpair_in
   rh2o   = rh2o_in
   rhoh2o = rhoh2o_in
   mwh2o  = mwh2o_in
   tmelt  = tmelt_in
   pi     = pi_in
   iulog  = iulog_in

end subroutine hetfrz_classnuc_init

!===================================================================================================

subroutine hetfrz_classnuc_calc( &
   deltat, t, p, supersatice,                 &
   fn,                                        &
   r3lx, icnlx,                               &
   frzbcimm, frzduimm,                        &
   frzbccnt, frzducnt,                        &
   frzbcdep, frzdudep,                        &
   hetraer, awcam, awfacm, dstcoat,                   &
   total_aer_num, coated_aer_num, uncoated_aer_num,  &
   total_interstitial_aer_num, total_cloudborne_aer_num, errstring)

   real(r8), intent(in) :: deltat            ! timestep [s]
   real(r8), intent(in) :: t                 ! temperature [K]
   real(r8), intent(in) :: p                 ! pressure [Pa]
   real(r8), intent(in) :: supersatice       ! supersaturation ratio wrt ice at 100%rh over water [ ]
   real(r8), intent(in) :: r3lx              ! volume mean drop radius [m]
   real(r8), intent(in) :: icnlx             ! in-cloud droplet concentration [cm-3]

   real(r8), intent(in) :: fn(3)               ! fraction activated [ ] for cloud borne aerosol number
                                               ! index values are 1:bc, 2:dust_a1, 3:dust_a3
   real(r8), intent(in) :: hetraer(3)          ! bc and dust mass mean radius [m]
   real(r8), intent(in) :: awcam(3)            ! modal added mass [mug m-3]
   real(r8), intent(in) :: awfacm(3)           ! (OC+BC)/(OC+BC+SO4)
   real(r8), intent(in) :: dstcoat(3)          ! coated fraction
   real(r8), intent(in) :: total_aer_num(3)    ! total bc and dust number concentration(interstitial+cloudborne) [#/cm^3]
   real(r8), intent(in) :: coated_aer_num(3)   ! coated bc and dust number concentration(interstitial)
   real(r8), intent(in) :: uncoated_aer_num(3) ! uncoated bc and dust number concentration(interstitial)
   real(r8), intent(in) :: total_interstitial_aer_num(3) ! total bc and dust concentration(interstitial)
   real(r8), intent(in) :: total_cloudborne_aer_num(3)   ! total bc and dust concentration(cloudborne)

   real(r8), intent(out) :: frzbcimm           ! het. frz by BC immersion nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzduimm           ! het. frz by dust immersion nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzbccnt           ! het. frz by BC contact nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzducnt           ! het. frz by dust contact nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzbcdep           ! het. frz by BC deposition nucleation [cm-3 s-1]
   real(r8), intent(out) :: frzdudep           ! het. frz by dust deposition nucleation [cm-3 s-1]

   character(len=*), intent(out) :: errstring

   ! local variables

   real(r8) :: aw(3)                           ! water activity [ ]
   real(r8) :: molal(3)                        ! molality [moles/kg]
   real(r8), parameter :: Mso4 = 96.06_r8

   integer, parameter :: id_bc   = 1
   integer, parameter :: id_dst1 = 2
   integer, parameter :: id_dst3 = 3
   logical :: do_bc, do_dst1, do_dst3

   real(r8), parameter :: n1 = 1.e19_r8           ! number of water molecules in contact with unit area of substrate [m-2]
   real(r8), parameter :: kboltz = 1.38e-23_r8    
   real(r8), parameter :: hplanck = 6.63e-34_r8
   real(r8), parameter :: rhplanck = 1._r8/hplanck
   real(r8), parameter :: amu = 1.66053886e-27_r8 
   real(r8), parameter :: nus = 1.e13_r8         ! frequ. of vibration [s-1] higher freq. (as in P&K, consistent with Anupam's data) 
   real(r8), parameter :: taufrz = 195.435_r8     ! time constant for falloff of freezing rate [s]
   real(r8), parameter :: rhwincloud = 0.98_r8    ! 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
   real(r8), parameter :: limfacbc = 0.01_r8      ! max. ice nucleating fraction soot
   real(r8), parameter :: pi = 4._r8*atan(1.0_r8)   
   real(r8) :: tc   
   real(r8) :: vwice   
   real(r8) :: rhoice   
   real(r8) :: sigma_iw                        ! [J/m2]   
   real(r8) :: sigma_iv                        ! [J/m2]   
   real(r8) :: esice                           ! [Pa]   
   real(r8) :: eswtr                           ! [Pa]   
   real(r8) :: rgimm   
   real(r8) :: rgdep   
   real(r8) :: dg0dep   
   real(r8) :: Adep   
   real(r8) :: dg0cnt   
   real(r8) :: Acnt   
   real(r8) :: rgimm_bc   
   real(r8) :: rgimm_dust_a1, rgimm_dust_a3   
   real(r8) :: dg0imm_bc  
   real(r8) :: dg0imm_dust_a1, dg0imm_dust_a3  
   real(r8) :: Aimm_bc
   real(r8) :: Aimm_dust_a1, Aimm_dust_a3
   real(r8) :: q, m, phi   
   real(r8) :: r_bc                            ! model radii of BC modes [m]   
   real(r8) :: r_dust_a1, r_dust_a3            ! model radii of dust modes [m]   
   real(r8) :: f_imm_bc 
   real(r8) :: f_imm_dust_a1, f_imm_dust_a3 
   real(r8) :: Jimm_bc
   real(r8) :: Jimm_dust_a1, Jimm_dust_a3
   real(r8) :: f_dep_bc 
   real(r8) :: f_dep_dust_a1, f_dep_dust_a3  
   real(r8) :: Jdep_bc 
   real(r8) :: Jdep_dust_a1, Jdep_dust_a3 
   real(r8) :: f_cnt_bc   
   real(r8) :: f_cnt_dust_a1,f_cnt_dust_a3   
   real(r8) :: Jcnt_bc
   real(r8) :: Jcnt_dust_a1,Jcnt_dust_a3 
   integer :: i

   !********************************************************
   ! Hoose et al., 2010 fitting parameters
   !********************************************************
   !freezing parameters for immersion freezing
   !real(r8),parameter :: theta_imm_bc = 40.17         ! contact angle [deg], converted to rad later
   !real(r8),parameter :: dga_imm_bc = 14.4E-20        ! activation energy [J]
   !real(r8),parameter :: theta_imm_dust = 30.98       ! contact angle [deg], converted to rad later
   !real(r8),parameter :: dga_imm_dust = 15.7E-20      ! activation energy [J]
   !freezing parameters for deposition nucleation
   !real(r8),parameter :: theta_dep_dust = 12.7       ! contact angle [deg], converted to rad later !Zimmermann et al (2008), illite
   !real(r8),parameter :: dga_dep_dust = -6.21E-21    ! activation energy [J]
   !real(r8),parameter :: theta_dep_bc = 28.          ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
   !real(r8),parameter :: dga_dep_bc = -2.E-19        ! activation energy [J]
   !********************************************************
   ! Wang et al., 2014 fitting parameters
   !********************************************************
   ! freezing parameters for immersion freezing
   real(r8),parameter :: theta_imm_bc = 48.0_r8            ! contact angle [deg], converted to rad later !DeMott et al (1990)
   real(r8),parameter :: dga_imm_bc = 14.15E-20_r8         ! activation energy [J]
   real(r8),parameter :: theta_imm_dust = 46.0_r8          ! contact angle [deg], converted to rad later !DeMott et al (2011) SD
   real(r8),parameter :: dga_imm_dust = 14.75E-20_r8       ! activation energy [J]
   ! freezing parameters for deposition nucleation
   real(r8),parameter :: theta_dep_dust = 20.0_r8          ! contact angle [deg], converted to rad later !Koehler et al (2010) SD
   real(r8),parameter :: dga_dep_dust = -8.1E-21_r8        ! activation energy [J]
   real(r8),parameter :: theta_dep_bc = 28._r8             ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
   real(r8),parameter :: dga_dep_bc = -2.E-19_r8           ! activation energy [J]

   real(r8) :: Kcoll_bc                                    ! collision kernel [cm3 s-1]
   real(r8) :: Kcoll_dust_a1                               ! collision kernel [cm3 s-1]
   real(r8) :: Kcoll_dust_a3                               ! collision kernel [cm3 s-1]

   logical :: tot_in = .false.

   !*****************************************************************************
   !                PDF theta model 
   !*****************************************************************************
   ! some variables for PDF theta model
   ! immersion freezing
   !
   ! With the original value of pdf_n_theta set to 101 the dust activation
   ! fraction between -15 and 0 C could be overestimated.  This problem was
   ! eliminated by increasing pdf_n_theta to 301.  To reduce the expense of
   ! computing the dust activation fraction the integral is only evaluated
   ! where dim_theta is non-zero.  This was determined to be between
   ! dim_theta index values of 53 through 113.  These loop bounds are
   ! hardcoded in the variables i1 and i2.
   !
   real(r8),parameter :: theta_min = 1._r8/180._r8*pi
   real(r8),parameter :: theta_max = 179._r8/180._r8*pi
   real(r8) :: x1_imm
   real(r8) :: x2_imm
   real(r8) :: norm_theta_imm
   real(r8),parameter :: imm_dust_mean_theta = 46.0_r8/180.0_r8*pi 
   real(r8),parameter :: imm_dust_var_theta = 0.01_r8
   real(r8) :: pdf_d_theta
   integer,parameter :: pdf_n_theta = 301
   integer,parameter :: i1 = 53
   integer,parameter :: i2 = 113
   real(r8) :: dim_theta(pdf_n_theta)
   real(r8) :: dim_f_imm_dust_a1(pdf_n_theta), dim_f_imm_dust_a3(pdf_n_theta)
   real(r8) :: dim_Jimm_dust_a1(pdf_n_theta), dim_Jimm_dust_a3(pdf_n_theta)
   real(r8) :: pdf_imm_theta(pdf_n_theta)
   real(r8) :: sum_imm_dust_a1, sum_imm_dust_a3
   logical :: pdf_imm_in = .true.
   !------------------------------------------------------------------------------------------------

   errstring = ' '

   if (pdf_imm_in) then
      pdf_d_theta = (179._r8-1._r8)/180._r8*pi/(pdf_n_theta-1)
      ! calculate the integral in the denominator
      x1_imm = (LOG(theta_min)-LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
      x2_imm = (LOG(theta_max)-LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
      norm_theta_imm = (ERF(x2_imm)-ERF(x1_imm))*0.5_r8
      dim_theta = 0.0_r8
      pdf_imm_theta = 0.0_r8
      do i = i1,i2
         dim_theta(i) = 1._r8/180._r8*pi+(i-1)*pdf_d_theta
         pdf_imm_theta(i) = exp(-((LOG(dim_theta(i))-LOG(imm_dust_mean_theta))**2._r8)/(2._r8*imm_dust_var_theta**2._r8))/ &
                                (dim_theta(i)*imm_dust_var_theta*SQRT(2*pi))/norm_theta_imm
      end do
   end if

   ! get saturation vapor pressures
   eswtr = svp_water(t)  ! 0 for liquid
   esice = svp_ice(t)  ! 1 for ice

   tc = t - tmelt
   rhoice = 916.7_r8-0.175_r8*tc-5.e-4_r8*tc**2
   vwice = mwh2o*amu/rhoice
   sigma_iw = (28.5_r8+0.25_r8*tc)*1E-3_r8
   sigma_iv = (76.1_r8-0.155_r8*tc + 28.5_r8+0.25_r8*tc)*1E-3_r8

   ! get mass mean radius
   r_bc = hetraer(1)    
   r_dust_a1 = hetraer(2)    
   r_dust_a3 = hetraer(3)    

   ! calculate collision kernels as a function of environmental parameters and aerosol/droplet sizes
   call collkernel(t, p, eswtr, rhwincloud, r3lx,         &
                   r_bc,                                  &  ! BC modes
                   r_dust_a1, r_dust_a3,                  &  ! dust modes
                   Kcoll_bc,                              &  ! collision kernel [cm3 s-1]
                   Kcoll_dust_a1, Kcoll_dust_a3)
        
   !*****************************************************************************
   !                take water activity into account 
   !*****************************************************************************
   !   solute effect
   aw(:) = 1._r8
   molal(:) = 0._r8

   ! The heterogeneous ice freezing temperatures of all IN generally decrease with
   ! increasing total solute mole fraction. Therefore, the large solution concentration
   ! will cause the freezing point depression and the ice freezing temperatures of all
   ! IN will get close to the homogeneous ice freezing temperatures. Since we take into
   ! account water activity for three heterogeneous freezing modes(immersion, deposition, 
   ! and contact), we utilize interstitial aerosols(not cloudborne aerosols) to calculate 
   ! water activity. 
   ! If the index of IN is 0, it means three freezing modes of this aerosol are depressed.

   do i = 1, 3
      !calculate molality
      if ( total_interstitial_aer_num(i) > 0._r8 ) then
         molal(i) = (1.e-6_r8*awcam(i)*(1._r8-awfacm(i))/(Mso4*total_interstitial_aer_num(i)*1.e6_r8))/ &
            (4*pi/3*rhoh2o*(MAX(r3lx,4.e-6_r8))**3)
         aw(i) = 1._r8/(1._r8+2.9244948e-2_r8*molal(i)+2.3141243e-3_r8*molal(i)**2+7.8184854e-7_r8*molal(i)**3)
      end if
   end do

   !*****************************************************************************
   !                immersion freezing begin 
   !*****************************************************************************    

   frzbcimm = 0._r8
   frzduimm = 0._r8
   frzbccnt = 0._r8
   frzducnt = 0._r8
   frzbcdep = 0._r8
   frzdudep = 0._r8

   ! critical germ size
   rgimm = 2*vwice*sigma_iw/(kboltz*t*LOG(supersatice))
   ! take solute effect into account
   rgimm_bc = rgimm
   rgimm_dust_a1 = rgimm
   rgimm_dust_a3 = rgimm

   ! if aw*Si<=1, the freezing point depression is strong enough to prevent freezing

   if (aw(id_bc)*supersatice > 1._r8 ) then
      do_bc   = .true.
      rgimm_bc = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_bc)*supersatice))
   else
      do_bc = .false.
   end if

   if (aw(id_dst1)*supersatice > 1._r8 ) then
      do_dst1 = .true.
      rgimm_dust_a1 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst1)*supersatice))
   else
      do_dst1 = .false.
   end if

   if (aw(id_dst3)*supersatice > 1._r8 ) then
      do_dst3 = .true.
      rgimm_dust_a3 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst3)*supersatice))
   else
      do_dst3 = .false.
   end if
    
   ! form factor
   ! only consider flat surfaces due to uncertainty of curved surfaces

   m = COS(theta_imm_bc*pi/180._r8)
   f_imm_bc = (2+m)*(1-m)**2/4._r8
   if (.not. pdf_imm_in) then
      m = COS(theta_imm_dust*pi/180._r8)
      f_imm_dust_a1 = (2+m)*(1-m)**2/4._r8

      m = COS(theta_imm_dust*pi/180._r8)
      f_imm_dust_a3 = (2+m)*(1-m)**2/4._r8
   else 
      dim_f_imm_dust_a1 = 0.0_r8
      dim_f_imm_dust_a3 = 0.0_r8
      do i = i1,i2
         m = cos(dim_theta(i))
         dim_f_imm_dust_a1(i) = (2+m)*(1-m)**2/4._r8

         m = cos(dim_theta(i))
         dim_f_imm_dust_a3(i) = (2+m)*(1-m)**2/4._r8
      end do
   end if

   ! homogeneous energy of germ formation
   dg0imm_bc = 4*pi/3._r8*sigma_iw*rgimm_bc**2
   dg0imm_dust_a1 = 4*pi/3._r8*sigma_iw*rgimm_dust_a1**2
   dg0imm_dust_a3 = 4*pi/3._r8*sigma_iw*rgimm_dust_a3**2

   ! prefactor
   Aimm_bc = n1*((vwice*rhplanck)/(rgimm_bc**3)*SQRT(3._r8/pi*kboltz*T*dg0imm_bc))
   Aimm_dust_a1 = n1*((vwice*rhplanck)/(rgimm_dust_a1**3)*SQRT(3._r8/pi*kboltz*T*dg0imm_dust_a1))
   Aimm_dust_a3 = n1*((vwice*rhplanck)/(rgimm_dust_a3**3)*SQRT(3._r8/pi*kboltz*T*dg0imm_dust_a3))

   ! nucleation rate per particle

   Jimm_bc = Aimm_bc*r_bc**2/SQRT(f_imm_bc)*EXP((-dga_imm_bc-f_imm_bc*dg0imm_bc)/(kboltz*T))
   if (.not. pdf_imm_in) then
      ! 1/sqrt(f)
      ! the expression of Chen et al. (sqrt(f)) may however lead to unphysical
      ! behavior as it implies J->0 when f->0 (i.e. ice nucleation would be
      ! more difficult on easily wettable materials). 
      Jimm_dust_a1 = Aimm_dust_a1*r_dust_a1**2/SQRT(f_imm_dust_a1)*EXP((-dga_imm_dust-f_imm_dust_a1*dg0imm_dust_a1)/(kboltz*T))
      Jimm_dust_a3 = Aimm_dust_a3*r_dust_a3**2/SQRT(f_imm_dust_a3)*EXP((-dga_imm_dust-f_imm_dust_a3*dg0imm_dust_a3)/(kboltz*T))
   end if

   if (pdf_imm_in) then
      dim_Jimm_dust_a1 = 0.0_r8
      dim_Jimm_dust_a3 = 0.0_r8
      do i = i1,i2
         ! 1/sqrt(f)
         dim_Jimm_dust_a1(i) = Aimm_dust_a1*r_dust_a1**2/SQRT(dim_f_imm_dust_a1(i))*EXP((-dga_imm_dust-dim_f_imm_dust_a1(i)* &
            dg0imm_dust_a1)/(kboltz*T))
         dim_Jimm_dust_a1(i) = max(dim_Jimm_dust_a1(i), 0._r8)

         dim_Jimm_dust_a3(i) = Aimm_dust_a3*r_dust_a3**2/SQRT(dim_f_imm_dust_a3(i))*EXP((-dga_imm_dust-dim_f_imm_dust_a3(i)* &
            dg0imm_dust_a3)/(kboltz*T))
         dim_Jimm_dust_a3(i) = max(dim_Jimm_dust_a3(i), 0._r8)
      end do
   end if

   ! Limit to 1% of available potential IN (for BC), no limit for dust 
   if (pdf_imm_in) then
      sum_imm_dust_a1 = 0._r8
      sum_imm_dust_a3 = 0._r8
      do i = i1,i2-1
         sum_imm_dust_a1 = sum_imm_dust_a1+0.5_r8*((pdf_imm_theta(i)*exp(-dim_Jimm_dust_a1(i)*deltat)+ &
            pdf_imm_theta(i+1)*exp(-dim_Jimm_dust_a1(i+1)*deltat)))*pdf_d_theta
         sum_imm_dust_a3 = sum_imm_dust_a3+0.5_r8*((pdf_imm_theta(i)*exp(-dim_Jimm_dust_a3(i)*deltat)+ &
            pdf_imm_theta(i+1)*exp(-dim_Jimm_dust_a3(i+1)*deltat)))*pdf_d_theta
      end do
      do i = i1,i2
      	   if (sum_imm_dust_a1 > 0.99_r8) then
               sum_imm_dust_a1 = 1.0_r8
           end if
           if (sum_imm_dust_a3 > 0.99_r8) then
               sum_imm_dust_a3 = 1.0_r8
           end if
      end do
      
   end if

   if (.not.tot_in) then
      if (do_bc) frzbcimm = frzbcimm+MIN(limfacbc*total_cloudborne_aer_num(id_bc)/deltat, &
                            total_cloudborne_aer_num(id_bc)/deltat*(1._r8-exp(-Jimm_bc*deltat))) 

      if (.not. pdf_imm_in) then
         if (do_dst1) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat, & 
                                 total_cloudborne_aer_num(id_dst1)/deltat*(1._r8-exp(-Jimm_dust_a1*deltat)))
         if (do_dst3) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat, &
                                 total_cloudborne_aer_num(id_dst3)/deltat*(1._r8-exp(-Jimm_dust_a3*deltat)))
      else
         if (do_dst1) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat,        &
                                 total_cloudborne_aer_num(id_dst1)/deltat*(1._r8-sum_imm_dust_a1))
         if (do_dst3) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat,        &
                                 total_cloudborne_aer_num(id_dst3)/deltat*(1._r8-sum_imm_dust_a3))
      end if

   else
      if (do_bc) frzbcimm = frzbcimm+MIN(limfacbc*fn(id_bc)*total_aer_num(id_bc)/deltat, & 
                            fn(id_bc)*total_aer_num(id_bc)/deltat*(1._r8-exp(-Jimm_bc*deltat))) 

      if (.not. pdf_imm_in) then
         if (do_dst1) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat, &
                                 fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1._r8-exp(-Jimm_dust_a1*deltat)))
         if (do_dst3) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat, &
                                 fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1._r8-exp(-Jimm_dust_a3*deltat)))
      else
         if (do_dst1) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat,        &
                                 fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1._r8-sum_imm_dust_a1))
         if (do_dst3) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat,        &
                                 fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1._r8-sum_imm_dust_a3))
      end if
   end if

   if (t > 263.15_r8) then
      frzduimm = 0._r8
      frzbcimm = 0._r8
   end if

   !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   Deposition nucleation
   !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! critical germ size
   ! assume 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
   rgdep=2*vwice*sigma_iv/(kboltz*t*LOG(rhwincloud*supersatice)) 

   ! form factor
   m = COS(theta_dep_bc*pi/180._r8)
   f_dep_bc = (2+m)*(1-m)**2/4._r8

   m = COS(theta_dep_dust*pi/180._r8)
   f_dep_dust_a1 = (2+m)*(1-m)**2/4._r8

   m = COS(theta_dep_dust*pi/180._r8)
   f_dep_dust_a3 = (2+m)*(1-m)**2/4._r8

   ! homogeneous energy of germ formation
   dg0dep = 4*pi/3._r8*sigma_iv*rgdep**2

   ! prefactor
   ! attention: division of small numbers
   Adep = (rhwincloud*eswtr)**2*(vwice/(mwh2o*amu))/(kboltz*T*nus)*SQRT(sigma_iv/(kboltz*T))

   ! nucleation rate per particle
   if (rgdep > 0) then
      Jdep_bc = Adep*r_bc**2/SQRT(f_dep_bc)*EXP((-dga_dep_bc-f_dep_bc*dg0dep)/(kboltz*T))
      Jdep_dust_a1 = Adep*r_dust_a1**2/SQRT(f_dep_dust_a1)*EXP((-dga_dep_dust-f_dep_dust_a1*dg0dep)/(kboltz*T))
      Jdep_dust_a3 = Adep*r_dust_a3**2/SQRT(f_dep_dust_a3)*EXP((-dga_dep_dust-f_dep_dust_a3*dg0dep)/(kboltz*T))
   else
      Jdep_bc = 0._r8
      Jdep_dust_a1 = 0._r8
      Jdep_dust_a3 = 0._r8
   end if

   ! Limit to 1% of available potential IN (for BC), no limit for dust 
   if (.not.tot_in) then
      if (do_bc) frzbcdep = frzbcdep+MIN(limfacbc*uncoated_aer_num(id_bc)/deltat, &
                                         uncoated_aer_num(id_bc)/deltat &
                                         *(1._r8-exp(-Jdep_bc*deltat))) 
      if (do_dst1) frzdudep = frzdudep+MIN(uncoated_aer_num(id_dst1)/deltat, &
                                           uncoated_aer_num(id_dst1)/deltat &
                                           *(1._r8-exp(-Jdep_dust_a1*deltat)))
      if (do_dst3) frzdudep = frzdudep+MIN(uncoated_aer_num(id_dst3)/deltat, &
                                           uncoated_aer_num(id_dst3)/deltat &
                                           *(1._r8-exp(-Jdep_dust_a3*deltat)))
   else
      if (do_bc) frzbcdep = frzbcdep+MIN(limfacbc*(1._r8-fn(id_bc)) &
                                         *(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat, &
                                          (1._r8-fn(id_bc))*(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat &
                                          *(1._r8-exp(-Jdep_bc*deltat))) 
      if (do_dst1) frzdudep = frzdudep+MIN((1._r8-fn(id_dst1)) &
                                           *(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat, &
                                           (1._r8-fn(id_dst1))*(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat &
                                           *(1._r8-exp(-Jdep_dust_a1*deltat))) 
      if (do_dst3) frzdudep = frzdudep+MIN((1._r8-fn(id_dst3)) &
                                           *(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat, &
                                           (1._r8-fn(id_dst3))*(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat &
                                           *(1._r8-exp(-Jdep_dust_a3*deltat))) 
   end if

   !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! contact nucleation
   !!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! form factor
   m = COS(theta_dep_bc*pi/180._r8)
   f_cnt_bc = (2+m)*(1-m)**2/4._r8

   m = COS(theta_dep_dust*pi/180._r8)
   f_cnt_dust_a1 = (2+m)*(1-m)**2/4._r8

   m = COS(theta_dep_dust*pi/180._r8)
   f_cnt_dust_a3 = (2+m)*(1-m)**2/4._r8

   ! homogeneous energy of germ formation
   dg0cnt = 4*pi/3._r8*sigma_iv*rgimm**2

   ! prefactor
   ! attention: division of small numbers
   Acnt = rhwincloud*eswtr*4*pi/(nus*SQRT(2*pi*mwh2o*amu*kboltz*T))

   ! nucleation rate per particle
   Jcnt_bc = Acnt*r_bc**2*EXP((-dga_dep_bc-f_cnt_bc*dg0cnt)/(kboltz*T))*Kcoll_bc*icnlx
   Jcnt_dust_a1 = Acnt*r_dust_a1**2*EXP((-dga_dep_dust-f_cnt_dust_a1*dg0cnt)/(kboltz*T))*Kcoll_dust_a1*icnlx
   Jcnt_dust_a3 = Acnt*r_dust_a3**2*EXP((-dga_dep_dust-f_cnt_dust_a3*dg0cnt)/(kboltz*T))*Kcoll_dust_a3*icnlx

   ! Limit to 1% of available potential IN (for BC), no limit for dust 
   if (.not.tot_in) then
      if (do_bc) frzbccnt = frzbccnt+MIN(limfacbc*uncoated_aer_num(id_bc)/deltat, &
                                         uncoated_aer_num(id_bc)/deltat &
                                         *(1._r8-exp(-Jcnt_bc*deltat))) 
      if (do_dst1) frzducnt = frzducnt+MIN(uncoated_aer_num(id_dst1)/deltat, &
                                           uncoated_aer_num(id_dst1)/deltat &
                                           *(1._r8-exp(-Jcnt_dust_a1*deltat)))
      if (do_dst3) frzducnt = frzducnt+MIN(uncoated_aer_num(id_dst3)/deltat, &
                                           uncoated_aer_num(id_dst3)/deltat &
                                           *(1._r8-exp(-Jcnt_dust_a3*deltat)))
   else
      if (do_bc) frzbccnt = frzbccnt+MIN(limfacbc*(1._r8-fn(id_bc))*(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat, &
                                         (1._r8-fn(id_bc))*(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat &
                                         *(1._r8-exp(-Jcnt_bc*deltat))) 
      if (do_dst1) frzducnt = frzducnt+MIN((1._r8-fn(id_dst1))*(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat, &
                                           (1._r8-fn(id_dst1))*(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat &
                                           *(1._r8-exp(-Jcnt_dust_a1*deltat)))
      if (do_dst3) frzducnt = frzducnt+MIN((1._r8-fn(id_dst3))*(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat, &
                                           (1._r8-fn(id_dst3))*(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat &
                                           *(1._r8-exp(-Jcnt_dust_a3*deltat)))
   end if
    
   if (frzducnt <= -1._r8) then
      write(iulog,*) 'hetfrz_classnuc_calc: frzducnt', frzducnt, Jcnt_dust_a1,Jcnt_dust_a3, &    
                                    Kcoll_dust_a1, Kcoll_dust_a3
      errstring = 'ERROR in hetfrz_classnuc_calc::frzducnt'
      return
   end if

end subroutine  hetfrz_classnuc_calc

!===================================================================================================

!-----------------------------------------------------------------------
!
! Purpose: calculate collision kernels as a function of environmental parameters and aerosol/droplet sizes
!
! Author: Corinna Hoose, UiO, October 2009
!
! Modifications: Yong Wang and Xiaohong Liu, UWyo, 12/2012
!-----------------------------------------------------------------------

subroutine collkernel( &
   t, pres, eswtr, rhwincloud, r3lx,       &
   r_bc,                                   &  ! BC modes
   r_dust_a1, r_dust_a3,                   &  ! dust modes
   Kcoll_bc,                               &  ! collision kernel [cm3 s-1]
   Kcoll_dust_a1, Kcoll_dust_a3)

   real(r8), intent(in) :: t                ! temperature [K]
   real(r8), intent(in) :: pres             ! pressure [Pa]
   real(r8), intent(in) :: eswtr            ! saturation vapor pressure of water [Pa]
   real(r8), intent(in) :: r3lx             ! volume mean drop radius [m]
   real(r8), intent(in) :: rhwincloud       ! in-cloud relative humidity over water [ ]
   real(r8), intent(in) :: r_bc             ! model radii of BC modes [m]
   real(r8), intent(in) :: r_dust_a1        ! model radii of dust modes [m]
   real(r8), intent(in) :: r_dust_a3        ! model radii of dust modes [m]

   real(r8), intent(out) :: Kcoll_bc        ! collision kernel [cm3 s-1]
   real(r8), intent(out) :: Kcoll_dust_a1
   real(r8), intent(out) :: Kcoll_dust_a3

   ! local variables
   real(r8) :: a, b, c, a_f, b_f, c_f, f
   real(r8) :: tc          ! temperature [deg C]
   real(r8) :: rho_air     ! air density [kg m-3]
   real(r8) :: viscos_air  ! dynamic viscosity of air [kg m-1 s-1]
   real(r8) :: Ktherm_air  ! thermal conductivity of air [J/(m s K)]
   real(r8) :: lambda      ! mean free path [m]
   real(r8) :: Kn          ! Knudsen number [ ]
   real(r8) :: Re          ! Reynolds number [ ]
   real(r8) :: Pr          ! Prandtl number [ ]
   real(r8) :: Sc          ! Schmidt number [ ]
   real(r8) :: vterm       ! terminal velocity [m s-1]
   real(r8) :: Ktherm      ! thermal conductivity of aerosol [J/(m s K)]
   real(r8) :: Dvap        ! water vapor diffusivity [m2 s-1]
   real(r8) :: Daer        ! aerosol diffusivity [m2 s-1]
   real(r8) :: latvap      ! latent heat of vaporization [J kg-1]
   real(r8) :: kboltz      ! Boltzmann constant [J K-1]
   real(r8) :: G           ! thermodynamic function in Cotton et al. [kg m-1 s-1]
   real(r8) :: r_a         ! aerosol radius [m]
   real(r8) :: f_t         ! factor by Waldmann & Schmidt [ ]
   real(r8) :: Q_heat      ! heat flux [J m-2 s-1]
   real(r8) :: Tdiff_cotton ! temperature difference between droplet and environment [K]
   real(r8) :: K_brownian,K_thermo_cotton,K_diffusio_cotton   ! collision kernels [m3 s-1]
   real(r8) :: K_total     ! total collision kernel [cm3 s-1]
   integer  :: i
   !------------------------------------------------------------------------------------------------
        
   Kcoll_bc      = 0._r8
   Kcoll_dust_a1 = 0._r8
   Kcoll_dust_a3 = 0._r8

   tc     = t - tmelt
   kboltz = 1.38065e-23_r8

   ! air viscosity for tc<0, from depvel_part.F90
   viscos_air = (1.718_r8+0.0049_r8*tc-1.2e-5_r8*tc*tc)*1.e-5_r8
   ! air density
   rho_air = pres/(rair*t)
   ! mean free path: Seinfeld & Pandis 8.6
   lambda = 2*viscos_air/(pres*SQRT(8/(pi*rair*t)))
   ! latent heat of vaporization, varies with T
   latvap = 1000*(-0.0000614342_r8*tc**3 + 0.00158927_r8*tc**2 - 2.36418_r8*tc + 2500.79_r8)
   ! droplet terminal velocity after Chen & Liu, QJRMS 2004
   a = 8.8462e2_r8
   b = 9.7593e7_r8
   c = -3.4249e-11_r8
   a_f = 3.1250e-1_r8
   b_f = 1.0552e-3_r8
   c_f = -2.4023_r8
   f = EXP(EXP(a_f + b_f*(LOG(r3lx))**3 + c_f*rho_air**1.5_r8))
   vterm = (a+ (b + c*r3lx)*r3lx)*r3lx*f

   ! Reynolds number
   Re = 2*vterm*r3lx*rho_air/viscos_air
   ! thermal conductivity of air: Seinfeld & Pandis eq. 15.75
   Ktherm_air = 1.e-3_r8*(4.39_r8+0.071_r8*t)  !J/(m s K)
   ! Prandtl number
   Pr = viscos_air*cpair/Ktherm_air
   ! water vapor diffusivity: Pruppacher & Klett 13-3
   Dvap = 0.211e-4_r8*(t/273.15_r8)*(101325._r8/pres) 
   ! G-factor = rhoh2o*Xi in Rogers & Yau, p. 104
   G = rhoh2o/((latvap/(rh2o*t) - 1)*latvap*rhoh2o/(Ktherm_air*t) &
       + rhoh2o*rh2o*t/(Dvap*eswtr))
      
   ! variables depending on aerosol radius
   ! loop over 3 aerosol modes
   do i = 1, 3
      if (i == 1) r_a = r_bc
      if (i == 2) r_a = r_dust_a1
      if (i == 3) r_a = r_dust_a3
      ! Knudsen number (Seinfeld & Pandis 8.1)
      Kn = lambda/r_a
      ! aerosol diffusivity
      Daer = kboltz*t*(1 + Kn)/(6*pi*r_a*viscos_air)
      ! Schmidt number
      Sc = viscos_air/(Daer*rho_air)

      ! Young (1974) first equ. on page 771
      K_brownian = 4*pi*r3lx*Daer*(1 + 0.3_r8*Re**0.5_r8*Sc**0.33_r8)

      ! thermal conductivities from Seinfeld & Pandis, Table 8.6
      if (i == 1) Ktherm = 4.2_r8 ! Carbon
      if (i == 2 .or. i == 3) Ktherm = 0.72_r8 ! clay
      ! form factor
      f_t = 0.4_r8*(1._r8 + 1.45_r8*Kn + 0.4_r8*Kn*EXP(-1._r8/Kn))      &
                  *(Ktherm_air + 2.5_r8*Kn*Ktherm)                      &
               /((1._r8 + 3._r8*Kn)*(2._r8*Ktherm_air + 5._r8*Kn*Ktherm+Ktherm))
      ! calculate T-Tc as in Cotton et al.
      Tdiff_cotton = -G*(rhwincloud - 1._r8)*latvap/Ktherm_air
      Q_heat = Ktherm_air/r3lx*(1._r8 + 0.3_r8*Re**0.5_r8*Pr**0.33_r8)*Tdiff_cotton
      K_thermo_cotton = 4._r8*pi*r3lx*r3lx*f_t*Q_heat/pres
      K_diffusio_cotton = -(1._r8/f_t)*(rh2o*t/latvap)*K_thermo_cotton
      K_total = 1.e6_r8*(K_brownian + K_thermo_cotton + K_diffusio_cotton)  ! convert m3/s -> cm3/s
      ! set K to 0 if negative
      if (K_total .lt. 0._r8) K_total = 0._r8

      if (i == 1) Kcoll_bc = K_total
      if (i == 2) Kcoll_dust_a1 = K_total
      if (i == 3) Kcoll_dust_a3 = K_total
        
   end do
      
end subroutine collkernel

!===================================================================================================


end module hetfrz_classnuc
