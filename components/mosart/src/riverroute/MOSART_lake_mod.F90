!
MODULE MOSART_lake_mod
! Description: core code of lake module in MOSART framework
! A simple wier formula (Kindsvater and Carter (1959))is used to estimate outflow from lakes 
! Developed by Wondie Yigzaw, July 2020. 
! REVISION HISTORY:
! 
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TPara, rtmCTL
    use rof_cpl_indices, only : nt_nliq, nt_nice
	use RtmVar         , only : iulog, ngeom, nlayers, rstraflag, lakeflag
	use MOSART_heat_mod
	use RtmTimeManager
    implicit none
	
    public mosart_lake
	
! !PUBLIC MEMBER FUNCTIONS:
    contains  
		
	subroutine mosart_lake(iunit,nt,localDeltaT,temp_Tr)
	! !DESCRIPTION: calculate the water temperature of reservoir.
	
		use shr_sys_mod , only : shr_sys_flush
	
		implicit none
		integer,  intent(in) :: iunit, nt
		real(r8), intent(in) :: localDeltaT,temp_Tr
		character(len=*),parameter :: subname = '(mosart_lake)'
        integer :: j,i,w,k,ii,l,ww,m,n,nn,mm,jmax,jmin							! indices

		real(r8) :: c_w = 4.188e3												! Water specific heat capacity  (J/kg/k)
		real(r8) :: F = 0.8_r8      											! dimensionless factor for wind sheltering by riparian vegetation, Wu et al, 2012
		real(r8) :: grav   = 9.8062_r8	      									! gravity constant (m/s2)
		real(r8) :: le= 2.501e6													! latent heat of vaporaization (J/kg)
		real(r8) :: rho_w = 1.e3												! Water density  (kg/m3)
		real(r8) :: rho_a = 1.177_r8											! Air density  (kg/m3)
		real(r8) :: st_bl = 5.67e-8												! Stefan-Boltzmann constant ~ W/m^2/K^4
		real(r8) :: t_frz   = 273.15_r8											! freezing temperature (K)
		real(r8) :: sar   	= 1.0_r8											! Surface area ratio
		integer  :: s_dtime														! number of sub hourly time step
		integer  :: d_n_n														! Adjusted layer number
		integer  :: outlet														! Outlet layer location
		integer  :: mix1,mix2													! used in convective mixing
		real(r8) :: dtime														! time step (sec)
		real(r8) :: es,ea														! Satuated, atmospheric vapor pressure(mb)
		real(r8) :: evap														! evaporation rate (mm/d)
		real(r8) :: kl															! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
		real(r8) :: d_res_sub													! Reservoir depth, taken as 0.95*(dam height) in GRanD database, (m)
		real(r8) :: df_eff(nlayers+1)  			          						! Effective diffusivity (molecular + eddy) [m^2/s]
		real(r8) :: rho_z(nlayers),rho_r     		   							! Depth based water density, density of inflow water   (kg/m^3)
		real(r8) :: drhodz(nlayers)												! d [rhow] /dz (kg/m**4)
		real(r8) :: m1(nlayers),m2(nlayers),m3(nlayers)							! used in tridiagonal matrix
		real(r8) :: Fx(nlayers)													! used in tridiagonal matrix
		real(r8) :: a(nlayers),b(nlayers),c(nlayers)    						! columns of tridiagonal matrix
		real(r8) :: r(nlayers)    												! right  diagonal tridiagonal matrix
		real(r8) :: ri                      									! Richardson number
		real(r8) :: q_in,t_in,q_ou			   									! inflow (m^3/s),inflow temperature (k), and outflow (m^3/s)
		real(r8) :: dv_in(nlayers),dv_ou(nlayers)		   						! volume increment/decrease at layer due to inflow/outflow(m^3)
		real(r8) :: d_evap,v_evap                        						! Evaporated depth (m) and volume (m^3)
		real(r8) :: delta_z                        								! depth change to calculate corresponding area/volume(m)
		real(r8) :: top_d                        								! top layer depth after mixing
		real(r8) :: v_zt0(nlayers+1)                        					! Total reservoir volume at depth z from surface(m^3)
		real(r8) :: v_mix                        								! Total volume of mixed layer(m^3)
		real(r8) :: t_s,t_out            										! reservoir surface and outflow temperature (k)
		real(r8) :: t_z(nlayers)            									! lake layer temperature
		real(r8) :: t_z_old(nlayers)            								! previous time lake layer temperature
		real(r8) :: bv_f                      									! brunt-vaisala frequency (/s**2)
		real(r8) :: eta      													! light extinction coefficient
		real(r8) :: sh_net        												! net short wave radiation
		real(r8) :: sh_mix        												! net short wave radiation in mixed layer
		real(r8) :: beta                          								! shortwave absorbtion factor
		real(r8) :: lw_abr  													! net atmospheric longwave absorbtion (W/m^2)
		real(r8) :: phi_o  														! net surface radiation (W/m^2)
		real(r8) :: phi_z(nlayers)  											! radiation absorbed by layer (W)
		real(r8) :: phi_x(nlayers+1)											! radiation absorbed by mixed layer (W/m^2)
		real(r8) :: lw_ems,lt_heat,sn_heat       								! longwave emission,latent heat,sensible heat (W/m^2)
		real(r8) :: alb_s                 										! surface albedo shortwave radiation
		real(r8) :: k_m 														! molecular diffusivity
		real(r8) :: enr_in(nlayers),enr_ou(nlayers) 							! Layer Energy from inflow,outflow (J/s,W)
		real(r8) :: enr_0(nlayers),enr_1(nlayers),enr_2(nlayers)				! Initial inner energy,Inner energy after advection,Inner energy after stratification(J/s,W)
		real(r8) :: enr_phi(nlayers) 											! Inner energy from solar/atmospheric radiation(J/s,W)
		real(r8) :: enr_err1,enr_err2 											! Energy error (w) before stratification and after triadiagonal solution
	    real(r8) :: th_en(nlayers)												! Layer thermal energy (j/s)
		real(r8) :: e_a,e_b,e_ab												! used in layer merging/split, factors
		real(r8) :: k_ew,k_ad(nlayers)											! Effective,wind,convection,advective kinetic energy (kg.m^2/s^2)
		real(r8) :: c_d															! Drag coefficient
		real(r8) :: cfa,cfw														! used in diffusion calculation
		real(r8) :: Fr(nlayers)  												! Froude number squared and inverted for diffusion coeff. calculation
		real(r8) :: dis_ad(nlayers)         									! rate of dissipation-inflow/outflow
		real(r8) :: dis_w                   									! rate of dissipation-wind
		real(r8) :: l_vel, s_vel												! used in Froude number calculation
		real(r8) :: tau															! shear stress
		real(r8) :: q_adv(nlayers)												! Layer flow rate (m^3/s)
		real(r8) :: denmix,tmix,tsum											! used in convective mixing
		real(r8) :: mixvol1,mixvol2,sumvol,num,dem 								! used in convective mixing
		real(r8) :: ta,tb,tab,dv_nt(nlayers),dd_za,dd_zb,d_va,dv_oua,delta_a	! used in layer merging/split
		real(r8) :: dv_inb,dv_ina,dv_oub,dv_ouab,dv_inab,dd_zab,d_vab,d_vb		! used in layer merging/split	
		real(r8) :: ddz_min,ddz_max      										! (nd) Minimum and Maximum layer thickness limit for merge/split

!**************************************************************************************************************************************************
		
		if (TRunoff%lake_flg(iunit) ==1) then	! Lake module active if there is natural lake
			! return
			do i=2,nlayers
               if (TRunoff%v_zt(iunit,i)<TRunoff%v_zt(iunit,i-1) .or. TRunoff%d_v(iunit,j)<0.0_r8)return
			end do
			do j = 2, ngeom+1	 
               if (TRunoff%a_di(iunit,j)<0._r8 .or. TRunoff%v_zti(iunit,j)<0._r8) return
			end do
			
		!	Calculate sub-time step for numerical stability			
			s_dtime = 12  			
			dtime   = localDeltaT/s_dtime		!	Sub-timestep for numerical purpose
			t_in    = temp_Tr 
			if(t_in<273.15_r8) t_in=273.15_r8
			q_in = TRunoff%lake_inflow(iunit)
			q_ou = -TRunoff%erout(iunit,nt)			
			
			!	Initialize for sub-timestep	
			! !***************************************************************************************************************		
			do ww = 1,s_dtime	!	Start calculation for each sub-timestep	************************************************			
				!	Initialize for sub-timestep
				d_n_n  		= TRunoff%d_ns(iunit)					! carry number of layers if there is merging/split
				d_res_sub 	= TRunoff%d_lake(iunit)
				do j = 1, nlayers+1	
					TRunoff%v_zt0(iunit,j) = TRunoff%v_zt(iunit,j)
					TRunoff%a_d0(iunit,j) = TRunoff%a_d(iunit,j)
					TRunoff%d_z0(iunit,j) = TRunoff%d_z(iunit,j)
				end do
				phi_o  		= 0._r8 !
				sh_net 		= 0._r8			
				do j = 1,nlayers   
					dv_nt(j)	= 0._r8
					dv_in(j)	= 0._r8
					dv_ou(j)	= 0._r8
					enr_0(j)	= 0._r8
					enr_in(j)	= 0._r8
					enr_ou(j)	= 0._r8
					enr_1(j)	= 0._r8
					enr_2(j)	= 0._r8
					enr_phi(j)	= 0._r8
					phi_z(j)	= 0._r8
					rho_z(j)	= 0._r8
					a(j)		= 0._r8
					b(j)		= 0._r8
					c(j)		= 0._r8
					r(j)		= 0._r8
					Fx(j)		= 0._r8
				end do
							
				! 	Intitialize/reassign layer storage
				!	Assign layer temperature and calculate reservoir density at depth z and incoming flow 	
				! 	Allocate layer temperature as old for assigning counter during averaging sub-timestep result		
				rho_r = den(t_in)				
				do j = 1, nlayers	
					if (j<=TRunoff%d_ns(iunit)) then
						TRunoff%temp_lake(iunit,j) = TRunoff%temp_lake(iunit,j)
						rho_z(j) = den(TRunoff%temp_lake(iunit,j))
						t_z_old(j) = TRunoff%temp_lake(iunit,j)
					else 
						TRunoff%temp_lake(iunit,j) = 0._r8
						rho_z(j) = 0._r8
						t_z_old(j) = 0._r8
					end if
				end do	
				
			! 	Calculation of Surface fluxes and heat source 
				!	Net shortwave radiation (w/m^2)
				if (THeat%coszen(iunit)>0._r8) then 
					alb_s = 0.05_r8/(THeat%coszen(iunit) + 0.15_r8)!THeat%albedo(iunit)!
				else
					alb_s = 0.06_r8	!0.06_r8
				end if
				! alb_s = 0.03_r8	!0.06_r8	
				sh_net 	= max(THeat%forc_solar(iunit)*(1._r8 - alb_s),0._r8)!			
				t_s    	= TRunoff%temp_lake(iunit,TRunoff%d_ns(iunit)) 
				lw_abr 	= (1._r8 - 0.03_r8)*THeat%forc_lwrad(iunit) !
				call QSat_str(t_s, THeat%forc_pbot(iunit),es)				
				lw_ems  = 0.97_r8*st_bl*t_s**4_r8    
				sn_heat = 1.5701_r8*THeat%forc_wind(iunit)* (t_s - THeat%forc_t(iunit))
				le 		= 1000._r8*(2499.64_r8 - 2.51_r8 * (t_s-273.15_r8)) 
				kl    	= 0.211_r8 + 0.103_r8 * THeat%forc_wind(iunit)* F
				
				evap  	= max(kl*(es - THeat%forc_vp(iunit))/100._r8,0._r8)  ! in mm/d
				
				if (t_s > t_frz) then
					lt_heat = rho_w * evap* le/(86.4e6)	!latent heat (w/m^2)
				else 
					lt_heat = 0.0_r8
				end if
				! 
				! 	Net surface heat flux  
				phi_o = sh_net + lw_abr - (lt_heat + lw_ems + sn_heat)!
				
		! 	Calculation of equivalent evaporated depth (m), and volume(m^3)
				if (q_in==0._r8 .or. q_in==q_ou)evap=0._r8 	!Avoid continuous water abstraction if there is no net inflow. Assumption is precipitation balances evaporation on longterm 
				d_evap 	= evap*dtime/(86.4e6)
				if(TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1) > 0.10_r8*(TRunoff%v_zti(iunit,ngeom+1)).and. q_ou < q_in) then 
					v_evap 	= max(d_evap*sar*TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1),0._r8)
				else
					v_evap	= 0._r8
				end if				
			
		! 	Calculation of flow contibution due to inflow/outflow 				
				! 	Calculate total mass (kg), layer energy (w)
				do j = 1, nlayers	
					if (j<=TRunoff%d_ns(iunit))  then
						TRunoff%v_zo(iunit,j) = TRunoff%d_v(iunit,j)
						enr_0(j) = TRunoff%temp_lake(iunit,j)*TRunoff%d_v(iunit,j)*rho_z(j)*c_w/dtime
					else
						TRunoff%v_zo(iunit,j) = 0._r8		
						enr_0(j) = 0._r8
					end if
				end do
				
				if (abs(q_in - q_ou)*dtime <= 0.5_r8*TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1) .and. TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1) > 0.6_r8*TRunoff%v_zti(iunit,ngeom+1) .and. TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1) < TRunoff%v_zti(iunit,ngeom+1)) then  ! Redistribute mass flux only if there is net flow advection and enough storage
				
					call flowdist(TRunoff%d_ns(iunit),q_in,t_in,q_ou,dtime,iunit,dv_in,dv_ou) 
					
					!	Resize layer thickness and numbers based on inflow/outflow contribution							
					do j = 1, nlayers
						if (j<=TRunoff%d_ns(iunit)-1 .and. TRunoff%d_ns(iunit)>1) then
							dv_nt(j)	= (dv_in(j)-dv_ou(j)-v_evap/(TRunoff%d_ns(iunit)-1))
						elseif(j==TRunoff%d_ns(iunit) .and. TRunoff%d_ns(iunit)>1) then
							dv_nt(j)	= (dv_in(j)-dv_ou(j))
						elseif(j==TRunoff%d_ns(iunit) .and. TRunoff%d_ns(iunit)==1) then
							dv_nt(j)	= (dv_in(j)-dv_ou(j)-v_evap)
						else
							dv_nt(j)	= 0._r8
						end if
					end do
					
					! 	Calculate layer volume (m3)
					do j = 1, nlayers	
						if (j<=TRunoff%d_ns(iunit))  then   
							TRunoff%v_zn(iunit,j) = ((TRunoff%d_v(iunit,j)+dv_nt(j)))
						else
							TRunoff%v_zn(iunit,j) = 0._r8
						end if
					end do
								
				! Calculate new layer thickness, depth, area/volume after net advection		
					call geometry(iunit,dv_nt)					
					
				! Merge/split layers					
					if (TRunoff%d_lake(iunit) <= 2.0_r8 .or. TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1)<=TRunoff%v_zti(iunit,1) .or. TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1)>TRunoff%v_zti(iunit,ngeom+1)) then ! skip if storage is too small/large
						call over_under_flow(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)						
					else
						ddz_min = 1.5_r8
						ddz_max = max(3.0_r8*TRunoff%ddz_local(iunit),5.0_r8)
						!	Check if layers are too small   
						if (TRunoff%d_ns(iunit) >=2) then
							do i=1,TRunoff%d_ns(iunit)-1
								if (TRunoff%dd_z(iunit,i) < ddz_min .and. TRunoff%d_ns(iunit)>1) then
									call layer_merge(iunit,i,dv_in,dv_ou,enr_0)
									if (TRunoff%d_ns(iunit)<1) then
										exit
									end if
								elseif (TRunoff%d_ns(iunit)==1) then
									call over_under_flow(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
									exit
								end if								
							end do
						end if						
						! 	Check if layers are too big
						if (TRunoff%d_ns(iunit) <= nlayers-1) then !
							do i=1,TRunoff%d_ns(iunit)
								if(TRunoff%dd_z(iunit,i) > ddz_max .and. TRunoff%d_ns(iunit)<nlayers) then 
									call layer_split(iunit,i,dv_in,dv_ou,enr_0)
								elseif (TRunoff%dd_z(iunit,i) > ddz_max .and. i<TRunoff%d_ns(iunit) .and. TRunoff%d_ns(iunit)==nlayers) then  ! allow further merging by increasing ddz_min
									ddz_min = max(0.5*TRunoff%ddz_local(iunit),3.0_r8)
									call layer_merge(iunit,i,dv_in,dv_ou,enr_0)
									if (TRunoff%d_ns(iunit)<1) then
										exit
									end if
									call layer_split(iunit,i,dv_in,dv_ou,enr_0)
									if (TRunoff%d_ns(iunit)==nlayers) then
										exit
									end if
								elseif(i==TRunoff%d_ns(iunit) .and. TRunoff%d_ns(iunit)==nlayers) then
									call over_under_flow(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
									exit
								end if								
							end do
						end if					
					end if		!layer merge/split
					
					!	Recalculare reservoir depth
					TRunoff%d_lake(iunit)=0._r8 
					do i=1,TRunoff%d_ns(iunit)
					   TRunoff%d_lake(iunit) = TRunoff%d_lake(iunit)+TRunoff%dd_z(iunit,i)
					end do
					! 	Recalculate density after layer change
					do j = 1, nlayers	
						if (j<=TRunoff%d_ns(iunit))  then
							rho_z(j) = den(TRunoff%temp_lake(iunit,j))
						else
							rho_z(j) = 0._r8
						end if
					end do
					
					! 	Calculate layer internal energy (w) due to inflow/outflow		
					do j = 1, nlayers	
						if (j<=TRunoff%d_ns(iunit))  then   
							enr_in(j) = (dv_in(j)*t_in*rho_r*c_w/dtime)		!Energy from inflow
							enr_ou(j) = (dv_ou(j)*TRunoff%temp_lake(iunit,j)*rho_z(j)*c_w/dtime)	!Energy loss due to outflow
						else
							enr_in(j) = 0._r8
							enr_ou(j) = 0._r8
						end if
					end do	
					
					if(TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1) <= 0.050_r8*TRunoff%v_zti(iunit,ngeom+1) .or. TRunoff%d_v(iunit,1) <= 0._r8 .or. (TRunoff%d_ns(iunit) ==1 .and. TRunoff%d_lake(iunit)<=2._r8)) then	! Reservoir virtually empty/shallow
						do j = 1, nlayers	
							if (j<=TRunoff%d_ns(iunit)) then
								TRunoff%temp_lake(iunit,j) = THeat%Tr(iunit)
							else 
								TRunoff%temp_lake(iunit,j) = 0._r8
							end if
						end do
					end if 				
					
				end if ! layer recalculation
	! *********************************************************************************************************************************
			! 	Calculate layer energy (w)
				do j = 1, nlayers	
					if (j<=TRunoff%d_ns(iunit))  then   
						enr_1(j) = enr_0(j)+ enr_in(j) - enr_ou(j)
					else
						enr_1(j) = 0._r8
					end if
				end do	
				
				! 	Adjust layer temperature for mixing due to inflow/outflow
				do j = 1, nlayers	
					if (j<=TRunoff%d_ns(iunit))  then   
						TRunoff%temp_lake(iunit,j) = enr_1(j)*dtime/(TRunoff%d_v(iunit,j)*rho_z(j)*c_w)
					else
						TRunoff%temp_lake(iunit,j) = 0._r8
					end if
				end do
				
			!	Check energy balance (w) after advective mixing	
				enr_err1 = (sum(enr_1) - (sum(enr_0)+ sum(enr_in) - sum(enr_ou)))
				
			!******************************************************************************
			! 	Calculate solar energy absorbed at each layer		
				eta  = 1.1925_r8*(TRunoff%d_lake(iunit))**(-0.424_r8) ! as used in Subin et al, 2011 (citing Hakanson, 1995) but modified for actual reservoir depth 
				beta = 0.175_r8 !
		
				do j=1,nlayers+1
					phi_x(j) = 0._r8
				end do
				
				if (sh_net > 0._r8 .and. TRunoff%d_ns(iunit) >1) then
					k=0
					top_d=TRunoff%d_lake(iunit)-TRunoff%d_z(iunit,TRunoff%d_ns(iunit)-k)
					if(top_d < 0.61_r8) then
						k=k+1
					else
						k=k+2
					end if	
					v_mix=(TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1)-TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)-k))	
					!	Solar radiation energy absorbed at the mixed zone
					sh_mix=sh_net*sar*(TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)-(1._r8-beta)*TRunoff%a_d(iunit,TRunoff%d_ns(iunit)-k))
					j=TRunoff%d_ns(iunit)-k
					do i=j,TRunoff%d_ns(iunit)
					   phi_z(i)=sh_mix*TRunoff%d_v(iunit,i)/v_mix
					end do
					phi_x(TRunoff%d_ns(iunit)+1)=sh_net
					phi_x(TRunoff%d_ns(iunit)-k)=(1._r8-beta)*sh_net
					if(k>0) then
						do i=1,k
						   ii=TRunoff%d_ns(iunit)-i+1
						   phi_x(ii)=(TRunoff%a_d(iunit,ii+1)*phi_x(ii+1)-phi_z(ii))/(TRunoff%a_d(iunit,ii))
						end do
					end if
					
					!	Solar radiation energy absorbed at sub layers
					l=TRunoff%d_ns(iunit)-k-1
					do j=1,l
					   i=(TRunoff%d_ns(iunit)-k)-j+1
					   phi_x(i-1)=phi_x(i)*exp(-eta*(TRunoff%d_z(iunit,i)-TRunoff%d_z(iunit,i-1)))
					end do
				
					!	Solar radiation energy absorbed in each layer
					j=TRunoff%d_ns(iunit)-k-1
					do i=1,j
					   phi_z(i)=(sar*TRunoff%a_d(iunit,i+1)*phi_x(i+1)-sar*TRunoff%a_d(iunit,i)*phi_x(i))
					end do
				elseif (sh_net > 0._r8 .and. TRunoff%d_ns(iunit)==1) then
					phi_z(TRunoff%d_ns(iunit))=sh_net*sar*TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)
				else
					do j=1,TRunoff%d_ns(iunit)
						phi_z(j) = 0._r8
					end do	
				end if		    
				
			! *********************************************************************************************************************************
			! 	Calculation of effective diffusion coefficient 	, Herb and Stefan (2004) citing Wu (1971)
				if (THeat%forc_wind(iunit) >= 15._r8) then
					c_d = 2.6e-3
				else
					c_d = 5.e-4*sqrt(THeat%forc_wind(iunit))
				end if
				tau = rho_a*c_d*THeat%forc_wind(iunit)**2._r8 ! Shear stress at surface
				s_vel = sqrt(tau/rho_w) ! Shear velocity at surface
				k_ew=tau*s_vel*sar*TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)*dtime	! Wind driven kinetic energy at surface
				dis_w = k_ew/(rho_w*TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1)*dtime)  ! rate of dissipation-wind
				cfw = 1.e-02!
				cfa = 1.e-05!
				k_m = 0.57_r8/(c_w*rho_w) 						  !molecular diffusivity
				df_eff(1)= 0._r8			!bottom interface
				df_eff(TRunoff%d_ns(iunit)+1)= 0._r8		!air interface
				do j = 2,TRunoff%d_ns(iunit) 
					q_adv(j) = max((dv_in(j)+dv_ou(j)),0._r8)
					k_ad(j)=0.5_r8*rho_w*q_adv(j)*dtime*(q_adv(j)/(TRunoff%Width_r(iunit)*TRunoff%dd_z(iunit,j)))**2._r8 ! Advection driven kinetic energy 
					dis_ad(j)= k_ad(j)/(rho_w*TRunoff%v_zt(iunit,j)*dtime)	! rate of dissipation-inflow/outflow					
				! Calculate Richardson number
					drhodz(j) = (rho_z(j-1)-rho_z(j))/0.5_r8*(TRunoff%dd_z(iunit,j)+TRunoff%dd_z(iunit,j-1))
					bv_f = max((grav/rho_w)*drhodz(j),0._r8)
					if (s_vel <= 0._r8) ri = 0._r8
					ri = bv_f/((s_vel/(0.4_r8*TRunoff%d_z(iunit,j)))**2._r8)				
				! Calculate Froude number
					l_vel = q_adv(j)*TRunoff%Length_r(iunit)/(sar*TRunoff%a_d(iunit,j)*TRunoff%dd_z(iunit,j))
					if (q_adv(j) <= 0._r8 .or. drhodz(j) <= 0._r8) Fr(j) = 0._r8
					Fr(j)= (grav*TRunoff%dd_z(iunit,j)*drhodz(j)/rho_w)/l_vel**2._r8					
				! Calculate diffusion coefficients				
					df_eff(j)=min(max(dtime**2._r8*((cfw*dis_w/(1+ri))+(0.5_r8*cfa*(dis_ad(j)+dis_ad(j-1))/(1+Fr(j)))),k_m),5.56e-03) !5.56e-03!
				end do
			
			!*****************************************************************
			! Calculate matrix elements
				do j = 1,TRunoff%d_ns(iunit)
					if (j == 1 .and. TRunoff%d_ns(iunit)>1) then
						m1(j) = 2_r8*dtime/(0.5_r8*(sar*TRunoff%a_d(iunit,j)+sar*TRunoff%a_d(iunit,j+1))*TRunoff%dd_z(iunit,j))
						m2(j) = m1(j)*sar*TRunoff%a_d(iunit,j+1)*df_eff(j+1)/(TRunoff%dd_z(iunit,j)+TRunoff%dd_z(iunit,j+1))
						m3(j) = 0._r8						
						Fx(j) = dtime*phi_z(j)/(TRunoff%d_v(iunit,j)*c_w*rho_z(j))
						a(j) = - (m2(j))
						b(j) = 1._r8 + (m2(j) + m3(j)) 
						c(j) = 0._r8 
						r(j) = TRunoff%temp_lake(iunit,j) + Fx(j) ! bottom boundary condition 
					elseif (j <= TRunoff%d_ns(iunit)-1 .and. TRunoff%d_ns(iunit)>2) then
						m1(j) = 2_r8*dtime/(0.5_r8*(sar*TRunoff%a_d(iunit,j)+sar*TRunoff%a_d(iunit,j+1))*TRunoff%dd_z(iunit,j))
						m2(j) = m1(j)*sar*TRunoff%a_d(iunit,j+1)*df_eff(j+1)/(TRunoff%dd_z(iunit,j)+TRunoff%dd_z(iunit,j+1))
						m3(j) = m1(j)*sar*TRunoff%a_d(iunit,j)*df_eff(j)/(TRunoff%dd_z(iunit,j)+TRunoff%dd_z(iunit,j-1))						
						Fx(j) = dtime*phi_z(j)/(TRunoff%d_v(iunit,j)*c_w*rho_z(j))
						a(j) = - m2(j)
						b(j) = 1._r8 + m2(j) + m3(j) 
						c(j) = - m3(j)
						r(j) = TRunoff%temp_lake(iunit,j) + Fx(j)
					elseif (j == TRunoff%d_ns(iunit)) then!top layer
						m1(j) = 2_r8*dtime/(0.5_r8*(sar*TRunoff%a_d(iunit,j)+sar*TRunoff%a_d(iunit,j+1))*TRunoff%dd_z(iunit,j))
						m2(j) = 0._r8
						m3(j) = m1(j)*sar*TRunoff%a_d(iunit,j)*df_eff(j)/(TRunoff%dd_z(iunit,j)+TRunoff%dd_z(iunit,j-1))
						Fx(j) = dtime*((phi_o-sh_net)*sar*TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)+phi_z(j))/(TRunoff%d_v(iunit,TRunoff%d_ns(iunit))*c_w*rho_z(j)) ! 
						a(j) = 0._r8
						b(j) = 1._r8 + (m2(j) + m3(j)) 
						c(j) = - (m3(j))
						r(j) = TRunoff%temp_lake(iunit,j) + Fx(j) ! top boundary condition 						
					end if
				end do	
				
				
			!	Solve for temperature
				call solve(a,b,c,r,iunit,TRunoff%d_ns(iunit))
				
            ! 	Check numerical instability from reservoir geometry data 	
				do i = 1,TRunoff%d_ns(iunit)
					if (isnan(TRunoff%temp_lake(iunit,i)) .or. TRunoff%temp_lake(iunit,i)>=huge(1._r8)) then
						TRunoff%temp_lake(iunit,i) = THeat%Tr(iunit)
						! To do: Check individual lake geometry/inflow/outflow 
						! write(iulog,*) subname,'Check geometry',iunit					
						return
					end if	
				end do
			
				! Calculate layer intermediate internal energy (w)
				do j = 1, nlayers	
					if (j<TRunoff%d_ns(iunit))  then   
						enr_2(j)   = TRunoff%temp_lake(iunit,j)*TRunoff%d_v(iunit,j)*rho_z(j)*c_w/dtime
						enr_phi(j) = phi_z(j)
					elseif (j==TRunoff%d_ns(iunit))  then
						enr_2(j)   = TRunoff%temp_lake(iunit,j)*TRunoff%d_v(iunit,j)*rho_z(j)*c_w/dtime
						enr_phi(j) = (phi_o-sh_net)*sar*TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)		!Enegry from net surface flux (w)
					else
						enr_2(j) = 0.
					end if
				end do				
			
		!	Check energy balance error (w) after triadiagonal matrix solution	
				enr_err2 = (sum(enr_2) - (sum(enr_1) + sum(enr_phi)))
				
		!***********************************************************************************************************************		
		! !	Solve convective mixing
			! Recalculate layer density 
				do j = 1,TRunoff%d_ns(iunit)   
					rho_z(j) = den(TRunoff%temp_lake(iunit,j))
				end do
				
			! Check if instability exists
				do j = 1,TRunoff%d_ns(iunit)-1 
					if(rho_z(j) < rho_z(j+1)) then
						! Start mixing layers
						mix1=j
						mix2=mix1
						sumvol=TRunoff%d_v(iunit,mix2)*1._r8
						tsum=TRunoff%temp_lake(iunit,mix2)*TRunoff%d_v(iunit,mix2)	
						
						mix2=mix2+1
						mixvol2=TRunoff%d_v(iunit,mix2)*1._r8
						sumvol=sumvol + mixvol2
						tsum=tsum+TRunoff%temp_lake(iunit,mix2)*mixvol2
						tmix=tsum/sumvol
						
						! Calculate density of mixed layer	
						denmix = den(tmix)
						
						! Check if instability exists below mixed layer	and mix layers	
						if(rho_z(mix1-1) < denmix .and. mix1 >= 2) then
							mix1=mix1-1	
						
							! Calculate temperature of mixed layer	
							mixvol1=TRunoff%d_v(iunit,mix1)*1._r8
							sumvol=sumvol + mixvol1
							tsum=tsum+TRunoff%temp_lake(iunit,mix1)*mixvol1
							tmix= tsum/sumvol
						end if
						
						! Calculate density of mixed layer	
						denmix = den(tmix)					
				  
						! Set new layer temperature and density 
						do i=mix1,mix2
							rho_z(i)=denmix
							TRunoff%temp_lake(iunit,i)=tmix
						end do
					end if
				end do
	! *********************************************************************************************************************************
			
			! 	Sum sub-timestep variables			
				d_res_sub = (d_res_sub + TRunoff%d_lake(iunit))/2.0_r8
				if (d_n_n == TRunoff%d_ns(iunit)) then
					do j = 1,nlayers 
						TRunoff%temp_lake(iunit,j) = (t_z_old(j) + TRunoff%temp_lake(iunit,j))/2.0_r8
					end do 
				else 
					do j = 1,nlayers
						if (j<=TRunoff%d_ns(iunit))	then		
							TRunoff%temp_lake(iunit,j) = (t_z_old(j) + TRunoff%temp_lake(iunit,j))/2.0_r8
						else
							TRunoff%temp_lake(iunit,j) = 0._r8
						end if
					end do 
				end if
				
				TRunoff%d_lake(iunit) = d_res_sub
				
			end do ! sub-timestep
						
			t_out = TRunoff%temp_lake(iunit,TRunoff%d_ns(iunit))
			if (t_out < 273.15_r8) t_out = 273.15_r8 
			
			THeat%Tr(iunit)= t_out
					
		end if ! stratification
    end subroutine mosart_lake

	subroutine flowdist(d_n,in_f,in_t,ou_f,dtime,iunit,dv_in,dv_ou)
!*******************************************************************************************************
! 	Calculation inflow/outflow contribution modified to distribute uniformly
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use RtmVar         , only : iulog, ngeom, nlayers
		implicit none
		integer, intent(in)  :: d_n,iunit
		real(r8) :: c_w = 4.188e3     		   			! Water specific heat capacity  (J/kg/k)
		real(r8), intent(in) :: dtime 
		real(r8), intent(in)  :: in_f,in_t,ou_f
		real(r8),dimension(nlayers), intent(out) :: dv_in,dv_ou    ! layer inflow/outflow (m3/s) 
		real(r8) :: rho_r,in_v,ou_v,m_in,denom    !
		integer :: j,jmax,jmin							! indices
	
		rho_r = 1000._r8*( 1.0_r8 - 1.9549e-05*(abs(in_t-277._r8))**1.68_r8)
		in_v = in_f*dtime
		ou_v = ou_f*dtime
		m_in = in_v*rho_r
		
	!	Initialize
		do j=1,nlayers
			dv_in(j)=0._r8
			dv_ou(j)=0._r8
		end do
				
!   Layer inflow and energy contribution, keep top layer constant by avoiding inflow/outflow
		jmin=1
		if(d_n>3)jmax=d_n-2 
		if(d_n<=3)jmax=d_n-1
		
		if(d_n==1) then ! Single layer
			do j=1,nlayers
				dv_in(j) = 0._r8
				dv_ou(j) = 0._r8
			end do
			dv_in(d_n) = in_v
			dv_ou(d_n) = ou_v
		else			
			do j=jmin,jmax
				denom=max((TRunoff%v_zt(iunit,jmax+1)-TRunoff%v_zt(iunit,jmin)),TRunoff%d_v(iunit,j))
				dv_in(j) = in_v*(TRunoff%d_v(iunit,j)/(denom))
				dv_ou(j) = ou_v*(TRunoff%d_v(iunit,j)/(denom))
			end do
		end if
		
	end subroutine flowdist


	subroutine geometry(iunit,dv_nt)
!*******************************************************************************************************
! 	Calculate new layer thickness, depth, area/volume after net advection
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: iunit              
		real(r8),dimension(nlayers), intent(in) :: dv_nt    ! layer net inflow/outflow (m3/s)
		real(r8) :: num,dem,delta_z,delta_a    !
		integer :: i,j,k,mm							! indices
	
		TRunoff%d_z(iunit,1)=0._r8
		TRunoff%a_d(iunit,1)=TRunoff%a_di(iunit,1)
		TRunoff%v_zt(iunit,1)=TRunoff%v_zti(iunit,1)								
		do i=1,nlayers 			! check layers for available volume to satisfy net outflow
			if (i<=TRunoff%d_ns(iunit))  then
				if(-dv_nt(i) > TRunoff%d_v(iunit,i))then !current layer collapses, hence remaining volume taken from next upper/lower layer
					TRunoff%v_zn(iunit,i)=0._r8
					if (i<TRunoff%d_ns(iunit)-1 .and. TRunoff%d_ns(iunit)>1) then
						TRunoff%v_zn(iunit,i+1)=((TRunoff%v_zn(iunit,i+1) - (-dv_nt(i)-TRunoff%d_v(iunit,i))))			
						if (TRunoff%v_zn(iunit,i+1)<0._r8) TRunoff%v_zn(iunit,i+1) = 0._r8								
					elseif (i==TRunoff%d_ns(iunit)-1 .and. TRunoff%d_ns(iunit)>2) then	!avoid top layer collapse, hence remaining mass taken from next lower layer
						TRunoff%v_zn(iunit,i-1)=((TRunoff%v_zn(iunit,i-1) - (-dv_nt(i)-TRunoff%d_v(iunit,i))))
						if (TRunoff%v_zn(iunit,i-1)<0._r8) TRunoff%v_zn(iunit,i-1) = 0._r8								
						mm=i-1
						do k=mm,TRunoff%d_ns(iunit)	
							TRunoff%v_zt(iunit,k+1)=(TRunoff%v_zt(iunit,k)+(TRunoff%v_zn(iunit,k)))
							do j=2,ngeom+1
								if (TRunoff%v_zt(iunit,k+1)>TRunoff%v_zti(iunit,j-1).and.TRunoff%v_zt(iunit,k+1)<=TRunoff%v_zti(iunit,j))then
									num = (TRunoff%v_zt(iunit,k+1)-TRunoff%v_zti(iunit,j-1))
									dem = (TRunoff%v_zti(iunit,j)-TRunoff%v_zti(iunit,j-1))
									delta_z  = (TRunoff%d_zi(iunit,j)-TRunoff%d_zi(iunit,j-1))*num/dem
									delta_a  = (TRunoff%a_di(iunit,j)-TRunoff%a_di(iunit,j-1))*num/dem
									TRunoff%d_z(iunit,k+1) = TRunoff%d_zi(iunit,j-1) + delta_z 
									TRunoff%a_d(iunit,k+1) = ((TRunoff%a_di(iunit,j-1) + max((delta_a),0.0_r8)))
								elseif (TRunoff%v_zt(iunit,k+1)>TRunoff%v_zti(iunit,ngeom+1))then ! inflow causes maximum storage
									num = (TRunoff%v_zt(iunit,k+1)-TRunoff%v_zti(iunit,ngeom))
									dem = (TRunoff%v_zti(iunit,ngeom+1)-TRunoff%v_zti(iunit,ngeom))
									delta_z  = (TRunoff%d_zi(iunit,ngeom+1)-TRunoff%d_zi(iunit,ngeom))*num/dem
									delta_a  = (TRunoff%a_di(iunit,TRunoff%d_ns(iunit)+1) -TRunoff%a_di(iunit,TRunoff%d_ns(iunit)))*num/dem
									TRunoff%d_z(iunit,k+1) = TRunoff%d_zi(iunit,ngeom) + delta_z
									TRunoff%a_d(iunit,k+1) = ((TRunoff%a_di(iunit,ngeom) + max((delta_a),0.0_r8)))
								end if
							end do
							TRunoff%d_z(iunit,k)=TRunoff%d_z(iunit,k+1)
							TRunoff%dd_z(iunit,k)=TRunoff%d_z(iunit,k+1)-TRunoff%d_z(iunit,k)
						end do		
					elseif ((i==TRunoff%d_ns(iunit)-1 .or. i==TRunoff%d_ns(iunit)) .and. TRunoff%d_ns(iunit)<=2) then	! Top layer collapses, skip outflow							
						if (TRunoff%v_zn(iunit,i)<0._r8) TRunoff%v_zn(iunit,i) = TRunoff%d_v(iunit,i)								
					end if
					TRunoff%v_zt(iunit,i+1)=(TRunoff%v_zt(iunit,i)+(TRunoff%v_zn(iunit,i)))
					do j=2,ngeom+1
						if (TRunoff%v_zt(iunit,i+1)>TRunoff%v_zti(iunit,j-1).and.TRunoff%v_zt(iunit,i+1)<=TRunoff%v_zti(iunit,j))then
							num = (TRunoff%v_zt(iunit,i+1)-TRunoff%v_zti(iunit,j-1))
							dem = (TRunoff%v_zti(iunit,j)-TRunoff%v_zti(iunit,j-1))
							delta_z  = (TRunoff%d_zi(iunit,j)-TRunoff%d_zi(iunit,j-1))*num/dem
							delta_a  = (TRunoff%a_di(iunit,j)-TRunoff%a_di(iunit,j-1))*num/dem
							TRunoff%d_z(iunit,i+1) = TRunoff%d_zi(iunit,j-1) + delta_z 
							TRunoff%a_d(iunit,i+1) = ((TRunoff%a_di(iunit,j-1) + max((delta_a),0.0_r8)))
						elseif (TRunoff%v_zt(iunit,i+1)>TRunoff%v_zti(iunit,ngeom+1))then! inflow causes maximum storage
							num = (TRunoff%v_zt(iunit,i+1)-TRunoff%v_zti(iunit,ngeom))
							dem = (TRunoff%v_zti(iunit,ngeom+1)-TRunoff%v_zti(iunit,ngeom))
							delta_z  = (TRunoff%d_zi(iunit,ngeom+1)-TRunoff%d_zi(iunit,ngeom))*num/dem
							delta_a  = (TRunoff%a_di(iunit,TRunoff%d_ns(iunit)+1)*1._r8 -TRunoff%a_di(iunit,TRunoff%d_ns(iunit)))*num/dem
							TRunoff%d_z(iunit,i+1) = TRunoff%d_zi(iunit,ngeom) + delta_z
							TRunoff%a_d(iunit,i+1) = ((TRunoff%a_di(iunit,ngeom) + max((delta_a),0.0_r8)))
						end if
					end do
					TRunoff%d_z(iunit,i)=TRunoff%d_z(iunit,i+1)
					TRunoff%dd_z(iunit,i)=max(TRunoff%d_z(iunit,i+1)-TRunoff%d_z(iunit,i),0.0_r8)
				else !enough volume, layers don't collapses
					TRunoff%v_zt(iunit,i+1)=(TRunoff%v_zt(iunit,i) + TRunoff%v_zn(iunit,i))!dv_nt(i)
					if (TRunoff%v_zt(iunit,i+1) < 1.05_r8*TRunoff%v_zti(iunit,ngeom+1)) then ! Allowable maximum storage should include only 5% (freeboard provision)
						do j=2,ngeom+1
							if (TRunoff%v_zt(iunit,i+1)>TRunoff%v_zti(iunit,j-1).and.TRunoff%v_zt(iunit,i+1)<=TRunoff%v_zti(iunit,j))then
								num = (TRunoff%v_zt(iunit,i+1)-TRunoff%v_zti(iunit,j-1))
								dem = (TRunoff%v_zti(iunit,j)-TRunoff%v_zti(iunit,j-1))
								delta_z  = (TRunoff%d_zi(iunit,j)-TRunoff%d_zi(iunit,j-1))*num/dem
								delta_a  = (TRunoff%a_di(iunit,j)-TRunoff%a_di(iunit,j-1))*num/dem
								TRunoff%d_z(iunit,i+1) = TRunoff%d_zi(iunit,j-1) + delta_z 
								TRunoff%a_d(iunit,i+1) = ((TRunoff%a_di(iunit,j-1) + max((delta_a),0.0_r8)))
							elseif (TRunoff%v_zt(iunit,i+1)>TRunoff%v_zti(iunit,ngeom+1))then ! inflow causes maximum storage
								num = (TRunoff%v_zt(iunit,i+1)-TRunoff%v_zti(iunit,ngeom))
								dem = (TRunoff%v_zti(iunit,ngeom+1)-TRunoff%v_zti(iunit,ngeom))
								delta_z  = (TRunoff%d_zi(iunit,ngeom+1)-TRunoff%d_zi(iunit,ngeom))*num/dem
								delta_a  = (TRunoff%a_di(iunit,TRunoff%d_ns(iunit)+1) -TRunoff%a_di(iunit,TRunoff%d_ns(iunit)))*num/dem
								TRunoff%d_z(iunit,i+1) = TRunoff%d_zi(iunit,ngeom) + delta_z
								TRunoff%a_d(iunit,i+1) = ((TRunoff%a_di(iunit,ngeom) + max((delta_a),0.0_r8)))
							end if
						end do
					else
						TRunoff%d_z(iunit,i+1)  = 1.05_r8*TRunoff%d_zi(iunit,ngeom+1)
						TRunoff%a_d(iunit,i+1)  = 1.05_r8*TRunoff%a_di(iunit,ngeom+1)
						TRunoff%v_zt(iunit,i+1) = 1.05_r8*TRunoff%v_zti(iunit,ngeom+1)
					end if
				end if 
				TRunoff%dd_z(iunit,i)=max(TRunoff%d_z(iunit,i+1)-TRunoff%d_z(iunit,i),0.0_r8)
			elseif (i==TRunoff%d_ns(iunit)+1)then
				TRunoff%d_z(iunit,i)  = TRunoff%d_z(iunit,i)
				TRunoff%a_d(iunit,i)  = TRunoff%a_d(iunit,i)
				TRunoff%v_zt(iunit,i) = TRunoff%v_zt(iunit,i)
				TRunoff%dd_z(iunit,i) = 0._r8
			else
				TRunoff%d_z(iunit,i)  = 0._r8
				TRunoff%a_d(iunit,i)  = 0._r8
				TRunoff%v_zt(iunit,i) = 0._r8
				TRunoff%dd_z(iunit,i) = 0._r8
			end if
		end do
		
		! 	Recalculate layer thickness and volume	
		do j = 1,nlayers
			if (j<=TRunoff%d_ns(iunit))  then   
				TRunoff%d_v(iunit,j) 	= max((TRunoff%v_zt(iunit,j+1) - TRunoff%v_zt(iunit,j)),0.0_r8)
				TRunoff%v_zn(iunit,j) 	= TRunoff%d_v(iunit,j)
			else
				TRunoff%d_v(iunit,j)	= 0._r8
				TRunoff%v_zn(iunit,j) 	= 0._r8
			end if
		end do
		
		!	Recalculare reservoir depth
		TRunoff%d_lake(iunit)=0._r8 
		do i=1,TRunoff%d_ns(iunit)
			TRunoff%d_lake(iunit) = TRunoff%d_lake(iunit)+TRunoff%dd_z(iunit,i)
		end do	
		
	end subroutine geometry

	subroutine layer_merge(iunit,i,dv_in,dv_ou,enr_0)
!*******************************************************************************************************
! 	Merge layer if the thickness is less than threshold ddz_min 
!*******************************************************************************************************   
		use shr_sys_mod , only : shr_sys_flush
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: iunit,i              
		real(r8),dimension(nlayers), intent(inout) :: dv_in,dv_ou,enr_0    ! layer inflow,outflow(m3/s), and energy(w) 
		real(r8) :: ta,tb,dd_za,dd_zb,e_a,e_b,d_va,d_vb
		real(r8) :: dv_oua,dv_oub,dv_ina,dv_inb    
		integer :: j,k,m							! indices				
		
		if(i<TRunoff%d_ns(iunit)-1 .or. (i==1 .and. TRunoff%d_ns(iunit)==2)) then
			m=i+1
		else 
			m=i-1 !Avoid merging to top layer
		end if
		ta 		= TRunoff%temp_lake(iunit,i)
		tb 		= TRunoff%temp_lake(iunit,m)
		dd_za	= TRunoff%dd_z(iunit,i)
		dd_zb	= TRunoff%dd_z(iunit,m)
		e_a 	= enr_0(i)
		e_b 	= enr_0(m)
		d_va	= TRunoff%d_v(iunit,i)
		d_vb	= TRunoff%d_v(iunit,m)
		dv_oua	= dv_ou(i)
		dv_oub	= dv_ou(m)
		dv_ina	= dv_in(i)
		dv_inb	= dv_in(m)
		
	!	Merge layers
		if(d_va+d_vb > 0._r8)then
			TRunoff%temp_lake(iunit,i)	= (d_va*ta+d_vb*tb)/(d_va+d_vb)
		else
		end if
	
	! 	Adjust new layer volume, thickness, mass and inflow/outflow 
		TRunoff%dd_z(iunit,i)	= dd_za+dd_zb
		TRunoff%d_v(iunit,i)	= d_va+d_vb
		dv_ou(i)				= dv_oua+dv_oub
		dv_in(i)				= dv_ina+dv_inb
		enr_0(i)				= e_a+e_b

	!	Re-number layers before collapse
		do j=i,nlayers-1
			if (j==i .and. i<(TRunoff%d_ns(iunit)-1)) then ! 	Identify lower (small and to be merged) and upper (larger) layer 
				dv_ou(j)	= dv_ou(i)
				dv_in(j)	= dv_in(i)
				enr_0(j)	= enr_0(i)
				TRunoff%temp_lake(iunit,j)	= TRunoff%temp_lake(iunit,i)
				TRunoff%dd_z(iunit,j)		= TRunoff%dd_z(iunit,i)
				TRunoff%d_v(iunit,j)		= TRunoff%d_v(iunit,i)
				TRunoff%d_z(iunit,j+1)		= TRunoff%d_z(iunit,i+2)
				TRunoff%v_zt(iunit,j+1)		= TRunoff%v_zt(iunit,i+2)
				TRunoff%a_d(iunit,j+1)		= TRunoff%a_d(iunit,i+2)		
			elseif (j<TRunoff%d_ns(iunit) .and. i<(TRunoff%d_ns(iunit)-1)) then
				dv_ou(j)	= dv_ou(j+1)
				dv_in(j)	= dv_in(j+1)
				enr_0(j)	= enr_0(j+1)
				TRunoff%temp_lake(iunit,j)	= TRunoff%temp_lake(iunit,j+1)
				TRunoff%dd_z(iunit,j)		= TRunoff%dd_z(iunit,j+1)
				TRunoff%d_v(iunit,j)		= TRunoff%d_v(iunit,j+1)
				TRunoff%d_z(iunit,j+1)		= TRunoff%d_z(iunit,j+2)
				TRunoff%v_zt(iunit,j+1)		= TRunoff%v_zt(iunit,j+2)
				TRunoff%a_d(iunit,j+1)		= TRunoff%a_d(iunit,j+2)									
			elseif (j==i .and. i==(TRunoff%d_ns(iunit)-1)) then
				dv_ou(j-1)	= dv_ou(i)
				dv_in(j-1)	= dv_in(i)
				enr_0(j-1)	= enr_0(i)
				TRunoff%temp_lake(iunit,j-1)	= TRunoff%temp_lake(iunit,i)
				TRunoff%dd_z(iunit,j-1)			= TRunoff%dd_z(iunit,i)
				TRunoff%d_v(iunit,j-1)			= TRunoff%d_v(iunit,i)
				TRunoff%d_z(iunit,j)			= TRunoff%d_z(iunit,j+1)
				TRunoff%v_zt(iunit,j)			= TRunoff%v_zt(iunit,j+1)
				TRunoff%a_d(iunit,j)			= TRunoff%a_d(iunit,j+1)				
				TRunoff%temp_lake(iunit,TRunoff%d_ns(iunit)-1)=TRunoff%temp_lake(iunit,TRunoff%d_ns(iunit))
				dv_ou(TRunoff%d_ns(iunit)-1)	= dv_ou(TRunoff%d_ns(iunit))
				dv_in(TRunoff%d_ns(iunit)-1)	= dv_in(TRunoff%d_ns(iunit))
				TRunoff%dd_z(iunit,TRunoff%d_ns(iunit)-1)	= TRunoff%dd_z(iunit,TRunoff%d_ns(iunit))
				TRunoff%d_v(iunit,TRunoff%d_ns(iunit)-1)	= TRunoff%d_v(iunit,TRunoff%d_ns(iunit))
				enr_0(TRunoff%d_ns(iunit)-1)			= enr_0(TRunoff%d_ns(iunit))
				TRunoff%d_z(iunit,TRunoff%d_ns(iunit))	= TRunoff%d_z(iunit,TRunoff%d_ns(iunit)+1)
				TRunoff%v_zt(iunit,TRunoff%d_ns(iunit))	= TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1)
				TRunoff%a_d(iunit,TRunoff%d_ns(iunit))	= TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)				
			elseif (j>=TRunoff%d_ns(iunit)) then
				dv_ou(j)	= 0._r8
				dv_in(j)	= 0._r8
				enr_0(j)	= 0._r8
				TRunoff%temp_lake(iunit,j)	= 0._r8
				TRunoff%dd_z(iunit,j)		= 0._r8
				TRunoff%d_v(iunit,j)		= 0._r8
				TRunoff%d_z(iunit,j+1)		= 0._r8
				TRunoff%v_zt(iunit,j+1)		= 0._r8
				TRunoff%a_d(iunit,j+1)		= 0._r8				
			end if
		end do		
		TRunoff%d_ns(iunit)	= TRunoff%d_ns(iunit)-1		

	end subroutine layer_merge
	
	subroutine layer_split(iunit,i,dv_in,dv_ou,enr_0)
!*******************************************************************************************************
! 	Split layer if the thickness is more than threshold ddz_max
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: iunit,i             
		real(r8),dimension(nlayers), intent(inout) :: dv_in,dv_ou,enr_0    ! layer inflow,outflow(m3/s), and energy(w) 
		real(r8) :: tab,e_ab,dd_zab,d_vab,dv_ouab,dv_inab,delta_z    
		integer :: j,k,m							! indices
				
		! Calculate layer geometric properties to be halved
		dd_zab 	= 0.5_r8*TRunoff%dd_z(iunit,i)
		e_ab   	= 0.5_r8*enr_0(i)
		d_vab	= TRunoff%d_v(iunit,i)
		tab		= TRunoff%temp_lake(iunit,i)
		dv_ouab	= dv_ou(i)
		dv_inab	= dv_in(i)
		!	Re-number layers before dividing layer 
		TRunoff%d_z(iunit,TRunoff%d_ns(iunit)+2)	= TRunoff%d_z(iunit,TRunoff%d_ns(iunit)+1)
		TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+2)	= TRunoff%v_zt(iunit,TRunoff%d_ns(iunit)+1)
		TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+2)	= TRunoff%a_d(iunit,TRunoff%d_ns(iunit)+1)
		m=i+1
		do j=m,TRunoff%d_ns(iunit)
			k=TRunoff%d_ns(iunit)-j+i+1
			TRunoff%temp_lake(iunit,k+1)	= TRunoff%temp_lake(iunit,k)
			dv_ou(k+1)	= dv_ou(k)
			dv_in(k+1)	= dv_in(k)
			TRunoff%d_z(iunit,k+1)	= TRunoff%d_z(iunit,k)
			TRunoff%a_d(iunit,k+1)	= TRunoff%a_d(iunit,k)
			TRunoff%v_zt(iunit,k+1)	= TRunoff%v_zt(iunit,k)
			TRunoff%d_v(iunit,k+1)	= TRunoff%d_v(iunit,k)
			TRunoff%dd_z(iunit,k+1)	= TRunoff%dd_z(iunit,k)
			enr_0(k+1)= enr_0(k)
		end do
	!	Divide layer in half and calculate corresponding properties
		TRunoff%dd_z(iunit,i)	= dd_zab
		TRunoff%dd_z(iunit,i+1)	= dd_zab
		TRunoff%d_z(iunit,i+1)	= TRunoff%d_z(iunit,i)+TRunoff%dd_z(iunit,i)
		do j=2,ngeom+1
			if (TRunoff%d_z(iunit,i+1)>TRunoff%d_zi(iunit,j-1).and.TRunoff%d_z(iunit,i+1)<=TRunoff%d_zi(iunit,j))then
				delta_z   = (TRunoff%d_z(iunit,i+1)-TRunoff%d_zi(iunit,j-1))/(TRunoff%d_zi(iunit,j)-TRunoff%d_zi(iunit,j-1))
				TRunoff%a_d(iunit,i+1)  =(delta_z*(TRunoff%a_di(iunit,j)-TRunoff%a_di(iunit,j-1))+TRunoff%a_di(iunit,j-1))
				TRunoff%v_zt(iunit,i+1) = (delta_z*(TRunoff%v_zti(iunit,j)-TRunoff%v_zti(iunit,j-1))+TRunoff%v_zti(iunit,j-1))
			elseif (TRunoff%d_z(iunit,i+1)> TRunoff%d_zi(iunit,ngeom+1))then
				delta_z   = (TRunoff%d_z(iunit,i+1)-TRunoff%d_zi(iunit,ngeom))/(TRunoff%d_zi(iunit,ngeom+1)-TRunoff%d_zi(iunit,ngeom))
				TRunoff%a_d(iunit,i+1)  = (delta_z*(TRunoff%a_di(iunit,ngeom+1)-TRunoff%a_di(iunit,ngeom))+TRunoff%a_di(iunit,ngeom))
				TRunoff%v_zt(iunit,i+1) = (delta_z*(TRunoff%v_zti(iunit,ngeom+1)-TRunoff%v_zti(iunit,ngeom))+TRunoff%v_zti(iunit,ngeom))
			end if
		end do
		TRunoff%d_v(iunit,i+1)	= d_vab-(TRunoff%v_zt(iunit,i+1)-TRunoff%v_zt(iunit,i))
		TRunoff%d_v(iunit,i)	= (TRunoff%v_zt(iunit,i+1)-TRunoff%v_zt(iunit,i))
		dv_ou(i+1)	= TRunoff%d_v(iunit,i+1)*dv_ouab/(TRunoff%d_v(iunit,i)+TRunoff%d_v(iunit,i+1))
		dv_ou(i)	= TRunoff%d_v(iunit,i)*dv_ouab/(TRunoff%d_v(iunit,i)+TRunoff%d_v(iunit,i+1))
		TRunoff%temp_lake(iunit,i)		= tab
		TRunoff%temp_lake(iunit,i+1)	= tab
		enr_0(i)	= e_ab
		enr_0(i+1)	= e_ab 
		dv_in(i+1)	= TRunoff%d_v(iunit,i+1)*dv_inab/(TRunoff%d_v(iunit,i)+TRunoff%d_v(iunit,i+1))
		dv_in(i)	= TRunoff%d_v(iunit,i)*dv_inab/(TRunoff%d_v(iunit,i)+TRunoff%d_v(iunit,i+1))
		TRunoff%d_ns(iunit)	= TRunoff%d_ns(iunit)+1	
		
	end subroutine layer_split

	subroutine over_under_flow(iunit,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
!*******************************************************************************************************
! 	Split layer if the thickness is more than threshold ddz_max
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: iunit,d_n_n              
		real(r8),intent(in) :: d_res_sub              
		real(r8),dimension(nlayers), intent(in) :: t_z_old 
		real(r8),dimension(nlayers), intent(out) :: dv_in,dv_ou    ! layer inflow,outflow(m3/s)
		integer :: j							! indices
		
		TRunoff%d_ns(iunit)  	= d_n_n				
		TRunoff%d_lake(iunit)	= d_res_sub 
		do j = 1, nlayers+1	
			if (j<=nlayers)  then   
				TRunoff%d_v(iunit,j) 		= TRunoff%v_zo(iunit,j)
				TRunoff%temp_lake(iunit,j) = t_z_old(j)
				TRunoff%v_zt(iunit,j) 	= TRunoff%v_zt0(iunit,j)
				TRunoff%a_d(iunit,j) 	= TRunoff%a_d0(iunit,j)
				TRunoff%d_z(iunit,j) 	= TRunoff%d_z0(iunit,j)
				TRunoff%dd_z(iunit,j)	= TRunoff%d_z0(iunit,j+1)-TRunoff%d_z0(iunit,j)
				dv_in(j)	= 0._r8
				dv_ou(j)	= 0._r8
			else
				TRunoff%v_zt(iunit,j) 	= TRunoff%v_zt0(iunit,j)
				TRunoff%a_d(iunit,j) 	= TRunoff%a_d0(iunit,j)
				TRunoff%d_z(iunit,j) 	= TRunoff%d_z0(iunit,j)
			end if
		end do

	end subroutine over_under_flow
	
	function den(t_z) result(rho)
	! calculate density from temperature 
		implicit none
		real(r8), intent(in) :: t_z! Temperature (k) 
		real(r8) :: rho! ! density (kg/m3)
			
		rho = 1000._r8*( 1._r8 - 1.9549e-05*(abs(t_z-277._r8))**1.68_r8) ! modified from Subin et al, 2011 with lake ice fraction = 0
			
    end function den 
	
    subroutine solve(a,b,c,r,iunit,n) !
		use shr_sys_mod , only : shr_sys_flush
	
		implicit none
!	 a:left coefficient , b:middle coefficient, c:right coefficient, r:right hand side (known), n - number of layers
		integer,intent(in) :: n,iunit
        real(r8),dimension(n),intent(in) :: a,b,c,r
        real(r8),dimension(n) :: bp,rp
        real(r8) :: m,tt 
        integer i		
		
! 	initialize c-prime and d-prime
		bp(1)=b(1)
		rp(1)=r(1)
		do i=2,n 
			tt=a(i)/b(i-1)
			bp(i) = b(i)-c(i-1)*tt
			rp(i) = r(i)-rp(i-1)*tt
		end do
		bp(n)=b(n)		
! 	initialize u
        TRunoff%temp_lake(iunit,n) = max(rp(n)/bp(n),273.15_r8)	! to be modified after including ice/snow effect
! 	Back substitution
        do i = n-1, 1, -1
			TRunoff%temp_lake(iunit,i) = max((rp(i)-c(i)*TRunoff%temp_lake(iunit,i+1))/bp(i),273.15_r8)
        end do
		
	end subroutine solve 
	
	subroutine QSat_str (T, p, es)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
		use shr_kind_mod , only: r8 => shr_kind_r8
		use shr_const_mod, only: SHR_CONST_TKFRZ
	!
	! !ARGUMENTS:
		implicit none
		real(r8), intent(in)  :: T        ! temperature (K)
		real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
		real(r8), intent(out) :: es       ! vapor pressure (pa)
!! !LOCAL VARIABLES:
!EOP
!
		real(r8) :: T_limit
		real(r8) :: td,vp,vp1,vp2
	!
	! For water vapor (temperature range 0C-100C)
	!
		real(r8), parameter :: a0 =  6.11213476_r8
		real(r8), parameter :: a1 =  0.444007856_r8
		real(r8), parameter :: a2 =  0.143064234e-01_r8
		real(r8), parameter :: a3 =  0.264461437e-03_r8
		real(r8), parameter :: a4 =  0.305903558e-05_r8
		real(r8), parameter :: a5 =  0.196237241e-07_r8
		real(r8), parameter :: a6 =  0.892344772e-10_r8
		real(r8), parameter :: a7 = -0.373208410e-12_r8
		real(r8), parameter :: a8 =  0.209339997e-15_r8
	!
	! For ice (temperature range -75C-0C)
	!
		real(r8), parameter :: c0 =  6.11123516_r8
		real(r8), parameter :: c1 =  0.503109514_r8
		real(r8), parameter :: c2 =  0.188369801e-01_r8
		real(r8), parameter :: c3 =  0.420547422e-03_r8
		real(r8), parameter :: c4 =  0.614396778e-05_r8
		real(r8), parameter :: c5 =  0.602780717e-07_r8
		real(r8), parameter :: c6 =  0.387940929e-09_r8
		real(r8), parameter :: c7 =  0.149436277e-11_r8
		real(r8), parameter :: c8 =  0.262655803e-14_r8
	!-----------------------------------------------------------------------

		T_limit = T - SHR_CONST_TKFRZ
		if (T_limit > 100.0_r8) T_limit=100.0_r8
		if (T_limit < -75.0_r8) T_limit=-75.0_r8

		td       = T_limit
		if (td >= 0.0_r8) then
		   es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
				+ td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
		else
		   es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
				+ td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
		endif

		es    = es    * 100._r8            ! pa
	
	end subroutine QSat_str
  
end MODULE MOSART_lake_mod