!
MODULE MOSART_stra_mod
! Description: core code of MOSART-heat reservoir stratification. 
! 
! Developed by Wondie Yigzaw, July 2019. 
! REVISION HISTORY:
! 
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TPara, rtmCTL
    use rof_cpl_indices, only : nt_nliq, nt_nice
    use RtmVar         , only : iulog, ngeom, nlayers
    use WRM_type_mod  , only :  ctlSubwWRM, WRMUnit, StorWater
    use MOSART_heat_mod
    use RtmTimeManager
    implicit none
	
    public stratification
	
! !PUBLIC MEMBER FUNCTIONS:
    contains  
		
	subroutine stratification(iunit, localDeltaT,nt)
	! !DESCRIPTION: calculate the water temperature of reservoir.
	
		use shr_sys_mod , only : shr_sys_flush
	
		implicit none
		integer,  intent(in) :: iunit, nt
		real(r8), intent(in) :: localDeltaT
		integer :: damID              
		character(len=*),parameter :: subname = '(stratification)'
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
		damID = WRMUnit%INVicell(iunit) 
			
		if (damID <= ctlSubwWRM%LocalNumDam .and. damID >0 .and. WRMUnit%MeanMthFlow(damID,13) > 0.01_r8) then	! Reservoir stratification only if regulation
			if(WRMUnit%d_resrv(damID) <= 10._r8 .or. WRMUnit%d_ns(damID) < 1 .or. WRMUnit%Depth(damID) <= 0.0_r8 .or. WRMUnit%Height(damID) <= 0.0_r8 .or. WRMUnit%geometry(damID) < 1.0_r8) then 
				THeat%Tr(iunit) = THeat%Tr(iunit)
				WRMUnit%resrv_surf(iunit)= THeat%Tr(iunit)
				return
			end if
			
		!	Calculate sub-time step for numerical stability			
			s_dtime = 12  			
			dtime   = localDeltaT/s_dtime		!	Sub-timestep for numerical purpose
			if(THeat%Tr(iunit)<273.15_r8)THeat%Tr(iunit)=273.15_r8
			t_in    = THeat%Tr(iunit) 
			q_in    = TRunoff%erin(iunit,nt)
			q_ou    = -TRunoff%erout(iunit,nt)
			
			!	Initialize for sub-timestep	
			! !***************************************************************************************************************		
			do ww = 1,s_dtime	!	Start calculation for each sub-timestep	************************************************			
				!	Initialize for sub-timestep
				d_n_n  		= WRMUnit%d_ns(damID)					! carry number of layers if there is merging/split
				d_res_sub 	= WRMUnit%d_resrv(damID)
				do j = 1, nlayers+1	
					WRMUnit%v_zt0(damID,j) = WRMUnit%v_zt(damID,j)
					WRMUnit%a_d0(damID,j) = WRMUnit%a_d(damID,j)
					WRMUnit%d_z0(damID,j) = WRMUnit%d_z(damID,j)
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
					if (j<=WRMUnit%d_ns(damID)) then
						WRMUnit%temp_resrv(damID,j) = WRMUnit%temp_resrv(damID,j)
						rho_z(j) = den(WRMUnit%temp_resrv(damID,j))
						t_z_old(j) = WRMUnit%temp_resrv(damID,j)
					else 
						WRMUnit%temp_resrv(damID,j) = 0._r8
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
				t_s    	= WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)) 
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
				if(WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1) > 0.10_r8*(WRMUnit%v_zti(damID,ngeom+1)).and. q_ou < q_in) then 
					v_evap 	= max(d_evap*sar*WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1),0._r8)
				else
					v_evap	= 0._r8
				end if				
			
		! 	Calculation of flow contibution due to inflow/outflow 				
				! 	Calculate total mass (kg), layer energy (w)
				do j = 1, nlayers	
					if (j<=WRMUnit%d_ns(damID))  then
						WRMUnit%v_zo(damID,j) = WRMUnit%d_v(damID,j)
						enr_0(j) = WRMUnit%temp_resrv(damID,j)*WRMUnit%d_v(damID,j)*rho_z(j)*c_w/dtime
					else
						WRMUnit%v_zo(damID,j) = 0._r8		
						enr_0(j) = 0._r8
					end if
				end do
				
				if (abs(q_in - q_ou)*dtime <= 0.5_r8*WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1) .and. WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1) > 0.6_r8*WRMUnit%v_zti(damID,ngeom+1) .and. WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1) < WRMUnit%v_zti(damID,ngeom+1)) then  ! Redistribute mass flux only if there is net flow advection and enough storage
				
					call flowdist(WRMUnit%d_ns(damID),q_in,t_in,q_ou,dtime,damID,dv_in,dv_ou) 
					
					!	Resize layer thickness and numbers based on inflow/outflow contribution							
					do j = 1, nlayers
						if (j<=WRMUnit%d_ns(damID)-1 .and. WRMUnit%d_ns(damID)>1) then
							dv_nt(j)	= (dv_in(j)-dv_ou(j)-v_evap/(WRMUnit%d_ns(damID)-1))
						elseif(j==WRMUnit%d_ns(damID) .and. WRMUnit%d_ns(damID)>1) then
							dv_nt(j)	= (dv_in(j)-dv_ou(j))
						elseif(j==WRMUnit%d_ns(damID) .and. WRMUnit%d_ns(damID)==1) then
							dv_nt(j)	= (dv_in(j)-dv_ou(j)-v_evap)
						else
							dv_nt(j)	= 0._r8
						end if
					end do
					
					! 	Calculate layer volume (m3)
					do j = 1, nlayers	
						if (j<=WRMUnit%d_ns(damID))  then   
							WRMUnit%v_zn(damID,j) = ((WRMUnit%d_v(damID,j)+dv_nt(j)))
						else
							WRMUnit%v_zn(damID,j) = 0._r8
						end if
					end do
								
				! Calculate new layer thickness, depth, area/volume after net advection		
					call geometry(damID,dv_nt)					
					
				! Merge/split layers					
					if (WRMUnit%d_resrv(damID) <= 2.0_r8 .or. WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1)<=WRMUnit%v_zti(damID,1) .or. WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1)>WRMUnit%v_zti(damID,ngeom+1)) then ! skip if storage is too small/large
						call over_under_flow(damID,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)						
					else
						ddz_min = 1.5_r8
						ddz_max = max(3.0_r8*WRMUnit%ddz_local(damID),5.0_r8)
						!	Check if layers are too small   
						if (WRMUnit%d_ns(damID) >=2) then
							do i=1,WRMUnit%d_ns(damID)-1
								if (WRMUnit%dd_z(damID,i) < ddz_min .and. WRMUnit%d_ns(damID)>1) then
									call layer_merge(damID,i,dv_in,dv_ou,enr_0)
									if (WRMUnit%d_ns(damID)<1) then
										exit
									end if
								elseif (WRMUnit%d_ns(damID)==1) then
									call over_under_flow(damID,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
									exit
								end if								
							end do
						end if						
						! 	Check if layers are too big
						if (WRMUnit%d_ns(damID) <= nlayers-1) then !
							do i=1,WRMUnit%d_ns(damID)
								if(WRMUnit%dd_z(damID,i) > ddz_max .and. WRMUnit%d_ns(damID)<nlayers) then 
									call layer_split(damID,i,dv_in,dv_ou,enr_0)
								elseif (WRMUnit%dd_z(damID,i) > ddz_max .and. i<WRMUnit%d_ns(damID) .and. WRMUnit%d_ns(damID)==nlayers) then  ! allow further merging by increasing ddz_min
									ddz_min = max(0.5*WRMUnit%ddz_local(damID),3.0_r8)
									call layer_merge(damID,i,dv_in,dv_ou,enr_0)
									if (WRMUnit%d_ns(damID)<1) then
										exit
									end if
									call layer_split(damID,i,dv_in,dv_ou,enr_0)
									if (WRMUnit%d_ns(damID)==nlayers) then
										exit
									end if
								elseif(i==WRMUnit%d_ns(damID) .and. WRMUnit%d_ns(damID)==nlayers) then
									call over_under_flow(damID,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
									exit
								end if								
							end do
						end if					
					end if		!layer merge/split
					
					!	Recalculare reservoir depth
					WRMUnit%d_resrv(damID)=0._r8 
					do i=1,WRMUnit%d_ns(damID)
					   WRMUnit%d_resrv(damID) = WRMUnit%d_resrv(damID)+WRMUnit%dd_z(damID,i)
					end do
					
					! 	Recalculate density after layer change
					do j = 1, nlayers	
						if (j<=WRMUnit%d_ns(damID))  then
							rho_z(j) = den(WRMUnit%temp_resrv(damID,j))
						else
							rho_z(j) = 0._r8
						end if
					end do
					
					! 	Calculate layer internal energy (w) due to inflow/outflow		
					do j = 1, nlayers	
						if (j<=WRMUnit%d_ns(damID))  then   
							enr_in(j) = (dv_in(j)*t_in*rho_r*c_w/dtime)		!Energy from inflow
							enr_ou(j) = (dv_ou(j)*WRMUnit%temp_resrv(damID,j)*rho_z(j)*c_w/dtime)	!Energy loss due to outflow
						else
							enr_in(j) = 0._r8
							enr_ou(j) = 0._r8
						end if
					end do	
					
					if(WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1) <= 0.050_r8*WRMUnit%v_zti(damID,ngeom+1) .or. WRMUnit%d_v(damID,1) <= 0._r8 .or. (WRMUnit%d_ns(damID) ==1 .and. WRMUnit%d_resrv(damID)<=2._r8)) then	! Reservoir virtually empty/shallow
						! write(iulog,*) 'Reservoir empty/shallow',WRMUnit%grandid(damID)
						do j = 1, nlayers	
							if (j<=WRMUnit%d_ns(damID)) then
								WRMUnit%temp_resrv(damID,j) = THeat%Tr(iunit)
							else 
								WRMUnit%temp_resrv(damID,j) = 0._r8
							end if
						end do
					end if 				
					
				end if ! layer recalculation
				
			! 	Calculate layer energy (w)
				do j = 1, nlayers	
					if (j<=WRMUnit%d_ns(damID))  then   
						enr_1(j) = enr_0(j)+ enr_in(j) - enr_ou(j)
					else
						enr_1(j) = 0._r8
					end if
				end do	
				
				! 	Adjust layer temperature for mixing due to inflow/outflow
				do j = 1, nlayers	
					if (j<=WRMUnit%d_ns(damID))  then   
						WRMUnit%temp_resrv(damID,j) = enr_1(j)*dtime/(WRMUnit%d_v(damID,j)*rho_z(j)*c_w)
					else
						WRMUnit%temp_resrv(damID,j) = 0._r8
					end if
				end do
				
			!	Check energy balance (w) after advective mixing	
				enr_err1 = (sum(enr_1) - (sum(enr_0)+ sum(enr_in) - sum(enr_ou)))
				
			!******************************************************************************
			! 	Calculate solar energy absorbed at each layer		
				eta  = 1.1925_r8*(WRMUnit%d_resrv(damID))**(-0.424_r8) ! as used in Subin et al, 2011 (citing Hakanson, 1995) but modified for actual reservoir depth 
				beta = 0.175_r8 !
		
				do j=1,nlayers+1
					phi_x(j) = 0._r8
				end do
				
				if (sh_net > 0._r8 .and. WRMUnit%d_ns(damID) >1) then
					k=0
					top_d=WRMUnit%d_resrv(damID)-WRMUnit%d_z(damID,WRMUnit%d_ns(damID)-k)
					if(top_d < 0.61_r8) then
						k=k+1
					else
						k=k+2
					end if	
					v_mix=(WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1)-WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)-k))	
					!	Solar radiation energy absorbed at the mixed zone
					sh_mix=sh_net*sar*(WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)-(1._r8-beta)*WRMUnit%a_d(damID,WRMUnit%d_ns(damID)-k))
					j=WRMUnit%d_ns(damID)-k
					do i=j,WRMUnit%d_ns(damID)
					   phi_z(i)=sh_mix*WRMUnit%d_v(damID,i)/v_mix
					end do
					phi_x(WRMUnit%d_ns(damID)+1)=sh_net
					phi_x(WRMUnit%d_ns(damID)-k)=(1._r8-beta)*sh_net
					if(k>0) then
						do i=1,k
						   ii=WRMUnit%d_ns(damID)-i+1
						   phi_x(ii)=(WRMUnit%a_d(damID,ii+1)*phi_x(ii+1)-phi_z(ii))/(WRMUnit%a_d(damID,ii))
						end do
					end if
					
					!	Solar radiation energy absorbed at sub layers
					l=WRMUnit%d_ns(damID)-k-1
					do j=1,l
					   i=(WRMUnit%d_ns(damID)-k)-j+1
					   phi_x(i-1)=phi_x(i)*exp(-eta*(WRMUnit%d_z(damID,i)-WRMUnit%d_z(damID,i-1)))
					end do
				
					!	Solar radiation energy absorbed in each layer
					j=WRMUnit%d_ns(damID)-k-1
					do i=1,j
					   phi_z(i)=(sar*WRMUnit%a_d(damID,i+1)*phi_x(i+1)-sar*WRMUnit%a_d(damID,i)*phi_x(i))
					   
					end do
				elseif (sh_net > 0._r8 .and. WRMUnit%d_ns(damID)==1) then
					phi_z(WRMUnit%d_ns(damID))=sh_net*sar*WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)
				else
					do j=1,WRMUnit%d_ns(damID)
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
				k_ew=tau*s_vel*sar*WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)*dtime	! Wind driven kinetic energy at surface
				dis_w = k_ew/(rho_w*WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1)*dtime)  ! rate of dissipation-wind
				cfw = 1.e-02!
				cfa = 1.e-05!
				k_m = 0.57_r8/(c_w*rho_w) 						  !molecular diffusivity
				df_eff(1)= 0._r8			!bottom interface
				df_eff(WRMUnit%d_ns(damID)+1)= 0._r8		!air interface
				do j = 2,WRMUnit%d_ns(damID) 
					q_adv(j) = max((dv_in(j)+dv_ou(j)),0._r8)
					k_ad(j)=0.5_r8*rho_w*q_adv(j)*dtime*(q_adv(j)/(WRMUnit%Width_r(damID)*WRMUnit%dd_z(damID,j)))**2._r8 ! Advection driven kinetic energy 
					dis_ad(j)= k_ad(j)/(rho_w*WRMUnit%v_zt(damID,j)*dtime)	! rate of dissipation-inflow/outflow					
				! Calculate Richardson number
					drhodz(j) = (rho_z(j-1)-rho_z(j))/0.5_r8*(WRMUnit%dd_z(damID,j)+WRMUnit%dd_z(damID,j-1))
					bv_f = max((grav/rho_w)*drhodz(j),0._r8)
					if (s_vel <= 0._r8) ri = 0._r8
					ri = bv_f/((s_vel/(0.4_r8*WRMUnit%d_z(damID,j)))**2._r8)				
				! Calculate Froude number
					l_vel = q_adv(j)*WRMUnit%Length_r(damID)/(sar*WRMUnit%a_d(damID,j)*WRMUnit%dd_z(damID,j))
					if (q_adv(j) <= 0._r8 .or. drhodz(j) <= 0._r8) Fr(j) = 0._r8
					Fr(j)= (grav*WRMUnit%dd_z(damID,j)*drhodz(j)/rho_w)/l_vel**2._r8					
				! Calculate diffusion coefficients				
					df_eff(j)=min(max(dtime**2._r8*((cfw*dis_w/(1+ri))+(0.5_r8*cfa*(dis_ad(j)+dis_ad(j-1))/(1+Fr(j)))),k_m),5.56e-03) !5.56e-03!
					
				end do
			
			!*****************************************************************
			! Calculate matrix elements
				do j = 1,WRMUnit%d_ns(damID)
					if (j == 1 .and. WRMUnit%d_ns(damID)>1) then
						m1(j) = 2_r8*dtime/(0.5_r8*(sar*WRMUnit%a_d(damID,j)+sar*WRMUnit%a_d(damID,j+1))*WRMUnit%dd_z(damID,j))
						m2(j) = m1(j)*sar*WRMUnit%a_d(damID,j+1)*df_eff(j+1)/(WRMUnit%dd_z(damID,j)+WRMUnit%dd_z(damID,j+1))
						m3(j) = 0._r8						
						Fx(j) = dtime*phi_z(j)/(WRMUnit%d_v(damID,j)*c_w*rho_z(j))
						a(j) = - (m2(j))
						b(j) = 1._r8 + (m2(j) + m3(j)) 
						c(j) = 0._r8 
						r(j) = WRMUnit%temp_resrv(damID,j) + Fx(j) ! bottom boundary condition 
					elseif (j <= WRMUnit%d_ns(damID)-1 .and. WRMUnit%d_ns(damID)>2) then
						m1(j) = 2_r8*dtime/(0.5_r8*(sar*WRMUnit%a_d(damID,j)+sar*WRMUnit%a_d(damID,j+1))*WRMUnit%dd_z(damID,j))
						m2(j) = m1(j)*sar*WRMUnit%a_d(damID,j+1)*df_eff(j+1)/(WRMUnit%dd_z(damID,j)+WRMUnit%dd_z(damID,j+1))
						m3(j) = m1(j)*sar*WRMUnit%a_d(damID,j)*df_eff(j)/(WRMUnit%dd_z(damID,j)+WRMUnit%dd_z(damID,j-1))						
						Fx(j) = dtime*phi_z(j)/(WRMUnit%d_v(damID,j)*c_w*rho_z(j))
						a(j) = - m2(j)
						b(j) = 1._r8 + m2(j) + m3(j) 
						c(j) = - m3(j)
						r(j) = WRMUnit%temp_resrv(damID,j) + Fx(j)
					elseif (j == WRMUnit%d_ns(damID)) then!top layer
						m1(j) = 2_r8*dtime/(0.5_r8*(sar*WRMUnit%a_d(damID,j)+sar*WRMUnit%a_d(damID,j+1))*WRMUnit%dd_z(damID,j))
						m2(j) = 0._r8
						m3(j) = m1(j)*sar*WRMUnit%a_d(damID,j)*df_eff(j)/(WRMUnit%dd_z(damID,j)+WRMUnit%dd_z(damID,j-1))
						Fx(j) = dtime*((phi_o-sh_net)*sar*WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)+phi_z(j))/(WRMUnit%d_v(damID,WRMUnit%d_ns(damID))*c_w*rho_z(j)) ! 
						a(j) = 0._r8
						b(j) = 1._r8 + (m2(j) + m3(j)) 
						c(j) = - (m3(j))
						r(j) = WRMUnit%temp_resrv(damID,j) + Fx(j) ! top boundary condition 						
					end if
				end do	
				
				
			!	Solve for temperature
				call solve(a,b,c,r,damID,WRMUnit%d_ns(damID))
				
            ! 	Check numerical instability from reservoir geometry data 	
				do i = 1,WRMUnit%d_ns(damID)
					if (isnan(WRMUnit%temp_resrv(damID,i)) .or. WRMUnit%temp_resrv(damID,i)>=huge(1._r8)) then
						! WRMUnit%temp_resrv(damID,i) = THeat%Tr(iunit)
						write(iulog,*) subname,'Check geometry',WRMUnit%grandid(damID) 					
						return
					end if	
				end do
			
				! Calculate layer intermediate internal energy (w)
				do j = 1, nlayers	
					if (j<WRMUnit%d_ns(damID))  then   
						enr_2(j)   = WRMUnit%temp_resrv(damID,j)*WRMUnit%d_v(damID,j)*rho_z(j)*c_w/dtime
						enr_phi(j) = phi_z(j)
					elseif (j==WRMUnit%d_ns(damID))  then
						enr_2(j)   = WRMUnit%temp_resrv(damID,j)*WRMUnit%d_v(damID,j)*rho_z(j)*c_w/dtime
						enr_phi(j) = (phi_o-sh_net)*sar*WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)		!Enegry from net surface flux (w)
					else
						enr_2(j) = 0.
					end if
				end do				
			
		!	Check energy balance error (w) after triadiagonal matrix solution	
				enr_err2 = (sum(enr_2) - (sum(enr_1) + sum(enr_phi)))
				
		!***********************************************************************************************************************		
		! !	Solve convective mixing
			! Recalculate layer density 
				do j = 1,WRMUnit%d_ns(damID)   
					rho_z(j) = den(WRMUnit%temp_resrv(damID,j))
				end do
				
			! Check if instability exists
				do j = 1,WRMUnit%d_ns(damID)-1 
					if(rho_z(j) < rho_z(j+1)) then
						! Start mixing layers
						mix1=j
						mix2=mix1
						sumvol=WRMUnit%d_v(damID,mix2)*1._r8
						tsum=WRMUnit%temp_resrv(damID,mix2)*WRMUnit%d_v(damID,mix2)	
						
						mix2=mix2+1
						mixvol2=WRMUnit%d_v(damID,mix2)*1._r8
						sumvol=sumvol + mixvol2
						tsum=tsum+WRMUnit%temp_resrv(damID,mix2)*mixvol2
						tmix=tsum/sumvol
						
						! Calculate density of mixed layer	
						denmix = den(tmix)
						
						! Check if instability exists below mixed layer	and mix layers	
						if(rho_z(mix1-1) < denmix .and. mix1 >= 2) then
							mix1=mix1-1	
						
							! Calculate temperature of mixed layer	
							mixvol1=WRMUnit%d_v(damID,mix1)*1._r8
							sumvol=sumvol + mixvol1
							tsum=tsum+WRMUnit%temp_resrv(damID,mix1)*mixvol1
							tmix= tsum/sumvol
						end if
						
						! Calculate density of mixed layer	
						denmix = den(tmix)					
				  
						! Set new layer temperature and density 
						do i=mix1,mix2
							rho_z(i)=denmix
							WRMUnit%temp_resrv(damID,i)=tmix
						end do
					end if
				end do
		
			! 	Sum sub-timestep variables			
				d_res_sub = (d_res_sub + WRMUnit%d_resrv(damID))/2.0_r8
				if (d_n_n == WRMUnit%d_ns(damID)) then
					do j = 1,nlayers 
						WRMUnit%temp_resrv(damID,j) = (t_z_old(j) + WRMUnit%temp_resrv(damID,j))/2.0_r8
					end do 
				else 
					do j = 1,nlayers
						if (j<=WRMUnit%d_ns(damID))	then		
							WRMUnit%temp_resrv(damID,j) = (t_z_old(j) + WRMUnit%temp_resrv(damID,j))/2.0_r8
						else
							WRMUnit%temp_resrv(damID,j) = 0._r8
						end if
					end do 
				end if
				
				WRMUnit%d_resrv(damID) = d_res_sub
				
			end do ! sub-timestep
						
			!********* determine outlet location based on reservoir purpose
			!********* 1=Water supply;2=Recreation;3=Other;4=Navigation;5=Irrigation;6=Hydroelectricity;7=Flood control;8=Fisheries; 9=Large Reservoir (depth > 95m); 10=Modified outlet
			if (WRMUnit%purpose(damID)== 1 .or. WRMUnit%purpose(damID)== 8) then
				WRMUnit%out_lc(damID) = 0.35_r8
			elseif (WRMUnit%purpose(damID)== 2 .or. WRMUnit%purpose(damID)== 3 .or. WRMUnit%purpose(damID)== 4 .or. WRMUnit%purpose(damID)== 5) then
				WRMUnit%out_lc(damID) = 0.15_r8
			elseif (WRMUnit%purpose(damID)== 6 .or. WRMUnit%purpose(damID)== 7) then
				WRMUnit%out_lc(damID) = 0.65_r8
			elseif (WRMUnit%purpose(damID)== 9) then 
				WRMUnit%out_lc(damID) = 0.25_r8
			elseif (WRMUnit%purpose(damID)== 10) then 
				WRMUnit%out_lc(damID) = WRMUnit%out_lc(damID)
			end if
			
			
			!********* outflow taken from layer at ~60% of the reservoir depth			
			! WRMUnit%out_lc(damID) = 0.3_r8
			if (WRMUnit%out_lc(damID) >= 0._r8) then 
				outlet = max(int(WRMUnit%d_ns(damID) - WRMUnit%out_lc(damID)*WRMUnit%d_ns(damID)),1)
				t_out = WRMUnit%temp_resrv(damID,outlet)
				if (t_out < 273.15_r8) t_out = 273.15_r8 
			else 
				if (WRMUnit%d_ns(damID)<= int(0.1_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID))
				elseif (WRMUnit%d_ns(damID)> int(0.1_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.2_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.075_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.2_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.3_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.15_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.3_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.4_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.225_r8*nlayers))	
				elseif (WRMUnit%d_ns(damID)> int(0.4_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.5_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.3_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.5_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.6_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.375_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.6_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.7_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.45_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.7_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.8_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.525_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.8_r8*nlayers) .and. WRMUnit%d_ns(damID)<= int(0.9_r8*nlayers)) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.6_r8*nlayers))
				elseif (WRMUnit%d_ns(damID)> int(0.9_r8*nlayers) .and. WRMUnit%d_ns(damID)<= nlayers) then
					t_out = WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-int(0.675_r8*nlayers))
				end if
				if (t_out < 273.15_r8) t_out = 273.15_r8 
			end if 	
			
			WRMUnit%resrv_surf(iunit)= WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID))
			THeat%Tr(iunit)= t_out
					
		end if ! stratification
    end subroutine stratification

	subroutine flowdist(d_n,in_f,in_t,ou_f,dtime,damID,dv_in,dv_ou)
!*******************************************************************************************************
! 	Calculation inflow/outflow contribution modified to distribute uniformly
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use WRM_type_mod  , only :  WRMUnit
		use RtmVar         , only : iulog, ngeom, nlayers
		implicit none
		integer, intent(in)  :: d_n,damID
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
				denom=max((WRMUnit%v_zt(damID,jmax+1)-WRMUnit%v_zt(damID,jmin)),WRMUnit%d_v(damID,j))
				dv_in(j) = in_v*(WRMUnit%d_v(damID,j)/(denom))
				dv_ou(j) = ou_v*(WRMUnit%d_v(damID,j)/(denom))
			end do
		end if
		
	end subroutine flowdist


	subroutine geometry(damID,dv_nt)
!*******************************************************************************************************
! 	Calculate new layer thickness, depth, area/volume after net advection
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use WRM_type_mod  , only :  WRMUnit
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: damID              
		real(r8),dimension(nlayers), intent(in) :: dv_nt    ! layer net inflow/outflow (m3/s)
		real(r8) :: num,dem,delta_z,delta_a    !
		integer :: i,j,k,mm							! indices
	
		WRMUnit%d_z(damID,1)=0._r8
		WRMUnit%a_d(damID,1)=WRMUnit%a_di(damID,1)
		WRMUnit%v_zt(damID,1)=WRMUnit%v_zti(damID,1)								
		do i=1,nlayers 			! check layers for available volume to satisfy net outflow
			if (i<=WRMUnit%d_ns(damID))  then
				if(-dv_nt(i) > WRMUnit%d_v(damID,i))then !current layer collapses, hence remaining volume taken from next upper/lower layer
					WRMUnit%v_zn(damID,i)=0._r8
					if (i<WRMUnit%d_ns(damID)-1 .and. WRMUnit%d_ns(damID)>1) then
						WRMUnit%v_zn(damID,i+1)=((WRMUnit%v_zn(damID,i+1) - (-dv_nt(i)-WRMUnit%d_v(damID,i))))			
						if (WRMUnit%v_zn(damID,i+1)<0._r8) WRMUnit%v_zn(damID,i+1) = 0._r8								
					elseif (i==WRMUnit%d_ns(damID)-1 .and. WRMUnit%d_ns(damID)>2) then	!avoid top layer collapse, hence remaining mass taken from next lower layer
						WRMUnit%v_zn(damID,i-1)=((WRMUnit%v_zn(damID,i-1) - (-dv_nt(i)-WRMUnit%d_v(damID,i))))
						if (WRMUnit%v_zn(damID,i-1)<0._r8) WRMUnit%v_zn(damID,i-1) = 0._r8								
						mm=i-1
						do k=mm,WRMUnit%d_ns(damID)	
							WRMUnit%v_zt(damID,k+1)=(WRMUnit%v_zt(damID,k)+(WRMUnit%v_zn(damID,k)))
							do j=2,ngeom+1
								if (WRMUnit%v_zt(damID,k+1)>WRMUnit%v_zti(damID,j-1).and.WRMUnit%v_zt(damID,k+1)<=WRMUnit%v_zti(damID,j))then
									num = (WRMUnit%v_zt(damID,k+1)-WRMUnit%v_zti(damID,j-1))
									dem = (WRMUnit%v_zti(damID,j)-WRMUnit%v_zti(damID,j-1))
									delta_z  = (WRMUnit%d_zi(damID,j)-WRMUnit%d_zi(damID,j-1))*num/dem
									delta_a  = (WRMUnit%a_di(damID,j)-WRMUnit%a_di(damID,j-1))*num/dem
									WRMUnit%d_z(damID,k+1) = WRMUnit%d_zi(damID,j-1) + delta_z 
									WRMUnit%a_d(damID,k+1) = ((WRMUnit%a_di(damID,j-1) + max((delta_a),0.0_r8)))
								elseif (WRMUnit%v_zt(damID,k+1)>WRMUnit%v_zti(damID,ngeom+1))then ! inflow causes maximum storage
									num = (WRMUnit%v_zt(damID,k+1)-WRMUnit%v_zti(damID,ngeom))
									dem = (WRMUnit%v_zti(damID,ngeom+1)-WRMUnit%v_zti(damID,ngeom))
									delta_z  = (WRMUnit%d_zi(damID,ngeom+1)-WRMUnit%d_zi(damID,ngeom))*num/dem
									delta_a  = (WRMUnit%a_di(damID,WRMUnit%d_ns(damID)+1) -WRMUnit%a_di(damID,WRMUnit%d_ns(damID)))*num/dem
									WRMUnit%d_z(damID,k+1) = WRMUnit%d_zi(damID,ngeom) + delta_z
									WRMUnit%a_d(damID,k+1) = ((WRMUnit%a_di(damID,ngeom) + max((delta_a),0.0_r8)))
								end if
							end do
							WRMUnit%d_z(damID,k)=WRMUnit%d_z(damID,k+1)
							WRMUnit%dd_z(damID,k)=WRMUnit%d_z(damID,k+1)-WRMUnit%d_z(damID,k)
						end do		
					elseif ((i==WRMUnit%d_ns(damID)-1 .or. i==WRMUnit%d_ns(damID)) .and. WRMUnit%d_ns(damID)<=2) then	! Top layer collapses, skip outflow							
						if (WRMUnit%v_zn(damID,i)<0._r8) WRMUnit%v_zn(damID,i) = WRMUnit%d_v(damID,i)								
					end if
					WRMUnit%v_zt(damID,i+1)=(WRMUnit%v_zt(damID,i)+(WRMUnit%v_zn(damID,i)))
					do j=2,ngeom+1
						if (WRMUnit%v_zt(damID,i+1)>WRMUnit%v_zti(damID,j-1).and.WRMUnit%v_zt(damID,i+1)<=WRMUnit%v_zti(damID,j))then
							num = (WRMUnit%v_zt(damID,i+1)-WRMUnit%v_zti(damID,j-1))
							dem = (WRMUnit%v_zti(damID,j)-WRMUnit%v_zti(damID,j-1))
							delta_z  = (WRMUnit%d_zi(damID,j)-WRMUnit%d_zi(damID,j-1))*num/dem
							delta_a  = (WRMUnit%a_di(damID,j)-WRMUnit%a_di(damID,j-1))*num/dem
							WRMUnit%d_z(damID,i+1) = WRMUnit%d_zi(damID,j-1) + delta_z 
							WRMUnit%a_d(damID,i+1) = ((WRMUnit%a_di(damID,j-1) + max((delta_a),0.0_r8)))
						elseif (WRMUnit%v_zt(damID,i+1)>WRMUnit%v_zti(damID,ngeom+1))then! inflow causes maximum storage
							num = (WRMUnit%v_zt(damID,i+1)-WRMUnit%v_zti(damID,ngeom))
							dem = (WRMUnit%v_zti(damID,ngeom+1)-WRMUnit%v_zti(damID,ngeom))
							delta_z  = (WRMUnit%d_zi(damID,ngeom+1)-WRMUnit%d_zi(damID,ngeom))*num/dem
							delta_a  = (WRMUnit%a_di(damID,WRMUnit%d_ns(damID)+1)*1._r8 -WRMUnit%a_di(damID,WRMUnit%d_ns(damID)))*num/dem
							WRMUnit%d_z(damID,i+1) = WRMUnit%d_zi(damID,ngeom) + delta_z
							WRMUnit%a_d(damID,i+1) = ((WRMUnit%a_di(damID,ngeom) + max((delta_a),0.0_r8)))
						end if
					end do
					WRMUnit%d_z(damID,i)=WRMUnit%d_z(damID,i+1)
					WRMUnit%dd_z(damID,i)=max(WRMUnit%d_z(damID,i+1)-WRMUnit%d_z(damID,i),0.0_r8)
				else !enough volume, layers don't collapses
					WRMUnit%v_zt(damID,i+1)=(WRMUnit%v_zt(damID,i) + WRMUnit%v_zn(damID,i))!dv_nt(i)
					if (WRMUnit%v_zt(damID,i+1) < 1.05_r8*WRMUnit%v_zti(damID,ngeom+1)) then ! Allowable maximum storage should include only 5% (freeboard provision)
						do j=2,ngeom+1
							if (WRMUnit%v_zt(damID,i+1)>WRMUnit%v_zti(damID,j-1).and.WRMUnit%v_zt(damID,i+1)<=WRMUnit%v_zti(damID,j))then
								num = (WRMUnit%v_zt(damID,i+1)-WRMUnit%v_zti(damID,j-1))
								dem = (WRMUnit%v_zti(damID,j)-WRMUnit%v_zti(damID,j-1))
								delta_z  = (WRMUnit%d_zi(damID,j)-WRMUnit%d_zi(damID,j-1))*num/dem
								delta_a  = (WRMUnit%a_di(damID,j)-WRMUnit%a_di(damID,j-1))*num/dem
								WRMUnit%d_z(damID,i+1) = WRMUnit%d_zi(damID,j-1) + delta_z 
								WRMUnit%a_d(damID,i+1) = ((WRMUnit%a_di(damID,j-1) + max((delta_a),0.0_r8)))
							elseif (WRMUnit%v_zt(damID,i+1)>WRMUnit%v_zti(damID,ngeom+1))then ! inflow causes maximum storage
								num = (WRMUnit%v_zt(damID,i+1)-WRMUnit%v_zti(damID,ngeom))
								dem = (WRMUnit%v_zti(damID,ngeom+1)-WRMUnit%v_zti(damID,ngeom))
								delta_z  = (WRMUnit%d_zi(damID,ngeom+1)-WRMUnit%d_zi(damID,ngeom))*num/dem
								delta_a  = (WRMUnit%a_di(damID,WRMUnit%d_ns(damID)+1) -WRMUnit%a_di(damID,WRMUnit%d_ns(damID)))*num/dem
								WRMUnit%d_z(damID,i+1) = WRMUnit%d_zi(damID,ngeom) + delta_z
								WRMUnit%a_d(damID,i+1) = ((WRMUnit%a_di(damID,ngeom) + max((delta_a),0.0_r8)))
							end if
						end do
					else
						WRMUnit%d_z(damID,i+1)  = 1.05_r8*WRMUnit%d_zi(damID,ngeom+1)
						WRMUnit%a_d(damID,i+1)  = 1.05_r8*WRMUnit%a_di(damID,ngeom+1)
						WRMUnit%v_zt(damID,i+1) = 1.05_r8*WRMUnit%v_zti(damID,ngeom+1)
					end if
				end if 
				WRMUnit%dd_z(damID,i)=max(WRMUnit%d_z(damID,i+1)-WRMUnit%d_z(damID,i),0.0_r8)
			elseif (i==WRMUnit%d_ns(damID)+1)then
				WRMUnit%d_z(damID,i)  = WRMUnit%d_z(damID,i)
				WRMUnit%a_d(damID,i)  = WRMUnit%a_d(damID,i)
				WRMUnit%v_zt(damID,i) = WRMUnit%v_zt(damID,i)
				WRMUnit%dd_z(damID,i) = 0._r8
			else
				WRMUnit%d_z(damID,i)  = 0._r8
				WRMUnit%a_d(damID,i)  = 0._r8
				WRMUnit%v_zt(damID,i) = 0._r8
				WRMUnit%dd_z(damID,i) = 0._r8
			end if
		end do
		
		! 	Recalculate layer thickness and volume	
		do j = 1,nlayers
			if (j<=WRMUnit%d_ns(damID))  then   
				WRMUnit%d_v(damID,j) 	= max((WRMUnit%v_zt(damID,j+1) - WRMUnit%v_zt(damID,j)),0.0_r8)
				WRMUnit%v_zn(damID,j) 	= WRMUnit%d_v(damID,j)
			else
				WRMUnit%d_v(damID,j)	= 0._r8
				WRMUnit%v_zn(damID,j) 	= 0._r8
			end if
		end do
		
		!	Recalculare reservoir depth
		WRMUnit%d_resrv(damID)=0._r8 
		do i=1,WRMUnit%d_ns(damID)
			WRMUnit%d_resrv(damID) = WRMUnit%d_resrv(damID)+WRMUnit%dd_z(damID,i)
		end do	
		
	end subroutine geometry

	subroutine layer_merge(damID,i,dv_in,dv_ou,enr_0)
!*******************************************************************************************************
! 	Merge layer if the thickness is less than threshold ddz_min 
!*******************************************************************************************************   
		use shr_sys_mod , only : shr_sys_flush
		use WRM_type_mod  , only :  WRMUnit
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: damID,i              
		real(r8),dimension(nlayers), intent(inout) :: dv_in,dv_ou,enr_0    ! layer inflow,outflow(m3/s), and energy(w) 
		real(r8) :: ta,tb,dd_za,dd_zb,e_a,e_b,d_va,d_vb
		real(r8) :: dv_oua,dv_oub,dv_ina,dv_inb    
		integer :: j,k,m							! indices				
		
		if(i<WRMUnit%d_ns(damID)-1 .or. (i==1 .and. WRMUnit%d_ns(damID)==2)) then
			m=i+1
		else 
			m=i-1 !Avoid merging to top layer
		end if
		ta 		= WRMUnit%temp_resrv(damID,i)
		tb 		= WRMUnit%temp_resrv(damID,m)
		dd_za	= WRMUnit%dd_z(damID,i)
		dd_zb	= WRMUnit%dd_z(damID,m)
		e_a 	= enr_0(i)
		e_b 	= enr_0(m)
		d_va	= WRMUnit%d_v(damID,i)
		d_vb	= WRMUnit%d_v(damID,m)
		dv_oua	= dv_ou(i)
		dv_oub	= dv_ou(m)
		dv_ina	= dv_in(i)
		dv_inb	= dv_in(m)
		
	!	Merge layers
		if(d_va+d_vb > 0._r8)then
			WRMUnit%temp_resrv(damID,i)	= (d_va*ta+d_vb*tb)/(d_va+d_vb)
		else
		end if
	
	! 	Adjust new layer volume, thickness, mass and inflow/outflow 
		WRMUnit%dd_z(damID,i)	= dd_za+dd_zb
		WRMUnit%d_v(damID,i)	= d_va+d_vb
		dv_ou(i)				= dv_oua+dv_oub
		dv_in(i)				= dv_ina+dv_inb
		enr_0(i)				= e_a+e_b

	!	Re-number layers before collapse
		do j=i,nlayers-1
			if (j==i .and. i<(WRMUnit%d_ns(damID)-1)) then ! 	Identify lower (small and to be merged) and upper (larger) layer 
				dv_ou(j)	= dv_ou(i)
				dv_in(j)	= dv_in(i)
				enr_0(j)	= enr_0(i)
				WRMUnit%temp_resrv(damID,j)	= WRMUnit%temp_resrv(damID,i)
				WRMUnit%dd_z(damID,j)		= WRMUnit%dd_z(damID,i)
				WRMUnit%d_v(damID,j)		= WRMUnit%d_v(damID,i)
				WRMUnit%d_z(damID,j+1)		= WRMUnit%d_z(damID,i+2)
				WRMUnit%v_zt(damID,j+1)		= WRMUnit%v_zt(damID,i+2)
				WRMUnit%a_d(damID,j+1)		= WRMUnit%a_d(damID,i+2)		
			elseif (j<WRMUnit%d_ns(damID) .and. i<(WRMUnit%d_ns(damID)-1)) then
				dv_ou(j)	= dv_ou(j+1)
				dv_in(j)	= dv_in(j+1)
				enr_0(j)	= enr_0(j+1)
				WRMUnit%temp_resrv(damID,j)	= WRMUnit%temp_resrv(damID,j+1)
				WRMUnit%dd_z(damID,j)		= WRMUnit%dd_z(damID,j+1)
				WRMUnit%d_v(damID,j)		= WRMUnit%d_v(damID,j+1)
				WRMUnit%d_z(damID,j+1)		= WRMUnit%d_z(damID,j+2)
				WRMUnit%v_zt(damID,j+1)		= WRMUnit%v_zt(damID,j+2)
				WRMUnit%a_d(damID,j+1)		= WRMUnit%a_d(damID,j+2)									
			elseif (j==i .and. i==(WRMUnit%d_ns(damID)-1)) then
				dv_ou(j-1)	= dv_ou(i)
				dv_in(j-1)	= dv_in(i)
				enr_0(j-1)	= enr_0(i)
				WRMUnit%temp_resrv(damID,j-1)	= WRMUnit%temp_resrv(damID,i)
				WRMUnit%dd_z(damID,j-1)			= WRMUnit%dd_z(damID,i)
				WRMUnit%d_v(damID,j-1)			= WRMUnit%d_v(damID,i)
				WRMUnit%d_z(damID,j)			= WRMUnit%d_z(damID,j+1)
				WRMUnit%v_zt(damID,j)			= WRMUnit%v_zt(damID,j+1)
				WRMUnit%a_d(damID,j)			= WRMUnit%a_d(damID,j+1)				
				WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID)-1)=WRMUnit%temp_resrv(damID,WRMUnit%d_ns(damID))
				dv_ou(WRMUnit%d_ns(damID)-1)	= dv_ou(WRMUnit%d_ns(damID))
				dv_in(WRMUnit%d_ns(damID)-1)	= dv_in(WRMUnit%d_ns(damID))
				WRMUnit%dd_z(damID,WRMUnit%d_ns(damID)-1)	= WRMUnit%dd_z(damID,WRMUnit%d_ns(damID))
				WRMUnit%d_v(damID,WRMUnit%d_ns(damID)-1)	= WRMUnit%d_v(damID,WRMUnit%d_ns(damID))
				enr_0(WRMUnit%d_ns(damID)-1)			= enr_0(WRMUnit%d_ns(damID))
				WRMUnit%d_z(damID,WRMUnit%d_ns(damID))	= WRMUnit%d_z(damID,WRMUnit%d_ns(damID)+1)
				WRMUnit%v_zt(damID,WRMUnit%d_ns(damID))	= WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1)
				WRMUnit%a_d(damID,WRMUnit%d_ns(damID))	= WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)				
			elseif (j>=WRMUnit%d_ns(damID)) then
				dv_ou(j)	= 0._r8
				dv_in(j)	= 0._r8
				enr_0(j)	= 0._r8
				WRMUnit%temp_resrv(damID,j)	= 0._r8
				WRMUnit%dd_z(damID,j)		= 0._r8
				WRMUnit%d_v(damID,j)		= 0._r8
				WRMUnit%d_z(damID,j+1)		= 0._r8
				WRMUnit%v_zt(damID,j+1)		= 0._r8
				WRMUnit%a_d(damID,j+1)		= 0._r8				
			end if
		end do		
		WRMUnit%d_ns(damID)	= WRMUnit%d_ns(damID)-1		

	end subroutine layer_merge
	
	subroutine layer_split(damID,i,dv_in,dv_ou,enr_0)
!*******************************************************************************************************
! 	Split layer if the thickness is more than threshold ddz_max
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use WRM_type_mod  , only :  WRMUnit
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: damID,i             
		real(r8),dimension(nlayers), intent(inout) :: dv_in,dv_ou,enr_0    ! layer inflow,outflow(m3/s), and energy(w) 
		real(r8) :: tab,e_ab,dd_zab,d_vab,dv_ouab,dv_inab,delta_z    
		integer :: j,k,m							! indices
				
		! Calculate layer geometric properties to be halved
		dd_zab 	= 0.5_r8*WRMUnit%dd_z(damID,i)
		e_ab   	= 0.5_r8*enr_0(i)
		d_vab	= WRMUnit%d_v(damID,i)
		tab		= WRMUnit%temp_resrv(damID,i)
		dv_ouab	= dv_ou(i)
		dv_inab	= dv_in(i)
		!	Re-number layers before dividing layer 
		WRMUnit%d_z(damID,WRMUnit%d_ns(damID)+2)	= WRMUnit%d_z(damID,WRMUnit%d_ns(damID)+1)
		WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+2)	= WRMUnit%v_zt(damID,WRMUnit%d_ns(damID)+1)
		WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+2)	= WRMUnit%a_d(damID,WRMUnit%d_ns(damID)+1)
		m=i+1
		do j=m,WRMUnit%d_ns(damID)
			k=WRMUnit%d_ns(damID)-j+i+1
			WRMUnit%temp_resrv(damID,k+1)	= WRMUnit%temp_resrv(damID,k)
			dv_ou(k+1)	= dv_ou(k)
			dv_in(k+1)	= dv_in(k)
			WRMUnit%d_z(damID,k+1)	= WRMUnit%d_z(damID,k)
			WRMUnit%a_d(damID,k+1)	= WRMUnit%a_d(damID,k)
			WRMUnit%v_zt(damID,k+1)	= WRMUnit%v_zt(damID,k)
			WRMUnit%d_v(damID,k+1)	= WRMUnit%d_v(damID,k)
			WRMUnit%dd_z(damID,k+1)	= WRMUnit%dd_z(damID,k)
			enr_0(k+1)= enr_0(k)
		end do
	!	Divide layer in half and calculate corresponding properties
		WRMUnit%dd_z(damID,i)	= dd_zab
		WRMUnit%dd_z(damID,i+1)	= dd_zab
		WRMUnit%d_z(damID,i+1)	= WRMUnit%d_z(damID,i)+WRMUnit%dd_z(damID,i)
		do j=2,ngeom+1
			if (WRMUnit%d_z(damID,i+1)>WRMUnit%d_zi(damID,j-1).and.WRMUnit%d_z(damID,i+1)<=WRMUnit%d_zi(damID,j))then
				delta_z   = (WRMUnit%d_z(damID,i+1)-WRMUnit%d_zi(damID,j-1))/(WRMUnit%d_zi(damID,j)-WRMUnit%d_zi(damID,j-1))
				WRMUnit%a_d(damID,i+1)  =(delta_z*(WRMUnit%a_di(damID,j)-WRMUnit%a_di(damID,j-1))+WRMUnit%a_di(damID,j-1))
				WRMUnit%v_zt(damID,i+1) = (delta_z*(WRMUnit%v_zti(damID,j)-WRMUnit%v_zti(damID,j-1))+WRMUnit%v_zti(damID,j-1))
			elseif (WRMUnit%d_z(damID,i+1)> WRMUnit%d_zi(damID,ngeom+1))then
				delta_z   = (WRMUnit%d_z(damID,i+1)-WRMUnit%d_zi(damID,ngeom))/(WRMUnit%d_zi(damID,ngeom+1)-WRMUnit%d_zi(damID,ngeom))
				WRMUnit%a_d(damID,i+1)  = (delta_z*(WRMUnit%a_di(damID,ngeom+1)-WRMUnit%a_di(damID,ngeom))+WRMUnit%a_di(damID,ngeom))
				WRMUnit%v_zt(damID,i+1) = (delta_z*(WRMUnit%v_zti(damID,ngeom+1)-WRMUnit%v_zti(damID,ngeom))+WRMUnit%v_zti(damID,ngeom))
			end if
		end do
		WRMUnit%d_v(damID,i+1)	= d_vab-(WRMUnit%v_zt(damID,i+1)-WRMUnit%v_zt(damID,i))
		WRMUnit%d_v(damID,i)	= (WRMUnit%v_zt(damID,i+1)-WRMUnit%v_zt(damID,i))
		dv_ou(i+1)	= WRMUnit%d_v(damID,i+1)*dv_ouab/(WRMUnit%d_v(damID,i)+WRMUnit%d_v(damID,i+1))
		dv_ou(i)	= WRMUnit%d_v(damID,i)*dv_ouab/(WRMUnit%d_v(damID,i)+WRMUnit%d_v(damID,i+1))
		WRMUnit%temp_resrv(damID,i)		= tab
		WRMUnit%temp_resrv(damID,i+1)	= tab
		enr_0(i)	= e_ab
		enr_0(i+1)	= e_ab 
		dv_in(i+1)	= WRMUnit%d_v(damID,i+1)*dv_inab/(WRMUnit%d_v(damID,i)+WRMUnit%d_v(damID,i+1))
		dv_in(i)	= WRMUnit%d_v(damID,i)*dv_inab/(WRMUnit%d_v(damID,i)+WRMUnit%d_v(damID,i+1))
		WRMUnit%d_ns(damID)	= WRMUnit%d_ns(damID)+1	
		
	end subroutine layer_split

	subroutine over_under_flow(damID,d_n_n,d_res_sub,t_z_old,dv_in,dv_ou)
!*******************************************************************************************************
! 	Split layer if the thickness is more than threshold ddz_max
!*******************************************************************************************************    
		use shr_sys_mod , only : shr_sys_flush
		use WRM_type_mod  , only :  WRMUnit
		use RtmVar         , only : iulog, ngeom, nlayers
		
		implicit none
		integer, intent(in) :: damID,d_n_n              
		real(r8),intent(in) :: d_res_sub              
		real(r8),dimension(nlayers), intent(in) :: t_z_old 
		real(r8),dimension(nlayers), intent(out) :: dv_in,dv_ou    ! layer inflow,outflow(m3/s)
		integer :: j							! indices
		
		WRMUnit%d_ns(damID)  	= d_n_n				
		WRMUnit%d_resrv(damID)	= d_res_sub 
		do j = 1, nlayers+1	
			if (j<=nlayers)  then   
				WRMUnit%d_v(damID,j) 		= WRMUnit%v_zo(damID,j)
				WRMUnit%temp_resrv(damID,j) = t_z_old(j)
				WRMUnit%v_zt(damID,j) 	= WRMUnit%v_zt0(damID,j)
				WRMUnit%a_d(damID,j) 	= WRMUnit%a_d0(damID,j)
				WRMUnit%d_z(damID,j) 	= WRMUnit%d_z0(damID,j)
				WRMUnit%dd_z(damID,j)	= WRMUnit%d_z0(damID,j+1)-WRMUnit%d_z0(damID,j)
				dv_in(j)	= 0._r8
				dv_ou(j)	= 0._r8
			else
				WRMUnit%v_zt(damID,j) 	= WRMUnit%v_zt0(damID,j)
				WRMUnit%a_d(damID,j) 	= WRMUnit%a_d0(damID,j)
				WRMUnit%d_z(damID,j) 	= WRMUnit%d_z0(damID,j)
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
	
    subroutine solve(a,b,c,r,damID,n) !
		use shr_sys_mod , only : shr_sys_flush
	
		implicit none
!	 a:left coefficient , b:middle coefficient, c:right coefficient, r:right hand side (known), n - number of layers
		integer,intent(in) :: n,damID
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
        WRMUnit%temp_resrv(damID,n) = max(rp(n)/bp(n),273.15_r8)	! to be modified after including ice/snow effect
		
! 	Back substitution
        do i = n-1, 1, -1
			WRMUnit%temp_resrv(damID,i) = max((rp(i)-c(i)*WRMUnit%temp_resrv(damID,i+1))/bp(i),273.15_r8)
			
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
  
end MODULE MOSART_stra_mod
