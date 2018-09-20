module five_intr

!-------------------------------------------
! Module to interface FIVE with E3SM

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,           only : pcols, pver, pverp, begchunk, endchunk
  use physics_types, only: physics_state
  use physics_buffer, only: physics_buffer_desc
  use constituents, only: pcnst
  use physconst, only: rair, gravit

  implicit none
  
  ! This is the number of layers to add between E3SM levels
  integer, parameter :: five_add_nlevels = 1 
  
  ! The bottom layer to which we will add layers to (set
  !   to a very large value to add all the way to surface)
  real, parameter :: five_bot_toadd = 90000._r8
  
  ! The top layer to which we will add layers to
  real, parameter :: five_top_toadd = 80000._r8
  
  ! These shouldn't be hardcoded +PAB figure this out!
  integer, parameter :: pver_five = 97
  integer, parameter :: pverp_five = 98
  
  ! define physics buffer indicies here
  integer :: t_five_idx, &
             q_five_idx, &
	     u_five_idx, &
	     v_five_idx, &
	     pmid_five_idx, &
	     pint_five_idx
  
!  real(r8) :: pmid_five(pver_five)
!  real(r8) :: pint_five(pverp_five) 

!  data pint_five /010.0000001490116_r8, 014.7650822997093_r8, 021.8007653951645_r8, 
!    032.1890085935593_r8, 047.5273340940475_r8, 070.1744973659515_r8, 
!    103.613221645355_r8, 152.985763549805_r8, 225.88472366333_r8, 333.520650863647_r8, 
!    492.445993423462_r8, 701.243877410889_r8, 974.236965179443_r8, 1320.52049636841_r8, 
!    1746.26712799072_r8, 2253.00045013428_r8, 2835.93883514404_r8, 3482.71141052246_r8, 
!    4190.54527282715_r8, 4943.69430541992_r8, 5718.21784973145_r8, 6484.81826782227_r8, 
!    7210.45989990234_r8, 7860.6071472168_r8, 8528.64761352539_r8, 9253.46145629883_r8, 
!    10039.8735046387_r8, 10893.1198120117_r8, 11818.8804626465_r8, 12823.3184814453_r8, 
!    13913.117980957_r8, 15095.5352783203_r8, 16378.4408569336_r8, 17770.3750610352_r8, 
!    19280.6045532227_r8, 20919.1833496094_r8, 22697.0169067383_r8, 24625.9414672852_r8, 
!    26718.7957763672_r8, 28989.5141601562_r8, 31453.2135009766_r8, 34126.2878417969_r8, 
!    37026.5380859375_r8, 40173.2666015625_r8, 43587.4237060547_r8, 47291.7388916016_r8, 
!    51201.9775390625_r8, 55125.927734375_r8, 58999.0539550781_r8, 62729.7058105469_r8, 
!    66334.2895507812_r8, 69765.3198242188_r8, 71370.5017089844_r8, 72975.6103515625_r8, 
!    74447.4975585938_r8, 75919.3542480469_r8, 77236.2976074219_r8, 78553.2104492188_r8, 
!    79695.3002929688_r8, 80837.3413085938_r8, 81786.9018554688_r8, 82736.42578125_r8, 
!    83509.4970703125_r8, 84282.6110839844_r8, 84966.1010742188_r8, 85649.6398925781_r8, 
!    86317.6025390625_r8, 86985.64453125_r8, 87637.0971679688_r8, 88288.4826660156_r8, 
!    88922.3022460938_r8, 89556.0607910156_r8, 90171.1975097656_r8, 90786.3037109375_r8, 
!    91381.7993164062_r8, 91977.197265625_r8, 92552.001953125_r8, 93126.7517089844_r8, 
!    93679.9011230469_r8, 94233.0444335938_r8, 94763.5986328125_r8, 95294.189453125_r8, 
!    95801.3000488281_r8, 96308.3679199219_r8, 96791.1010742188_r8, 97273.8525390625_r8, 
!    97731.4025878906_r8, 98188.96484375_r8, 98620.5017089844_r8, 99052.1057128906_r8, 
!    99375.7019042969_r8, 99699.2858886719_r8, 99849.5971679688_r8, 100000._r8/
  
!  data pmid_five /012.3825412243605_r8, 018.2829238474369_r8, 026.9948869943619_r8, &
!    039.8581713438034_r8, 058.8509157299995_r8, 086.8938595056534_r8, 128.29949259758_r8, &
!    189.435243606567_r8, 279.702687263489_r8, 412.983322143555_r8, 596.844935417175_r8, &
!    837.740421295166_r8, 1147.37873077393_r8, 1533.39381217957_r8, 1999.6337890625_r8, &
!    2544.46964263916_r8, 3159.32512283325_r8, 3836.6283416748_r8, 4567.11978912354_r8, &
!    5330.95607757568_r8, 6101.51805877686_r8, 6847.6390838623_r8, 7535.53352355957_r8, &
!    8194.62738037109_r8, 8891.05453491211_r8, 9646.66748046875_r8, 10466.4966583252_r8, &
!    11356.0001373291_r8, 12321.0994720459_r8, 13368.2182312012_r8, 14504.3266296387_r8, &
!    15736.988067627_r8, 17074.4079589844_r8, 18525.4898071289_r8, 20099.893951416_r8, &
!    21808.1001281738_r8, 23661.4791870117_r8, 25672.3686218262_r8, 27854.1549682617_r8, &
!    30221.3638305664_r8, 32789.7506713867_r8, 35576.4129638672_r8, 38599.90234375_r8, &
!    41880.3451538086_r8, 45439.5812988281_r8, 49246.858215332_r8, 53163.9526367188_r8, &
!    57062.4908447266_r8, 60864.3798828125_r8, 64531.9976806641_r8, 68049.8046875_r8, &
!    70567.9107666016_r8, 72173.0560302734_r8, 73711.5539550781_r8, 75183.4259033203_r8, & 
!    76577.8259277344_r8, 77894.7540283203_r8, 79124.2553710938_r8, 80266.3208007812_r8, &
!    81312.1215820312_r8, 82261.6638183594_r8, 83122.9614257812_r8, 83896.0540771484_r8, &
!    84624.3560791016_r8, 85307.8704833984_r8, 85983.6212158203_r8, 86651.6235351562_r8, &
!    87311.3708496094_r8, 87962.7899169922_r8, 88605.3924560547_r8, 89239.1815185547_r8, &
!    89863.6291503906_r8, 90478.7506103516_r8, 91084.0515136719_r8, 91679.4982910156_r8, &
!    92264.599609375_r8, 92839.3768310547_r8, 93403.3264160156_r8, 93956.4727783203_r8, &
!    94498.3215332031_r8, 95028.8940429688_r8, 95547.7447509766_r8, 96054.833984375_r8, &
!    96549.7344970703_r8, 97032.4768066406_r8, 97502.6275634766_r8, 97960.1837158203_r8, &
!    98404.7332763672_r8, 98836.3037109375_r8, 99213.9038085938_r8, 99537.4938964844_r8, &
!    99774.4415283203_r8, 99924.7985839844_r8, 100375.329450090_r8, 100759.95623136699_r8, &
!    101144.583012644_r8, 101626.967611472_r8/  
    

    
!  data pint_five /0.100000001490116, 0.147650822997093, 0.218007653951645, &
!    0.321890085935593, 0.475273340940475, 0.701744973659515, &
!    1.03613221645355, 1.52985763549805, 2.2588472366333, 3.33520650863647, &
!    4.92445993423462, 7.01243877410889, 9.74236965179443, 13.2052049636841, &
!    17.4626712799072, 22.5300045013428, 28.3593883514404, 34.8271141052246, &
!    41.9054527282715, 49.4369430541992, 57.1821784973145, 64.8481826782227, &
!    72.1045989990234, 78.606071472168, 85.2864761352539, 92.5346145629883, &
!    100.398735046387, 108.931198120117, 118.188804626465, 128.233184814453, &
!    139.13117980957, 150.955352783203, 163.784408569336, 177.703750610352, &
!    192.806045532227, 209.191833496094, 226.970169067383, 246.259414672852, &
!    267.187957763672, 289.895141601562, 314.532135009766, 341.262878417969, &
!    370.265380859375, 401.732666015625, 435.874237060547, 472.917388916016, !&
!    512.019775390625, 551.25927734375, 589.990539550781, 627.297058105469, &
!    663.342895507812, 697.653198242188, 713.705017089844, 729.756103515625, &
!    744.474975585938, 759.193542480469, 772.362976074219, 785.532104492188, &
!    796.953002929688, 808.373413085938, 817.869018554688, 827.3642578125, &
!    835.094970703125, 842.826110839844, 849.661010742188, 856.496398925781, & 
!    863.176025390625, 869.8564453125, 876.370971679688, 882.884826660156, &
!    889.223022460938, 895.560607910156, 901.711975097656, 907.863037109375, & 
!    913.817993164062, 919.77197265625, 925.52001953125, 931.267517089844, &
!    936.799011230469, 942.330444335938, 947.635986328125, 952.94189453125, &
!    958.013000488281, 963.083679199219, 967.911010742188, 972.738525390625, &
!    977.314025878906, 981.8896484375, 986.205017089844, 990.521057128906, &
!    993.757019042969, 996.992858886719, 998.495971679688, 1000/   
  
  contains 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !    
  
  subroutine five_readnl(nlfile)
  
    use spmd_utils,    only: masterproc  
    use units,         only: getunit, freeunit
    use namelist_utils,  only: find_group_name
    use cam_abortutils,         only: endrun
     
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
    integer :: iunit, read_status
    
    namelist /five_nl/ five_add_nlevels
    
    !----- Begin Code -----
    
    !  Initialize some FIVE parameters
    five_add_nlevels = 1 
    
    !  Read namelist to determine FIVE parameters
    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(nlfile), status='old' ) 
      
      call find_group_name(iunit, 'five_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=five_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('five_readnl:  error reading namelist')
         end if
      end if        
      
      close(unit=iunit)
      call freeunit(iunit)
    end if
    
#ifdef SPMD
! Broadcast namelist variables
      call mpibcast(five_add_nlevels,            1,   mpiint,   0, mpicom)
#endif           
  
  end subroutine five_readnl
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine five_register_e3sm
  
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols
    
    ! Define PBUF for prognostics 
    call pbuf_add_field('T_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), t_five_idx)
    call pbuf_add_field('Q_FIVE',        'global', dtype_r8, (/pcols,pver_five,pcnst,dyn_time_lvls/), q_five_idx)
    call pbuf_add_field('U_FIVE',         'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), u_five_idx)
    call pbuf_add_field('V_FIVE',         'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), v_five_idx)
    ! Not prognostic, but save PMID and PINT for FIVE
    call pbuf_add_field('PMID_FIVE', 'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), pmid_five_idx)
    call pbuf_add_field('PINT_FIVE', 'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), pint_five_idx)
  
  end subroutine five_register_e3sm 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine init_five_profiles(phys_state, pbuf2d)
  
    ! Purpose is to initialize FIVE profiles.  
    ! Here we initialize using the state variables (on E3SM grid),
    !   then interpolate onto the FIVE grid
    use time_manager, only: is_first_step
    use physics_buffer, only: pbuf_set_field, pbuf_get_chunk
  
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf2d_chunk(:)
    
    real(r8) :: t_host(pver)
    real(r8) :: q_host(pver)
    real(r8) :: u_host(pver)
    real(r8) :: v_host(pver)
    
    real(r8) :: t_five(pcols,pver_five)
    real(r8) :: q_five(pcols,pver_five,pcnst)
    real(r8) :: u_five(pcols,pver_five)
    real(r8) :: v_five(pcols,pver_five)
    
    integer :: ncol, i, p, lchnk, k, kh
    
    real(r8) :: pint_host(pcols,pverp)
    real(r8) :: pint_five(pcols,pverp_five) ! interface on 
    real(r8) :: pmid_five(pcols,pver_five)   
    
    do lchnk = begchunk,endchunk
    
      ncol = phys_state(lchnk)%ncol
      pint_host = phys_state(lchnk)%pint
      pbuf2d_chunk => pbuf_get_chunk(pbuf2d,lchnk)
    
      ! First define the FIVE grid.  
      !  Add layers to the interface grid, then interpolate
      !  onto the pmid grid
      do i=1,ncol       
        kh=1
        do k=1,pver
          
	  pint_five(i,kh) = pint_host(i,k)
	  kh=kh+1
	  if (pint_host(i,k) .le. five_bot_toadd .and. &
	      pint_host(i,k) .ge. five_top_toadd) then
	  
	    incr=(pint_host(i,k)+pint_host(i,k+1))/(five_add_nlevels+1)
	    
	    do ki=1,five_add_nlevels
	      pint_five(i,kh) = pint_five(i,kh-1)+incr
	      kh=kh+1
	    enddo
	    
	  endif
	  
	enddo
      enddo
      
      ! Now define five_mid layers
      do i=1,ncol
        do k=1,pver_five
	  pmid_five(i,k) = (pint_five(i,k)+pint_five(i,k+1))/2.0_r8
	enddo
      enddo
    
      do i=1,ncol
    
        t_host(:) = phys_state(lchnk)%t(i,:)
        u_host(:) = phys_state(lchnk)%u(i,:)
        v_host(:) = phys_state(lchnk)%v(i,:)
      
        call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),t_host,t_five(i,:),pver,pver_five)
        call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),u_host,u_five(i,:),pver,pver_five)
        call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),v_host,v_five(i,:),pver,pver_five)
      
        do p=1,pcnst
          q_host(:) = phys_state(lchnk)%q(i,:,p)
          call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),q_host,q_five(i,:,p),pver,pver_five)
        enddo
    
      enddo
    
      call pbuf_set_field(pbuf2d_chunk, t_five_idx, t_five(:,:))
      call pbuf_set_field(pbuf2d_chunk, q_five_idx, q_five(:,:,:))
      call pbuf_set_field(pbuf2d_chunk, u_five_idx, u_five(:,:))
      call pbuf_set_field(pbuf2d_chunk, v_five_idx, v_five(:,:))
      ! Also save pmid values to be used by parameterizations
      call pbuf_set_field(pbuf2d_chunk, pmid_five_idx, pmid_five(:,:))
      call pbuf_set_field(pbuf2d_chunk, pint_five_idx, pint_five(:,:))
    
    enddo
  
  end subroutine init_five_profiles
  
  ! ======================================== !
  !                                          !
  ! ======================================== !  
  
  subroutine five_syncronize_e3sm(&
                state, dtime, &
		pint_five,pmid_five, &
		t_five, u_five, v_five, q_five)

    ! Subroutine will syncronize FIVE variables
    !  with the E3SM state.  
		
    type(physics_state), intent(in) :: state
    real(r8) :: dtime
    real(r8) :: pmid_five(pcols,pver_five)
    real(r8) :: pint_five(pcols,pverp_five)
    real(r8) :: t_five(pcols,pver_five)
    real(r8) :: u_five(pcols,pver_five)
    real(r8) :: v_five(pcols,pver_five)
    real(r8) :: q_five(pcols,pver_five,pcnst)    
		
    integer :: k, i, p, ncol
    
    real(r8) :: t_five_low(pcols,pver)
    real(r8) :: u_five_low(pcols,pver)
    real(r8) :: v_five_low(pcols,pver)
    real(r8) :: q_five_low(pcols,pver,pcnst)
    
    real(r8) :: t_five_tend_low(pcols,pver)
    real(r8) :: u_five_tend_low(pcols,pver)
    real(r8) :: v_five_tend_low(pcols,pver)
    real(r8) :: q_five_tend_low(pcols,pver,pcnst)
    
    real(r8) :: t_five_tend(pcols,pver_five)
    real(r8) :: u_five_tend(pcols,pver_five)
    real(r8) :: v_five_tend(pcols,pver_five)
    real(r8) :: q_five_tend(pcols,pver_five,pcnst)
    
    real(r8) :: pdel_five(pcols,pver)
        
    ncol = state%ncol
  
    ! Compute pressure differences on FIVE grid  
    do k=1,pver
      do i=1,ncol
        pdel_five(i,k) = pint(i,k+1)-pint(i,k)
      enddo
    enddo
  
    ! First compute a mass weighted average of the five variables onto the 
    !   the lower resolution E3SM grid   
    do i=1,ncol
      ! Mass weight for temperature
      call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    t_five(i,:),state%t(i,:))
			    
      ! Mass weight for u wind
      call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    u_five(i,:),state%u(i,:))
			    
      ! Mass weight for v wind
      call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    v_five(i,:),state%v(i,:))			    
      
      do p=1,pcnst
        call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    q_five(i,:,q),state%q(i,:,p))
      enddo
    enddo  
   
    ! Next compute the tendency of FIVE variables from the state
    do k=1,pver
      do i=1,ncol
        t_five_tend_low(i,k) = (state%t(i,k) - t_five_low(i,k))/dtime 
        u_five_tend_low(i,k) = (state%u(i,k) - u_five_low(i,k))/dtime 
        v_five_tend_low(i,k) = (state%v(i,k) - v_five_low(i,k))/dtime
       
        do p=1,pcnst
          q_five_tend_low(i,k,p) = (state%q(i,k,p) - q_five_low(i,k,p))/dtime
        enddo
        
      enddo
    enddo

    ! Now interpolate this tendency to the FIVE grid   
    do i=1,ncol
      call linear_interp(state%pmid(i,:),pmid_five,t_five_tend_low(i,1:pver),&
                      t_five_tend(i,1:pver_five),pver,pver_five)
      call linear_interp(state%pmid(i,:),pmid_five,u_five_tend_low(i,1:pver),&
                      u_five_tend(i,1:pver_five),pver,pver_five)	
      call linear_interp(state%pmid(i,:),pmid_five,v_five_tend_low(i,1:pver),&
                      v_five_tend(i,1:pver_five),pver,pver_five)	
      do p=1,pcnst
        call linear_interp(state%pmid(i,:),pmid_five,q_five_tend_low(i,1:pver,p), &
                      q_five_tend(i,1:pver_five,p),pver,pver_five)
      enddo	
    enddo      	      

    ! Now update FIVE prognostic variables based on this tendency
    do k=1,pver_five
      do i=1,ncol
        t_five(i,k) = t_five(i,k) + dtime * t_five_tend(i,k)
        u_five(i,k) = u_five(i,k) + dtime * u_five_tend(i,k)
        v_five(i,k) = v_five(i,k) + dtime * v_five_tend(i,k)
       
        do p=1,pcnst
          q_five(i,k,p) = q_five(i,k,p) + dtime * q_five_tend(i,k,p)
        enddo
       
      enddo
    enddo    
		
    return
		
   end subroutine five_syncronize_e3sm
   
  ! ======================================== !
  !                                          !
  ! ======================================== !  
  
  subroutine masswgt_vert_avg(&
               t_host,t_high,pdel_host,pdel_high,&
	       pint_host,pmid_high,pmid_host,&
	       var_high,var_host)
	       
    ! Purpose is to compute the mass weighted vertical 
    !  average from the FIVE high resolution grid onto the
    !  E3SM grid
    implicit none
   
    real(r8), intent(in) :: pmid_high(pver_five)
    real(r8), intent(in) :: pint_host(pver+1)
    real(r8), intent(in) :: pmid_host(pver)
    real(r8), intent(in) :: var_high(pver_five)
    real(r8), intent(in) :: t_host(pver)
    real(r8), intent(in) :: t_high(pver_five)
    real(r8), intent(in) :: pdel_host(pver)
    real(r8), intent(in) :: pdel_high(pver_five)
    real(r8), intent(out) :: var_host(pver)
    
    integer :: i, k, kh
    real(r8) :: rho_high(pver_five)
    real(r8) :: rho_host(pver)
    real(r8) :: dz_high(pver_five)
    real(r8) :: dz_host(pver)
    
    ! Define the density on the host and 
    !  high resolution grids    
    do k=1,pver
      rho_host(k) = pmid_host(k)/(rair*t_host(k))
    enddo
    
    do k=1,pver_five
      rho_high(k) = pmid_high(k)/(rair*t_high(k))
    enddo
    
    ! define dz_host and dz_high
    do k=1,pver
      dz_host(k) = pdel_host(k)/(rho_host(k)*gravit)
    enddo
    
    do k=1,pver_five
      dz_high(k) = pdel_high(k)/(rho_high(k)*gravit)
    enddo
    
    ! Initialize host variable
    var_host(:) = 0._r8
    
    kh=1
    do k=1,pver
      
      do while ( pmid_high(kh) .lt. pint_host(k+1) .and. &
                 pmid_high(kh) .gt. pint_host(k))
      
        var_host(k) = var_host(k) + rho_high(kh) * &
	  var_high(kh) * dz_high(kh)
	  
	kh = kh + 1 ! increase high res model by one layer
	
      end do ! end while loop for kh
      
      var_host(k) = var_host(k)/(rho_host(k)*dz_host(k))
      
    enddo
    
    return
  
  end subroutine masswgt_vert_avg     
  
  ! ======================================== !
  !                                          !
  ! ======================================== !  
  
  subroutine linear_interp(x1,x2,y1,y2,km1,km2)
    implicit none

    integer :: km1, km2, k1, k2
    real(KIND=8) :: x1(km1), y1(km1)
    real(KIND=8) :: x2(km2), y2(km2)

    do k2=1,km2

      if( x2(k2) <= x1(1) ) then
        y2(k2) = y1(1) + (y1(2)-y1(1))*(x2(k2)-x1(1))/(x1(2)-x1(1))
      elseif( x2(k2) >= x1(km1) ) then
        y2(k2) = y1(km1) + (y1(km1)-y1(km1-1))*(x2(k2)-x1(km1))/(x1(km1)-x1(km1-1))    
      else
        do k1 = 2,km1
          if( (x2(k2)>=x1(k1-1)).and.(x2(k2)<x1(k1)) ) then
            y2(k2) = y1(k1-1) + (y1(k1)-y1(k1-1))*(x2(k2)-x1(k1-1))/(x1(k1)-x1(k1-1))
          endif
        enddo
      endif
     
    enddo

  end subroutine linear_interp  

end module five_intr
