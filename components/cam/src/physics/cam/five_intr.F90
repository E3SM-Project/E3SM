module five_intr

!-------------------------------------------
! Module to interface FIVE With E3SM

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,           only : pcols, pver, pverp, begchunk, endchunk

  implicit none
  
  integer :: five_add_nlevels
  integer, parameter :: pver_five = 97
  integer, parameter :: pverp_five = 98
  
  ! define physics buffer indicies here
  integer :: t_five_idx, &
             q_five_idx, &
	     u_five_idx, &
	     v_five_idx
  
  real(r8) :: pmid_five(pver_five)
  real(r8) :: pint_five(pverp_five) 
  
  data pmid_five /012.3825412243605_r8, 018.2829238474369_r8, 026.9948869943619_r8, &
    039.8581713438034_r8, 058.8509157299995_r8, 086.8938595056534_r8, 128.29949259758_r8, &
    189.435243606567_r8, 279.702687263489_r8, 412.983322143555_r8, 596.844935417175_r8, &
    837.740421295166_r8, 1147.37873077393_r8, 1533.39381217957_r8, 1999.6337890625_r8, &
    2544.46964263916_r8, 3159.32512283325_r8, 3836.6283416748_r8, 4567.11978912354_r8, &
    5330.95607757568_r8, 6101.51805877686_r8, 6847.6390838623_r8, 7535.53352355957_r8, &
    8194.62738037109_r8, 8891.05453491211_r8, 9646.66748046875_r8, 10466.4966583252_r8, &
    11356.0001373291_r8, 12321.0994720459_r8, 13368.2182312012_r8, 14504.3266296387_r8, &
    15736.988067627_r8, 17074.4079589844_r8, 18525.4898071289_r8, 20099.893951416_r8, &
    21808.1001281738_r8, 23661.4791870117_r8, 25672.3686218262_r8, 27854.1549682617_r8, &
    30221.3638305664_r8, 32789.7506713867_r8, 35576.4129638672_r8, 38599.90234375_r8, &
    41880.3451538086_r8, 45439.5812988281_r8, 49246.858215332_r8, 53163.9526367188_r8, &
    57062.4908447266_r8, 60864.3798828125_r8, 64531.9976806641_r8, 68049.8046875_r8, &
    70567.9107666016_r8, 72173.0560302734_r8, 73711.5539550781_r8, 75183.4259033203_r8, & 
    76577.8259277344_r8, 77894.7540283203_r8, 79124.2553710938_r8, 80266.3208007812_r8, &
    81312.1215820312_r8, 82261.6638183594_r8, 83122.9614257812_r8, 83896.0540771484_r8, &
    84624.3560791016_r8, 85307.8704833984_r8, 85983.6212158203_r8, 86651.6235351562_r8, &
    87311.3708496094_r8, 87962.7899169922_r8, 88605.3924560547_r8, 89239.1815185547_r8, &
    89863.6291503906_r8, 90478.7506103516_r8, 91084.0515136719_r8, 91679.4982910156_r8, &
    92264.599609375_r8, 92839.3768310547_r8, 93403.3264160156_r8, 93956.4727783203_r8, &
    94498.3215332031_r8, 95028.8940429688_r8, 95547.7447509766_r8, 96054.833984375_r8, &
    96549.7344970703_r8, 97032.4768066406_r8, 97502.6275634766_r8, 97960.1837158203_r8, &
    98404.7332763672_r8, 98836.3037109375_r8, 99213.9038085938_r8, 99537.4938964844_r8, &
    99774.4415283203_r8, 99924.7985839844_r8, 100375.329450090_r8, 100759.95623136699_r8, &
    101144.583012644_r8, 101626.967611472_r8/  
    

    
  data pint_five /0.100000001490116, 0.147650822997093, 0.218007653951645, &
    0.321890085935593, 0.475273340940475, 0.701744973659515, &
    1.03613221645355, 1.52985763549805, 2.2588472366333, 3.33520650863647, &
    4.92445993423462, 7.01243877410889, 9.74236965179443, 13.2052049636841, &
    17.4626712799072, 22.5300045013428, 28.3593883514404, 34.8271141052246, &
    41.9054527282715, 49.4369430541992, 57.1821784973145, 64.8481826782227, &
    72.1045989990234, 78.606071472168, 85.2864761352539, 92.5346145629883, &
    100.398735046387, 108.931198120117, 118.188804626465, 128.233184814453, &
    139.13117980957, 150.955352783203, 163.784408569336, 177.703750610352, &
    192.806045532227, 209.191833496094, 226.970169067383, 246.259414672852, &
    267.187957763672, 289.895141601562, 314.532135009766, 341.262878417969, &
    370.265380859375, 401.732666015625, 435.874237060547, 472.917388916016, &
    512.019775390625, 551.25927734375, 589.990539550781, 627.297058105469, &
    663.342895507812, 697.653198242188, 713.705017089844, 729.756103515625, &
    744.474975585938, 759.193542480469, 772.362976074219, 785.532104492188, &
    796.953002929688, 808.373413085938, 817.869018554688, 827.3642578125, &
    835.094970703125, 842.826110839844, 849.661010742188, 856.496398925781, & 
    863.176025390625, 869.8564453125, 876.370971679688, 882.884826660156, &
    889.223022460938, 895.560607910156, 901.711975097656, 907.863037109375, & 
    913.817993164062, 919.77197265625, 925.52001953125, 931.267517089844, &
    936.799011230469, 942.330444335938, 947.635986328125, 952.94189453125, &
    958.013000488281, 963.083679199219, 967.911010742188, 972.738525390625, &
    977.314025878906, 981.8896484375, 986.205017089844, 990.521057128906, &
    993.757019042969, 996.992858886719, 998.495971679688, 1000/   
  
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
    call pbuf_add_field('T_FIVE',       'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), t_five_idx)
    call pbuf_add_field('Q_FIVE',        'global', dtype_r8, (/pcols,pverp_five,pcnst,dyn_time_lvls/), q_five_idx)
    call pbuf_add_field('U_FIVE',         'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), u_five_idx)
    call pbuf_add_field('V_FIVE',         'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), v_five_idx)
  
  end subroutine five_register_e3sm
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine five_init_e3sm(phys_state, pbuf2d)
  
    use time_manager, only: is_first_step
  
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    
    real(r8) :: t_host(pver)
    real(r8) :: q_host(pver)
    real(r8) :: u_host(pver)
    real(r8) :: v_host(pver)
    
    real(r8) :: t_five(pcols,pver_five)
    real(r8) :: q_five(pcols,pver_five,pcnst)
    real(r8) :: u_five(pcols,pver_five)
    real(r8) :: v_five(pcols,pver_five)
    
    integer :: ncol, i, p
    
    ncol = phys_state%ncol
    
    do i=1,ncol
    
      t_host(:) = phys_state%t(i,:)
      u_host(:) = phys_state%u(i,:)
      v_host(:) = phys_state%v(i,:)
      
      call linear_interp(phys_state%pmid(i,:),pmid_five,t_host,t_five(i,:),pver,pver_five)
      call linear_interp(phys_state%pmid(i,:),pmid_five,u_host,u_five(i,:),pver,pver_five)
      call linear_interp(phys_state%pmid(i,:),pmid_five,v_host,v_five(i,:),pver,pver_five)
      
      do p=1,pcnst
        q_host(:) = phys_state%q(i,:,p)
        call linear_interp(phys_state%pmid(i,:),pmid_five,q_host,q_five(i,:,p),pver,pver_five)
      enddo
    
    enddo
    
    call pbuf_set_field(pbuf2d, t_five_idx, t_five)
    call pbuf_set_field(pbuf2d, q_five_idx, q_five)
    call pbuf_set_field(pbuf2d, u_five_idx, u_five)
    call pbuf_set_field(pbuf2d, v_five_idx, v_five)

  ! Here we initialize the FIVE grid
  
  ! For now, for testing purposes, let's just
  !  hardcode the levels in.  Future implementaitons
  !  will involve algorithm to generate the needed
  !  levels
  
  end subroutine five_init_e3sm
  
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
