module five_intr

!-------------------------------------------
! Module to interface FIVE with E3SM physics packages

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,           only : pcols, pver, pverp, begchunk, endchunk
  use physics_types, only: physics_state
  use physics_buffer, only: physics_buffer_desc
  use constituents, only: pcnst
  use physconst, only: rair, gravit
  use cam_logfile, only: iulog

  implicit none
  
  ! This is the number of layers to add between E3SM levels
  integer :: five_add_nlevels = 1 
  
  ! The bottom layer to which we will add layers to (set
  !   to a very large value to add all the way to surface)
  real, parameter :: five_bot_toadd = 100000._r8
  
  ! The top layer to which we will add layers to
  real, parameter :: five_top_toadd = 80000._r8
  ! TASK: the parameters above should probably be moved to the
  !   E3SM atm_in namelist
  
  ! These shouldn't be hardcoded, but currently they are
  !   because PBUF needs these to be initialized.
  !   TASK: figure out a better way so these can be computed
  !   on the fly and not user specified!
  integer, parameter :: pver_five = 87
  integer, parameter :: pverp_five = 88
  
  ! define physics buffer indicies here for the FIVE
  !  variables added to the PBUF
  integer :: t_five_idx, &
             q_five_idx, &
	     u_five_idx, &
	     v_five_idx, &
	     pmid_five_idx, &
	     pint_five_idx 
  
  contains 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !    
  
  subroutine five_readnl(nlfile)
  
    ! Currently this subroutine is not called,
    !   here as a placeholder once a more flexible
    !   framework is put in to define layers added
    !   at the namelist level
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
  
    ! Register FIVE fields to the physics buffer
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols
    
    ! Define PBUF for prognostics 
    call pbuf_add_field('T_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), t_five_idx)
    call pbuf_add_field('Q_FIVE',       'global', dtype_r8, (/pcols,pver_five,pcnst,dyn_time_lvls/), q_five_idx)
    call pbuf_add_field('U_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), u_five_idx)
    call pbuf_add_field('V_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), v_five_idx)
    
    ! Define PBUF for non-prognostic variables
    ! Probably also need to save things like cloud fraction
    call pbuf_add_field('PMID_FIVE', 'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), pmid_five_idx)
    call pbuf_add_field('PINT_FIVE', 'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), pint_five_idx)
  
  end subroutine five_register_e3sm 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine init_five_profiles(phys_state, pbuf2d)
  
    ! Purpose is to initialize FIVE profiles.
    ! 1)  The interface and midlayers are defined for FIVE
    !     based on the user input
    ! 2)  We initialize the FIVE variables from the E3SM state,
    !     then interpolate onto the FIVE grid

    use time_manager, only: is_first_step
    use physics_buffer, only: pbuf_set_field, pbuf_get_chunk, &
                        dyn_time_lvls
  
    ! Input/output variables
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! Local Variables
    type(physics_buffer_desc), pointer :: pbuf2d_chunk(:)
    
    ! Prognostics on the E3SM grid
    real(r8) :: t_host(pver) ! temperature
    real(r8) :: q_host(pver) ! tracers
    real(r8) :: u_host(pver) ! U-wind
    real(r8) :: v_host(pver) ! V-wind
    
    ! Prognostics on the FIVE grid
    real(r8) :: t_five(pcols,pver_five)
    real(r8) :: q_five(pcols,pver_five,pcnst)
    real(r8) :: u_five(pcols,pver_five)
    real(r8) :: v_five(pcols,pver_five)
    
    integer :: ncol, i, p, lchnk, k, kh, ki, n
    
    real(r8) :: pint_host(pcols,pverp)
    real(r8) :: pint_five(pcols,pverp_five) ! interface on 
    real(r8) :: pmid_five(pcols,pver_five)
    real :: incr   
    
    ! Loop over all the physics "chunks" in E3SM
    do lchnk = begchunk,endchunk
    
      ncol = phys_state(lchnk)%ncol ! number of columns in this chunk
      pint_host = phys_state(lchnk)%pint ! E3SM interface pressures
      pbuf2d_chunk => pbuf_get_chunk(pbuf2d,lchnk)
    
      ! First define the FIVE grid.
      ! This is done based on the parameter "five_add_nlevels".
      ! If this value is set to 1, then one layer will be added in between
      ! the default E3SM grid, but only in the area of the grid 
      ! that we twll it to.  This is based on "five_bot_toadd" and
      ! "five_top_toadd" parameters.    
      !
      ! Add layers to the interface grid, then interpolate
      ! onto the pmid grid
      do i=1,ncol       
        kh=1 ! kh is layer index on FIVE grid
        do k=1,pverp ! k is layer index on E3SM grid
          
	  pint_five(i,kh) = pint_host(i,k) ! Copy pre-existing level 
	                                   ! to FIVE grid
	  kh=kh+1
	  ! Test to see if this is within the bounds of the grid
	  ! we want to add layers to
	  if (pint_host(i,k) .le. five_bot_toadd .and. &
	      pint_host(i,k) .ge. five_top_toadd) then
	  
	    ! If we are inside the portion of grid we want to add layers
	    ! to then compute the pressure increment (resolution)
	    incr=(pint_host(i,k+1)-pint_host(i,k))/(five_add_nlevels+1)
	    
	    ! Define the new layer
	    do ki=1,five_add_nlevels
	      pint_five(i,kh) = pint_five(i,kh-1)+incr
	      kh=kh+1
	    enddo
	    
	  endif
	  
	enddo
	write(iulog,*) 'NUMBER OF FIVE LEVELS', kh-1

      enddo
      
      ! Now define five_mid layers
      do i=1,ncol
        do k=1,pver_five
	  pmid_five(i,k) = (pint_five(i,k)+pint_five(i,k+1))/2.0_r8
	enddo
      enddo
    
      ! Now we want to initialize stuff on the FIVE grid.  
      ! Here we initialize the prognostics (temperature, tracers, u, v)
      ! Also store the FIVE pressure levels on the PBUF 
      do i=1,ncol
    
        ! Copy variables from state to FIVE, on the E3SM grid
        t_host(:) = phys_state(lchnk)%t(i,:)
        u_host(:) = phys_state(lchnk)%u(i,:)
        v_host(:) = phys_state(lchnk)%v(i,:)
      
        ! Now interpolate onto the FIVE grid
        call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),t_host,t_five(i,:),pver,pver_five)
        call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),u_host,u_five(i,:),pver,pver_five)
        call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),v_host,v_five(i,:),pver,pver_five)
      
        ! For Q constituents 
        do p=1,pcnst
          q_host(:) = phys_state(lchnk)%q(i,:,p)
          call linear_interp(phys_state(lchnk)%pmid(i,:),pmid_five(i,:),q_host,q_five(i,:,p),pver,pver_five)
        enddo
    
      enddo

      ! Store all needed variables onto the PBUF, so they can be used by the 
      !  parameterizations.  Loop over all the dynamics time levels, which is needed
      !  when Eulerian Dynamical core is used since that value is 2 (otherwise the second
      !  dynamics time level will be initialized with NaNs
      do n=1,dyn_time_lvls
        call pbuf_set_field(pbuf2d_chunk, t_five_idx, t_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, q_five_idx, q_five,(/1,1,1,n/),(/pcols,pver_five,pcnst,1/))
        call pbuf_set_field(pbuf2d_chunk, u_five_idx, u_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, v_five_idx, v_five,(/1,1,n/),(/pcols,pver_five,1/))
        ! Also save pmid values to be used by parameterizations
        call pbuf_set_field(pbuf2d_chunk, pmid_five_idx, pmid_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, pint_five_idx, pint_five,(/1,1,n/),(/pcols,pverp_five,1/))
      enddo
    
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
    !   with the E3SM state. 
    ! This subroutine should be called BEFORE
    !   a parameterization is called that uses FIVE
    
    ! 1) A mass weighted vertical average will
    !    be called to get FIVE values of T, q, u, v
    !    onto E3SM grid.
    ! 2) The tendency will be computed on the 
    !    the E3SM grid
    ! 3) This tendency will be interpolated onto
    !    the high resolution FIVE grid
    ! 4) Finally, this tendency will be applied to the 
    !    FIVE prognostics, for updated values to be 
    !    used by parameterizations 
	
    ! Input/Output variables		
    type(physics_state), intent(in) :: state
    real(r8), intent(in) :: dtime
    real(r8), intent(in) :: pmid_five(pcols,pver_five)
    real(r8), intent(in) :: pint_five(pcols,pverp_five)
    real(r8), intent(inout) :: t_five(pcols,pver_five)
    real(r8), intent(inout) :: u_five(pcols,pver_five)
    real(r8), intent(inout) :: v_five(pcols,pver_five)
    real(r8), intent(inout) :: q_five(pcols,pver_five,pcnst)    
	
    ! Local variables		
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
    
    real(r8) :: pdel_five(pcols,pver_five)
        
    ncol = state%ncol
  
    ! Compute pressure differences on FIVE grid  
    do k=1,pver_five
      do i=1,ncol
        pdel_five(i,k) = pint_five(i,k+1)-pint_five(i,k)
      enddo
    enddo
  
    ! First compute a mass weighted average of the five variables onto the 
    !   the lower resolution E3SM grid   
    do i=1,ncol
      ! Mass weighted vertical average for temperature
      call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    t_five(i,:),t_five_low(i,:))
			    
      ! Mass weighted vertical average for u wind
      call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    u_five(i,:),u_five_low(i,:))
			    
      ! Mass weighted vertical average for v wind
      call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    v_five(i,:),v_five_low(i,:))			    
      
      ! Mass weighted vertical average for tracers 
      do p=1,pcnst
        call masswgt_vert_avg(state%t(i,:),t_five(i,:),state%pdel(i,:),pdel_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    q_five(i,:,p),q_five_low(i,:,p))
      enddo
      
    enddo  
   
    ! Next compute the tendency of FIVE variables from the state, this is 
    !   done on the E3SM grid
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

    ! Now interpolate this tendency to the higher resolution FIVE grid 
    !  TASK: interpolation scheme needs to be updated  
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

    ! Finally, update FIVE prognostic variables based on this tendency
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
    real(r8), intent(in) :: pint_host(pverp)
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
    real(r8) :: rho_host_avg(pver)
    real(r8) :: dz_high(pver_five)
    real(r8) :: dz_host(pver)
    
    ! Define the density on the E3SM grid   
    do k=1,pver
      rho_host(k) = pmid_host(k)/(rair*t_host(k))
    enddo
    
    ! Compute density on FIVE grid
    do k=1,pver_five
      rho_high(k) = pmid_high(k)/(rair*t_high(k))
    enddo
    
    ! define dz on E3SM grid
    do k=1,pver
      dz_host(k) = pdel_host(k)/(rho_host(k)*gravit)
    enddo
    
    ! define dz on the FIVE grid
    do k=1,pver_five
      dz_high(k) = pdel_high(k)/(rho_high(k)*gravit)
    enddo
    
    ! Initialize host variable
    var_host(:) = 0._r8
    rho_host_avg(:) = 0._r8
    
    kh=1 ! Vertical index for FIVE 
    do k=1,pver ! Vertical index for E3SM
       
      ! Check to see how many FIVE layers there are within
      !   an E3SM layer 
      do while ( pmid_high(kh) .lt. pint_host(k+1) .and. &
                 pmid_high(kh) .gt. pint_host(k))
      
        var_host(k) = var_host(k) + rho_high(kh) * &
	  var_high(kh) * dz_high(kh)
	rho_host_avg(k) = rho_host_avg(k) + rho_high(kh) * &
	  dz_high(kh)
	  
	kh = kh + 1 ! increase high res model by one layer
	if (kh .gt. pver_five) goto 10
	
      end do ! end while loop for kh
10 continue
      
      ! Compute rho on host grid
      rho_host_avg(k) = rho_host_avg(k)/dz_host(k)
      
      ! Compute the mass weighted vertical average
!      var_host(k) = var_host(k)/(rho_host(k)*dz_host(k))
      var_host(k) = var_host(k)/(rho_host_avg(k)*dz_host(k))
      
    enddo
    
    return
  
  end subroutine masswgt_vert_avg     
  
  ! ======================================== !
  !                                          !
  ! ======================================== !  
  
  subroutine linear_interp(x1,x2,y1,y2,km1,km2)
    implicit none

    ! Simple linear interpolation routine
    ! WILL LIKELY BE REPLACED
    
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
