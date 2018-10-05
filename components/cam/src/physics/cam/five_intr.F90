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
  integer :: five_add_nlevels = 2 
  
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
  integer, parameter :: pver_five = 102
  integer, parameter :: pverp_five = 103
  
  ! define physics buffer indicies here for the FIVE
  !  variables added to the PBUF
  integer :: t_five_idx, &
             q_five_idx, &
	     u_five_idx, &
	     v_five_idx, &
	     pmid_five_idx, &
	     pint_five_idx 
	     
  integer :: five_bot_k, five_top_k ! indicees where
                                    ! levels are added
  
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
    
    real(r8) :: dz_e3sm(pcols,pver)
    real(r8) :: dz_five(pcols,pver_five)
    
    real(r8) :: rho_e3sm(pcols,pver)
    real(r8) :: rho_five(pcols,pver_five)    
    
    real(r8) :: pdel_five(pcols,pver_five)
        
    ncol = state%ncol
  
    ! Compute pressure differences on FIVE grid  
    do k=1,pver_five
      do i=1,ncol
        pdel_five(i,k) = pint_five(i,k+1)-pint_five(i,k)
      enddo
    enddo
    
    ! Compute grid stuff needed for vertical averaging and tendencies
    do i=1,ncol
      call compute_five_grids(state%t(i,:),state%pdel(i,:),state%pmid(i,:),pver,&
             dz_e3sm(i,:),rho_e3sm(i,:))
      call compute_five_grids(t_five(i,:),pdel_five(i,:),pmid_five(i,:),pver_five,&
             dz_five(i,:),rho_five(i,:))
    enddo
  
    ! First compute a mass weighted average of the five variables onto the 
    !   the lower resolution E3SM grid   
    do i=1,ncol
      ! Mass weighted vertical average for temperature
      call masswgt_vert_avg(rho_e3sm(i,:),rho_five(i,:),dz_e3sm(i,:),dz_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    t_five(i,:),t_five_low(i,:))
			    
      ! Mass weighted vertical average for u wind
      call masswgt_vert_avg(rho_e3sm(i,:),rho_five(i,:),dz_e3sm(i,:),dz_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    u_five(i,:),u_five_low(i,:))
			    
      ! Mass weighted vertical average for v wind
      call masswgt_vert_avg(rho_e3sm(i,:),rho_five(i,:),dz_e3sm(i,:),dz_five(i,:),&
                            state%pint(i,:),pmid_five(i,:),state%pmid(i,:),&
			    v_five(i,:),v_five_low(i,:))			    
      
      ! Mass weighted vertical average for tracers 
      do p=1,pcnst
        call masswgt_vert_avg(rho_e3sm(i,:),rho_five(i,:),dz_e3sm(i,:),dz_five(i,:),&
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
  
  subroutine compute_five_grids(&
               temp,pdel,pmid,nlev,&
	       dz,rho)
	       
    implicit none

    ! Compute thickness and density on the E3SM
    !   and FIVE grid	
    real(r8), intent(in) :: temp(nlev)
    real(r8), intent(in) :: pdel(nlev)
    real(r8), intent(in) :: pmid(nlev)  
    integer, intent(in) :: nlev  
    
    real(r8), intent(out) :: dz(nlev)
    real(r8), intent(out) :: rho(nlev)
    
    integer :: k
    
    ! Compute density
    do k=1,nlev
      rho(k) = pmid(k)/(rair*temp(k))
    enddo          
	      
    ! Compute dz	     
    do k=1,nlev
      dz(k) = pdel(k)/(rho(k)*gravit)	      
    enddo	      
	      
    return
	       
  end subroutine compute_five_grids	 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !   
  
  subroutine masswgt_vert_avg(&
               rho_host,rho_high,dz_host,dz_high,&
	       pint_host,pmid_high,pmid_host,&
	       var_high,var_host)
	       
    ! Purpose is to compute the mass weighted vertical 
    !  average from the FIVE high resolution grid onto the
    !  E3SM grid
    implicit none

    real(r8), intent(in) :: rho_high(pver_five)
    real(r8), intent(in) :: rho_host(pver)   
    real(r8), intent(in) :: dz_high(pver_five)
    real(r8), intent(in) :: dz_host(pver)
    real(r8), intent(in) :: pmid_high(pver_five)
    real(r8), intent(in) :: pint_host(pverp)
    real(r8), intent(in) :: pmid_host(pver)
    real(r8), intent(in) :: var_high(pver_five)
    real(r8), intent(out) :: var_host(pver)
    
    integer :: i, k, kh

    real(r8) :: rho_host_avg(pver)
    
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
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine tendency_low_to_high(&
               zm_in, zi_in, &
	       zm_five_in, &
	       rho_low_in, rho_five_in, &
	       ten_low,ten_high)
    implicit none
    
    ! Input variables
    real(r8), intent(in) :: zm_in(pver) ! midpoint levels E3SM
    real(r8), intent(in) :: zi_in(pverp) ! interface levels E3SM
    real(r8), intent(in) :: zm_five_in(pver_five) ! midpoint levels for FIVE
    real(r8), intent(in) :: ten_low(pver)
    real(r8), intent(in) :: rho_low_in(pver)
    real(r8), intent(in) :: rho_five_in(pver_five)

    ! Output variables
    real(r8), intent(out) :: ten_high(pver_five)
    
    ! Local variables
    integer, parameter :: ml=1, mu=1, lda=2*ml+mu+1
    real(r8) :: zm(pver) ! flipped version of zm_in
    real(r8) :: zi(pver) ! flipped version of zi_in
    real(r8) :: df_z(pver) 
    real(r8) :: rho_low(pver)
    real(r8) :: rho_zs(pver_five)
    real(r8) :: zm_five(pver) ! flipped
    real(r8) :: dz ! distance of lowest level
    real(r8) :: alpha
    real(r8) :: adz_dn(pver)
    real(r8) :: adz_up(pver)
    real(r8) :: a(lda,pver)
    real(r8) :: adz(pver) ! ratio of the thickness of scalar levels to dz
    real(r8) :: adzw(pver) ! ratio of the thickness of w levels to dz
    real(r8) :: ipiv(pver)
    integer :: info
    real(r8) :: weight(pver_five)
    real(r8) :: rdf_z(pver)
    real(r8) :: rdf_zs_ml(pver) ! mid-layer target value
    real(r8) :: df_zs(pver_five)
    real(r8), dimension(pver) :: rdf_zm, rdf_zm_dn, rdf_zm_up
    real(r8), dimension(pver) :: c0, c1, c2, c3, ic1, ic2, ic3
    real(r8) :: b(pver)
    
    logical, dimension(pver) :: spurious
    logical :: cnd1, cnd2, cnd3, cnd4, cnd5
    
    integer :: i, i1, i2, i3, i4, i5
    integer :: k, km3, km2, km1, k00
    integer :: kp1, kp2, ierr
    integer :: zi1, zi2
    integer :: i5zi1
    
    ! Construct tendency profile from E3SM grid to FIVE grid
    
    ! The constructed tendency profile satisfies layer mean average
    
    ! If no levels are added via FIVE just return a copy
    if (five_add_nlevels .eq. 0) then
      ten_high(1:pver) = ten_low(1:pver)
      return
    endif
    
    ! First let's "flip" things so lowest level index = 1
    do k=1,pver
      zm(k) = zm_in(pver-k+1)
      rho_low(k) = rho_low_in(pver-k+1)
      df_z(k) = ten_low(pver-k+1)
    enddo
    
    do k=1,pverp
      zi(k) = zi_in(pverp-k+1)
    enddo
    
    do k=1,pver_five
      zm_five(k) = zm_five_in(pver_five-k+1)
      rho_zs(k) = rho_five_in(pver_five-k+1)
    enddo  
    
    zi1 = pver-five_bot_k+1
    zi2 = pver-five_top_k+1
    
    dz=zi(2)
    
    ! define adz and adzw
    do k=2,pver
      adzw(k) = (zm(k)-zm(k-1))/dz
    enddo
    adzw(1) = 1._r8
    
    adz(1) = 1._r8
    do k=2,pver-1
      adz(k) = 0.5*(zm(k+1)-zm(k-1))/dz
    enddo
    adz(pver) = adzw(pver)
    
    ! Prepare coefficients

    ! If the location of the mid-layer point is optionally specified then following variables
    ! are required to be computed every time this function is called.
			
    ! Distance from the mid-layer level to the host model lower/upper interface value divided
    ! by dz. Here the mid-layer level is z(k).
    
    do k=1,pver
      adz_dn(k) = (zm(k) - zi(k))/dz
    enddo
    
    do k=1,pver
      adz_up(k) = (zi(k+1) - zm(k))/dz
    enddo
    
    ! For solving system of equations
    ! The j-th column of the matrix A is stored in the j-th column of the array a as follows:
    !    a(ml+mu+1+i-j,j) = A(i,j) for max(1,j-mu)<=i<=min(nzm,j+ml)
    ! Set up superdiagonal, diagonal, subdiagonal component of the matrix
    a(:,:) = 0._r8
    ! Superdiagonal 
    do k=2,pver
      a(2,k)=adz_up(k-1)**2/(adz(k)*adzw(k-1))*0.5_r8
    enddo
    ! Diagonal
    k = 1
    a(3,1) = (adz_dn(1) + adzw(1) + adz_up(1) * adz_dn(2) / adz(2))/adzw(1) * 0.5_r8
    do k=2,pver
      kp1=min(k+1,pver)
      a(3,k) = (adz_up(k-1) * adz_dn(k) / adz(k) + adzw(k) & 
               + adz_up(k) * adz_dn(k+1) / adz(kp1) ) / adzw(k) * 0.5
    enddo
    ! Subdiagonal
    do k=1, pver
      a(4,k) = adz_dn(k)**2 / (adz(k) * adzw(k))*0.5
    enddo
    
    ! Factor the matrix with LAPACK, BLAS
    call dgbtrf( pver, pver, ml, mu, a, lda, ipiv, info )    
    
    ! For interpolation
    i5=0
    do k=1,zi1-1
      i5 = i5 + 1
    enddo
    
    i5zi1 = i5
    do k = zi1, zi2
      i1 = i5 + 1
      i2 = i1 + five_add_nlevels / 2 - 1 
      i3 = i1 + five_add_nlevels / 2
      i4 = i1 + five_add_nlevels / 2 + 1
      i5 = i1 + five_add_nlevels 
      ! weight for linear interpolation
      weight(i1:i2) = ( zm_five(i1:i2) - zi(k) ) / ( zm(k) - zi(k) )
      weight(i3) = 1.0_r8
      weight(i4:i5) = ( zm_five(i4:i5) - zi(k+1) ) / ( zm(k) - zi(k+1) ) 
      ! c1, c2, c3
      c0(k) = 2.0 * (zi(k+1) - zi(k))
      c1(k) = zi(k+1) - zi(k)
      c2(k) = zm_five(i3) - zi(k)
      c3(k) = zi(k+1) - zm_five(i3)
    enddo
    ic1(:) = 1.0_r8 / c1(:)
    ic2(:) = 1.0_r8 / c2(:) 
    ic3(:) = 1.0_r8 / c3(:)
    
    ! add flag computed?
    
    ! Mass weight inout value
    rdf_z(:) = rho_low(:) * ten_low(:)
    
    ! Solve system of equations to get mid-layer target value
    ! Solve the linear system with LAPACK, BLAS
    ! input b wil be solution when output
    b(:) = rdf_z(:)
    call dgbtrs('n', pver, ml, mu, 1, a, lda, ipiv, b, pver, info)
    rdf_zs_ml(:) = b(:)
    
    ! Interface target value
    rdf_zm(1) = rdf_zs_ml(1)
    do k = 2, pver-1
      rdf_zm(k) = adz_dn(k) / adz(k) * rdf_zs_ml(k-1) & 
        + adz_up(k-1) / adz(k) * rdf_zs_ml(k)
    enddo
    rdf_zm(pver) = 0.0 !domain top tendency
    
    ! Detection and correction of grid-scale violation for df_zm 
    !  Zerroukat et al. (2005 QJRMS)
    spurious(:) = .false.
    do k = 1, pver
      km3 = MAX( k-3, 1 )
      km2 = MAX( k-2, 1 )
      km1 = MAX( k-1, 1 )
      k00 = MIN( k, pver )
      kp1 = MIN( k+1, pver )
      kp2 = MIN( k+2, pver )
      cnd1 = ( rdf_zm(k) - rdf_z(km1) ) * ( rdf_z(k00) - rdf_zm(k) ) < 0.0
      cnd2 = ( rdf_z(km1) - rdf_z(km2) ) * ( rdf_z(kp1) - rdf_z(k00) ) >= 0.0
      cnd3 = ( rdf_z(km1) - rdf_z(km2) ) * ( rdf_z(km2) - rdf_z(km3) ) <= 0.0
      cnd4 = ( rdf_z(kp1) - rdf_z(k00) ) * ( rdf_z(kp2) - rdf_z(kp1) ) <= 0.0
      cnd5 = ( rdf_zm(k) - rdf_z(km1) ) * ( rdf_z(km1) - rdf_z(km2) ) <= 0.0
      if ( cnd1 .and. ( cnd2 .or. cnd3 .or. cnd4 .or. cnd5 ) ) then
        spurious(k) = .false.
	alpha = ABS( rdf_zm(k) - rdf_z(k00) ) - ABS( rdf_zm(k) - rdf_z(km1) )
	alpha = SIGN( 0.5_r8, alpha ) + 0.5_r8 ! Heaviside step function, alpha = 0 or 1
        rdf_zm(k) = alpha * rdf_z(km1) + ( 1.0_r8 - alpha ) * rdf_z(k00)
      endif
    enddo
    
    ! Store rdf_zm into rdf_zm_up and rdf_zm_dn
    rdf_zm_up(:) = rdf_zm(1:pver)
    rdf_zm_dn(:) = rdf_zm(1:pver)  
    
    ! Detection and correction of grid-scale violation for rdf_zs_ml for monotonic layer
    ! - Recompute rdf_zs_zl with updated rdf_zm
    ! - For monotonic layer, check if it is bounded between rdf_zm(k) and rdf_zm(k+1)
    !   If not bounded, assign rdf_zs_ml to closest rdf_zm, and compute other interface value
    !   and store it into either rdf_zm_up or rdf_zm_dn. That level will have two interface
    !   values for upper and lower layer.    
    spurious(:) = .FALSE.
    do k = 1, pver
      km1 = MAX( k-1, 1 )
      kp1 = MIN( k+1, pver )
      rdf_zs_ml(k) = ic1(k) * ( c0(k)*rdf_z(k) - c2(k)*rdf_zm(k) - c3(k)*rdf_zm(k+1) )
      cnd1 = ( rdf_z(k) - rdf_z(km1) ) * ( rdf_z(kp1) - rdf_z(k) ) >= 0.0_r8
      cnd2 = ( rdf_zs_ml(k) - rdf_zm(k) ) * ( rdf_zm(k+1) - rdf_zs_ml(k) ) < 0.0_r8
      if ( cnd1 .AND. cnd2 ) then
        ! Inflection within a monotonic layer
	spurious(k) = .TRUE.
	alpha = ABS( rdf_zs_ml(k) - rdf_zm(k) ) - ABS( rdf_zs_ml(k) - rdf_zm(k+1) )
	alpha = SIGN( 0.5_r8, alpha ) + 0.5_r8 ! alpha = 0 or 1
	rdf_zs_ml(k) = alpha * rdf_zm(k+1) + ( 1.0_r8 - alpha ) * rdf_zm(k)
	if ( alpha < 0.5_r8 ) then
	  rdf_zm_up(k) = ic3(k) * ( c0(k)*rdf_z(k) - c1(k)*rdf_zs_ml(k) - c2(k)*rdf_zm(k) )
	else
	  rdf_zm_dn(k) = ic2(k) * (c0(k)*rdf_z(k) - c1(k)*rdf_zs_ml(k) - c3(k)*rdf_zm(k+1))
	endif
      endif
    enddo
    
    ! Remove discountinuity at interface level as many as possible
    ! - For monotonic layer, set rdf_z_dn(k) = rdf_z_up(k-1) and rdf_z_up(k) = rdf_z_dn(k+1),
    !   then re-compute rdf_zs_ml.
    ! - Check if new rdf_zs_ml is bounded between rdf_zm_dn and rdf_zm_up, and if not, set
    !   rdf_zs_ml to the closer interface value and compute the other interfaec value.              
    do k = 1, pver
      km1 = MAX( k-1, 1 )
      kp1 = MIN( k+1, pver )
      cnd1 = ( rdf_zs_ml(k) - rdf_zm_dn(k) ) * ( rdf_zm_up(k) - rdf_zs_ml(k) ) >= 0.0_r8
      if ( cnd1 ) then
        ! Monotonic layer
					
	! Re-set df_zm_dn and df_zm_up
	rdf_zm_dn(k) = rdf_zm_up(km1)
	rdf_zm_up(k) = rdf_zm_dn(kp1)
	! Re-compute df_zs
	rdf_zs_ml(k) = ic1(k) * ( c0(k)*rdf_z(k) - c2(k)*rdf_zm_dn(k) - c3(k)*rdf_zm_up(k) )
	!
	cnd2 = ( rdf_zs_ml(k) - rdf_zm_dn(k) ) * ( rdf_zm_up(k) - rdf_zs_ml(k) ) < 0.0_r8

	if ( cnd2 ) then
	  ! Non-monotonic profile
	  spurious(k) = .TRUE.
	  alpha = ABS( rdf_zs_ml(k) - rdf_zm_dn(k) ) - ABS( rdf_zs_ml(k) - rdf_zm_up(k) )
	  alpha = SIGN( 0.5_r8, alpha ) + 0.5_r8
	  rdf_zs_ml(k) = alpha * rdf_zm_up(k) + ( 1.0_r8 - alpha ) * rdf_zm_dn(k)
	  if ( alpha < 0.5_r8 ) then
	    rdf_zm_up(k) = ic3(k) * (c0(k)*rdf_z(k)-c1(k)*rdf_zs_ml(k)-c2(k)*rdf_zm_dn(k))
	  else
	    rdf_zm_dn(k) = ic2(k) * (c0(k)*rdf_z(k)-c1(k)*rdf_zs_ml(k)-c3(k)*rdf_zm_up(k))
	  endif
	endif
	
      endif
    enddo
    
    ! Construct the tendency profile
    i5 = 0
    ! Below zi1
    do k=1,zi1-1
      i5=i5+1
      df_zs(i5) = df_z(k)
    enddo

    ! between zi1 and zi2
    do k = zi1, zi2
      i1 = i5 + 1
      i2 = i1 + five_add_nlevels / 2 - 1 ! i2 = i3-1 or i2 = i1
      i3 = i1 + five_add_nlevels / 2     ! zs(i3) = z(k)
      i4 = i1 + five_add_nlevels / 2 + 1 ! i4 = i3+1 or i4 = i5
      i5 = i1 + five_add_nlevels
      ! Compute df_zs for all zs levels
      ! zi(k) < zs(i1:i2) < z(k)
      df_zs(i1:i2) = ( 1.0_r8 - weight(i1:i2) ) * rdf_zm_dn(k) + weight(i1:i2) * rdf_zs_ml(k)
      ! zs(i3)
      df_zs(i3) = rdf_zs_ml(k)
      ! z(k) < zs(i4:i5) < zi(k+1)
      df_zs(i4:i5) = ( 1.0_r8 - weight(i4:i5) ) * rdf_zm_up(k) + weight(i4:i5) * rdf_zs_ml(k)
      ! Non-mass weighted
      df_zs(i1:i5) = df_zs(i1:i5) / rho_zs(i1:i5)
    enddo
    
    ! Above zi2
    do k = zi2+1,pver
      i5 = i5 + 1
      df_zs(i5) = df_z(k)
    enddo    
    
    ! Finally, "unflip" the tendency
    do k = 1,pver_five
      ten_high(k) = df_zs(pver_five-k+1)
    enddo
  
  end subroutine tendency_low_to_high  

end module five_intr
