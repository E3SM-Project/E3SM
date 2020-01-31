module five_intr

!-------------------------------------------
! Module to interface FIVE with E3SM physics packages
!

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,           only : pcols, pver, pverp, begchunk, endchunk
  use physics_types, only: physics_state
  use physics_buffer, only: physics_buffer_desc
  use constituents, only: pcnst
  use physconst, only: rair, gravit, cpair, zvir
  use cam_logfile, only: iulog

  implicit none
! HHLEE
  public:: five_readnl ! read namelist from file 

  ! This is the number of layers to add between E3SM levels
  !  NOTE: This must be an EVEN number, due to limitations
  !  in the tendency interpolation scheme
  integer :: five_add_nlevels_0
  integer :: five_add_nlevels_1

  ! Determine which layers you want to add levels to.  NOTE:
  !  this refers to the base or reference pressure profiles
  !  (as pressure levels can change from time step to time step
  !  in an E3SM simulation).
  
  ! The bottom layer to which we will add layers to (set
  !   to a very large value to add all the way to surface, though
  !   currently setting this value to surface results in model
  !   crashes in SCM, so need to investigate)
  
  integer  :: five_num_patches
  
  real(r8) :: five_bot_toadd_0 
  real(r8) :: five_bot_toadd_1 
  
  ! The top layer to which we will add layers to
  real(r8) :: five_top_toadd_0 
  real(r8) :: five_top_toadd_1
 
  ! Number of five levels.  These are computed on the fly
  ! using five_bot_toadd_0, five_top_toadd_0, five_add_nlevels_0
  ! using five_bot_toadd_1, five_top_toadd_1, five_add_nlevels_1

  integer :: pver_five, pverp_five

  ! Pressure values as initialized in hycoef.F90
  real(r8), allocatable :: alev_five(:) ! midpoint FIVE pressures (pascals)
  real(r8), allocatable :: ailev_five(:)! interface FIVE pressures (pascals)
  real(r8), allocatable :: hyai_five(:)
  real(r8), allocatable :: hyam_five(:)
  real(r8), allocatable :: hybi_five(:)
  real(r8), allocatable :: hybm_five(:)  
  

  real(r8), parameter :: ps0 = 1.0e5_r8
  
  ! define physics buffer indicies here for the FIVE
  !  variables added to the PBUF
  integer :: t_five_idx       ,&
             q_five_idx       ,&
	     u_five_idx       ,&
	     v_five_idx       ,&
	     pmid_five_idx    ,&
	     pint_five_idx    

 ! add sharing FIVE variables
  integer :: cld_five_idx     ,&
             concld_five_idx  ,&
             ast_five_idx     ,&
             dei_five_idx     ,&
             des_five_idx     ,&
             mu_five_idx      ,&
             lambdac_five_idx ,&
             iciwp_five_idx   ,&
             iclwp_five_idx   ,&
             icswp_five_idx   ,&
             cldfsnow_five_idx,&
	     wp2_five_idx     ,&
	     radf_idx         ,&
	     wp2_idx          ,&
	     wp3_idx          ,&
	     wpthlp_idx       ,&
	     wprtp_idx        ,&
	     rtpthlp_idx      ,&
	     rtp2_idx         ,&
	     thlp2_idx        ,&
	     up2_idx          ,&
	     vp2_idx          ,&
	     upwp_idx         ,&
	     vpwp_idx 
  

     
  integer :: five_bot_k_0, five_top_k_0  ! indicees where
                                         ! levels are added

  integer :: five_bot_k_1, five_top_k_1  ! indicees where
                                         ! levels are added

  contains 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !    
  
  subroutine five_readnl(nlfile)
    
    use spmd_utils,    only: masterproc  
    use units,         only: getunit, freeunit
    use namelist_utils,  only: find_group_name
    use cam_abortutils,         only: endrun
    use mpishorthand    
 
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
    integer :: iunit, read_status
    
    namelist /five_nl/ five_num_patches, five_add_nlevels_0, five_add_nlevels_1, five_bot_toadd_0, five_bot_toadd_1, five_top_toadd_0, five_top_toadd_1
    
    !----- Begin Code -----
     write(*,*)"Entering function: five_readnl"
    !  Initialize some FIVE parameters
    !five_add_nlevels = 1 
    
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
    write (*,*)"five_num_patches = ",five_num_patches
#ifdef SPMD
! Broadcast namelist variables
      call mpibcast(five_num_patches  ,            1,   mpir8,    0, mpicom)
      call mpibcast(five_add_nlevels_0,            1,   mpiint,   0, mpicom)

      call mpibcast(five_bot_toadd_0  ,            1,   mpir8,    0, mpicom)
      call mpibcast(five_top_toadd_0  ,            1,   mpir8,    0, mpicom)

      call mpibcast(five_add_nlevels_1,            1,   mpiint,   0, mpicom)
      call mpibcast(five_bot_toadd_1  ,            1,   mpir8,    0, mpicom)
      call mpibcast(five_top_toadd_1  ,            1,   mpir8,    0, mpicom)
#endif 

    ! Now, based on five_add_nlevels_0, five_bot_toadd_0, and five_top_toadd_0
    ! Now, based on five_add_nlevels_1, five_bot_toadd_1, and five_top_toadd_1
    !   we will determine the number of FIVE levels we will need
    !   for our simulation.            
  write(*,*)"Leaving function: five_readnl"
  write(*,*)"  "
  end subroutine five_readnl
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine five_register_e3sm
  
    ! Register FIVE fields to the physics buffer
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols
    
      write(*,*)"Entering function: five_register_e3sm"
    ! Define PBUF for prognostics
    call pbuf_add_field('T_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls      /), t_five_idx)
    call pbuf_add_field('Q_FIVE',       'global', dtype_r8, (/pcols,pver_five,pcnst,dyn_time_lvls/), q_five_idx) 
    call pbuf_add_field('U_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls      /), u_five_idx)
    call pbuf_add_field('V_FIVE',       'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls      /), v_five_idx)
    
    ! Define PBUF for non-prognostic variables
    ! Probably also need to save things like cloud fraction
    call pbuf_add_field('PMID_FIVE', 'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), pmid_five_idx)
    call pbuf_add_field('PINT_FIVE', 'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), pint_five_idx) 
    
    ! sharing FIVE variables
    call pbuf_add_field('CLD_FIVE',     'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), cld_five_idx) 
    call pbuf_add_field('CONCLD_FIVE',  'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), concld_five_idx) 
    call pbuf_add_field('AST_FIVE',     'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), ast_five_idx) 
    call pbuf_add_field('DEI_FIVE',     'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), dei_five_idx) 
    call pbuf_add_field('DES_FIVE',     'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), des_five_idx) 
    call pbuf_add_field('MU_FIVE',      'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), mu_five_idx) 
    call pbuf_add_field('LAMBDAC_FIVE', 'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), lambdac_five_idx) 
    call pbuf_add_field('ICIWP_FIVE',   'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), iciwp_five_idx) 
    call pbuf_add_field('ICLWP_FIVE',   'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), iclwp_five_idx) 
    call pbuf_add_field('ICSWP_FIVE',   'global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), icswp_five_idx) 
    call pbuf_add_field('CLDFSNOW_FIVE','global', dtype_r8, (/pcols,pver_five,dyn_time_lvls/), cldfsnow_five_idx) 
    
    call pbuf_add_field('RAD_CLUBB',   'global', dtype_r8, (/pcols,pver_five               /), radf_idx)    
    call pbuf_add_field('WP2_nadv',    'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), wp2_idx)
    call pbuf_add_field('WP3_nadv',    'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), wp3_idx)
    call pbuf_add_field('WPTHLP_nadv', 'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), wpthlp_idx)
    call pbuf_add_field('WPRTP_nadv',  'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), wprtp_idx)
    call pbuf_add_field('RTPTHLP_nadv','global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), rtpthlp_idx)
    call pbuf_add_field('RTP2_nadv',   'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), rtp2_idx)
    call pbuf_add_field('THLP2_nadv',  'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), thlp2_idx)
    call pbuf_add_field('UP2_nadv',    'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), up2_idx)
    call pbuf_add_field('VP2_nadv',    'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), vp2_idx)    
    call pbuf_add_field('UPWP',        'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), upwp_idx)
    call pbuf_add_field('VPWP',        'global', dtype_r8, (/pcols,pverp_five,dyn_time_lvls/), vpwp_idx)    
  
  write(*,*)"Leaving function: five_register_e3sm"
  write(*,*)"  "
  end subroutine five_register_e3sm 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !

  subroutine init_five_heights(ailev, &
                               alev,  &
                               hyai,  &
                               hyam,  &
                               hybi,  &
                               hybm)
  
    ! Purpose is to initialize FIVE heights.  This is called 
    !   dyn_comp.F90.  We need to do this here to initialize
    !   height coordinates for output purposes, so output
    !   can be written to CAM history tapes.  
    !   NOTE: that the pressure values determined here
    !   will not necessarily be the values used in the 
    !   simulation.  This is because pressure levels can 
    !   vary from timestep to timestep based on surface pressure.
    !   Here we will find the hybrid points
    !   and the actual pressure values will be computed in
    !   init_five_profiles.
    
    use cam_history_support, only: add_vert_coord
    
    ! Input variables
    real(r8), intent(in) :: alev(pver)   ! midpoint pressure from 
                                         ! E3SM host model
    real(r8), intent(in) :: ailev(pverp) ! interface pressure from
                                         ! E3SM host model
    real(r8), intent(in) :: hyai(pverp)  ! hybrid points from host
    real(r8), intent(in) :: hyam(pver)
    real(r8), intent(in) :: hybi(pverp)
    real(r8), intent(in) :: hybm(pver)
    
    ! Local variables
     real(r8),allocatable :: patchPressureLo(:)
     real(r8),allocatable :: patchPressureHi(:)  
     
    integer,allocatable :: patchIndexLo(:)
    integer,allocatable :: patchIndexHi(:)    
    
    integer ,allocatable :: fiveAddNLevels(:)  
    integer :: iPatch 
    
    integer :: low_ind, high_ind

    integer :: kh, ki, k, i
    real(r8) :: incr, incr2, incr3
    write(*,*)"Entering function: init_five_heights"  
    write(*,*)"five_num_patches = ",five_num_patches
  
    allocate(patchPressureLo (five_num_patches ))
    allocate(patchPressureHi (five_num_patches ))

    allocate(patchIndexLo    (five_num_patches ))
    allocate(patchIndexHi    (five_num_patches ))
      
    allocate(fiveAddNLevels  (five_num_patches ))
    
    ! put five_bot_toadd_0 and five_bot_toadd_1 into patchPressureLo
    patchPressureLo(1) = five_bot_toadd_0
!    patchPressureLo(2) = five_bot_toadd_1

    ! put five_top_toadd_0 and five_top_toadd_1 into patchPressureHi
    patchPressureHi(1) = five_top_toadd_0
 !   patchPressureHi(2) = five_top_toadd_1

    !put  five_add_nlevels_0  and  five_add_nlevels_1 into fiveAddNLevels 
    fiveAddNLevels(1) = five_add_nlevels_0 
  !  fiveAddNLevels(2) = five_add_nlevels_1

    write(*,*)  "patchPressureLo(1) =",patchPressureLo(1)
    write(*,*)  "patchPressureHi(1) =",patchPressureHi(1)
    write(*,*) ""
   ! write(*,*)  "patchPressureLo(2) =",patchPressureLo(2)
   ! write(*,*)  "patchPressureHi(2) =",patchPressureHi(2)
    write(*,*) ""
    write(*,*)  "fiveAddNLevels(1)  =",fiveAddNLevels(1)
    !write(*,*)  "fiveAddNLevels(2)  =",fiveAddNLevels(2)

    
    ! First we need to determine how many pver_five and
    !   and pverp_five levels there are.  This may seem repetitive 
    !   with the code below it but is needed to appease E3SM 
    !   order of operations and to allow us to add levels
    !   on the fly

   
    kh = 0
    !calculate pverp_five
      do k=1,pverp 
        kh=kh+1
        
        !if (ailev(k) .le. five_bot_toadd .and. ailev(k) .ge. five_top_toadd) then
       do iPatch = 1, five_num_patches
         if ((k .le. pverp) .and. ailev(k) .le. patchPressureLo(iPatch) .and. ailev(k) .ge. patchPressureHi(iPatch)) then
	   kh = kh + fiveAddNLevels(iPatch) 
	 endif
       enddo 
      enddo
    
    ! pverp_five = pverp + sum over patches( fiveAddNLevels(iPatch)*number tagged(iPatch))
    pverp_five = kh
    pver_five  = pverp_five-1

    Write(*,*) "pver        = ", pver
    Write(*,*) "pverp       = ", pverp
    Write(*,*) "kh          = ", kh
    Write(*,*) "pver_five   = ", pver_five
    Write(*,*) "pverp_five  = ", pverp_five
    
    allocate(alev_five (pver_five )) ! midpoint FIVE pressures (pascals)
    allocate(ailev_five(pverp_five)) ! interface FIVE pressures (pascals)
    allocate(hyai_five (pverp_five))
    allocate(hyam_five (pver_five ))
    allocate(hybi_five (pverp_five))
    allocate(hybm_five (pver_five ))               
    
    ! Note that here we are actually adding "layers" to the
    !   hybrid coordinates.  The reason is that the 
    !   pressure levels of E3SM can change from time step
    !   to timestep based on the surface pressure.  
    !   The hybrid coordinates are the only thing related to
    !   the height profile that remains constant. 

     !kh is an index of the five grid
     kh       = 1 
     
     !initialize integer counters for patch boundaries
     do iPatch = 1,five_num_patches   
       
      ! these integers record the bottom and top of patches
      patchIndexLo(iPatch) = pver
      patchIndexHi(iPatch) = 1  
     enddo

     do k=1,pverp ! k is layer index on E3SM grid

        ! Copy preexisting points to FIVE grid
        hyai_five(kh) = hyai(k)
        hybi_five(kh) = hybi(k)

        write(*,*)"hyai_five(kh) hybi_five(kh) kh,ailev(k) = ",hyai_five(kh),hybi_five(kh),kh,ailev(k)
        
      
        kh=kh+1
       
        !loop over patches 
      do iPatch = 1,five_num_patches   
       
	if ((k .le. pverp) .and. ailev(k) .le.  patchPressureLo(iPatch) .and. ailev(k) .ge.  patchPressureHi(iPatch)) then

     	 ! record the index 
   	 if (patchIndexHi(iPatch) .eq. 1) patchIndexHi(iPatch)= k
	 patchIndexLo(iPatch) = k
	  
	 ! If we are inside the portion of grid we want to add layers
	 ! to then compute the pressure increment (resolution)
	 incr2=(hyai(k+1)-hyai(k))/(fiveAddNlevels(iPatch)+1)
	 incr3=(hybi(k+1)-hybi(k))/(fiveAddNlevels(iPatch)+1)
	    
	 ! Define new data
	 do ki=1,fiveAddNlevels(iPatch)
	  hyai_five(kh) = hyai_five(kh-1)+incr2
	  hybi_five(kh) = hybi_five(kh-1)+incr3
         
          write(*,*)"hyai_five(kh) hybi_five(kh) kh,ailev(k) = **",hyai_five(kh),hybi_five(kh),kh,ailev(k)
	  kh=kh+1
	 enddo
	    
        endif
         
     ! write patch info
     write(iulog,*) 'iPatch = :',iPatch
     write(iulog,*) 'Number of FIVE levels in this patch: ', kh-1    
    
     ! write bottom and top indices
     write(iulog,*) 'Index of bottom of lower layer of this patch to add FIVE levels to:',patchIndexLo(iPatch)
     write(iulog,*) 'Index of top    of lower layer of this patch to add FIVE levels to:',patchIndexHi(iPatch)
    
      enddo
     enddo

     !temporary work-around until global variables are arrays
     five_bot_k_0 = patchIndexLo(1)
     five_top_k_0 = patchIndexHi(1)

!     five_bot_k_1 = patchIndexLo(2)
 !    five_top_k_1 = patchIndexHi(2)
         
     write(*,*)"five_bot_k_0 = ",five_bot_k_0
     write(*,*)"five_top_k_0 = ",five_top_k_0
     write(*,*) ""
     write(*,*)"five_bot_k_1 = ",five_bot_k_1
     write(*,*)"five_top_k_1 = ",five_top_k_1

    !  the reference  = base pressure: ps0 = 1.0e5_r8
    ailev_five(:pverp_five) = 0.01_r8*ps0*(hyai_five(:pverp_five) + hybi_five(:pverp_five)) 
    
    ! Define five_mid layers
    do k=1,pver_five
      alev_five(k) = (ailev_five(k)+ailev_five(k+1))/2.0_r8
    enddo      
	
    ! Add vertical coordinate to the history output field
    call add_vert_coord('lev_five', pver_five,  'FIVE hybrid level at midpoints (1000*(A+B))' , 'hPa', alev_five, positive='down',&
                         standard_name='atmosphere_hybrid_sigma_pressure_coordinate')

    call add_vert_coord('ilev_five', pverp_five,                               &
         'FIVE hybrid level at interfaces (1000*(A+B))', 'hPa', ailev_five, &
         positive='down',                                                    &
         standard_name='atmosphere_hybrid_sigma_pressure_coordinate')	
        
   
  write(*,*)"Leaving function: init_five_heights"
  write(*,*)"  "  
  end subroutine init_five_heights 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !  
  
  subroutine init_five_profiles(phys_state, pbuf2d)
  
    ! Purpose is to initialize FIVE profiles.
    ! 1)  The interface and midlayers are defined for FIVE
    !     based on the coordinate values and reference and base
    !     pressures
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
    
    ! Add sharing FIVE variables
    real(r8) :: cld_five(pcols,pver_five)
    real(r8) :: concld_five(pcols,pver_five)
    real(r8) :: ast_five(pcols,pver_five)
    real(r8) :: dei_five(pcols,pver_five)
    real(r8) :: des_five(pcols,pver_five)
    real(r8) :: mu_five(pcols,pver_five)
    real(r8) :: lambdac_five(pcols,pver_five)
    real(r8) :: iciwp_five(pcols,pver_five)
    real(r8) :: iclwp_five(pcols,pver_five)
    real(r8) :: icswp_five(pcols,pver_five)
    real(r8) :: cldfsnow_five(pcols,pver_five)
    
    integer :: ncol, i, p, lchnk, k, kh, ki, n
    integer :: low_ind, high_ind
    
    real(r8) :: pint_host(pcols,pverp)
    real(r8) :: pmid_host(pcols,pver)
    real(r8) :: pint_five(pcols,pverp_five) 
    real(r8) :: pmid_five(pcols,pver_five)
    real(r8) :: incr   
    real(r8) :: psr
    
    ! Really most of this code below may not be needed
    !   At the first time FIVE is syncronized with E3SM
    !   the FIVE states should be naturally updated.
    !   All FIVE variables could be initialized to zero.  
    !   However, will keep for now.  Can experiment in the
    !   future with removing.  
    
    ! Loop over all the physics "chunks" in E3SM
    write(*,*)"Entering function: init_five_profiles"  
    do lchnk = begchunk,endchunk
    
      ncol = phys_state(lchnk)%ncol ! number of columns in this chunk
      pint_host = phys_state(lchnk)%pint ! E3SM interface pressures
      pmid_host = phys_state(lchnk)%pmid ! E3SM midpoint pressures
      pbuf2d_chunk => pbuf_get_chunk(pbuf2d,lchnk)
      
      ! Now define five interface layers
      !  Note that here the base and reference pressures could
      !  be different
      do i=1,ncol
	psr = pint_host(i,pverp)
	pint_five(i,:pverp_five) = ps0*hyai_five(:pverp_five) + psr*hybi_five(:pverp_five)
      enddo
      
      ! Now define five_mid layers
      do i=1,ncol
        do k = 1,pver_five
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
      
      ! Initialize sharing FIVE variables
      cld_five(:,:) = 0._r8
      concld_five(:,:) = 0._r8
      ast_five(:,:) = 0._r8
      dei_five(:,:) = 0._r8
      des_five(:,:) = 0._r8
      mu_five(:,:) = 0._r8
      lambdac_five(:,:) = 0._r8
      iciwp_five(:,:) = 0._r8
      iclwp_five(:,:) = 0._r8
      icswp_five(:,:) = 0._r8
      cldfsnow_five(:,:) = 0._r8

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
        ! new sharing FIVE avariables
        call pbuf_set_field(pbuf2d_chunk, cld_five_idx,      cld_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, concld_five_idx,   concld_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, ast_five_idx,      ast_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, dei_five_idx,      dei_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, des_five_idx,      des_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, mu_five_idx,       mu_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, lambdac_five_idx,  lambdac_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, iciwp_five_idx,    iciwp_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, iclwp_five_idx,    iclwp_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, icswp_five_idx,    icswp_five,(/1,1,n/),(/pcols,pver_five,1/))
        call pbuf_set_field(pbuf2d_chunk, cldfsnow_five_idx, cldfsnow_five,(/1,1,n/),(/pcols,pver_five,1/))
      enddo
    
    enddo
  
  write(*,*)"Leaving function: init_five_profiles" 
  write(*,*)"  " 
  end subroutine init_five_profiles
  
  ! ======================================== !
  !                                          !
  ! ======================================== !  
  
  subroutine five_syncronize_e3sm(&
                state, dtime, p0, &
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
    real(r8), intent(in) :: p0 ! reference pressure
    real(r8), intent(inout) :: pmid_five(pcols,pver_five)
    real(r8), intent(inout) :: pint_five(pcols,pverp_five)
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
    real(r8) :: zm_five(pcols,pver_five)
    real(r8) :: zi_five(pcols,pverp_five)
    
    real(r8) :: psr
! HHLEE 20190117
    real(r8) :: five_water, e3sm_water, ratio, diff
    write(*,*)"Entering function: five_syncronize_e3sm"  
    ncol = state%ncol
    
    ! First update the FIVE pressure levels using
    !   the hybrid coordinates and the host model 
    !   surface pressure
    do i=1,ncol
      psr = state%pint(i,pverp)
      pint_five(i,:pverp_five) = ps0*hyai_five(:pverp_five) + psr*hybi_five(:pverp_five)
    enddo
      
    ! Now define five_mid layers
    do i=1,ncol
      do k=1,pver_five
        pmid_five(i,k) = (pint_five(i,k)+pint_five(i,k+1))/2.0_r8
      enddo
    enddo 
  
    ! Compute pressure differences on FIVE grid  
    do k = 1,pver_five
      do i=1,ncol
        pdel_five(i,k) = pint_five(i,k+1)-pint_five(i,k)
      enddo
    enddo
    
    ! Compute grid stuff needed for vertical averaging and tendencies
    do i=1,ncol
      write(*,*)"First Call to compute_five_grids"   
      call compute_five_grids(state%t(i,:),state%pdel(i,:),state%pmid(i,:),pver,&
             dz_e3sm(i,:),rho_e3sm(i,:))
      write(*,*)"Second Call to compute_five_grids"   
      call compute_five_grids(t_five(i,:),pdel_five(i,:),pmid_five(i,:),pver_five,&
             dz_five(i,:),rho_five(i,:))
    enddo
    
    ! Compute height arrays needed for mass weighted averaging and interpolation
    do i=1,ncol
      call compute_five_heights(pmid_five(i,:),pint_five(i,:),t_five(i,:),&
             q_five(i,:,1),q_five(i,:,2),pdel_five(i,:),pver_five,p0,&
	     zm_five(i,:),zi_five(i,:))
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

    ! Now interpolate this tendency to the higher resolution FIVE grid, 
    !   using the interpolation method of Sheng and Zwiers (1998), 
    !   as documented in Yamaguchi et al. (2017) Appendix B
    do i=1,ncol
    
      call tendency_low_to_high(state%zm(i,:),state%zi(i,:),zm_five(i,:),&
             rho_e3sm(i,:),rho_five(i,:),t_five_tend_low(i,:),t_five_tend(i,:)) 
	     
      call tendency_low_to_high(state%zm(i,:),state%zi(i,:),zm_five(i,:),&
             rho_e3sm(i,:),rho_five(i,:),u_five_tend_low(i,:),u_five_tend(i,:))
	     
      call tendency_low_to_high(state%zm(i,:),state%zi(i,:),zm_five(i,:),&
             rho_e3sm(i,:),rho_five(i,:),v_five_tend_low(i,:),v_five_tend(i,:))
	     
      do p=1,pcnst
        call tendency_low_to_high(state%zm(i,:),state%zi(i,:),zm_five(i,:),&
             rho_e3sm(i,:),rho_five(i,:),q_five_tend_low(i,:,p),q_five_tend(i,:,p))     
      enddo	     	             	
      
    enddo      	     

    ! Finally, update FIVE prognostic variables based on this tendency, so 
    !   complete syncronization with E3SM
    do k=1,pver_five
      do i=1,ncol
        t_five(i,k) = t_five(i,k) + dtime * t_five_tend(i,k)
        u_five(i,k) = u_five(i,k) + dtime * u_five_tend(i,k)
        v_five(i,k) = v_five(i,k) + dtime * v_five_tend(i,k)
       
        do p=1,pcnst
          q_five(i,k,p) = q_five(i,k,p) + dtime * q_five_tend(i,k,p)
! HHLEE 20190117
          q_five(i,k,p) = max(q_five(i,k,p),0._r8)
        enddo
       
      enddo
    enddo    
! HHLEE 20190117 
! water adjustment. This adjustment can solve the conservation issue in clubb.
    do i = 1, ncol 
       do p = 1, pcnst
            five_water = 0._r8
            e3sm_water = 0._r8 
            five_water = sum(q_five(i,:,p)*pdel_five(i,:)/gravit)
            e3sm_water = sum(state%q(i,:,p)*state%pdel(i,:)/gravit)

            if (five_water .gt. 0._r8) then
               ratio = e3sm_water/five_water
               q_five(i,:,p) = ratio * q_five(i,:,p)
            end if

            five_water = sum(q_five(i,:,p)*pdel_five(i,:)/gravit)
            if (five_water .ne. e3sm_water .and. five_water .gt. 0._r8) then
               diff = five_water - e3sm_water
               do k = 1, pver_five
                  ratio = (q_five(i,k,p)*pdel_five(i,k)/gravit)/five_water
                  q_five(i,k,p) = q_five(i,k,p) - diff*ratio
               end do
            end if

       end do
    end do

  

   write(*,*)"Leaving function: five_syncronize_e3sm"
   write(*,*)"  "  
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
     write(*,*)"Entering function: compute_five_grids"  	
    ! Compute density

    write(*,*)"nlev = ",nlev
    do k = 1,nlev
     write(*,*)"temp(k) pmid(k) k = ",temp(k),pmid(k),k
         
     rho(k) = pmid(k)/(rair*temp(k))
    enddo          
	      
    ! Compute dz	     
    do k = 1,nlev
     write(*,*)"rho(k) pmid(k) k = ",rho(k),pmid(k),k
          
     dz(k) = pdel(k)/(rho(k)*gravit)	      
    enddo	      
	      
   

   write(*,*)"Leaving function: compute_five_grids"
   write(*,*)"  "  	
  return       
  end subroutine compute_five_grids	 
  
  ! ======================================== !
  !                                          !
  ! ======================================== ! 
  
  subroutine compute_five_heights(&
               pmid, pint, t, qv, ql, &
	       pdel, nlev, p0, &
	       zm, zi)
	       
    ! Compute heights for the FIVE grid
    implicit none
    
    ! Input variables
    real(r8), intent(in) :: pmid(nlev)	 
    real(r8), intent(in) :: pint(nlev+1)
    real(r8), intent(in) :: t(nlev)
    real(r8), intent(in) :: qv(nlev)
    real(r8), intent(in) :: ql(nlev)
    real(r8), intent(in) :: pdel(nlev)
    real(r8), intent(in) :: p0
    
    integer, intent(in) :: nlev 
    
    ! Output variables
    real(r8), intent(out) :: zm(nlev)
    real(r8), intent(out) :: zi(nlev+1)

    ! Local variables
    real(r8) :: exner(nlev), tv(nlev), thv(nlev)
    real(r8) :: hkl, hkk
    
    integer :: i, k
    write(*,*)"Entering function: compute_five_heights"  
    ! Compute the Exner values using pressure
    do k = 1,nlev
      exner(k) = 1._r8/(pmid(k)/p0)**(rair/cpair)
    enddo   
    
    ! Define virtual temperature
    do k = 1,pver_five
      tv(k) = t(k)*(1._r8+zvir*qv(k)-ql(k)) 
      thv(k) = exner(k)*tv(k)
    enddo 
    
    ! First compute the heights (in [m]) on the FIVE grid
    ! Do this in a consistent manner how it is done in geopotential.F90
      
    zi(nlev+1) = 0.0_r8 ! Surface always zero by definition
    do k = nlev, 1, -1
      
      hkl = pdel(k)/pmid(k)
      hkk = 0.5_r8 * hkl
	
      zm(k) = zi(k+1) + (rair/gravit)*tv(k)*hkk
      zi(k) = zi(k+1) + (rair/gravit)*tv(k)*hkl
      
    enddo 
  
  
    write(*,*)"Leaving function: compute_five_heights" 
    write(*,*)"  " 
  
    return
  end subroutine compute_five_heights
  
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
    write(*,*)"Entering function: masswgt_vert_avg"  
    ! Initialize host variable
    var_host(:) = 0._r8
    rho_host_avg(:) = 0._r8
    
    kh=1 ! Vertical index for FIVE 
    do k = 1,pver ! Vertical index for E3SM
       
     ! Check to see how many FIVE layers there are within
     !   an E3SM layer 
     do while ( pmid_high(kh) .lt. pint_host(k+1) .and. &
                pmid_high(kh) .gt. pint_host(k))
      
       var_host(k) = var_host(k) + rho_high(kh) * var_high(kh) * dz_high(kh)
	rho_host_avg(k) = rho_host_avg(k) + rho_high(kh) * &
	  dz_high(kh)
	  
	kh = kh + 1 ! increase high res model by one layer
	if (kh .gt. pver_five) goto 10
	
      end do ! end while loop for kh
10 continue
      
      ! Compute rho on host grid
      rho_host_avg(k) = rho_host_avg(k)/dz_host(k)
      
      var_host(k) = var_host(k)/(rho_host_avg(k)*dz_host(k))

    enddo
   
    write(*,*)"Leaving function: masswgt_vert_avg"  
    write(*,*)"  "
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
    write(*,*)"Entering function: linear_interp"  
    do k2 = 1,km2

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
  
  write(*,*)"Leaving function: linear_interp"
  write(*,*)"  "  
  end subroutine linear_interp  
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine tendency_low_to_high(&
               zm_in, zi_in, &
	       zm_five_in, &
	       rho_low_in, rho_five_in, &
	       ten_low,ten_high)

    ! Subroutine to compute the tendency from the low
    !   resolution to high resolution E3SM grid.  This code
    !   was adopted from the FIVE implementation into SAM,
    !   obtained by Peter Bogenschutz from Tak Yamaguchi.
    !   SAM code assumes that the surface is index 1, whereas
    !   E3SM assumes that the surface is index pver(p).  Thus,
    !   for simplicity, we need to "flip" the arrays.  
    ! This interpolation routine is based on that of 
    !   Sheng and Zwiers (1998) and implemented as 
    !   documented in Yamaguchi et al. (2017) appendix B.
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
    real(r8) :: zi(pverp) ! flipped version of zi_in
    real(r8) :: df_z(pver) 
    real(r8) :: rho_low(pver)
    real(r8) :: rho_zs(pver_five)
    real(r8) :: zm_five(pver_five) ! flipped
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
    real(r8), dimension(pverp) :: rdf_zm
    real(r8), dimension(pver) :: rdf_zm_dn, rdf_zm_up
    real(r8), dimension(pver) :: c0, c1, c2, c3, ic1, ic2, ic3
    real(r8) :: b(pver)
    
    logical, dimension(pver) :: spurious
    logical :: cnd1, cnd2, cnd3, cnd4, cnd5
    
    logical :: do_limit
    
    integer :: i, i1, i2, i3, i4, i5
    integer :: k, km3, km2, km1, k00
    integer :: kp1, kp2, ierr

    ! these variables are patch dependent
    integer :: zi1, zi2
   
    integer :: i5zi1
    integer :: iPatch, lastPatchIndexHi
    
    integer,allocatable :: patchIndexLo(:)
    integer,allocatable :: patchIndexHi(:)    
    integer ,allocatable :: fiveAddNLevels(:)  
    
    write(*,*)"Entering function: tendency_low_to_high"  
    ! Construct tendency profile from E3SM grid to FIVE grid
    
    allocate(patchIndexLo  (five_num_patches ))
    allocate(patchIndexHi  (five_num_patches ))
    allocate(fiveAddNLevels(five_num_patches ))
    

    !indices on the unrefined grid that mark the bottom and top of layers
    patchIndexLo(1) = five_bot_k_0
    !patchIndexLo(2) = five_bot_k_1

    patchIndexHi(1) = five_top_k_0
    !patchIndexHi(2) = five_top_k_1
 
    write(*,*)" patchIndexLo(1) = ", patchIndexLo(1)
    !write(*,*)" patchIndexLo(2) = ", patchIndexLo(2)

    write(*,*)" patchIndexHi(1) = ", patchIndexHi(1)       
    !write(*,*)" patchIndexHi(2) = ", patchIndexHi(2)
    
    !put  five_add_nlevels_0  and  five_add_nlevels_1 into fiveAddNLevels 
    fiveAddNLevels(1) = five_add_nlevels_0 
    !fiveAddNLevels(2) = five_add_nlevels_1

    ! The constructed tendency profile satisfies layer mean average
    
    ! If no levels are added via FIVE just return a copy !!!!!change this
    if (fiveAddNlevels(1) .eq. 0) then
      ten_high(1:pver) = ten_low(1:pver)
      return
    endif
    
    ! First let's "flip" things so lowest level index = 1
    do k = 1,pver
      zm(k) = zm_in(pver-k+1)
      rho_low(k) = rho_low_in(pver-k+1)
      df_z(k) = ten_low(pver-k+1)
    enddo
    
    do k = 1,pverp
      zi(k) = zi_in(pverp-k+1)
    enddo
    
    do k = 1,pver_five
      zm_five(k) = zm_five_in(pver_five-k+1)
      rho_zs(k) = rho_five_in(pver_five-k+1)
    enddo  
    
    dz=zi(2)
    
    ! define adz and adzw
    do k = 2,pver
      adzw(k) = (zm(k)-zm(k-1))/dz
    enddo
    adzw(1) = 1._r8
    
    adz(1) = 1._r8
    do k = 2,pver-1
      adz(k) = 0.5*(zm(k+1)-zm(k-1))/dz
    enddo
    adz(pver) = adzw(pver)
    
    ! Prepare coefficients

    ! If the location of the mid-layer point is optionally specified then following variables
    ! are required to be computed every time this function is called.
			
    ! Distance from the mid-layer level to the host model lower/upper interface value divided
    ! by dz. Here the mid-layer level is z(k).
    
    do k = 1,pver
      adz_dn(k) = (zm(k) - zi(k))/dz
    enddo
    
    do k = 1,pver
      adz_up(k) = (zi(k+1) - zm(k))/dz
    enddo
    
    ! For solving system of equations
    ! The j-th column of the matrix A is stored in the j-th column of the array a as follows:
    !    a(ml+mu+1+i-j,j) = A(i,j) for max(1,j-mu)<=i<=min(nzm,j+ml)
    ! Set up superdiagonal, diagonal, subdiagonal component of the matrix
    a(:,:) = 0._r8
    ! Superdiagonal 
    do k = 2,pver
      a(2,k)=adz_up(k-1)**2/(adz(k)*adzw(k-1))*0.5_r8
    enddo

    ! Diagonal
    a(3,1) = (adz_dn(1) + adzw(1) + adz_up(1) * adz_dn(2) / adz(2))/adzw(1) * 0.5_r8
    do k = 2,pver
      kp1=min(k+1,pver)
      a(3,k) = (adz_up(k-1) * adz_dn(k) / adz(k) + adzw(k) & 
               + adz_up(k) * adz_dn(kp1) / adz(kp1) ) / adzw(k) * 0.5
	       !+ adz_dn(kp1) used to be adz_dn(k)
    enddo
    ! Subdiagonal
    do k = 1, pver
      a(4,k) = adz_dn(k)**2 / (adz(k) * adzw(k))*0.5
    enddo
    
    ! Factor the matrix with LAPACK, BLAS
    call dgbtrf( pver, pver, ml, mu, a, lda, ipiv, info )    
    
    c0(:) = 0._r8
    c1(:) = 0._r8
    c2(:) = 0._r8
    c3(:) = 0._r8
    ic1(:) = 0._r8
    ic2(:) = 0._r8
    ic3(:) = 0._r8
    
    weight(:) = 0._r8
    
    !loop over patches
     do iPatch = 1,five_num_patches
     
     zi1 = pver - patchIndexLo(iPatch)+ 1
     zi2 = pver - patchIndexHi(iPatch)+ 1
  
     write(*,*)"zi1,iPatch = ",zi1,iPatch
     write(*,*)"zi2,iPatch = ",zi2,iPatch
     
     i5    = zi1 - 1
     i5zi1 = i5

     do k = zi1, zi2
      i1 = i5 + 1
      i2 = i1 + fiveAddNlevels(iPatch) / 2 - 1 
      i3 = i1 + fiveAddNlevels(iPatch) / 2
      i4 = i1 + fiveAddNlevels(iPatch) / 2 + 1
      i5 = i1 + fiveAddNlevels(iPatch) 

      ! weight for linear interpolation
      weight(i1:i2) = ( zm_five(i1:i2) - zi(k) ) / ( zm(k) - zi(k) )
      weight(i3   ) = 1.0_r8
      weight(i4:i5) = ( zm_five(i4:i5) - zi(k+1) ) / ( zm(k) - zi(k+1) ) 

      ! c1, c2, c3
      c0(k) = 2.0_r8 * (zi(k+1) - zi(k))
      c1(k) = zi(k+1) - zi(k)
      c2(k) = zm_five(i3) - zi(k)
      c3(k) = zi(k+1) - zm_five(i3)
      
      ! ic1, ic2, ic3
      ic1(k) = 1.0_r8 / c1(k)
      ic2(k) = 1.0_r8 / c2(k) 
      ic3(k) = 1.0_r8 / c3(k)            
    
    enddo
    
   enddo
   
    ! add flag computed?
    
    ! Mass weight inout value
    rdf_z(:) = rho_low(:) * df_z(:)
    
    ! Solve system of equations to get mid-layer target value
    ! Solve the linear system with LAPACK, BLAS
    ! input b wil be solution when output
    b(:) = rdf_z(:)
    call dgbtrs('n', pver, ml, mu, 1, a, lda, ipiv, b, pver, info)
    rdf_zs_ml(:) = b(:)
    
    ! Interface target value
    rdf_zm(1) = rdf_zs_ml(1)
    do k = 2, pver
      rdf_zm(k) = adz_dn(k) / adz(k) * rdf_zs_ml(k-1) + adz_up(k-1) / adz(k) * rdf_zs_ml(k)
    enddo

    rdf_zm(pverp) = 0.0_r8 !domain top tendency
    
    do_limit = .true.
    if (do_limit) then 
    
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
    rdf_zm_up(:) = rdf_zm(2:pverp)
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
    
      rdf_zs_ml(k) = ic1(k) * ( c0(k)*rdf_z(k) - c2(k)*rdf_zm(k) - c3(k)*rdf_zm(k+1) ) !+PAB change
      cnd1 = ( rdf_z(k) - rdf_z(km1) ) * ( rdf_z(kp1) - rdf_z(k) ) >= 0.0_r8
      cnd2 = ( rdf_zs_ml(k) - rdf_zm(k) ) * ( rdf_zm(k+1) - rdf_zs_ml(k) ) < 0.0_r8 !+PAB change
  
     if ( cnd1 .AND. cnd2 ) then
        ! Inflection within a monotonic layer
	spurious(k) = .TRUE.
	alpha = ABS( rdf_zs_ml(k) - rdf_zm(k) ) - ABS( rdf_zs_ml(k) - rdf_zm(k+1) ) !+PAB change
	alpha = SIGN( 0.5_r8, alpha ) + 0.5_r8 ! alpha = 0 or 1
	rdf_zs_ml(k) = alpha * rdf_zm(k+1) + ( 1.0_r8 - alpha ) * rdf_zm(k) !+PAB change
	if ( alpha < 0.5_r8 ) then
	  rdf_zm_up(k) = ic3(k) * ( c0(k)*rdf_z(k) - c1(k)*rdf_zs_ml(k) - c2(k)*rdf_zm(k) )
	else
	  rdf_zm_dn(k) = ic2(k) * (c0(k)*rdf_z(k) - c1(k)*rdf_zs_ml(k) - c3(k)*rdf_zm(k+1)) !+PAB change
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
    
    else
    
    rdf_zm_up(:) = rdf_zm(2:pverp)
    rdf_zm_dn(:) = rdf_zm(1:pver)
    
    endif
    
    ! Construct the tendency profile 
    
    !loop over patches
    lastPatchIndexHi = 0
    do iPatch = 1,five_num_patches
     
     zi1 = pver - patchIndexLo(iPatch)+ 1
     zi2 = pver - patchIndexHi(iPatch)+ 1
  
     write(*,*)"zi1,iPatch = **",zi1,iPatch
     write(*,*)"zi2,iPatch = **",zi2,ipatch


     i5 = zi1 - 1
     do k = lastPatchIndexHi + 1,zi1 - 1
       df_zs(k) = df_z(k)
     enddo

    ! between zi1_1 and zi2
    do k = zi1, zi2
   
     i1 = i5 + 1
     i2 = i1 + fiveAddNlevels(iPatch) / 2 - 1 ! i2 = i3-1 or i2 = i1
     i3 = i1 + fiveAddNlevels(iPatch) / 2     ! zs(i3) = z(k)
     i4 = i1 + fiveAddNlevels(iPatch) / 2 + 1 ! i4 = i3+1 or i4 = i5
     i5 = i1 + fiveAddNlevels(iPatch)
      
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

    lastPatchIndexHi = patchIndexHi(iPatch)
    !end loop over patches
   enddo 

    ! Above last Patch
    do k = lastPatchIndexHi + 1,pver
      i5 = i5 + 1
      df_zs(i5) = df_z(k)
    enddo    
    
    ! Finally, "unflip" the tendency
    do k = 1,pver_five
      ten_high(k) = df_zs(pver_five-k+1)
    enddo
  
   write(*,*)"Leaving function: tendency_low_to_high"
   write(*,*)"  "  
   deallocate(patchIndexLo)
   deallocate(patchIndexHi)
  
  end subroutine tendency_low_to_high 
  
  ! ======================================== !
  !                                          !
  ! ======================================== !
  
  subroutine find_level_match_index(&
           e3sm_pmid,five_pmid,five_pint,top_lev_e3sm,top_lev_five)
  
    ! Purpose is to find FIVE level indicee that matches
    !  the level of a given e3sm level indicee
    implicit none
  
    real(r8), intent(in) :: e3sm_pmid(pver)
    real(r8), intent(in) :: five_pmid(pver)
    real(r8), intent(in) :: five_pint(pverp)
    integer, intent(in) :: top_lev_e3sm
    integer, intent(out) :: top_lev_five
    
    real(r8) :: pmid_target
    integer :: i, k
    write(*,*)"Entering function: find_level_match_index" 
    pmid_target = e3sm_pmid(top_lev_e3sm)
    
    do k=1,pver_five
      ! if the pressures have not been modified, then just
      !   retun the indicee where they are equal
      if (five_pmid(k) .eq. pmid_target) then
        top_lev_five = k
	return
      endif     

      ! IF the best above didn't pass then see where
      !  the level lies between pressure interfaces      
      if (five_pint(k) .lt. pmid_target .and. five_pint(k+1) .gt. pmid_target) then
        top_lev_five = k
	return	  
      endif	      
    enddo
     
  write(*,*)"Leaving function: find_level_match_index"
  write(*,*)"  " 
  return     
  end subroutine find_level_match_index

end module five_intr
