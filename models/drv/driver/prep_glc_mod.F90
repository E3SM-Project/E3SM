module prep_glc_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8 
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_glc, num_inst_lnd, num_inst_frc, &
                              num_inst_ocn  
  use seq_comm_mct    , only: CPLID, GLCID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs 
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata  
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: glc, lnd, ocn

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_glc_init
  public :: prep_glc_mrg

  public :: prep_glc_accum
  public :: prep_glc_accum_avg

  public :: prep_glc_calc_l2x_gx
  public :: prep_glc_calc_o2x_gx  

  public :: prep_glc_get_l2x_gx
  public :: prep_glc_get_l2gacc_lx
  public :: prep_glc_get_l2gacc_lx_cnt
  
  public :: prep_glc_get_o2x_gx
  public :: prep_glc_get_x2gacc_gx
  public :: prep_glc_get_x2gacc_gx_cnt  
  
  public :: prep_glc_get_mapper_SFl2g
  public :: prep_glc_get_mapper_SFo2g 

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_SFl2g
  type(seq_map), pointer :: mapper_SFo2g

  ! attribute vectors 
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: o2x_gx(:) ! Ocn export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated
  
  type(mct_aVect), pointer :: x2gacc_gx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: x2gacc_gx_cnt ! x2gacc_gx: number of time samples accumulated  

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
  
  real(r8), allocatable ::  oceanTemperature(:)
  real(r8), allocatable ::  oceanSalinity(:)
  real(r8), allocatable ::  oceanHeatTransferVelocity(:)
  real(r8), allocatable ::  oceanSaltTransferVelocity(:)
  real(r8), allocatable ::  interfacePressure(:)
  real(r8), allocatable ::  iceTemperature(:)
  real(r8), allocatable ::  iceTemperatureDistance(:)
  real(r8), allocatable ::  outInterfaceSalinity(:)
  real(r8), allocatable ::  outInterfaceTemperature(:)
  real(r8), allocatable ::  outFreshwaterFlux(:)
  real(r8), allocatable ::  outOceanHeatFlux(:)
  real(r8), allocatable ::  outIceHeatFlux(:)
  
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_glc_init(infodata, lnd_c2_glc, ocn_c2_glc)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and mapping variables
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: lnd_c2_glc ! .true.  => lnd to glc coupling on
    logical                  , intent(in)    :: ocn_c2_glc ! .true.  => ocn to glc coupling on    
    !
    ! Local Variables
    integer                          :: eli, egi, eoi
    integer                          :: lsize_l
    integer                          :: lsize_g 
    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                          :: glc_present   ! .true. => glc is present
    logical                          :: samegrid_go   ! .true. => samegrid ocean and glc
    character(CL)                    :: lnd_gnam      ! lnd grid
    character(CL)                    :: glc_gnam      ! glc grid
    character(CL)                    :: ocn_gnam      ! ocn grid    
    type(mct_avect), pointer         :: l2x_lx  
    type(mct_avect), pointer         :: x2g_gx
    type(mct_avect), pointer         :: o2x_ox
    type(mct_avect), pointer         :: g2x_gx 
          
    character(*), parameter          :: subname = '(prep_glc_init)'
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         glc_present=glc_present       , &
         lnd_gnam=lnd_gnam             , &
         glc_gnam=glc_gnam             , &
	 ocn_gnam=ocn_gnam)

    allocate(mapper_SFl2g)
    allocate(mapper_SFo2g)
    
    if (glc_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)
	    
       g2x_gx => component_get_c2x_cx(glc(1))
       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(g2x_gx)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)

       allocate(l2x_gx(num_inst_lnd))
       allocate(l2gacc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2x_gx(eli) ,lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))
          
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2gacc_lx(eli), lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       enddo
       l2gacc_lx_cnt = 0

       allocate(o2x_gx(num_inst_ocn))
       allocate(x2gacc_gx(num_inst_glc))
       
       do eoi = 1,num_inst_ocn    
          call mct_aVect_init(o2x_gx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_g)
          call mct_aVect_zero(o2x_gx(eoi))
       enddo

       do egi = 1,num_inst_glc
          call mct_avect_init(x2gacc_gx(egi), x2g_gx, lsize_g)
          call mct_aVect_zero(x2gacc_gx(egi))
       end do           

       x2gacc_gx_cnt = 0            

       if (lnd_c2_glc) then
	 if (iamroot_CPLID) then
           write(logunit,*) ' '
           write(logunit,F00) 'Initializing mapper_SFl2g'
	 end if
	 call seq_map_init_rearrolap(mapper_SFl2g, lnd(1), glc(1), 'mapper_SFl2g')
	 call shr_sys_flush(logunit)
       end if

       if (ocn_c2_glc) then	       
	   if (iamroot_CPLID) then
              write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_SFo2g'
	   end if
           samegrid_go = .true.
	 if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_go = .false.

	 call seq_map_init_rcfile(mapper_SFo2g, ocn(1), glc(1), &
	 'seq_maps.rc','ocn2glc_fmapname:','ocn2glc_fmaptype:',samegrid_go, &
	 'mapper_SFo2g initialization',esmf_map_flag)
       end if

       !Initialize module-level arrays associated with compute_melt_fluxes
       allocate(oceanTemperature(lsize_g))
       allocate(oceanSalinity(lsize_g))
       allocate(oceanHeatTransferVelocity(lsize_g))
       allocate(oceanSaltTransferVelocity(lsize_g))
       allocate(interfacePressure(lsize_g))
       allocate(iceTemperature(lsize_g))
       allocate(iceTemperatureDistance(lsize_g))
       allocate(outInterfaceSalinity(lsize_g))
       allocate(outInterfaceTemperature(lsize_g))
       allocate(outFreshwaterFlux(lsize_g))
       allocate(outOceanHeatFlux(lsize_g))
       allocate(outIceHeatFlux(lsize_g))

       call shr_sys_flush(logunit)

    end if

  end subroutine prep_glc_init

  !================================================================================================

  subroutine prep_glc_accum(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate glc inputs
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli, egi
    type(mct_avect), pointer :: l2x_lx
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter :: subname = '(prep_glc_accum)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       l2x_lx => component_get_c2x_cx(lnd(eli))
       if (l2gacc_lx_cnt == 0) then
          call mct_avect_copy(l2x_lx, l2gacc_lx(eli))
       else
          call mct_avect_accum(l2x_lx, l2gacc_lx(eli))
       endif
    end do
    l2gacc_lx_cnt = l2gacc_lx_cnt + 1
   
    do egi = 1,num_inst_glc
       x2g_gx => component_get_x2c_cx(glc(egi))
       if (x2gacc_gx_cnt == 0) then
          call mct_avect_copy(x2g_gx, x2gacc_gx(egi))
       else
          call mct_avect_accum(x2g_gx, x2gacc_gx(egi))
       endif
    end do
    x2gacc_gx_cnt = x2gacc_gx_cnt + 1
    
    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_accum

  !================================================================================================

  subroutine prep_glc_accum_avg(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of glc inputs
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli, egi
    type(mct_avect), pointer :: x2g_gx  
    
    character(*), parameter :: subname = '(prep_glc_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (l2gacc_lx_cnt > 1) then
       do eli = 1,num_inst_lnd
          call mct_avect_avg(l2gacc_lx(eli), l2gacc_lx_cnt)
       end do
    end if
    l2gacc_lx_cnt = 0
    
    do egi = 1,num_inst_glc
       ! temporary formation of average
       if (x2gacc_gx_cnt > 1) then
          call mct_avect_avg(x2gacc_gx(egi), x2gacc_gx_cnt)
       end if
       
       ! ***NOTE***THE FOLLOWING ACTUALLY MODIFIES x2g_gx
       x2g_gx   => component_get_x2c_cx(glc(egi)) 
       call mct_avect_copy(x2gacc_gx(egi), x2g_gx)
    enddo
    x2gacc_gx_cnt = 0
    
    call t_drvstopf  (trim(timer))
    
  end subroutine prep_glc_accum_avg

  !================================================================================================
  
  subroutine prep_glc_mrg(infodata, timer_mrg) 

    !---------------------------------------------------------------
    ! Description
    ! Merge glc inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer :: egi, eli, eoi
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter  :: subname = '(prep_glc_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       eoi = mod((egi-1),num_inst_ocn) + 1

       x2g_gx => component_get_x2c_cx(glc(egi))
       
       call prep_glc_merge(l2x_gx(eli), o2x_gx(eoi), x2g_gx)
       
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg

  !================================================================================================

  subroutine prep_glc_merge( s2x_g, o2x_g, x2g_g )

    !----------------------------------------------------------------------- 
    ! Arguments
    type(mct_aVect), intent(inout)  :: s2x_g  ! input
    type(mct_aVect), intent(inout)  :: o2x_g  ! input
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !----------------------------------------------------------------------- 

    integer       :: ngflds,noflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    integer, save :: index_o2x_So_tglc
    integer, save :: index_x2g_So_tglc    
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: s2x_sharedindices
    type(mct_aVect_sharedindices),save :: o2x_sharedindices
        
    character(*), parameter   :: subname = '(prep_glc_merge) '

    !----------------------------------------------------------------------- 

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    ngflds = mct_aVect_nRattr(x2g_g)
    noflds = mct_aVect_nRattr(o2x_g)

    if (first_time) then
        
	index_o2x_So_tglc = mct_aVect_indexRA(x2g_g,'So_tglc')
	index_x2g_So_tglc = mct_aVect_indexRA(o2x_g,'So_tglc')
	
	!Jer: document fields being passed to GLC
       allocate(mrgstr(ngflds))
       do i = 1,ngflds
          field = mct_aVect_getRList2c(i, x2g_g)
          mrgstr(i) = subname//'x2g%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(s2x_g, x2g_g, s2x_SharedIndices)
       call mct_aVect_setSharedIndices(o2x_g, x2g_g, o2x_SharedIndices)       

       !--- document copy operations ---
       do i=1,s2x_SharedIndices%shared_real%num_indices
          i1=s2x_SharedIndices%shared_real%aVindices1(i)
          o1=s2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, s2x_g)
          mrgstr(o1) = trim(mrgstr(o1))//' = s2x%'//trim(field)
       enddo
       do i=1,o2x_SharedIndices%shared_real%num_indices
          i1=o2x_SharedIndices%shared_real%aVindices1(i)
          o1=o2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, o2x_g)
          mrgstr(o1) = trim(mrgstr(o1))//' = o2x%'//trim(field)
       enddo       
    endif

    ! Create input glc state directly from land snow output state
    call mct_aVect_copy(aVin=s2x_g, aVout=x2g_g, vector=mct_usevector, sharedIndices=s2x_SharedIndices)
    
    !Create input glc state directly from ocean model output state (once remapped to glc grid)
    call mct_aVect_copy(aVin=o2x_g, aVout=x2g_g, vector=mct_usevector, sharedIndices=o2x_SharedIndices)   

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,ngflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif

    first_time = .false.

  end subroutine prep_glc_merge

  !================================================================================================

  subroutine prep_glc_calc_o2x_gx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_gx
    
    ! Arguments
    character(len=*), intent(in) :: timer
    
    character(*), parameter :: subname = '(prep_glc_calc_o2x_gx)'
    ! Local Variables
    integer eoi
    type(mct_avect), pointer :: o2x_ox
    
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)   
    do eoi = 1,num_inst_ocn
      o2x_ox => component_get_c2x_cx(ocn(eoi))
      call seq_map_map(mapper_SFo2g, o2x_ox, o2x_gx(eoi), &
           fldlist='So_tglc',norm=.true.)
    enddo
    
    call t_drvstopf  (trim(timer))     
  end subroutine prep_glc_calc_o2x_gx    

  !================================================================================================

  subroutine prep_glc_calc_l2x_gx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_gx (note that l2x_gx is a local module variable)
    ! Also l2x_gx is really the accumulated l2xacc_lx mapped to l2x_gx
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       call seq_map_map(mapper_SFl2g, l2gacc_lx(eli), l2x_gx(eli), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  subroutine prep_glc_calc_calculate_subshelf_boundary_fluxes
  
    !---------------------------------------------------------------
    ! Description
    ! On the ice sheet grid, calculate shelf boundary fluxes

    use shr_const_mod , only: SHR_CONST_KAPPA_LAND_ICE

    ! Local Variables

    integer :: gsize, egi, n, err

    !---------------------------------------------------------------

    gsize = mct_aVect_lsize(x2g_gx)
    err = 0
    
    do egi = 1,num_inst_glc
      
      for n=1:gsize
        !Extract coupler fields used as input to compute_melt_fluxes to local arrays...
	
        oceanTemperature(n)	     = x2g_gx(egi)%rAttr(index_x2g_So_blt,n)
	oceanSalinity(n)	     = x2g_gx(egi)%rAttr(index_x2g_So_bls,n)
	oceanHeatTransferVelocity(n) = x2g_gx(egi)%rAttr(index_x2g_So_htv,n)
	oceanSaltTransferVelocity(n) = x2g_gx(egi)%rAttr(index_x2g_So_hsv,n)
	interfacePressure(n)	     = x2g_gx(egi)%rAttr(index_x2g_So_phieff,n)
	
	iceTemperature(n)	     = g2x_gx(egi)%rAttr(index_g2x_Sg_tbot,n)
	iceTemperatureDistance(n)    = g2x_gx(egi)%rAttr(index_g2x_Sg_dztbot,n)
	
	!...and initialize local compute_melt_fluxes output arrays.
	outInterfaceSalinity(n)     = 0.0_r8
        outInterfaceTemperature(n)  = 0.0_r8
        outFreshwaterFlux(n)	     = 0.0_r8
	outOceanHeatFlux(n)	     = 0.0_r8
	outIceHeatFlux(n)	     = 0.0_r8
      end
    
      call compute_melt_fluxes(oceanTemperature,&
                               oceanSalinity,&
                               oceanHeatTransferVelocity,&
                               oceanSaltTransferVelocity,&
                               interfacePressure,&
                               outInterfaceSalinity,&
                               outInterfaceTemperature,&
                               outFreshwaterFlux,&
                               outOceanHeatFlux,&
                               outIceHeatFlux,&
                               gsize,&
                               err,&
                               iceTemperature,&
                               iceTemperatureDistance,&
                               SHR_CONST_KAPPA_LAND_ICE &
                               )
      
      for n=1:gsize
      
        !Assign outputs from compute_melt_fluxes back into coupler attributes
	
	g2x_gx(egi)%rAttr(index_g2x_Sg_blis,n) = outInterfaceSalinity(n)      !to ocean
        g2x_gx(egi)%rAttr(index_g2x_Sg_blit,n) = outInterfaceTemperature(n)   !to ocean
	g2x_gx(egi)%rAttr(index_g2x_Fogx_qiceho,n) = outOceanHeatFlux(n)      !to ocean
        g2x_gx(egi)%rAttr(index_g2x_Fogx_qicelo,n) = outFreshwaterFlux(n)     !to ocean... need unit conversion?
	x2g_gx(egi)%rAttr(index_x2g_Fogx_qicehi,n) = outIceHeatFlux(n)        !to ice sheet
        x2g_gx(egi)%rAttr(index_x2g_Fogx_qiceli,n) = outFreshwaterFlux(n)     !to ice sheet... need unit conversion?
	
      end
      
    end
 
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

  end subroutine prep_glc_calc_calculate_subshelf_boundary_fluxes

  !================================================================================================

  function prep_glc_get_l2x_gx()
    type(mct_aVect), pointer :: prep_glc_get_l2x_gx(:)
    prep_glc_get_l2x_gx => l2x_gx(:)   
  end function prep_glc_get_l2x_gx

  function prep_glc_get_l2gacc_lx()
    type(mct_aVect), pointer :: prep_glc_get_l2gacc_lx(:)
    prep_glc_get_l2gacc_lx => l2gacc_lx(:)   
  end function prep_glc_get_l2gacc_lx

  function prep_glc_get_l2gacc_lx_cnt()
    integer, pointer :: prep_glc_get_l2gacc_lx_cnt
    prep_glc_get_l2gacc_lx_cnt => l2gacc_lx_cnt
  end function prep_glc_get_l2gacc_lx_cnt

  function prep_glc_get_o2x_gx()
    type(mct_aVect), pointer :: prep_glc_get_o2x_gx(:)
    prep_glc_get_o2x_gx => o2x_gx(:)   
  end function prep_glc_get_o2x_gx

  function prep_glc_get_x2gacc_gx()
    type(mct_aVect), pointer :: prep_glc_get_x2gacc_gx(:)
    prep_glc_get_x2gacc_gx => x2gacc_gx(:)   
  end function prep_glc_get_x2gacc_gx    

  function prep_glc_get_x2gacc_gx_cnt()
    integer, pointer :: prep_glc_get_x2gacc_gx_cnt
    prep_glc_get_x2gacc_gx_cnt => x2gacc_gx_cnt
  end function prep_glc_get_x2gacc_gx_cnt        

  function prep_glc_get_mapper_SFl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_SFl2g
    prep_glc_get_mapper_SFl2g => mapper_SFl2g  
  end function prep_glc_get_mapper_SFl2g
  
  function prep_glc_get_mapper_SFo2g()
    type(seq_map), pointer :: prep_glc_get_mapper_SFo2g
    prep_glc_get_mapper_SFo2g=> mapper_SFo2g  
  end function prep_glc_get_mapper_SFo2g 

end module prep_glc_mod
