module prep_glc_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8
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
  
  public :: prep_glc_get_mapper_Sl2g
  public :: prep_glc_get_mapper_Fl2g

  public :: prep_glc_get_mapper_So2g
  public :: prep_glc_get_mapper_Fo2g
  
  public :: prep_glc_calculate_subshelf_boundary_fluxes
  
  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_merge
  private :: prep_glc_map_one_field_lnd2glc

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sl2g
  type(seq_map), pointer :: mapper_Fl2g
  type(seq_map), pointer :: mapper_So2g
  type(seq_map), pointer :: mapper_Fo2g

  ! attribute vectors 
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: o2x_gx(:) ! Ocn export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated
  
  type(mct_aVect), pointer :: x2gacc_gx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: x2gacc_gx_cnt ! x2gacc_gx: number of time samples accumulated  

  type(mct_aVect), pointer :: o2gacc_ox(:) ! Ocn export, lnd grid, cpl pes - allocated in driver
  integer        , target :: o2gacc_ox_cnt ! number of time samples accumulated  

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

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
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
    integer                         ::  lsize_o
    logical                          :: samegrid_lg   ! samegrid land and glc
    logical                          :: samegrid_go   ! .true. => samegrid ocean and glc
    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    logical                          :: glc_present   ! .true. => glc is present
    
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

    allocate(mapper_Sl2g)
    allocate(mapper_Fl2g)
    allocate(mapper_So2g)
    allocate(mapper_Fo2g)   

    g2x_gx => component_get_c2x_cx(glc(1))
    x2g_gx => component_get_x2c_cx(glc(1))
    lsize_g = mct_aVect_lsize(g2x_gx)

    if (glc_present .and. lnd_c2_glc) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)
              
       allocate(l2x_gx(num_inst_lnd))
       allocate(l2gacc_lx(num_inst_lnd))
       
       do eli = 1,num_inst_lnd
          call mct_aVect_init(l2x_gx(eli), rList=seq_flds_x2g_fields, lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))

          call mct_aVect_init(l2gacc_lx(eli), rList=seq_flds_l2x_fields_to_glc, lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       enddo
       l2gacc_lx_cnt = 0

       samegrid_lg = .true.
       if (trim(lnd_gnam) /= trim(glc_gnam)) samegrid_lg = .false.

       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_Sl2g'
       end if
       call seq_map_init_rcfile(mapper_Sl2g, lnd(1), glc(1), &
            'seq_maps.rc', 'lnd2glc_smapname:', 'lnd2glc_smaptype:', samegrid_lg, &
            'mapper_Sl2g initialization', esmf_map_flag)
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_Fl2g'
       end if
       call seq_map_init_rcfile(mapper_Fl2g, lnd(1), glc(1), &
            'seq_maps.rc', 'lnd2glc_fmapname:', 'lnd2glc_fmaptype:', samegrid_lg, &
            'mapper_Fl2g initialization', esmf_map_flag)

       call shr_sys_flush(logunit)
          
    end if

    if (glc_present .and. ocn_c2_glc) then

       allocate(o2x_gx(num_inst_ocn))
       do eoi = 1,num_inst_ocn
	  call mct_aVect_init(o2x_gx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_g)
	  call mct_aVect_zero(o2x_gx(eoi))
       enddo

       allocate(x2gacc_gx(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(x2gacc_gx(egi), x2g_gx, lsize_g)
          call mct_aVect_zero(x2gacc_gx(egi))
       end do

       o2x_ox => component_get_c2x_cx(ocn(1))
       lsize_o = mct_aVect_lsize(o2x_ox)

       x2gacc_gx_cnt = 0
       samegrid_go = .true.
       if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_go = .false.
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_So2g'
       end if
       call seq_map_init_rcfile(mapper_So2g, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_smapname:','ocn2glc_smaptype:',samegrid_go, &
       'mapper_So2g initialization',esmf_map_flag)
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_Fo2g'
       end if
       call seq_map_init_rcfile(mapper_Fo2g, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_fmapname:','ocn2glc_fmaptype:',samegrid_go, &
       'mapper_Fo2g initialization',esmf_map_flag)

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
	  call mct_avect_accum(x2g_gx, x2gacc_gx(egi))  !This is accumulating into x2gacc_gx
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
       x2g_gx	=> component_get_x2c_cx(glc(egi)) 
       call mct_avect_copy(x2gacc_gx(egi), x2g_gx)
    enddo
    x2gacc_gx_cnt = 0
    

    
    call t_drvstopf  (trim(timer))
    
  end subroutine prep_glc_accum_avg

  !================================================================================================
  
  subroutine prep_glc_mrg(infodata, fractions_gx, timer_mrg) 

    !---------------------------------------------------------------
    ! Description
    ! Merge glc inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_gx(:)
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer :: egi, eli, efi, eoi
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter  :: subname = '(prep_glc_mrg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       efi = mod((egi-1),num_inst_frc) + 1
       eoi = mod((egi-1),num_inst_ocn) + 1

       x2g_gx => component_get_x2c_cx(glc(egi)) 
       call prep_glc_merge(l2x_gx(eli), fractions_gx(efi), x2g_gx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg

  !================================================================================================

  subroutine prep_glc_merge( l2x_g, fractions_g, x2g_g )

    !----------------------------------------------------------------------- 
    ! Description
    ! "Merge" land forcing for glc input.
    !
    ! State fields are copied directly, meaning that averages are taken just over the
    ! land-covered portion of the glc domain.
    !
    ! Flux fields are downweighted by landfrac, which effectively sends a 0 flux from the
    ! non-land-covered portion of the glc domain.
    !
    ! Arguments
    type(mct_aVect), intent(inout)  :: l2x_g  ! input
    type(mct_aVect), intent(in)     :: fractions_g
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !----------------------------------------------------------------------- 

    integer       :: num_land_to_glc_flux_fields
    integer       :: num_land_to_glc_state_fields
    integer       :: nflds
    integer       :: i,n
    integer       :: mrgstr_index
    integer       :: index_l2x
    integer       :: index_x2g
    integer       :: index_lfrac
    integer       :: lsize
    logical       :: iamroot
    real(r8)      :: lfrac
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    character(*), parameter   :: subname = '(prep_glc_merge) '

    !----------------------------------------------------------------------- 

    call seq_comm_getdata(CPLID, iamroot=iamroot)
    lsize = mct_aVect_lsize(x2g_g)
    
    num_land_to_glc_flux_fields = shr_string_listGetNum(trim(seq_flds_l2x_fluxes_to_glc))
    num_land_to_glc_state_fields = shr_string_listGetNum(trim(seq_flds_l2x_states_to_glc))
    
    if (first_time) then
       nflds = mct_aVect_nRattr(l2x_g)
       if (nflds /= (num_land_to_glc_flux_fields + num_land_to_glc_state_fields)) then
          write(logunit,*) subname,' ERROR: nflds /= num_land_to_glc_flux_fields + num_land_to_glc_state_fields: ', &
               nflds, num_land_to_glc_flux_fields, num_land_to_glc_state_fields
          call shr_sys_abort(subname//' ERROR: nflds /= num_land_to_glc_flux_fields + num_land_to_glc_state_fields')
       end if
          
       allocate(mrgstr(nflds))
    end if

    mrgstr_index = 1

    do i = 1, num_land_to_glc_state_fields
       call seq_flds_getField(field, i, seq_flds_l2x_states_to_glc)
       index_l2x = mct_aVect_indexRA(l2x_g, trim(field))
       index_x2g = mct_aVect_indexRA(x2g_g, trim(field))

       if (first_time) then
          mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
               ' = l2x%'//trim(field)
       end if

       do n = 1, lsize
          x2g_g%rAttr(index_x2g,n) = l2x_g%rAttr(index_l2x,n)
       end do

       mrgstr_index = mrgstr_index + 1
    enddo

    index_lfrac = mct_aVect_indexRA(fractions_g,"lfrac")
    do i = 1, num_land_to_glc_flux_fields
       call seq_flds_getField(field, i, seq_flds_x2g_fluxes)
       index_l2x = mct_aVect_indexRA(l2x_g, trim(field))
       index_x2g = mct_aVect_indexRA(x2g_g, trim(field))
       
       if (first_time) then
          mrgstr(mrgstr_index) = subname//'x2g%'//trim(field)//' =' // &
               ' = lfrac*l2x%'//trim(field)
       end if
       
       do n = 1, lsize
          lfrac = fractions_g%rAttr(index_lfrac,n)
          x2g_g%rAttr(index_x2g,n) = l2x_g%rAttr(index_l2x,n) * lfrac
       end do

       mrgstr_index = mrgstr_index + 1
    end do

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,nflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif

    first_time = .false.

  end subroutine prep_glc_merge


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
      call seq_map_map(mapper_So2g, o2x_ox, o2x_gx(eoi), &
	   fldlist='So_blt:So_bls:So_htv:So_stv:So_rhoeff',norm=.true.)
    enddo
    
    call t_drvstopf  (trim(timer))     
  end subroutine prep_glc_calc_o2x_gx	 

  !================================================================================================


  !================================================================================================

  subroutine prep_glc_calc_l2x_gx(fractions_lx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_gx (note that l2x_gx is a local module variable)
    ! Also l2x_gx is really the accumulated l2xacc_lx mapped to l2x_gx
    !
    use shr_string_mod, only : shr_string_listGetNum
    ! Arguments
    type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi, eli, efi
    integer :: num_flux_fields
    integer :: num_state_fields
    integer :: field_num
    character(len=cl) :: fieldname
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)

    num_flux_fields = shr_string_listGetNum(trim(seq_flds_x2g_fluxes))
    num_state_fields = shr_string_listGetNum(trim(seq_flds_x2g_states))
    
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       efi = mod((egi-1),num_inst_frc) + 1
       
       do field_num = 1, num_flux_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_fluxes)
          call prep_glc_map_one_field_lnd2glc(egi=egi, eli=eli, &
               fieldname = fieldname, &
               fractions_lx = fractions_lx(efi), &
               mapper = mapper_Fl2g)
       end do

       do field_num = 1, num_state_fields
          call seq_flds_getField(fieldname, field_num, seq_flds_x2g_states)
          call prep_glc_map_one_field_lnd2glc(egi=egi, eli=eli, &
               fieldname = fieldname, &
               fractions_lx = fractions_lx(efi), &
               mapper = mapper_Sl2g)
       end do
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  subroutine prep_glc_map_one_field_lnd2glc(egi, eli, fieldname, fractions_lx, mapper)
    ! Maps a single field from the land grid to the glc grid.
    !
    ! Note that we remap each field separately because each field needs its own
    ! vertical gradient calculator.

    use vertical_gradient_calculator_2nd_order, only : vertical_gradient_calculator_2nd_order_type
    use glc_elevclass_mod, only : glc_get_num_elevation_classes
    use map_lnd2glc_mod, only : map_lnd2glc

    ! Arguments
    integer, intent(in) :: egi  ! glc instance index
    integer, intent(in) :: eli  ! lnd instance index
    character(len=*), intent(in) :: fieldname  ! base name of field to map (without elevation class suffix)
    type(mct_aVect) , intent(in) :: fractions_lx  ! fractions on the land grid, for this frac instance
    type(seq_map), intent(inout) :: mapper
    !
    ! Local Variables
    type(mct_avect), pointer :: g2x_gx, l2x_lx
    type(vertical_gradient_calculator_2nd_order_type) :: gradient_calculator
    !---------------------------------------------------------------

    g2x_gx => component_get_c2x_cx(glc(egi))
    l2x_lx => component_get_c2x_cx(lnd(eli))

    gradient_calculator = vertical_gradient_calculator_2nd_order_type( &
         attr_vect = l2x_lx, &
         fieldname = fieldname, &
         toponame = 'Sl_topo', &
         min_elevation_class = 1, &
         max_elevation_class = glc_get_num_elevation_classes())
    call map_lnd2glc(l2x_l = l2x_lx, &
         landfrac_l = fractions_lx, &
         g2x_g = g2x_gx, &
         fieldname = fieldname, &
         gradient_calculator = gradient_calculator, &
         mapper = mapper, &
         l2x_g = l2x_gx(eli))

  end subroutine prep_glc_map_one_field_lnd2glc

  !================================================================================================

  subroutine prep_glc_calculate_subshelf_boundary_fluxes
  
    !---------------------------------------------------------------
    ! Description
    ! On the ice sheet grid, calculate shelf boundary fluxes

    use glc_cpl_indices
    use mpaso_cpl_indices
    
    use shr_const_mod , only: SHR_CONST_KAPPA_LAND_ICE

    ! Local Variables

    integer :: gsize, n
    type(mct_aVect), pointer :: o2x_ox ! Ocn export, ocn grid, cpl pes
    type(mct_aVect), pointer :: x2g_gx ! Glc import, glc grid, cpl pes
    type(mct_aVect), pointer :: g2x_gx ! Glc import, glc grid, cpl pes
    
    character(*), parameter :: subname = '(prep_glc_calculate_subshelf_boundary_fluxes)'
    !--------------------------------------------------------------- 

    gsize = mct_aVect_lsize(o2x_gx(1))

    o2x_ox => component_get_c2x_cx(ocn(1))
    g2x_gx => component_get_c2x_cx(glc(1))
    x2g_gx => component_get_x2c_cx(glc(1))

    !TO DO: CONSIDER HOW ENSEMBLES ARE GOING TO BE IMPLEMENTED (i.e. change hard-coded '1' index to e*i indices

    !Remap relevant ocean variables to ice sheet grid.
    !Done here instead of in glc-frequency mapping so it happens within ocean coupling interval.
    call seq_map_map(mapper_So2g, o2x_ox, o2x_gx(1), &
    fldlist='So_blt:So_bls:So_htv:So_stv:So_rhoeff',norm=.true.)

    do n=1,gsize
      !Extract coupler fields used as input to compute_melt_fluxes to local arrays...

      oceanTemperature(n) =               o2x_gx(1)%rAttr(index_o2x_So_blt,n)
      oceanSalinity(n) =                  o2x_gx(1)%rAttr(index_o2x_So_bls,n)
      oceanHeatTransferVelocity(n) =      o2x_gx(1)%rAttr(index_o2x_So_htv,n)
      oceanSaltTransferVelocity(n) =      o2x_gx(1)%rAttr(index_o2x_So_stv,n)
      interfacePressure(n) =              o2x_gx(1)%rAttr(index_o2x_So_rhoeff,n)

      iceTemperature(n) =                 g2x_gx%rAttr(index_g2x_Sg_tbot,n)
      iceTemperatureDistance(n) =         g2x_gx%rAttr(index_g2x_Sg_dztbot,n)

      !...and initialize local compute_melt_fluxes output arrays.
      outInterfaceSalinity(n)     =       0.0_r8
      outInterfaceTemperature(n)  =       0.0_r8
      outFreshwaterFlux(n)  =             0.0_r8
      outOceanHeatFlux(n) =               0.0_r8
      outIceHeatFlux(n) =                 0.0_r8
    end do

    call compute_melt_fluxes(oceanTemperature,&
                             oceanSalinity,&
                             oceanHeatTransferVelocity,&
                             oceanSaltTransferVelocity,&
                             interfacePressure,&
                             iceTemperature,&
                             iceTemperatureDistance, &
                             outInterfaceSalinity,&
                             outInterfaceTemperature,&
                             outFreshwaterFlux,&
                             outOceanHeatFlux,&
                             outIceHeatFlux,&
                             gsize)

    do n=1,gsize
      !Assign outputs from compute_melt_fluxes back into coupler attributes
      g2x_gx%rAttr(index_g2x_Sg_blis,n) =     outInterfaceSalinity(n)    !to ocean
      g2x_gx%rAttr(index_g2x_Sg_blit,n) =     outInterfaceTemperature(n) !to ocean
      g2x_gx%rAttr(index_g2x_Fogx_qiceho,n) = outOceanHeatFlux(n)        !to ocean
      g2x_gx%rAttr(index_g2x_Fogx_qicelo,n)=  outFreshwaterFlux(n)       !to ocean
      
      x2g_gx%rAttr(index_x2g_Fogx_qicehi,n) = outIceHeatFlux(n)          !to ice sheet 
      x2g_gx%rAttr(index_x2g_Fogx_qiceli,n) = outFreshwaterFlux(n)       !to ice sheet
    end do

    !Remap ocean-side outputs back onto ocean grid done in call to prep_ocn_shelf_calc_g2x_ox

  end subroutine prep_glc_calculate_subshelf_boundary_fluxes

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

  function prep_glc_get_mapper_Sl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Sl2g
    prep_glc_get_mapper_Sl2g => mapper_Sl2g  
  end function prep_glc_get_mapper_Sl2g

  function prep_glc_get_mapper_Fl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Fl2g
    prep_glc_get_mapper_Fl2g => mapper_Fl2g  
  end function prep_glc_get_mapper_Fl2g

  function prep_glc_get_mapper_So2g()
    type(seq_map), pointer :: prep_glc_get_mapper_So2g
    prep_glc_get_mapper_So2g=> mapper_So2g  
  end function prep_glc_get_mapper_So2g
  
  function prep_glc_get_mapper_Fo2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Fo2g
    prep_glc_get_mapper_Fo2g=> mapper_Fo2g  
  end function prep_glc_get_mapper_Fo2g  

!***********************************************************************
!
!  routine ocn_forcing_compute_melt_fluxes
!
!> \brief   Computes ocean and ice melt fluxes, etc.
!> \author  Xylar Asay-Davis
!> \date    3/27/2015
!>  This routine computes melt fluxes (melt rate, temperature fluxes
!>  into the ice and the ocean, and salt flux) as well as the interface
!>  temperature and salinity.  This routine expects an ice temperature
!>  in the bottom layer of ice and ocean temperature and salinity in
!>  the top ocean layer as well as the pressure at the ice/ocean interface.
!>
!>  The ocean heat and salt transfer velocities are determined based on
!>  observations of turbulent mixing rates in the under-ice boundary layer.
!>  They should be the product of the friction velocity and a (possibly
!>  spatially variable) non-dimenional transfer coefficient.
!>  
!>  The iceTemperatureDistance is the distance between the location
!>  where the iceTemperature is supplied and the ice-ocean interface,
!>  used to compute a temperature gradient.  The ice thermal conductivity,
!>  SHR_CONST_KAPPA_LAND_ICE, is zero for the freezing solution from Holland and Jenkins
!>  (1999) in which the ice is purely insulating.
!
!-----------------------------------------------------------------------

  subroutine compute_melt_fluxes( &
    oceanTemperature, &
    oceanSalinity, &
    oceanHeatTransferVelocity, &
    oceanSaltTransferVelocity, &
    interfacePressure, &
    iceTemperature, &
    iceTemperatureDistance, &
    outInterfaceSalinity, &
    outInterfaceTemperature, &
    outFreshwaterFlux, &
    outOceanHeatFlux, &
    outIceHeatFlux, &
    nCells) !{{{

    use shr_const_mod    , only: SHR_CONST_CPICE,  &
                                 SHR_CONST_CPSW,   &
                                 SHR_CONST_LATICE, &
                                 SHR_CONST_RHOICE, &
                                 SHR_CONST_RHOSW,  &
                                 SHR_CONST_DTF_DP, &
                                 SHR_CONST_DTF_DS, &
                                 SHR_CONST_TF0,    &
                                 SHR_CONST_KAPPA_LAND_ICE

    !-----------------------------------------------------------------
    !
    ! input variables
    !
    !-----------------------------------------------------------------

    real (kind=r8), dimension(:), intent(in) :: &
      oceanTemperature, &          !< Input: ocean temperature in top layer
      oceanSalinity, &             !< Input: ocean salinity in top layer
      oceanHeatTransferVelocity, & !< Input: ocean heat transfer velocity
      oceanSaltTransferVelocity, & !< Input: ocean salt transfer velocity
      interfacePressure, &         !< Input: pressure at the ice-ocean interface
      iceTemperature, &            !< Input: ice temperature in bottom layer
      iceTemperatureDistance       !< Input: distance to ice temperature from ice-ocean interface      

    integer, intent(in) :: nCells !< Input: number of cells in each array

    !-----------------------------------------------------------------
    !
    ! output variables
    !
    !-----------------------------------------------------------------

    real (kind=r8), dimension(:), intent(out) :: &
      outInterfaceSalinity, &    !< Output: ocean salinity at the interface
      outInterfaceTemperature, & !< Output: ice/ocean temperature at the interface
      outFreshwaterFlux, &   !< Output: ocean thickness flux (melt rate)
      outOceanHeatFlux, & !< Output: the temperature flux into the ocean
      outIceHeatFlux      !< Output: the temperature flux into the ice

    !-----------------------------------------------------------------
    !
    ! local variables
    !
    !-----------------------------------------------------------------

    real (kind=r8) :: T0, transferVelocityRatio, Tlatent, nu, a, b, c, eta, &
                         iceHeatFluxCoeff, iceDeltaT
    integer :: iCell
    character(*), parameter :: subname = '(compute_melt_fluxes)'
  
    Tlatent = SHR_CONST_LATICE/SHR_CONST_CPSW
    do iCell = 1, nCells

      iceHeatFluxCoeff = SHR_CONST_RHOICE*SHR_CONST_CPICE*SHR_CONST_KAPPA_LAND_ICE/iceTemperatureDistance(iCell)
      nu = iceHeatFluxCoeff/(SHR_CONST_RHOSW*SHR_CONST_CPSW*oceanHeatTransferVelocity(iCell))
      iceDeltaT = T0 - iceTemperature(iCell)

      T0 = SHR_CONST_TF0 + SHR_CONST_DTF_DP*interfacePressure(iCell)
      transferVelocityRatio = oceanSaltTransferVelocity(iCell)/oceanHeatTransferVelocity(iCell)

      a = -SHR_CONST_DTF_DS*(1.0_r8 + nu)
      b = transferVelocityRatio*Tlatent - nu*iceDeltaT + oceanTemperature(iCell) - T0
      c = -transferVelocityRatio*Tlatent

      ! a is strictly positive; c is strictly negative so we never get imaginary roots
      ! The positive root is the one we want (salinity is strictly positive)
      outInterfaceSalinity(iCell) = (-b + sqrt(b**2 - 4.0_r8*a*c*oceanSalinity(iCell)))/(2.0_r8*a)
      if (outInterfaceSalinity(iCell) .le. 0.0_r8) then
	write(logunit,*) subname,' ERROR: Negative interface salinity in subshelf boundary calculation: ', &
	                           'iCell,outInterfaceSalinity(iCell)=', iCell,outInterfaceSalinity(iCell)
        call shr_sys_abort(subname//' ERROR: Negative interface salinity in subshelf boundary calculation.')
      end if
      outInterfaceTemperature(iCell) = SHR_CONST_DTF_DS*outInterfaceSalinity(iCell)+T0

      outFreshwaterFlux(iCell) = SHR_CONST_RHOSW*oceanSaltTransferVelocity(iCell) &
        * (oceanSalinity(iCell)/outInterfaceSalinity(iCell) - 1.0_r8)

      ! According to Jenkins et al. (2001), the temperature fluxes into the ocean are:
      !   1. the advection of meltwater into the top layer (or removal for freezing)
      !   2. the turbulent transfer of heat across the boundary layer, based on the termal driving
      outOceanHeatFlux(iCell) = SHR_CONST_CPSW*(outFreshwaterFlux(iCell)*outInterfaceTemperature(iCell) &
        - SHR_CONST_RHOSW*oceanHeatTransferVelocity(iCell)*(oceanTemperature(iCell)-outInterfaceTemperature(iCell)))

      ! the temperature fluxes into the ice are:
      !   1. the advection of ice at the interface temperature out of the domain due to melting
      !      (or in due to freezing)
      !   2. the diffusion (if any) of heat into the ice, based on temperature difference between
      !      the reference point in the ice (either the surface or the middle of the bottom layer)
      !      and the interface
      outIceHeatFlux(iCell) = -SHR_CONST_CPICE*outFreshwaterFlux(iCell)*outInterfaceTemperature(iCell)

      outIceHeatFlux(iCell) = outIceHeatFlux(iCell) &
        - iceHeatFluxCoeff*(iceTemperature(iCell) - outInterfaceTemperature(iCell))

    end do

  !--------------------------------------------------------------------

  end subroutine compute_melt_fluxes !}}}
end module prep_glc_mod
