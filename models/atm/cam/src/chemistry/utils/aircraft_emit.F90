module aircraft_emit
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Manages reading and interpolation of aircraft aerosols
! 
! Authors: Chih-Chieh (Jack) Chen and Cheryl Craig -- February 2010
! 
!-----------------------------------------------------------------------
  use perf_mod,     only : t_startf, t_stopf

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog

  implicit none
  private
  save 

  public :: aircraft_emit_init
  public :: aircraft_emit_adv
  public :: aircraft_emit_register
  public :: aircraft_emit_readnl

  type :: forcing_air
     real(r8)              :: mw
     character(len=256) :: filelist
     character(len=256) :: filename
     real(r8), pointer     :: times(:)
     real(r8), pointer     :: levi(:)
     character(len=11)  :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld),pointer       :: fields(:)
     type(trfile)              :: file
  end type forcing_air

  type(forcing_air),allocatable :: forcings_air(:)

  integer, parameter :: N_AERO = 10 
  character(len=11)    :: aero_names(N_AERO) = (/'ac_HC      ','ac_NOX     ','ac_PMNV    ',&
                          'ac_PMSO    ','ac_PMFO    ','ac_FUELBURN','ac_CO2     ','ac_H2O     ',&
                          'ac_SOX     ','ac_CO      '/)

  real(r8), parameter :: molmass(N_AERO) = 1._r8

  logical :: advective_tracer(N_AERO) = (/.false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false.,.false./)
  character(len=3) :: mixtype(N_AERO) = (/'wet','wet','wet','wet','wet','wet','wet','wet','wet','wet'/)

  real(r8) :: cptmp = 666.0_r8
  real(r8) :: qmin = 0.0_r8
  logical :: readiv = .false.
  logical :: has_fixed_ubs = .false.
  logical :: cam_outfld = .false.

  integer            :: index_map(N_AERO)
  character(len=256) :: air_specifier(N_AERO)=''
  character(len=24)  :: air_type = 'CYCLICAL_LIST' ! 'CYCLICAL_LIST'



  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''

  logical            :: rmv_file = .false.

  integer :: number_flds

  integer :: aircraft_cnt = 0
  character(len=16) :: spc_name_list(N_AERO)
  character(len=256) :: spc_flist(N_AERO),spc_fname(N_AERO)

contains

  subroutine aircraft_emit_register()

!------------------------------------------------------------------
! **** Add the aircraft aerosol data to the physics buffer ****
!------------------------------------------------------------------
    use ppgrid,         only: pver,pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8
    use tracer_data,    only: incr_filename
    use constituents,   only: cnst_add

    integer :: i,idx, mm, ind, n
    character(len=16) :: spc_name
    character(len=256) :: filelist, curr_filename
    character(len=128) :: long_name
    logical            :: has_fixed_ubc=.false.
    logical            :: read_iv=.false.

    !------------------------------------------------------------------
    ! Return if air_specifier is blank (no aircraft data to process)
    !------------------------------------------------------------------
    if (air_specifier(1) == "") return

! count aircraft emission species used in the simulation
    count_emis: do n=1,N_AERO
	
	if( len_trim(air_specifier(n) ) == 0 ) then
	    exit count_emis
        endif

        i = scan(air_specifier(n),'->')
        spc_name = trim(adjustl(air_specifier(n)(:i-1)))
        filelist = trim(adjustl(air_specifier(n)(i+2:)))

        mm = get_aircraft_ndx(spc_name)
        if( mm < 1 ) then
	 call endrun('aircraft_emit_register: '//trim(spc_name)//' is not in the aircraft emission dataset')
        endif

        aircraft_cnt = aircraft_cnt + 1
        call pbuf_add_field(aero_names(mm),'physpkg',dtype_r8,(/pcols,pver/),idx)

        spc_flist(aircraft_cnt) = filelist
        spc_name_list(aircraft_cnt) = spc_name
        index_map(aircraft_cnt) = mm

        curr_filename=''	
        datapath=''
        spc_fname(aircraft_cnt) = incr_filename( curr_filename, filenames_list=spc_flist(aircraft_cnt), datapath=datapath)

        if( advective_tracer(mm) ) then
          long_name = 'aircraft_'//trim(spc_name)
          call cnst_add(aero_names(mm),molmass(mm),cptmp,qmin,ind,longname=long_name,readiv=read_iv, &
                        mixtype=mixtype(mm),cam_outfld=cam_outfld,fixed_ubc=has_fixed_ubc)
        endif

    enddo count_emis
! count aircraft emission species used in the simulation
  
  endsubroutine aircraft_emit_register

  subroutine aircraft_emit_init()
!-------------------------------------------------------------------
! **** Initialize the aircraft aerosol data handling ****
!-------------------------------------------------------------------
    use cam_history,    only: addfld, phys_decomp, add_default
    use ppgrid,         only: pver
    use tracer_data,    only: trcdata_init
    use physics_types,  only: physics_state
    use ppgrid,         only: begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    character(len=16)  :: spc_name

    integer :: ndx, istat, i, astat, m, n, mm, c
    
    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------
    if (aircraft_cnt == 0 ) return

    if (masterproc) write(iulog,*) ' '

    !-----------------------------------------------------------------------
    !       allocate forcings type array
    !-----------------------------------------------------------------------
    allocate( forcings_air(aircraft_cnt), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'aircraft_emit_init: failed to allocate forcings_air array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    !       setup the forcings_air type array
    !-----------------------------------------------------------------------
    species_loop : do m = 1,aircraft_cnt

          allocate( forcings_air(m)%sectors(1), stat=astat )
          if( astat/= 0 ) then
             write(iulog,*) 'aircraft_emit_init: failed to allocate forcings_air%sectors array; error = ',astat
             call endrun
          end if

          allocate( forcings_air(m)%fields(1), stat=astat )
          if( astat/= 0 ) then
             write(iulog,*) 'aircraft_emit_init: failed to allocate forcings_air%fields array; error = ',astat
             call endrun
          end if

          spc_name = spc_name_list(m)
          !-----------------------------------------------------------------------
          !         default settings
          !-----------------------------------------------------------------------
          forcings_air(m)%file%stepTime    = .true.  ! Aircraft data is not to be interpolated in time
          forcings_air(m)%file%cyclical_list    = .true.  ! Aircraft data cycles over the filename list
          forcings_air(m)%file%weight_by_lat     = .true.  ! Aircraft data -  interpolated with latitude weighting
          forcings_air(m)%file%conserve_column = .true. ! Aircraft data - vertically interpolated to conserve the total column
          forcings_air(m)%species          = spc_name
          forcings_air(m)%sectors          = spc_name ! Only one species per file for aircraft data
          forcings_air(m)%nsectors         = 1
          forcings_air(m)%filelist         = spc_flist(m)
!         forcings_air(m)%file%curr_filename    = spc_fname(m)
          forcings_air(m)%filename         = spc_fname(m)
          call addfld( trim(spc_name),  '1/s', pver, 'A', &
                       'aircraft emission '//trim(spc_name),   phys_decomp )
          call add_default( trim(spc_name), 1, ' ' )
    end do species_loop

    if (masterproc) then
       !-----------------------------------------------------------------------
       !            diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'aircraft_emit_init: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'aircraft_emit timing specs'
       write(iulog,*) 'type = ',air_type
       write(iulog,*) ' '
       write(iulog,*) 'there are ',aircraft_cnt,' species of aircraft emission'
       do m = 1,aircraft_cnt
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',m
          write(iulog,*) 'species = ',trim(forcings_air(m)%species)
          write(iulog,*) 'filelist= ',trim(forcings_air(m)%filelist)
       end do
       write(iulog,*) ' '
    endif

    !------------------------------------------------------------------
    !       Initialize the aircraft file processing
    !------------------------------------------------------------------
    do m=1,aircraft_cnt

       allocate (forcings_air(m)%file%in_pbuf(size(forcings_air(m)%sectors)))
       forcings_air(m)%file%in_pbuf(:) = .true.
       call trcdata_init( forcings_air(m)%sectors, forcings_air(m)%filename, forcings_air(m)%filelist, datapath, &
                          forcings_air(m)%fields, forcings_air(m)%file, rmv_file, 0, 0, 0, air_type)
        

       number_flds = 0
       if (associated(forcings_air(m)%fields)) number_flds = size( forcings_air(m)%fields )

       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(iulog,*) 'There are no aircraft aerosols'
             write(iulog,*) ' '
             call endrun
          endif
       end if
   end do


  end subroutine aircraft_emit_init



  subroutine aircraft_emit_adv( state, pbuf2d)
!-------------------------------------------------------------------
! **** Advance to the next aircraft data ****
!-------------------------------------------------------------------

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry       ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz                ! J/K/molecule
! C.-C. Chen
    use physconst,    only : gravit, rearth
    use phys_grid,    only : get_wght_all_p
    
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_field, pbuf_get_chunk
  
    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    integer :: ind,c,ncol,i,caseid,m,n
    real(r8) :: to_mmr(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)
   
! C.-C. Chen
    real(r8) :: wght(pcols)

   !------------------------------------------------------------------
   ! Return if aircraft_cnt is zero (no aircraft data to process)
   !------------------------------------------------------------------
    if (aircraft_cnt == 0 ) return
    call t_startf('All_aircraft_emit_adv')

   !-------------------------------------------------------------------
   !    For each field, read more data if needed and interpolate it to the current model time
   !-------------------------------------------------------------------
    do m = 1, aircraft_cnt
       call advance_trcdata( forcings_air(m)%fields, forcings_air(m)%file, state, pbuf2d)
    
   !-------------------------------------------------------------------
   !    set the tracer fields with the correct units
   !-------------------------------------------------------------------
       do i = 1,number_flds

! C.-C. Chen, adding case 4  for kg/sec
          select case ( to_lower(trim(forcings_air(m)%fields(i)%units(:GLC(forcings_air(m)%fields(i)%units)))) )
          case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
             caseid = 1
          case ('kg/kg','mmr')
             caseid = 2
          case ('mol/mol','mole/mole','vmr','fraction')
             caseid = 3
          case ('kg/kg/sec')
              caseid = 4
          case default
             print*, 'aircraft_emit_adv: units = ',trim(forcings_air(m)%fields(i)%units) ,' are not recognized'
             call endrun('aircraft_emit_adv: units are not recognized')
          end select

          ind = index_map(i)

!$OMP PARALLEL DO PRIVATE (C, NCOL, TO_MMR, tmpptr, pbuf_chnk)
          do c = begchunk,endchunk
             ncol = state(c)%ncol

! C.-C. Chen, turning emission data to mixing ratio
             call get_wght_all_p(c,ncol,wght(:ncol))

             if (caseid == 1) then
                to_mmr(:ncol,:) = (molmass(ind)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
             elseif (caseid == 2) then
                to_mmr(:ncol,:) = 1._r8
             elseif (caseid == 4) then
!                do n=1,ncol
!                  to_mmr(n,:) = 1.0_r8/(rearth*rearth*wght(n)*state(c)%pdel(n,:)/gravit)
!                end do
                to_mmr(:ncol,:) = 1.0_r8
             else
                to_mmr(:ncol,:) = molmass(ind)/mwdry
             endif
             pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
             call pbuf_get_field(pbuf_chnk, forcings_air(m)%fields(i)%pbuf_ndx, tmpptr )
   
             tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)

             call outfld( forcings_air(m)%fields(i)%fldnam, &
                  tmpptr, ncol, state(c)%lchnk )

          enddo
       enddo
    enddo

    call t_stopf('All_aircraft_emit_adv')
  end subroutine aircraft_emit_adv

  subroutine aircraft_emit_readnl(nlfile)
!-------------------------------------------------------------------
! **** Read in the aircraft_emit namelist *****
!-------------------------------------------------------------------
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aircraft_emit_readnl'

   character(len=256) :: aircraft_specifier(N_AERO)
   character(len=24)  :: aircraft_type 

   namelist /aircraft_emit_nl/  aircraft_specifier, aircraft_type
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   aircraft_specifier= air_specifier
   aircraft_type     = air_type

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aircraft_emit_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aircraft_emit_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(aircraft_specifier,len(aircraft_specifier(1))*N_AERO,     mpichar, 0, mpicom)
   call mpibcast(aircraft_type,     len(aircraft_type),                    mpichar, 0, mpicom)
#endif

   ! Update module variables with user settings.
   air_specifier  = aircraft_specifier
   air_type       = aircraft_type

 end subroutine aircraft_emit_readnl

 integer function get_aircraft_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_aircraft_ndx = 0
    do i = 1,N_AERO
      if ( trim(name) == trim(aero_names(i)) ) then
        get_aircraft_ndx = i
        return
      endif
    enddo

  end function get_aircraft_ndx

end module aircraft_emit
