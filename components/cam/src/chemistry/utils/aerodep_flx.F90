!-------------------------------------------------------------------
! Manages reading and interpolation of prescribed aerosol deposition 
! fluxes.  These are the deposition fluxes sent to the surface.
!
! Created by: Francis Vitt
!-------------------------------------------------------------------
module aerodep_flx

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog
  use ppgrid,       only : pcols, pver, begchunk, endchunk

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: aerodep_flx_init
  public :: aerodep_flx_adv
  public :: aerodep_flx_readnl
  public :: aerodep_flx_prescribed

  logical :: has_aerodep_flx = .false.
  integer, parameter, public :: N_BULK = 14
  integer, parameter, public :: N_MODAL = 22
  integer :: number_flds

  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=32)  :: specifier(N_MODAL) = ' '

  ! for bulk aerosol fluxes

  character(len=12), parameter :: bulk_names(N_BULK) = (/ &
       'BCDEPWET    ', 'BCPHODRY    ', 'BCPHIDRY    ',  &
       'OCDEPWET    ', 'OCPHODRY    ', 'OCPHIDRY    ',  &
       'DSTX01DD    ', 'DSTX02DD    ', 'DSTX03DD    ', 'DSTX04DD    ', &
       'DSTX01WD    ', 'DSTX02WD    ', 'DSTX03WD    ', 'DSTX04WD    ' /)

  integer :: index_bulk_map(N_BULK)

  integer :: ibcphiwet,ibcphidry,ibcphodry
  integer :: iocphiwet,iocphidry,iocphodry

  integer :: idstdry1,idstdry2,idstdry3,idstdry4
  integer :: idstwet1,idstwet2,idstwet3,idstwet4

  ! for modal aerosol fluxes

  character(len=12), parameter :: modal_names(N_MODAL) = (/ &
       'bc_a1DDF    ', 'bc_c1DDF    ', 'pom_a1DDF   ', 'pom_c1DDF   ',  &
       'soa_a1DDF   ', 'soa_c1DDF   ', 'soa_a2DDF   ', 'soa_c2DDF   ',  &
       'dst_a1DDF   ', 'dst_c1DDF   ', 'dst_a3DDF   ', 'dst_c3DDF   ',  &
       'bc_a1SFWET  ', 'bc_c1SFWET  ', 'pom_a1SFWET ', 'pom_c1SFWET ',  &
       'soa_a1SFWET ', 'soa_c1SFWET ', 'dst_a1SFWET ', 'dst_c1SFWET ',  &
       'dst_a3SFWET ', 'dst_c3SFWET ' /)

  integer :: index_modal_map(N_MODAL)

  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4
  integer, parameter :: idx_dst1 = 5
  integer, parameter :: idx_dst3 = 6
  integer, parameter :: idx_ncl3 = 7
  integer, parameter :: idx_so43 = 8
  
  integer, parameter :: nmodal_idxs = 8

  integer :: idx_bc1_dryis = -1
  integer :: idx_bc1_drycw = -1
  integer :: idx_pom1_dryis = -1
  integer :: idx_pom1_drycw = -1
  integer :: idx_soa1_dryis = -1
  integer :: idx_soa1_drycw = -1
  integer :: idx_soa2_dryis = -1
  integer :: idx_soa2_drycw = -1
  integer :: idx_dst1_dryis = -1
  integer :: idx_dst1_drycw = -1
  integer :: idx_dst3_dryis = -1
  integer :: idx_dst3_drycw = -1

  integer :: idx_bc1_wetis = -1
  integer :: idx_bc1_wetcw = -1
  integer :: idx_pom1_wetis = -1
  integer :: idx_pom1_wetcw = -1
  integer :: idx_soa1_wetis = -1
  integer :: idx_soa1_wetcw = -1
  integer :: idx_dst1_wetis = -1
  integer :: idx_dst1_wetcw = -1
  integer :: idx_dst3_wetis = -1
  integer :: idx_dst3_wetcw = -1

  logical :: modal_fluxes = .false.

contains

!-------------------------------------------------------------------
! parses the list of dep fluxes specified in aerodep_flx_specifier namelist
! variable and sets up index variables
!-------------------------------------------------------------------
  subroutine aerodep_flx_init()
    
    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, horiz_only
    use physics_buffer, only : physics_buffer_desc
    use modal_aero_deposition, only : modal_aero_deposition_init
    
    implicit none

    integer :: ndx, istat, i

    if ( has_aerodep_flx ) then
       if ( masterproc ) then
          write(iulog,*) 'aero dep fluxes are prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       has_aerodep_flx = .false.
       if (masterproc) then
          write(iulog,*) 'aerodep_flx_init: no aerosol deposition fluxes have been specified'
       endif
       return
    end if

    index_bulk_map(:) = -1
    index_modal_map(:) = -1

    do i = 1,number_flds

       ndx = get_ndx( fields(i)%fldnam, bulk_names )
       if (ndx >0) then
          index_bulk_map(ndx) = i
       else
          ndx = get_ndx( fields(i)%fldnam, modal_names )
          if (ndx >0) then
             index_modal_map(ndx) = i
          endif
       endif
       if (ndx>0) then
          call addfld(trim(fields(i)%fldnam)//'_D', horiz_only, 'A',fields(i)%units, 'prescribed aero dep' )
       else
          call endrun('aerodep_flx_init: aerosol flux name not recognized: '//trim(fields(i)%fldnam))
       endif
    enddo

    modal_fluxes = any(index_modal_map(:)>0)

    if (modal_fluxes) then

       idx_bc1_dryis  = index_modal_map(1)
       idx_bc1_drycw  = index_modal_map(2)
       idx_pom1_dryis = index_modal_map(3)
       idx_pom1_drycw = index_modal_map(4)
       idx_soa1_dryis = index_modal_map(5)
       idx_soa1_drycw = index_modal_map(6)
       idx_soa2_dryis = index_modal_map(7)
       idx_soa2_drycw = index_modal_map(8)
       idx_dst1_dryis = index_modal_map(9)
       idx_dst1_drycw = index_modal_map(10)
       idx_dst3_dryis = index_modal_map(11)
       idx_dst3_drycw = index_modal_map(12)

       idx_bc1_wetis  = index_modal_map(13)
       idx_bc1_wetcw  = index_modal_map(14)
       idx_pom1_wetis = index_modal_map(15)
       idx_pom1_wetcw = index_modal_map(16)
       idx_soa1_wetis = index_modal_map(17)
       idx_soa1_wetcw = index_modal_map(18)
       idx_dst1_wetis = index_modal_map(19)
       idx_dst1_wetcw = index_modal_map(20)
       idx_dst3_wetis = index_modal_map(21)
       idx_dst3_wetcw = index_modal_map(22)

       call modal_aero_deposition_init( bc1_ndx=idx_bc1,   pom1_ndx=idx_pom1, soa1_ndx=idx_soa1, &
                                        soa2_ndx=idx_soa2, dst1_ndx=idx_dst1, dst3_ndx=idx_dst3, &
                                        ncl3_ndx=idx_ncl3, so43_ndx=idx_so43 )
    else

       ibcphiwet = index_bulk_map(1)
       ibcphodry = index_bulk_map(2)
       ibcphidry = index_bulk_map(3)
       iocphiwet = index_bulk_map(4)
       iocphodry = index_bulk_map(5)
       iocphidry = index_bulk_map(6)
       idstdry1  = index_bulk_map(7)
       idstdry2  = index_bulk_map(8)
       idstdry3  = index_bulk_map(9)
       idstdry4  = index_bulk_map(10)
       idstwet1  = index_bulk_map(11)
       idstwet2  = index_bulk_map(12)
       idstwet3  = index_bulk_map(13)
       idstwet4  = index_bulk_map(14)

    endif

  end subroutine aerodep_flx_init

!-------------------------------------------------------------------
! sets namelist options
!-------------------------------------------------------------------
subroutine aerodep_flx_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aerodep_flx_readnl'

   character(len=32)  :: aerodep_flx_specifier(N_MODAL)
   character(len=256) :: aerodep_flx_file
   character(len=256) :: aerodep_flx_filelist
   character(len=256) :: aerodep_flx_datapath
   character(len=32)  :: aerodep_flx_type
   logical            :: aerodep_flx_rmfile
   integer            :: aerodep_flx_cycle_yr
   integer            :: aerodep_flx_fixed_ymd
   integer            :: aerodep_flx_fixed_tod

   namelist /aerodep_flx_nl/ &
      aerodep_flx_specifier, &
      aerodep_flx_file,      &
      aerodep_flx_filelist,  &
      aerodep_flx_datapath,  &
      aerodep_flx_type,      &
      aerodep_flx_rmfile,    &
      aerodep_flx_cycle_yr,  &
      aerodep_flx_fixed_ymd, &
      aerodep_flx_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   aerodep_flx_specifier= specifier
   aerodep_flx_file     = filename
   aerodep_flx_filelist = filelist
   aerodep_flx_datapath = datapath
   aerodep_flx_type     = datatype
   aerodep_flx_rmfile   = rmv_file
   aerodep_flx_cycle_yr = cycle_yr
   aerodep_flx_fixed_ymd= fixed_ymd
   aerodep_flx_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aerodep_flx_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aerodep_flx_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(aerodep_flx_specifier,len(aerodep_flx_specifier(1))*N_MODAL,     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_file,     len(aerodep_flx_file),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_filelist, len(aerodep_flx_filelist), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_datapath, len(aerodep_flx_datapath), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_type,     len(aerodep_flx_type),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(aerodep_flx_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = aerodep_flx_specifier
   filename   = aerodep_flx_file
   filelist   = aerodep_flx_filelist
   datapath   = aerodep_flx_datapath
   datatype   = aerodep_flx_type
   rmv_file   = aerodep_flx_rmfile
   cycle_yr   = aerodep_flx_cycle_yr
   fixed_ymd  = aerodep_flx_fixed_ymd
   fixed_tod  = aerodep_flx_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_aerodep_flx = .true.

end subroutine aerodep_flx_readnl

!-------------------------------------------------------------------
! sets the aerosol deposition fluxes in the cam_out structure 
! to be sent to the surface models
!-------------------------------------------------------------------
  subroutine aerodep_flx_set( cam_out, ncol, lchnk )
    use camsrfexch,       only : cam_out_t     

    type(cam_out_t),     intent(inout) :: cam_out
    integer,             intent(in)    :: ncol, lchnk
    
    if( .not. has_aerodep_flx ) return
    
    if (modal_fluxes) then
       call set_modal_fluxes( cam_out, ncol, lchnk )
    else
       call set_bulk_fluxes( cam_out, ncol, lchnk )
    endif

  end subroutine aerodep_flx_set

!-------------------------------------------------------------------
! advances the prescribed fluxes to the current time step
!-------------------------------------------------------------------
  subroutine aerodep_flx_adv( state, pbuf2d, cam_out )

    use tracer_data,      only : advance_trcdata
    use physics_types,    only : physics_state
    use camsrfexch,       only : cam_out_t
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: c, ncol
    
    if( .not. has_aerodep_flx ) return

    call advance_trcdata( fields, file, state, pbuf2d  )

!$OMP PARALLEL DO PRIVATE (C, NCOL)
    do c = begchunk, endchunk
       ncol = state(c)%ncol
       call aerodep_flx_set( cam_out(c), ncol, c )
    enddo

  end subroutine aerodep_flx_adv

!-------------------------------------------------------------------
! returns true if aerosol dep fluxes are prescribed from dataset
!-------------------------------------------------------------------
  function aerodep_flx_prescribed()
    logical :: aerodep_flx_prescribed
    aerodep_flx_prescribed = has_aerodep_flx
  endfunction aerodep_flx_prescribed

! private methods
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_bulk_fluxes( cam_out, ncol, lchnk )
    use camsrfexch,            only : cam_out_t     

    ! Arguments
    type(cam_out_t), intent(inout) :: cam_out
    integer,         intent(in)    :: ncol, lchnk

    call set_fluxes( cam_out%bcphiwet, ibcphiwet, ncol, lchnk )
    call set_fluxes( cam_out%bcphidry, ibcphidry, ncol, lchnk )
    call set_fluxes( cam_out%bcphodry, ibcphodry, ncol, lchnk )

    call set_fluxes( cam_out%ocphiwet, iocphiwet, ncol, lchnk )
    call set_fluxes( cam_out%ocphidry, iocphidry, ncol, lchnk )
    call set_fluxes( cam_out%ocphodry, iocphodry, ncol, lchnk )

    call set_fluxes( cam_out%dstdry1, idstdry1, ncol, lchnk )
    call set_fluxes( cam_out%dstdry2, idstdry2, ncol, lchnk )
    call set_fluxes( cam_out%dstdry3, idstdry3, ncol, lchnk )
    call set_fluxes( cam_out%dstdry4, idstdry4, ncol, lchnk )

    call set_fluxes( cam_out%dstwet1, idstwet1, ncol, lchnk )
    call set_fluxes( cam_out%dstwet2, idstwet2, ncol, lchnk )
    call set_fluxes( cam_out%dstwet3, idstwet3, ncol, lchnk )
    call set_fluxes( cam_out%dstwet4, idstwet4, ncol, lchnk )

  end subroutine set_bulk_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_modal_fluxes( cam_out, ncol, lchnk )
    use camsrfexch,            only : cam_out_t     
    use modal_aero_deposition, only : set_srf_drydep, set_srf_wetdep

    ! Arguments
    type(cam_out_t), intent(inout) :: cam_out
    integer,         intent(in)    :: ncol, lchnk

    ! local vars
    integer :: i
    real(r8) :: aerdepdryis(pcols,nmodal_idxs)
    real(r8) :: aerdepdrycw(pcols,nmodal_idxs)
    real(r8) :: aerdepwetis(pcols,nmodal_idxs)
    real(r8) :: aerdepwetcw(pcols,nmodal_idxs)

    ! bin the fluxes as using modal_aero_deposition...

    aerdepdryis(:,:) = 0._r8
    aerdepdrycw(:,:) = 0._r8
    aerdepwetis(:,:) = 0._r8
    aerdepwetcw(:,:) = 0._r8

    call set_fluxes( aerdepwetis(:ncol,idx_bc1 ), idx_bc1_wetis , ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_bc1 ), idx_bc1_wetcw , ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_pom1), idx_pom1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_pom1), idx_pom1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_soa1), idx_soa1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_soa1), idx_soa1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst1), idx_dst1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst1), idx_dst1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst3), idx_dst3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst3), idx_dst3_wetcw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_bc1 ), idx_bc1_dryis , ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_bc1 ), idx_bc1_drycw , ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_pom1), idx_pom1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_pom1), idx_pom1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_soa1), idx_soa1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_soa1), idx_soa1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_soa2), idx_soa2_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_soa2), idx_soa2_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst1), idx_dst1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst1), idx_dst1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst3), idx_dst3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst3), idx_dst3_drycw, ncol, lchnk )

    call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

  end  subroutine set_modal_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_fluxes( fluxes, fld_indx, ncol, lchnk )
    use cam_history,  only : outfld

    real(r8), intent(inout) :: fluxes(:)
    integer,  intent(in)    :: fld_indx, ncol, lchnk

    integer :: i

    if (fld_indx<1) return

    do i = 1,ncol
       ! modal aero wet dep history fields are negative
       fluxes(i) = fields(fld_indx)%data(i,1,lchnk)
    enddo

    call outfld(trim(fields(fld_indx)%fldnam)//'_D', fluxes(:ncol), ncol, lchnk )

  endsubroutine set_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  integer function get_ndx( name, list )

    implicit none
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: list(:)

    integer :: i
    integer :: maxnum

    maxnum = size(list)

    get_ndx = -1
    do i = 1, maxnum
      if ( trim(name) == trim(list(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module aerodep_flx
