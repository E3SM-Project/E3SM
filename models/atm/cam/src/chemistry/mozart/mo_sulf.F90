

      module mo_sulf
!---------------------------------------------------------------
!	... Annual cycle for sulfur
!---------------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8

      use abortutils,   only : endrun
      use cam_logfile,  only : iulog
      use tracer_data,  only : trfld,trfile
      use physics_types,only : physics_state
      use ppgrid,       only : begchunk, endchunk
      use physics_buffer, only : physics_buffer_desc
      use ppgrid, only : pcols, pver

      implicit none

      private
      public  :: sulf_inti, set_sulf_time, sulf_interp

      save

      type(trfld), pointer :: fields(:) => null()
      type(trfile) :: file

      logical :: read_sulf = .false.

      contains 

      subroutine sulf_inti( sulf_file )
!-----------------------------------------------------------------------
! 	... Open netCDF file containing annual sulfur data.  Initialize
!           arrays with the data to be interpolated to the current time.
!
!           It is assumed that the time coordinate is increasing
!           and represents calendar days; range = [1.,366.).
!-----------------------------------------------------------------------
      use spmd_utils,    only : masterproc
      use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx
      use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, interp_type
      use tracer_data,   only : trcdata_init
      use cam_history,   only : addfld, phys_decomp

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: sulf_file

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: ndxs(5), so4_ndx

      character(len=1), parameter :: filename = ' '
      character(len=1), parameter :: filelist = ' '
      character(len=1), parameter :: datapath = ' '
      character(len=8), parameter :: datatype = 'CYCLICAL'
      logical         , parameter :: rmv_file = .false.
      integer         , parameter :: cycle_yr  = 0
      integer         , parameter :: fixed_ymd = 0
      integer         , parameter :: fixed_tod = 0
      character(len=8), parameter :: fld_names(1) = (/'SULFATE '/)

      ndxs(1) = get_rxt_ndx( 'usr_N2O5_aer' )
      ndxs(2) = get_rxt_ndx( 'usr_NO3_aer' )
      ndxs(3) = get_rxt_ndx( 'usr_NO2_aer' )
      ndxs(4) = get_rxt_ndx( 'usr_HO2_aer' )
      ndxs(5) = get_rxt_ndx( 'het1' )
      so4_ndx = get_spc_ndx('SO4')

      read_sulf = any( ndxs > 0) .and. (so4_ndx < 0)

      if ( .not. read_sulf ) return

      allocate(file%in_pbuf(size(fld_names)))
      file%in_pbuf(:) = .false. 
      call trcdata_init( fld_names, sulf_file, filelist, datapath, fields, file, &
           rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

      call addfld('SULFATE','VMR', pver, 'I', 'sulfate data', phys_decomp )

      end subroutine sulf_inti

      subroutine set_sulf_time( pbuf2d, state )
!--------------------------------------------------------------------
!	... Check and set time interpolation indicies
!--------------------------------------------------------------------
      use tracer_data,  only : advance_trcdata

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
      type(physics_state), intent(in):: state(begchunk:endchunk)                 

      if ( .not. read_sulf ) return

      call advance_trcdata( fields, file, state, pbuf2d  )

      end subroutine set_sulf_time

      subroutine sulf_interp( ncol, lchnk, ccm_sulf )
!-----------------------------------------------------------------------
! 	... Time interpolate sulfatei to current time
!-----------------------------------------------------------------------
      use cam_history,  only : outfld

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)   :: ncol              ! columns in chunk
      integer, intent(in)   :: lchnk             ! chunk number
      real(r8), intent(out) :: ccm_sulf(:,:)     ! output sulfate

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------

      ccm_sulf(:,:) = 0._r8

      if ( .not. read_sulf ) return

      ccm_sulf(:ncol,:) = fields(1)%data(:ncol,:,lchnk)

      call outfld( 'SULFATE', ccm_sulf(:ncol,:), ncol, lchnk )

      end subroutine sulf_interp

      end module mo_sulf
