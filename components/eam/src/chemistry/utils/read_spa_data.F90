module read_spa_data

use tracer_data

implicit none

private

public :: read_spa_data_init, read_spa_data_adv, read_spa_data_register

type(trfld), pointer :: spa_fields(:)
type(trfile)         :: spa_file
logical              :: rmv_file = .false.

!fields to read from the file
character(len=16), parameter :: pbuf_names(1) = (/'AER_G_SW_band1  '/)

contains

  !-------------------------------------------------------------------
  ! registers spa fields to the phys buffer
  !-------------------------------------------------------------------
  subroutine read_spa_data_register()

    use ppgrid,         only: pver,pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: i,idx

    do i = 1,size(pbuf_names)
       call pbuf_add_field(pbuf_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
    enddo
  endsubroutine read_spa_data_register


  subroutine read_spa_data_init

    character(len=32)  :: specifier(1)

    specifier(1) = 'AER_G_SW_band1'

    !Tracer data routine init
    allocate (spa_file%in_pbuf(size(specifier)))
    spa_file%in_pbuf(:) = .true.

    call trcdata_init( specifier, 'bsingh_spa_file_lat_lon_11.nc', '', '/compyfs/sing201/lat_lon', spa_fields, spa_file, &
         rmv_file, 1, 0, 0, 'CYCLICAL')

  end subroutine read_spa_data_init

  subroutine read_spa_data_adv( state, pbuf2d )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk, &
         pbuf_get_index

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    call advance_trcdata( spa_fields, spa_file, state, pbuf2d )


  end subroutine read_spa_data_adv
end module read_spa_data
