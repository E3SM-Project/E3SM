module read_spa_data

!BSINGH - I didn't add any comments here due to time constrain but I will add comments later

use tracer_data

implicit none

private

public :: read_spa_data_init, read_spa_data_adv, read_spa_data_register

type(trfld), pointer :: spa_fields(:)
type(trfile)         :: spa_file
logical              :: rmv_file = .false.

!fields to read from the file
character(len=16), parameter :: pbuf_names(59) = &
     [ 'AER_G_SW_0   ', 'AER_G_SW_1   ', 'AER_G_SW_2   ', 'AER_G_SW_3   ', 'AER_G_SW_4   ', &
       'AER_G_SW_5   ', 'AER_G_SW_6   ', 'AER_G_SW_7   ', 'AER_G_SW_8   ', 'AER_G_SW_9   ', &
       'AER_G_SW_10  ', 'AER_G_SW_11  ', 'AER_G_SW_12  ', 'AER_G_SW_13  ', &
       'AER_SSA_SW_0 ', 'AER_SSA_SW_1 ', 'AER_SSA_SW_2 ', 'AER_SSA_SW_3 ', 'AER_SSA_SW_4 ', &
       'AER_SSA_SW_5 ', 'AER_SSA_SW_6 ', 'AER_SSA_SW_7 ', 'AER_SSA_SW_8 ', 'AER_SSA_SW_9 ', &
       'AER_SSA_SW_10', 'AER_SSA_SW_11', 'AER_SSA_SW_12', 'AER_SSA_SW_13', &
       'AER_TAU_LW_0 ', 'AER_TAU_LW_1 ', 'AER_TAU_LW_2 ', 'AER_TAU_LW_3 ', 'AER_TAU_LW_4 ', &
       'AER_TAU_LW_5 ', 'AER_TAU_LW_6 ', 'AER_TAU_LW_7 ', 'AER_TAU_LW_8 ', 'AER_TAU_LW_9 ', &
       'AER_TAU_LW_10', 'AER_TAU_LW_11', 'AER_TAU_LW_12', 'AER_TAU_LW_13', 'AER_TAU_LW_14', &
       'AER_TAU_LW_15', &
       'AER_TAU_SW_0 ', 'AER_TAU_SW_1 ', 'AER_TAU_SW_2 ', 'AER_TAU_SW_3 ', 'AER_TAU_SW_4 ', &
       'AER_TAU_SW_5 ', 'AER_TAU_SW_6 ', 'AER_TAU_SW_7 ', 'AER_TAU_SW_8 ', 'AER_TAU_SW_9 ', &
       'AER_TAU_SW_10', 'AER_TAU_SW_11', 'AER_TAU_SW_12','AER_TAU_SW_13', &
       'CCN3         ']

character(len=16), parameter :: specifier(59) = pbuf_names(59)


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

    !Tracer data routine init
    allocate (spa_file%in_pbuf(size(specifier)))
    spa_file%in_pbuf(:) = .true.

!    call trcdata_init( specifier, 'unfied_SPA_file_lat_lon.nc', '', '/compyfs/sing201/lat_lon', spa_fields, spa_file, &
!         rmv_file, 1, 0, 0, 'CYCLICAL')
    call trcdata_init( specifier, 'unfied_SPA_file_lat_lon.nc', '', '/global/cscratch1/sd/bsingh/users/hassan', spa_fields, spa_file, &
         rmv_file, 1, 0, 0, 'CYCLICAL')


  end subroutine read_spa_data_init


  subroutine read_spa_data_adv( state, pbuf2d )

    !advance fields in time and interpolate (space and time)
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
