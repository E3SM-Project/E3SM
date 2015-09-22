! The glimmer_to_dycore module contains the Fortran side of the Glimmer-DyCore
! interface.  It uses the routines in dycore_to_glim_extern.cpp to create one
! or more instances of a dynamic core ice sheet model.  The dycore_model_index is
! the only parameter needed by glimmer_to_dycore subroutines to interact with a 
! specific instance of a dynamic core model. DMR--5/24/10   

module glimmer_to_dycore
  !*FD glimmer_to_dycore contains Fortran routines to couple Glimmer to a
  ! dynamic core model.
  use glide_types
  !use mpi_mod
  use parallel
  !use simple_forcing

  contains

  subroutine gtd_init_dycore_interface()
    call dycore_init_registry()
!print *,"Past dycore_init_registry"
  end subroutine gtd_init_dycore_interface

  subroutine gtd_delete_dycore_interface()
    call dycore_reset_registry()
  end subroutine gtd_delete_dycore_interface

  subroutine gtd_init_dycore(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index

    integer*4 error_code
    integer*4 dycore_type  ! 0=BISICLES, 1=Ymir
    character(8),DIMENSION(2) :: dycore_names = (/"BISICLES","Ymir    "/)

print *,"In gtd_init, dycore_model_index = ",dycore_model_index

!print *,'in init -- topg ndims,shape = ',size(shape(model%geometry%topg)),shape(model%geometry%topg)
    dycore_type = model%options%external_dycore_type
!print *,"In gtd_init_dycore, calling get_new_model"
    call dycore_get_new_model(dycore_type,dycore_model_index,error_code)
!print *,"In gtd_init_dycore, calling get_set_var routines"
    call gtd_set_geometry_vars(model,dycore_model_index)
!print *,"In gtd_init_dycore, past set_geometry_vars"
    call gtd_set_velocity_vars(model,dycore_model_index)
    call gtd_set_numerics_vars(model,dycore_model_index)
    call gtd_set_temper_vars(model,dycore_model_index)
    call gtd_set_climate_vars(model,dycore_model_index)
    call gtd_set_mpi_vars(model,dycore_model_index)

    print *,"In gtd_init_dycore, dycore_type, dycore_index  =  " , &
             dycore_names(dycore_type+1),dycore_model_index
    call dycore_init_model(dycore_type,dycore_model_index, &
            trim(model%options%dycore_input_file)//char(0),error_code)

  end subroutine gtd_init_dycore

  subroutine gtd_run_dycore(dycore_model_index,cur_time,time_inc)
    integer*4 dycore_model_index
    real(sp) cur_time, time_inc

    call dycore_run_model(dycore_model_index,cur_time,time_inc)
  end subroutine gtd_run_dycore

  subroutine gtd_delete_dycore(dycore_model_index)
    integer*4 dycore_model_index

    call dycore_delete_model(dycore_model_index)
  end subroutine gtd_delete_dycore

  subroutine gtd_set_dim_info(shape,dim_info)
    integer, dimension(:), intent(in) :: shape
    integer*8, dimension(:), intent(inout) :: dim_info 

    dim_info = 0
    dim_info(1) = size(shape)
    dim_info(2:1+dim_info(1)) = shape
  end subroutine gtd_set_dim_info

  subroutine gtd_set_geometry_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index

    integer*4 shape2, rank
    character*20 var_name
    character*20 dtype_name
    integer*4 var_name_len, dtype_name_len

    integer*8 dim_info(11)
    integer*8 dim_info2(2)
    integer*8 ewlbl, ewubl, nslbl, nsubl

!    print *,"In gtd_set_geometry_vars, dycore_model_index = ",dycore_model_index
    
!    print *,'thck ndims,shape = ',size(shape(model%geometry%thck)),shape(model%geometry%thck)
!    print *,'topg ndims,shape = ',size(shape(model%geometry%topg)),shape(model%geometry%topg)

!    print *,'usrf ndims,shape = ',size(shape(model%geometry%usrf)),shape(model%geometry%usrf)

    dtype_name = 'geometry'//char(0)
     
    var_name = 'thck'//char(0)
    !call gtd_set_dim_info(shape(model%geometry%thck),dim_info)
    call dycore_set_ptr_double_var(model%geometry%thck,var_name,dtype_name,dycore_model_index)

    var_name = 'topg'//char(0)
    !call gtd_set_dim_info(shape(model%geometry%topg),dim_info)
    call dycore_set_ptr_double_var(model%geometry%topg,var_name,dtype_name,dycore_model_index)

    var_name = 'usrf'//char(0)
    !call gtd_set_dim_info(shape(model%geometry%usrf),dim_info)
    call dycore_set_ptr_double_var(model%geometry%usrf,var_name,dtype_name,dycore_model_index)


    print *,"this_rank, ewlb, ewub, nslb, nsub", this_rank,  ewlb, ewub, nslb, nsub
 
! (DFM -2/12/13) since ewlb, et al contain local grid info, use dim_info to 
! pass in global index space info
    dim_info(1) = 3
    dim_info(2) = model%general%upn
    dim_info(3) = global_ewn
    dim_info(4) = global_nsn 


!    dtype_name = 'geometry'
!    dtype_name_len = 8

    ! use age to get dim_info for now (only 3d var in geometry derived type)
!    call gtd_set_dim_info(shape(model%geometry%age),dim_info)

    print *, "dim_info = ", dim_info(1), dim_info(2), dim_info(3), dim_info(4)

    var_name = 'dimInfo'//char(0)
    dim_info2(1) = 1
    dim_info2(2) = dim_info(1) + 1
    call dycore_copy_in_long_var(dim_info,var_name,dtype_name,dim_info2, dycore_model_index)

    ewlbl = ewlb
    ewubl = ewub
    nslbl = nslb
    nsubl = nsub

    dim_info2(1) = 1
    dim_info2(2) = 1
    var_name = 'ewlb'//char(0)
    call dycore_copy_in_long_var(ewlbl,var_name,dtype_name,dim_info2, dycore_model_index)
    var_name = 'ewub'//char(0)
    call dycore_copy_in_long_var(ewubl,var_name,dtype_name,dim_info2, dycore_model_index)
    var_name = 'nslb'//char(0)
    call dycore_copy_in_long_var(nslbl,var_name,dtype_name,dim_info2, dycore_model_index)
    var_name = 'nsub'//char(0)
    call dycore_copy_in_long_var(nsubl,var_name,dtype_name,dim_info2, dycore_model_index)

!    print *,"leaving gtd_set_geometry_vars, dim_info =  ",dim_info
  end subroutine gtd_set_geometry_vars 


  subroutine gtd_set_velocity_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index
    
    character*20 var_name
    character*20 dtype_name
    integer*4 var_name_len, dtype_name_len

    integer*8 dim_info(11)
    integer*8 dim_info2(2)

!    print *,"In copy_velocity_vars, dycore_model_index = ",dycore_model_index

    dtype_name = 'velocity'//char(0)

    print *,'uvel ndims,shape = ',size(shape(model%velocity%uvel)),shape(model%velocity%uvel)

    print *,'vvel ndims,shape = ',size(shape(model%velocity%vvel)),shape(model%velocity%vvel)

    print *,'wvel ndims,shape = ',size(shape(model%velocity%wvel)),shape(model%velocity%wvel)


    var_name = 'uvel'//char(0)       
    call dycore_set_ptr_double_var(model%velocity%uvel,var_name, &
                              dtype_name,dycore_model_index);
    var_name = 'vvel'//char(0)
    call dycore_set_ptr_double_var(model%velocity%vvel,var_name, &
                              dtype_name,dycore_model_index);
    var_name = 'wvel'//char(0)
    call dycore_set_ptr_double_var(model%velocity%wvel,var_name, &
                              dtype_name,dycore_model_index);
    var_name = 'wgrd'//char(0)
    call dycore_set_ptr_double_var(model%velocity%wgrd,var_name, &
                              dtype_name,dycore_model_index);

!    print *,'beta ndims,shape = ',size(shape(model%velocity%beta)),shape(model%velocity%beta)

    var_name = 'btrc'//char(0)
    call dycore_set_ptr_double_var(model%velocity%beta,var_name, &
                                  dtype_name,dycore_model_index);

    call gtd_set_dim_info(shape(model%velocity%uvel),dim_info)

    var_name = 'dimInfo'//char(0)
    dim_info2(1) = 1
    dim_info2(2) = 4
    call dycore_copy_in_long_var(dim_info,var_name,dtype_name,dim_info2,dycore_model_index)
  end subroutine gtd_set_velocity_vars  

  subroutine gtd_set_numerics_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index
    
    character*20 var_name
    character*20 dtype_name
    integer*4 var_name_len, dtype_name_len
    integer*8 dim_info2(2)

    dtype_name = 'numerics'//char(0)

    dim_info2(1) = 1
    dim_info2(2) = 1

    var_name = 'dew'//char(0)    
    call dycore_copy_in_double_var(model%numerics%dew,var_name,dtype_name,dim_info2,dycore_model_index)
    var_name = 'dns'//char(0)
    call dycore_copy_in_double_var(model%numerics%dns,var_name,dtype_name,dim_info2,dycore_model_index)

  end subroutine gtd_set_numerics_vars

  subroutine gtd_set_temper_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index
    character*20 var_name
    character*20 dtype_name

    integer*8 dim_info(11), dim_info2(2)
    
    dtype_name = 'temper'//char(0)

    var_name = 'temp'//char(0)       
    call dycore_set_ptr_double_var(model%temper%temp,var_name,dtype_name,dycore_model_index)

    var_name = 'bheatflx'//char(0)       
    call dycore_set_ptr_double_var(model%temper%bheatflx,var_name,dtype_name,dycore_model_index)

    var_name = 'bmlt'//char(0)       
    call dycore_set_ptr_double_var(model%temper%bmlt,var_name,dtype_name,dycore_model_index)
      
    print *,'temp ndims,shape = ',size(shape(model%temper%temp)),shape(model%temper%temp)

    print *,'bheatflx ndims,shape = ',size(shape(model%temper%bheatflx)),shape(model%temper%bheatflx)

    print *,'bmlt ndims,shape = ',size(shape(model%temper%bmlt)),shape(model%temper%bmlt)

    call gtd_set_dim_info(shape(model%temper%temp),dim_info)

    var_name = 'dimInfo'//char(0)
    dim_info2(1) = 1
    dim_info2(2) = dim_info(1) + 1
    call dycore_copy_in_long_var(dim_info,var_name,dtype_name,dim_info2,dycore_model_index)    
  end subroutine gtd_set_temper_vars

  subroutine gtd_set_climate_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index
    character*20 var_name
    character*20 dtype_name

    integer*8 dim_info(11), dim_info2(2)
    
    dtype_name = 'climate'//char(0)

    var_name = 'acab'//char(0)       
    call dycore_set_ptr_double_var(model%climate%acab,var_name,dtype_name,dycore_model_index)
    var_name = 'acab_tavg'//char(0)       
    call dycore_set_ptr_double_var(model%climate%acab_tavg,var_name,dtype_name,dycore_model_index)
    var_name = 'calving'//char(0)       
    call dycore_set_ptr_double_var(model%climate%calving,var_name,dtype_name,dycore_model_index)
      
    call gtd_set_dim_info(shape(model%climate%acab),dim_info)
    ! print *,"In climate set, dim_info: ",dim_info
    var_name = 'dimInfo'//char(0)
    dim_info2(1) = 1
    dim_info2(2) = dim_info(1) + 1
    call dycore_copy_in_long_var(dim_info,var_name,dtype_name,dim_info2,dycore_model_index)    

    var_name = 'eus'
    dim_info2(1) = 1
    dim_info2(2) = 1
    ! eus parm isn't being set during initialization, so commented out here:
    !call dycore_copy_in_double_var(model%climate%eus,var_name,dtype_name,dim_info2,dycore_model_index)
    !print *,"eus: ",model%climate%eus

  end subroutine gtd_set_climate_vars
  
  subroutine gtd_set_mpi_vars(model,dycore_model_index)
    type(glide_global_type) :: model
    integer*4 dycore_model_index
    character*20 var_name
    character*20 dtype_name

    integer*8 dim_info(11), dim_info2(2)

    ! integer,save :: comm, tasks, this_rank -- from parallel_mpi.F90
    integer*8 communicator, process_count, my_rank


    communicator = comm
    process_count = tasks
    my_rank = this_rank
      
    dtype_name = 'mpi_vars'//char(0)

    dim_info2(1) = 1
    dim_info2(2) = 1
    var_name = 'communicator'//char(0)
    call dycore_copy_in_long_var(communicator,var_name,dtype_name,dim_info2, dycore_model_index)
    var_name = 'process_count'//char(0)
    call dycore_copy_in_long_var(process_count,var_name,dtype_name,dim_info2, dycore_model_index)
    var_name = 'my_rank'//char(0)
    call dycore_copy_in_long_var(my_rank,var_name,dtype_name,dim_info2, dycore_model_index)

  end subroutine gtd_set_mpi_vars

end module glimmer_to_dycore
