    ! main particle type (used by particle_framework)
    type mpas_particle_type
       ! IO block number
       !integer :: ioBlock = 0
       ! pool containing particle data needed for communication
       ! between computational cells
       type (mpas_pool_type), pointer :: haloDataPool => NULL()
       ! pool containing particle data not needing to be passed
       ! between processors
       type (mpas_pool_type), pointer :: nonhaloDataPool => NULL()
    end type mpas_particle_type

    ! just a list of particles (main unit that code operates on)
    type mpas_particle_list_type
       type (mpas_particle_type), pointer :: particle => NULL()
       ! next item in the list
       type (mpas_particle_list_type), pointer :: next => NULL()
       ! previous item in the list
       type (mpas_particle_list_type), pointer :: prev => NULL()
    end type mpas_particle_list_type

    type mpas_list_of_particle_list_type
       type (mpas_particle_list_type), pointer :: list => NULL()
    end type mpas_list_of_particle_list_type
