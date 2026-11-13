!
! topo smoother uses HOMME's laplace operator
! for topo smoothign, we run the iterations at ~ viscous CFL = 0.5
!
! for RRM grids (could also be used for cubed sphere grids) we use
! the tensor HV which has built in resolution scalings:
!   ne=0
!   meshfile = path_to_exodus_grid.g
!   smooth_phis_numcycle=12   ! recommended
!   hypervis_scalings=2       ! turn on Tensor with dx^2 scalings
!   smooth_phis_nudt=4e-16    ! valid for all resolutions
!
! For cubed-sphere grids:
!   ne=30
!   smooth_phis_numcycle=16  ! RECOMMENDED.  controls amount of smoothing
!   hypervis_scaling = 0     ! turn on constant coefficient HV
!   smooth_phis_nudt =  28e7 * ( 30/NE)**2
!
!   precomputed hypervis_scaling=0 values:
!     NE         smooth_phis_nudt
!     16         98e7
!     30         28e7
!     60         7e7
!     120        18e6
!     240        44e5
!
!
&ctl_nl
ne = 4
mesh_file = 'none'
smooth_phis_numcycle = 16
smooth_phis_nudt = 1.6e10    ! resolution dependent
hypervis_scaling = 0         ! 2 for RRM grids
/

&vert_nl
/

&analysis_nl
tool = 'topo_gll_to_smoothed'
infilenames = 'movies/phis-baroclinic1.nc'
output_dir = "movies/"
output_timeunits=1
output_frequency=1
output_varnames1='PHIS'  ! homme will output goes. rename afterwards
output_type='netcdf'
io_stride = 16
/
