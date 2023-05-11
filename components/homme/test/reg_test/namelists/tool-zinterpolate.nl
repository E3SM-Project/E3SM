&ctl_nl                                                                                                           
ne = 4
mesh_file='none'
/                                                                                                                 

&vert_nl
! for E3SM output, set this correctly:
!vanalytic=0
!vfile_mid     = "/home/mataylo/codes/homme/test/vcoord/acme-72m.ascii"
!vfile_int     = "/home/mataylo/codes/homme/test/vcoord/acme-72i.ascii"
/                                                                                                                 

&analysis_nl                                                                                                      
tool = 'interpolate_tool'
interpolate_analysis=.true.
infilenames = "movies/phis-smoothed1.nc"
output_dir = "./"                                                                                                 
output_timeunits=1                                                                                                
output_frequency=1                                                                                                
output_varnames1='all'
output_type='netcdf'                                                                                              
interp_gridtype=1   ! 1=equi-angle cap grid.  2=gauss.  3=equ-angle offset(no poles)
interp_type=1       ! 0=native SE, 1=bilinear
io_stride = 16                                                                                                    
/                                      