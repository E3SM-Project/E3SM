&ctl_nl
ne = 256
mesh_file='none'
/

&vert_nl
! for E3SM output, set this correctly:
!vanalytic=0
vfile_mid     = "/global/u2/t/taylorm/codes/acme/components/homme/test/vcoord/scream-128m.ascii"
vfile_int     = "/global/u2/t/taylorm/codes/acme/components/homme/test/vcoord/scream-128i.ascii"
/

&analysis_nl
tool = 'interpolate_tool'
interpolate_analysis=.true.
infilenames = "/global/cscratch1/sd/taylorm/e3sm_scratch/cori-knl/ape256b/run/ape256b.cam.h2.0001-02-02-00000.nc"
output_dir = ""
output_timeunits=1
output_frequency=1
output_varnames1='all'
output_type='pnetcdf64' 
interp_gridtype=1   ! 1=equi-angle cap grid.  2=gauss.  3=equ-angle offset(no poles)
interp_type=1       ! 0=native SE, 1=bilinear
!interp_nlat=3073   ! overide defaults if you want
!interp_nlon=6144
io_stride = 68
/
