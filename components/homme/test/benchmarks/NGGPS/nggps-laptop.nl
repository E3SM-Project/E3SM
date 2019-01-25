&ctl_nl
ndays = 1
NThreads=-1   ! use maximum threads available
partmethod    = 4
topology      = "cube"
test_case     = "asp_baroclinic"
u_perturb = 1
rotate_grid = 0
ne=4
qsize = 3
! tstep=40  360 steps = 4h  (Bechmark reports 2h time)
nmax = 360
statefreq=360
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep=40
rsplit=3
qsplit = 1
tstep_type = 5
integration   = "explicit"
nu=4.5e17
nu_div=4.5e17
nu_p=4.5e17
nu_q=4.5e17
nu_s=4.5e17
nu_top = 0e5
se_ftype     = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_order = 2
hypervis_subcycle=1
hypervis_subcycle_q=1
theta_hydrostatic_mode = .true.
/
&vert_nl
vform         = "ccm"
vfile_mid = './camm-26.ascii'
vfile_int = './cami-26.ascii'
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

&analysis_nl
! to compare with EUL ref solution:
! interp_nlat = 32
! interp_nlon = 64
! interp_gridtype=2
 
 output_timeunits=1,1  !0=timesteps, 1=days, 2=hours, 3=seconds
 output_frequency=1,0
 output_start_time=0,0
 output_end_time=30000,30000
 output_varnames1='ps','zeta','dp3d','T','u','v','omega','Q','Q2','Q3'
 output_varnames2=
 io_stride=8
 output_type = 'netcdf' 
/

