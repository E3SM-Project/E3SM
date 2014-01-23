&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb      = 1
rotate_grid = 0
ne=16
qsize         = 4
ndays          = 9
statefreq     = 144
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep=600
rsplit=3
qsplit = 1
tstep_type = 5
energy_fixer  = -1
integration   = "explicit"
smooth        = 0
nu=7e15
nu_div=7e15
nu_p=7e15
nu_q=7e15
nu_s=-1
nu_top = 2.5e5
se_ftype     = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_order = 2
hypervis_subcycle=3
hypervis_subcycle_q=1
/
&solver_nl
precon_method = "identity"
maxits        = 500
tol           = 1.e-9
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "./camm-30.ascii"
vfile_int     = "./cami-30.ascii"
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

&analysis_nl
! to compare with EUL ref solution:
! interp_nlat = 512
! interp_nlon = 1024
 interp_gridtype=2
 
 output_timeunits=1,1
 output_frequency=3,3
 output_start_time=0,0
 output_end_time=30000,30000
 output_varnames1='ps','zeta','dp3d'
 output_varnames2='Q','Q2','Q3','Q4','Q5'
 io_stride=8
 output_type = 'netcdf' 
/

