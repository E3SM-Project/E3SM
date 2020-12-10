&ctl_nl
partmethod        = 4
topology          = "cube"
test_case         = "jw_baroclinic"
u_perturb         = 1
rotate_grid       = 0
ne                = 2
qsize             = 4
ndays             = 1
statefreq         = 9999
mesh_file         = '/dev/null'
tstep             = 150
rsplit            = 2
qsplit            = 1 
integration       = "explicit"
smooth            = 0
nu                = 7e15
nu_div            = 1e15
nu_p              = 7e15
nu_q              = 7e15
nu_s              =-1
nu_top            = 2.5e5
se_ftype          = 0
limiter_option    = 9
vert_remap_q_alg  = 1
hypervis_scaling  = 0
hypervis_order    = 2
hypervis_subcycle = 2
hypervis_subcycle_tom = 0
theta_hydrostatic_mode = false
theta_advect_form = 0
tstep_type        = 10
/
&vert_nl
vform     = "ccm"
vfile_mid = './data/scream-128m.ascii'
vfile_int = './data/scream-128i.ascii'
/
