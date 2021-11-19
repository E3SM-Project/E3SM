&ctl_nl
disable_diagnostics    = .false.
partmethod             = 4
topology               = "cube"
mesh_file              = '/dev/null'
test_case              = "jw_baroclinic"
rotate_grid            = 0
ne                     = ${HOMME_TEST_NE}
tstep                  = ${HOMME_TEST_TIME_STEP}
u_perturb              = 1
statefreq              = 9999
tstep_type             = ${HOMME_TTYPE}
dt_remap_factor        = ${HOMME_TEST_REMAP_FACTOR}
dt_tracer_factor       = ${HOMME_TEST_TRACERS_FACTOR}
nu                     = ${HOMME_TEST_NU}
nu_div                 = ${HOMME_TEST_NU}
nu_p                   = ${HOMME_TEST_NU}
nu_q                   = ${HOMME_TEST_NU}
nu_s                   =-1
nu_top                 = ${HOMME_TEST_NUTOP}
se_ftype               = ${HOMME_SE_FTYPE}
limiter_option         = ${HOMME_TEST_LIM}
vert_remap_q_alg       = 1
hypervis_scaling       = ${HOMME_TEST_HVSCALING}
hypervis_order         = 2
hypervis_subcycle      = ${HOMME_TEST_HVS}
hypervis_subcycle_tom  = ${HOMME_TEST_HVS_TOM}
theta_hydrostatic_mode = ${HOMME_THETA_HY_MODE}
theta_advect_form      = ${HOMME_THETA_FORM}
moisture               = '${HOMME_TEST_MOISTURE}'
/
&vert_nl
vfile_mid = './vcoord/${HOMME_TEST_VCOORD_MID_FILE}'
vfile_int = './vcoord/${HOMME_TEST_VCOORD_INT_FILE}'
/
