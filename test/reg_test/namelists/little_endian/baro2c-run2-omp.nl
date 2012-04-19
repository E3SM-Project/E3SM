&ctl_nl
NThreads      = 4
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
ne            = 6
qsize = 1
tracer_advection_formulation=1
nmax          = 3
statefreq     = 1
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = -100
restartfile   = "./restart/R000002160"
runtype = 2
tstep         =  30
integration   = "explicit"
smooth        = 0.0000
energy_fixer  = 0
nu            = 4e16
nu_s          = 4e16
nu_p          = 0
nu_top = 0
hypervis_order = 2
hypervis_subcycle = 1
u_perturb      = 1
initial_total_mass = 0
limiter_option = 4
/
&solver_nl
precon_method = "identity"
maxits        = 50
tol           = 1.e-7
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.25D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/camm-26.fbin.littleendian"
vfile_int     = "vcoord/cami-26.fbin.littleendian"
/
&analysis_nl
 output_timeunits=0
 output_frequency=0
/

