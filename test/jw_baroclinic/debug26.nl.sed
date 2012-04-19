&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
ne            = NE
qsize         = 0
tracer_advection_formulation = 0
nmax          = NMAX
statefreq     = SFREQ
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = TSTEP
integration   = "explicit"
smooth        = 0.0000
energy_fixer  = 0
nu            = NU1
nu_s          = NU2
nu_p          = NUP
nu_top        = 0
hypervis_order = ORDER
hypervis_subcycle = SUBCYCLE
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

