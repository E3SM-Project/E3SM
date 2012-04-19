&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "baroclinic"
ne            = 9
ndays         = 12
statefreq     = 10
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = 600.0
integration   = "semi_imp"
smooth        = 0.05
nu            = 7.0e5
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
filter_mu     = 0.20D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "../vcoord/sabm-20.fbin"
vfile_int     = "../vcoord/sabi-20.fbin"
/
&analysis_nl
output_timeunits=1,
output_frequency=1,
/
