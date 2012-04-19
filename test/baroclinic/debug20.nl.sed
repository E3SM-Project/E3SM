&ctl_nl
NThreads      = 8
partmethod    = 4
topology      = "cube"
test_case     = "baroclinic"
ne            = NE
qsize          = 0
nmax          = NMAX
statefreq     = SFREQ
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = TSTEP
integration   = "semi_imp"
smooth        = 0.000
nu            = NU1
nu_s          = NU2
hypervis_order = ORDER
hypervis_subcycle = SUBCYCLE
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
filter_mu     = 0.04
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/sabm-20.fbin.littleendian"
vfile_int     = "vcoord/sabi-20.fbin.littleendian"
/
&analysis_nl
 output_timeunits=1,
 output_frequency=0,
 output_varnames1='T','u','v','zeta'
/
