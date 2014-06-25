&ctl_nl
NThreads                     = 1
partmethod                   = 4
topology                     = "cube"
test_case                    = "asp_baroclinic"
rotate_grid                  = 0
ne                           = 15 !120 !30 !30
qsize                        = 4
ntrac                        = 4
test_cfldep                  = .TRUE.
tstep_type                   = 1 
ndays                        = 9
statefreq                    = 240
runtype                      = 0
tstep                        = 180 ! 90 !20 !90 !50
qsplit                       = 4 !4
rk_stage_user                = 3
integration                  = "explicit"
smooth                       = 0
nu                           = 1e16 !1.1e13 !9.6e14 !2e14
nu_s                         = -1       ! use same value as nu
nu_q                         = 1e16 !1.1e13 !9.6e14 !2e14    
nu_p                         = 0
limiter_option               = 8 
energy_fixer                 = -1
hypervis_order               = 2
hypervis_subcycle            = 1
u_perturb                    = 1
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
vfile_mid     = "../vcoord/camm-26.fbin.littleendian"
vfile_int     = "../vcoord/cami-26.fbin.littleendian"
/

&prof_inparm
profile_outpe_num   = 100
profile_single_file	= .true.
profile_barrier = .true.
/

&analysis_nl
 output_timeunits  = 1,1
 output_frequency  = 9,9
 output_start_time = 0,0
 output_end_time   = 20,20
 output_varnames1  = 'ps', 'zeta' !, 'DIFFT'
!  output_varnames2  = 'Q', 'Q2', 'Q3', 'Q4', 'phys_lat','phys_lon'
 output_varnames2  = 'Q', 'Q2', 'Q3', 'Q4' , 'zeta', 'C', 'C2', 'C3', 'C4', 'phys_lat','phys_lon'
! output_varnames2  = 'C', 'C2', 'C3', 'C4', 'zeta', 'phys_lat','phys_lon'
!  output_varnames2  = 'Q3', 'Q4' ,'C', 'C3', 'zeta', 'phys_lat','phys_lon'
 io_stride         = 8
 interp_nlat       = 256
 interp_nlon       = 512
 output_dir ="./movies/"
 output_type       = 'netcdf' 
/

