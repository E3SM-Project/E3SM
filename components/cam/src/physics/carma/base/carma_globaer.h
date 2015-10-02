! CARMA Type aliases
! ---------------------
! This file containts shortcut names that map the variable names that
! were traditionally used in the common blocks by the Fortran 77 version
! of CARMA (globeaer.h), to the corresponding structure members in the
! Fortran 90 version of CARMA. This allows the older code to be
! converted to F90 with minimal changes, but without adding any
! processing overhead.
! ---------------------------------------------

! NOTE: Using macros causes some limitations:
!
!   1) You can not have another #define as a parameter to a macro. This causes
!      multiple expansions of the parameter. To prevent this, assign the parameter
!      to a varaible and use the variable in the macro.
!
!   2) You can not have comments on the same line as a macro. Put comments on the
!       line before the one with the macro.
!
!   3) Not all fortran preprocessors support the CPP's handling of recursion for
!      macro names. To work out of the box with the broadest number of fortran
!      compilers this requires making the field name different from the macro
!      or it will recursively try to replace the macro again (or report an
!      error message about recursion. Intel and IBM compilers handle it properly,
!      but Portland Group does not. To work around this problem, fields in the
!      cstate and carma structure are preceeded by f_ to make their names unique.

#define NZ            cstate%f_NZ
#define NZP1          cstate%f_NZP1
#define NGAS          carma%f_NGAS
#define NBIN          carma%f_NBIN
#define NGROUP        carma%f_NGROUP
#define NELEM         carma%f_NELEM
#define NSOLUTE       carma%f_NSOLUTE
#define NWAVE         carma%f_NWAVE

!  Model logical units for I/O
#define LUNOPRT       carma%f_LUNOPRT

!  Model startup control variables
#define do_print      carma%f_do_print

!  Gridding Information
#define igridv        cstate%f_igridv
#define igridh        cstate%f_igridh
#define xmet          cstate%f_xmet
#define ymet          cstate%f_ymet
#define zmet          cstate%f_zmet
#define zmetl         cstate%f_zmetl
#define xc            cstate%f_xc
#define yc            cstate%f_yc
#define zc            cstate%f_zc
#define dx            cstate%f_dx
#define dy            cstate%f_dy
#define dz            cstate%f_dz
#define zl            cstate%f_zl
#define lon           cstate%f_lon
#define lat           cstate%f_lat

! Element object
#define elemname(ielem)       carma%f_element(ielem)%f_name
#define rhoelem(ibin, ielem)  carma%f_element(ielem)%f_rho(ibin)
#define igelem(ielem)         carma%f_element(ielem)%f_igroup
#define itype(ielem)          carma%f_element(ielem)%f_itype
#define icomp(ielem)          carma%f_element(ielem)%f_icomposition
#define isolelem(ielem)       carma%f_element(ielem)%f_isolute

! Gas object
#define gasname(igas)         carma%f_gas(igas)%f_name
#define gwtmol(igas)          carma%f_gas(igas)%f_wtmol
#define ivaprtn(igas)         carma%f_gas(igas)%f_ivaprtn
#define igcomp(igas)          carma%f_gas(igas)%f_icomposition
#define dgc_threshold(igas)   carma%f_gas(igas)%f_dgc_threshold
#define ds_threshold(igas)    carma%f_gas(igas)%f_ds_threshold

! Group object
#define groupname(igroup)       carma%f_group(igroup)%f_name
#define nelemg(igroup)          carma%f_group(igroup)%f_nelem
#define ncore(igroup)           carma%f_group(igroup)%f_ncore
#define ishape(igroup)          carma%f_group(igroup)%f_ishape
#define ienconc(igroup)         carma%f_group(igroup)%f_ienconc
#define imomelem(igroup)        carma%f_group(igroup)%f_imomelem
#define solfac(igroup)          carma%f_group(igroup)%f_solface
#define scavcoef(igroup)        carma%f_group(igroup)%f_scavcoef
#define if_sec_mom(igroup)      carma%f_group(igroup)%f_if_sec_mom
#define is_grp_fractal(igroup)  carma%f_group(igroup)%f_is_fractal
#define is_grp_ice(igroup)      carma%f_group(igroup)%f_is_ice
#define is_grp_cloud(igroup)    carma%f_group(igroup)%f_is_cloud
#define is_grp_sulfate(igroup)  carma%f_group(igroup)%f_is_sulfate
#define grp_do_vtran(igroup)    carma%f_group(igroup)%f_grp_do_vtran
#define grp_do_drydep(igroup)   carma%f_group(igroup)%f_grp_do_drydep
#define irhswell(igroup)        carma%f_group(igroup)%f_irhswell
#define irhswcomp(igroup)       carma%f_group(igroup)%f_irhswcomp
#define rmrat(igroup)           carma%f_group(igroup)%f_rmrat
#define eshape(igroup)          carma%f_group(igroup)%f_eshape
#define r(ibin,igroup)          carma%f_group(igroup)%f_r(ibin)
#define rmass(ibin,igroup)      carma%f_group(igroup)%f_rmass(ibin) 
#define vol(ibin,igroup)        carma%f_group(igroup)%f_vol(ibin)
#define dr(ibin,igroup)         carma%f_group(igroup)%f_dr(ibin)
#define dm(ibin,igroup)         carma%f_group(igroup)%f_dm(ibin)
#define rmassup(ibin,igroup)    carma%f_group(igroup)%f_rmassup(ibin)
#define rmin(igroup)            carma%f_group(igroup)%f_rmin
#define rmassmin(igroup)        carma%f_group(igroup)%f_rmassmin
#define rup(ibin,igroup)        carma%f_group(igroup)%f_rup(ibin)
#define rlow(ibin,igroup)       carma%f_group(igroup)%f_rlow(ibin)
#define icorelem(icore,igroup)  carma%f_group(igroup)%f_icorelem(icore)
#define ifallrtn(igroup)        carma%f_group(igroup)%f_ifallrtn
#define arat(ibin,igroup)       carma%f_group(igroup)%f_arat(ibin) 
#define rrat(ibin,igroup)       carma%f_group(igroup)%f_rrat(ibin) 
#define rprat(ibin,igroup)      carma%f_group(igroup)%f_rprat(ibin) 
#define qext(iwave,ibin,igroup) carma%f_group(igroup)%f_qext(iwave,ibin)
#define ssa(iwave,ibin,igroup)  carma%f_group(igroup)%f_ssa(iwave,ibin)
#define do_mie(igroup)          carma%f_group(igroup)%f_do_mie
#define imiertn(igroup)         carma%f_group(igroup)%f_imiertn
#define dpc_threshold(igroup)   carma%f_group(igroup)%f_dpc_threshold
#define rmon(igroup)            carma%f_group(igroup)%f_rmon
#define df(ibin,igroup)         carma%f_group(igroup)%f_df(ibin)
#define nmon(ibin,igroup)       carma%f_group(igroup)%f_nmon(ibin)
#define falpha(igroup)          carma%f_group(igroup)%f_falpha
#define neutral_volfrc(igroup)  carma%f_group(igroup)%f_neutral_volfrc

! Solute object
#define solname(isolute)      carma%f_solute(isolute)%f_name
#define sol_ions(isolute)     carma%f_solute(isolute)%f_ions
#define solwtmol(isolute)     carma%f_solute(isolute)%f_wtmol
#define rhosol(isolute)       carma%f_solute(isolute)%f_rho

! Optical properties
#define wave          carma%f_wave
#define dwave         carma%f_dwave
#define do_wave_emit  carma%f_do_wave_emit

!  Model option & control variables
#define do_clearsky   carma%f_do_clearsky
#define do_cnst_rlh   carma%f_do_cnst_rlh
#define do_coag       carma%f_do_coag
#define do_detrain    carma%f_do_detrain
#define do_fixedinit  carma%f_do_fixedinit
#define do_grow       carma%f_do_grow
#define do_explised   carma%f_do_explised
#define do_incloud    carma%f_do_incloud
#define do_pheat      carma%f_do_pheat
#define do_pheatatm   carma%f_do_pheatatm
#define do_print_init carma%f_do_print_init
#define do_step       carma%f_do_step
#define do_substep    carma%f_do_substep
#define do_thermo     carma%f_do_thermo
#define do_vdiff      carma%f_do_vdiff
#define do_vtran      carma%f_do_vtran
#define do_drydep     carma%f_do_drydep
#define if_nuc        carma%f_if_nuc
#define time          cstate%f_time
#define dtime         cstate%f_dtime
#define dtime_orig    cstate%f_dtime_orig
#define nretries      cstate%f_nretries
#define dtmin         carma%f_dtmin
#define dtmax         carma%f_dtmax
#define conmax        carma%f_conmax
#define maxsubsteps   carma%f_maxsubsteps
#define minsubsteps   carma%f_minsubsteps
#define maxretries    carma%f_maxretries
#define ifall         carma%f_ifall
#define icoagop       carma%f_icoagop
#define icollec       carma%f_icollec
#define itbnd_pc      carma%f_itbnd_pc
#define ibbnd_pc      carma%f_ibbnd_pc
#define inucgas       carma%f_inucgas
#define igrowgas      carma%f_igrowgas
#define nnuc2elem     carma%f_nnuc2elem
#define ievp2elem     carma%f_ievp2elem
#define nnucelem      carma%f_nnucelem
#define inucproc      carma%f_inucproc
#define inuc2elem     carma%f_inuc2elem
#define inucelem      carma%f_inucelem
#define inuc2bin      carma%f_inuc2bin
#define ievp2bin      carma%f_ievp2bin
#define nnucbin       carma%f_nnucbin
#define inucbin       carma%f_inucbin
#define dt_threshold  carma%f_dt_threshold
#define igash2o       carma%f_igash2o
#define igash2so4     carma%f_igash2so4
#define igasso2       carma%f_igasso2
#define tstick        carma%f_tstick
#define gsticki       carma%f_gsticki
#define gstickl       carma%f_gstickl
#define cstick        carma%f_cstick

#define max_nsubstep  cstate%f_max_nsubstep
#define max_nretry    cstate%f_max_nretry
#define nstep         cstate%f_nstep
#define nsubstep      cstate%f_nsubstep
#define nretry        cstate%f_nretry
#define zsubsteps     cstate%f_zsubsteps
  
!  Particle grid structure
#define diffmass      carma%f_diffmass
#define rhop          cstate%f_rhop
#define r_wet         cstate%f_r_wet
#define rlow_wet      cstate%f_rlow_wet
#define rup_wet       cstate%f_rup_wet
#define rhop_wet      cstate%f_rhop_wet
#define r_ref         cstate%f_r_ref
#define rhop_ref      cstate%f_rhop_ref

!  Atmospheric structure
#define rhoa          cstate%f_rhoa
#define rhoa_wet      cstate%f_rhoa_wet
#define t             cstate%f_t
#define p             cstate%f_p
#define pl            cstate%f_pl
#define relhum        cstate%f_relhum
#define wtpct         cstate%f_wtpct
#define told          cstate%f_told
#define rmu           cstate%f_rmu
#define thcond        cstate%f_thcond
#define thcondnc      cstate%f_thcondnc
#define dkz           cstate%f_dkz

! Model primary vars
#define pc            cstate%f_pc
#define pcd           cstate%f_pcd
#define pc_surf       cstate%f_pc_surf
#define gc            cstate%f_gc
#define sedimentationflux       cstate%f_sedimentationflux
#define cldfrc        cstate%f_cldfrc
#define rhcrit        cstate%f_rhcrit

!  Model secondary variables
#define pcl           cstate%f_pcl
#define gcl           cstate%f_gcl
#define d_gc          cstate%f_d_gc
#define d_t           cstate%f_d_t
#define dpc_sed       cstate%f_dpc_sed
#define pconmax       cstate%f_pconmax
#define coaglg        cstate%f_coaglg
#define coagpe        cstate%f_coagpe
#define rnuclg        cstate%f_rnuclg
#define rnucpe        cstate%f_rnucpe
#define rhompe        cstate%f_rhompe
#define pc_nucl       cstate%f_pc_nucl
#define growpe        cstate%f_growpe
#define evappe        cstate%f_evappe
#define coreavg       cstate%f_coreavg
#define coresig       cstate%f_coresig
#define evdrop        cstate%f_evdrop
#define evcore        cstate%f_evcore
#define growlg        cstate%f_growlg
#define evaplg        cstate%f_evaplg
#define gasprod       cstate%f_gasprod
#define rlheat        cstate%f_rlheat
#define cmf           cstate%f_cmf
#define totevap       cstate%f_totevap
#define pc_topbnd     cstate%f_pc_topbnd
#define pc_botbnd     cstate%f_pc_botbnd
#define ftoppart      cstate%f_ftoppart
#define fbotpart      cstate%f_fbotpart
#define cmf           cstate%f_cmf
#define totevap       cstate%f_totevap
#define too_small     cstate%f_too_small
#define too_big       cstate%f_too_big
#define nuc_small     cstate%f_nuc_small
#define rlprod        cstate%f_rlprod
#define phprod        cstate%f_phprod

!  Coagulation kernels and bin pair mapping
#define ck0           carma%f_ck0
#define grav_e_coll0  carma%f_grav_e_coll0
#define icoag         carma%f_icoag
#define icoagelem     carma%f_icoagelem
#define icoagelem_cm  carma%f_icoagelem_cm
#define kbin          carma%f_kbin

#define ckernel       cstate%f_ckernel
#define pkernel       carma%f_pkernel

#define volx          carma%f_volx
#define ilow          carma%f_ilow
#define jlow          carma%f_jlow
#define iup           carma%f_iup
#define jup           carma%f_jup
#define npairl        carma%f_npairl
#define npairu        carma%f_npairu

!   Coagulation group pair mapping
#define iglow         carma%f_iglow
#define jglow         carma%f_jglow
#define igup          carma%f_igup
#define jgup          carma%f_jgup

!  Particle fall velocities, transport rates, and coagulation kernels
#define bpm           cstate%f_bpm
#define vf            cstate%f_vf
#define re            cstate%f_re
#define vf_const      carma%f_vf_const
#define vd            cstate%f_vd

!  Condensational growth parameters
#define diffus        cstate%f_diffus
#define rlhe          cstate%f_rlhe
#define rlhm          cstate%f_rlhm
#define pvapl         cstate%f_pvapl
#define pvapi         cstate%f_pvapi
#define surfctwa      cstate%f_surfctwa
#define surfctiw      cstate%f_surfctiw
#define surfctia      cstate%f_surfctia
#define akelvin       cstate%f_akelvin
#define akelvini      cstate%f_akelvini
#define ft            cstate%f_ft
#define gro           cstate%f_gro
#define gro1          cstate%f_gro1
#define gro2          cstate%f_gro2
#define supsatl       cstate%f_supsatl
#define supsati       cstate%f_supsati
#define supsatlold    cstate%f_supsatlold
#define supsatiold    cstate%f_supsatiold
#define scrit         cstate%f_scrit
#define rlh_nuc       carma%f_rlh_nuc
#define radint        cstate%f_radint
#define partheat      cstate%f_partheat
#define dtpart        cstate%f_dtpart
#define pratt         carma%f_pratt
#define prat          carma%f_prat
#define pden1         carma%f_pden1
#define palr          carma%f_palr
