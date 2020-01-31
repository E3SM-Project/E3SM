#define CLDFRAC

module microphysics

! main interface to Morrison microphysics.
! original implementation by Peter Blossey, UW

use params, only: lcond, lsub, fac_cond, fac_sub, ggr

use grid, only: nx,ny,nzm,nz, &  !grid dimensions; nzm = nz-1 # of scalar lvls
     dimx1_s,dimx2_s,dimy1_s,dimy2_s, & ! actual scalar-array dimensions in x,y
     dz, adz, dostatis, masterproc, &
     doSAMconditionals, dosatupdnconditionals

use vars, only: pres, rho, dt, dtn, w, t, tlatqi, condavg_mask, &
     ncondavg, condavgname, condavglongname
use vars, only: tke2, tk2
use params, only: doprecip, docloud, doclubb

use module_mp_GRAUPEL, only: GRAUPEL_INIT, M2005MICRO_GRAUPEL, &
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! use graupel
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
#if (defined CRM && defined MODAL_AERO)
      domodal_aero,     &   ! use modal aerosol from the CAM
#endif
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for aerosol_mode=1 (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for aerosol_mode=2
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed ! option to specify pgam (exponent of cloud water's gamma distn)

#ifdef CRM
  use abortutils,     only: endrun
#endif

implicit none

logical :: isallocatedMICRO = .false.

integer :: nmicro_fields ! total number of prognostic water vars

real, allocatable, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing
integer :: index_water_vapor ! separate water vapor index used by SAM

real, allocatable, dimension(:) :: lfac
integer, allocatable, dimension(:) :: flag_wmass, flag_precip, flag_number
integer, allocatable, dimension(:) :: flag_micro3Dout

integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

real, allocatable, dimension(:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
real, allocatable, dimension(:,:,:) :: reffc, reffi
real, allocatable, dimension(:,:,:) :: cloudliq

real, allocatable, dimension(:,:) :: & ! statistical arrays
     mkwle, & ! resolved vertical flux
     mkwsb, & ! SGS vertical flux
     mksed, & ! sedimentation vertical flux
     mkadv, & ! tendency due to vertical advection
     mkdiff, &! tendency due to vertical diffusion
     mklsadv, & ! tendency due to large-scale vertical advection
     mfrac, & ! fraction of domain with microphysical quantity > 1.e-6
     stend, & ! tendency due to sedimentation
     mtend, & ! tendency due to microphysical processes (other than sedimentation)
     mstor, & ! storage terms of microphysical variables 
     trtau    ! optical depths of various species

real, allocatable, dimension(:) :: tmtend

real :: sfcpcp, sfcicepcp

! arrays with names/units for microphysical outputs in statistics.
character*3, allocatable, dimension(:) :: mkname
character*80, allocatable, dimension(:) :: mklongname
character*10, allocatable, dimension(:) :: mkunits
real, allocatable, dimension(:) :: mkoutputscale
logical douse_reffc, douse_reffi

! You can also have some additional, diagnostic, arrays, for example, total
! nonprecipitating cloud water, etc:

!bloss: array which holds temperature tendency due to microphysics
real, allocatable, dimension(:,:,:), SAVE :: tmtend3d

#ifdef CRM
real, allocatable, dimension(:) ::  qpevp   !sink of precipitating water due to evaporation (set to zero here)
real, allocatable, dimension(:) ::  qpsrc   !source of precipitation microphysical processes (set to mtend)
#endif 

real, allocatable, dimension(:,:,:)  :: wvar  ! the vertical velocity variance from subgrid-scale motion,
                                              ! which is needed in droplet activation.
#ifdef CRM
! hm 7/26/11 new output
real, public, allocatable, dimension(:,:,:)  :: aut1  !
real, public, allocatable, dimension(:,:,:)  :: acc1  !
real, public, allocatable, dimension(:,:,:)  :: evpc1  !
real, public, allocatable, dimension(:,:,:)  :: evpr1  !
real, public, allocatable, dimension(:,:,:)  :: mlt1  !
real, public, allocatable, dimension(:,:,:)  :: sub1  !
real, public, allocatable, dimension(:,:,:)  :: dep1  !
real, public, allocatable, dimension(:,:,:)  :: con1  !

real, public, allocatable, dimension(:,:,:)  :: aut1a  !
real, public, allocatable, dimension(:,:,:)  :: acc1a  !
real, public, allocatable, dimension(:,:,:)  :: evpc1a  !
real, public, allocatable, dimension(:,:,:)  :: evpr1a  !
real, public, allocatable, dimension(:,:,:)  :: mlt1a  !
real, public, allocatable, dimension(:,:,:)  :: sub1a  !
real, public, allocatable, dimension(:,:,:)  :: dep1a  !
real, public, allocatable, dimension(:,:,:)  :: con1a  !
#endif

#ifdef CLUBB_CRM
logical :: doclubb_tb = .false.    ! use clubb as a turbulence scheme only +++mhwang
                                  ! so liquid water is diagnosed based on saturaiton adjustment
#endif /*CLUBB_CRM*/

CONTAINS

!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
subroutine micro_setparm()
  use vars
#ifdef CLUBB_CRM
  use module_mp_graupel, only: NNUCCD_REDUCE_COEF, NNUCCC_REDUCE_COEF
#endif
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder
  
   NAMELIST /MICRO_M2005/ &
#ifdef CLUBB_CRM
      NNUCCD_REDUCE_COEF, NNUCCC_REDUCE_COEF, &
#endif
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! graupel species has qualities of hail
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization in place of KK(2000)
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for aerosol_mode=1 (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for aerosol_mode=2
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed, & ! option to specify pgam (exponent of cloud water's gamma distn)
      douse_reffc, &        ! use computed effective radius in radiation computation
      douse_reffi           ! use computed effective ice size in radiation computation

   !bloss: Create dummy namelist, so that we can figure out error code
   !       for a mising namelist.  This lets us differentiate between
   !       missing namelists and those with an error within the namelist.
   NAMELIST /BNCUIODSBJCB/ place_holder

   ! define default values for namelist variables
   doicemicro = .true.        ! use ice
   dograupel = .true.         ! use graupel
   dohail = .false.           ! graupel species has properties of graupel
   dosb_warm_rain = .false.   ! use KK (2000) warm rain scheme by default
   dopredictNc = .true.       ! prognostic cloud droplet number
#if (defined CRM && defined MODAL_AERO)
   domodal_aero = .true.      ! use modal aerosol
#endif
#ifdef CLUBB_CRM
     dosubgridw = .true.        ! Use clubb's w'^2 for sgs w
     aerosol_mode = 2           ! use lognormal CCN relationship
     doarcticicenucl = .false.  ! use mid-latitude parameters
     docloudedgeactivation  = .true. ! activate droplets at cloud base, and edges
#else
     aerosol_mode = 2      
     dosubgridw = .true.     
     doarcticicenucl = .false.  ! use mid-latitude parameters
     docloudedgeactivation  = .true. 
#endif /*CLUBB_CRM*/
   douse_reffc = .false.  ! use computed effective radius in rad computations?
   douse_reffi = .false.  ! use computed effective radius in rad computations?

   Nc0 = 100. ! default droplet number concentration
   
   ccnconst = 120.            ! maritime value (/cm3), adapted from Rasmussen 
   ccnexpnt = 0.4             !   et al (2002) by Hugh Morrison et al.  Values
                              !   of 1000. and 0.5 suggested for continental
!   aer_rm1 = 0.052           ! two aerosol mode defaults from MPACE (from Hugh)
!   aer_sig1 = 2.04
!   aer_n1 = 72.2
!   aer_rm2 = 1.3
!   aer_sig2 = 2.5
!   aer_n2 = 1.8

   aer_rm1 = 0.052           ! two aerosol mode defaults (from mhwang for testing in global models)
   aer_sig1 = 2.04
   aer_n1 = 2500 
   aer_rm2 = 1.3
   aer_sig2 = 2.5
   aer_n2 = 1.8
   
   dofix_pgam = .false.
   pgam_fixed = 5. ! middle range value -- corresponds to radius dispersion ~ 0.4

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
!  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
  
  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
!  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
!  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_M2005 namelist
!  read (55,MICRO_M2005,IOSTAT=ios)

!  if (ios.ne.0) then
!     !namelist error checking
!     if(ios.ne.ios_missing_namelist) then
!        write(*,*) '****** ERROR: bad specification in MICRO_M2005 namelist'
!        call task_abort()
!     elseif(masterproc) then
!        write(*,*) '****************************************************'
!        write(*,*) '****** No MICRO_M2005 namelist in prm file *********'
!        write(*,*) '****************************************************'
!     end if
!  end if
!  close(55)

   if(.not.doicemicro) dograupel=.false.

   if(dohail.and..NOT.dograupel) then
      if(masterproc) write(*,*) 'dograupel must be .true. for dohail to be used.'
      call task_abort()
   end if

   ! write namelist values out to file for documentation
!   if(masterproc) then
!      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.options_namelist', form='formatted', position='append')    
!      write (unit=55,nml=MICRO_M2005,IOSTAT=ios)
!      write(55,*) ' '
!      close(unit=55)
!   end if

   ! scale values of parameters for m2005micro
   aer_rm1 = 1.e-6*aer_rm1 ! convert from um to m
   aer_rm2 = 1.e-6*aer_rm2 
   aer_n1 = 1.e6*aer_n1 ! convert from #/cm3 to #/m3
   aer_n2 = 1.e6*aer_n2
   
  nmicro_fields = 1 ! start with water vapor and cloud water mass mixing ratio
#ifdef CLUBB_CRM
  if(docloud.or.doclubb) then
#else
  if(docloud) then
#endif
!bloss/qt     nmicro_fields = nmicro_fields + 1 ! add cloud water mixing ratio
     if(dopredictNc) nmicro_fields = nmicro_fields + 1 ! add cloud water number concentration (if desired)
  end if
  if(doprecip)    nmicro_fields = nmicro_fields + 2 ! add rain mass and number (if desired)
  if(doicemicro)  nmicro_fields = nmicro_fields + 4 ! add snow and cloud ice number and mass (if desired)
  if(dograupel)   nmicro_fields = nmicro_fields + 2 ! add graupel mass and number (if desired)

  ! specify index of various quantities in micro_field array
  !  *** note that not all of these may be used if(.not.doicemicro) ***
  iqv = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
!bloss/qt  iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
  
!bloss/qt: cloud liquid water no longer prognosed
  if(dopredictNc) then
     incl = 2  ! cloud water number mixing ratio [#/kg dry air]
     iqr = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
     inr = 4   ! rain number mixing ratio [#/kg dry air]
     iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = 8   ! snow number mixing ratio [#/kg dry air]
     iqg = 9  ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = 10  ! graupel number mixing ratio [#/kg dry air]
  else
     iqr = 2   ! rain mass mixing ratio [kg H2O / kg dry air]
     inr = 3   ! rain number mixing ratio [#/kg dry air]
     iqci = 4  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = 5  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = 6   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = 7   ! snow number mixing ratio [#/kg dry air]
     iqg = 8   ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = 9  ! graupel number mixing ratio [#/kg dry air]
  end if

  ! stop if icemicro is specified without precip -- we don't support this right now.
  if((doicemicro).and.(.not.doprecip)) then
     if(masterproc) write(*,*) 'Morrison 2005 Microphysics does not support both doice and .not.doprecip'
     call task_abort()
  end if
  index_water_vapor = iqv ! set SAM water vapor flag

  if(.not.isallocatedMICRO) then
     ! allocate microphysical variables
     allocate(micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,nmicro_fields), &
          fluxbmk(nx,ny,nmicro_fields), fluxtmk(nx,ny,nmicro_fields), &
          reffc(nx,ny,nzm), reffi(nx,ny,nzm), &
          mkwle(nz,nmicro_fields), mkwsb(nz,nmicro_fields), &
          mkadv(nz,nmicro_fields), mkdiff(nz,nmicro_fields), &
          mklsadv(nz,nmicro_fields), &
          stend(nzm,nmicro_fields), mtend(nzm,nmicro_fields), &
          mfrac(nzm,nmicro_fields), trtau(nzm,nmicro_fields), &
          mksed(nzm,nmicro_fields), tmtend(nzm), &
          mstor(nzm,nmicro_fields),  &
          cloudliq(nx,ny,nzm), &
          tmtend3d(nx,ny,nzm), flag_micro3Dout(nmicro_fields), &
          flag_wmass(nmicro_fields), flag_precip(nmicro_fields), &
          flag_number(nmicro_fields), lfac(nmicro_fields), &
          mkname(nmicro_fields), mklongname(nmicro_fields), &
          mkunits(nmicro_fields), mkoutputscale(nmicro_fields), STAT=ierr)

#ifdef CRM
     allocate (qpevp(nz), qpsrc(nz), STAT=ierr)
#endif
     allocate (wvar(nx,ny,nzm), STAT=ierr)

#ifdef CRM
! hm 7/26/11, add new output
     allocate (aut1(nx,ny,nzm), STAT=ierr)
     allocate (acc1(nx,ny,nzm), STAT=ierr)
     allocate (evpc1(nx,ny,nzm), STAT=ierr)
     allocate (evpr1(nx,ny,nzm), STAT=ierr)
     allocate (mlt1(nx,ny,nzm), STAT=ierr)
     allocate (sub1(nx,ny,nzm), STAT=ierr)
     allocate (dep1(nx,ny,nzm), STAT=ierr)
     allocate (con1(nx,ny,nzm), STAT=ierr)

     allocate (aut1a(nx,ny,nzm), STAT=ierr)
     allocate (acc1a(nx,ny,nzm), STAT=ierr)
     allocate (evpc1a(nx,ny,nzm), STAT=ierr)
     allocate (evpr1a(nx,ny,nzm), STAT=ierr)
     allocate (mlt1a(nx,ny,nzm), STAT=ierr)
     allocate (sub1a(nx,ny,nzm), STAT=ierr)
     allocate (dep1a(nx,ny,nzm), STAT=ierr)
     allocate (con1a(nx,ny,nzm), STAT=ierr)
#endif

     if(ierr.ne.0) then
        write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
        call task_abort()
     else
        isallocatedMICRO = .true.
     end if

  ! zero out statistics variables associated with cloud ice sedimentation
  !   in Marat's default SAM microphysics
  tlatqi = 0.

  ! initialize these arrays
  micro_field = 0.
  cloudliq = 0. !bloss/qt: auxially cloud liquid water variable, analogous to qn in MICRO_SAM1MOM
  fluxbmk = 0.
  fluxtmk = 0.
  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.
  mklsadv = 0.
  mstor =0.

  wvar = 0.

#ifdef CRM
! hm 7/26/11, new output
  aut1 = 0.
  acc1 = 0.
  evpc1 = 0.
  evpr1 = 0.
  mlt1 = 0.
  sub1 = 0.
  dep1 = 0.
  con1 = 0.
  aut1a = 0.
  acc1a = 0.
  evpc1a = 0.
  evpr1a = 0.
  mlt1a = 0.
  sub1a = 0.
  dep1a = 0.
  con1a = 0.
#endif

  ! initialize flag arrays to all mass, no number, no precip
  flag_wmass = 1
  flag_number = 0
  flag_precip = 0
  flag_micro3Dout = 0

  end if

  compute_reffc = douse_reffc
  compute_reffi = douse_reffi

end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the 
!   beginning of each run, initial or restart:
subroutine micro_init()

  use vars
#if (defined CRM && defined MODAL_AERO)
  use drop_activation, only: drop_activation_init
#endif
  
  implicit none
  
  real, dimension(nzm) :: qc0, qi0

! Commented out by dschanen UWM 23 Nov 2009 to avoid a linking error
! real, external :: satadj_water 
  integer :: k

  ! initialize flag arrays
  if(dopredictNc) then
     ! Cloud droplet number concentration is a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass  = (/1,0,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,1,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns
           flag_wmass  = (/1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1/)
           flag_number = (/0,1,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, Nc, qr, Nr
           flag_wmass  = (/1,0,1,0/)
           flag_precip = (/0,0,1,1/)
           flag_number = (/0,1,0,1/)
        else
          !bloss/qt: qt, Nc
           flag_wmass  = (/1,0/)
           flag_precip = (/0,0/)
           flag_number = (/0,1/)
        end if
     end if
  else
     ! Cloud droplet number concentration is NOT a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass  = (/1,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns
           flag_wmass  = (/1,1,0,1,0,1,0/)
           flag_precip = (/0,1,1,0,0,1,1/)
           flag_number = (/0,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, qr, Nr
           flag_wmass  = (/1,1,0/)
           flag_precip = (/0,1,1/)
           flag_number = (/0,0,1/)
        else
          !bloss/qt: only total water variable is needed for no-precip, 
          !            fixed droplet number, warm cloud and no cloud simulations.
           flag_wmass  = (/1/)
           flag_precip = (/0/)
           flag_number = (/0/)
        end if
     end if
  end if

  ! output all microphysical fields to 3D output files if using more than
  !   just docloud.  Otherwise, rely on basic SAM outputs
#ifdef CLUBB_CRM
  if((docloud.OR.doclubb).AND.(doprecip.OR.dopredictNc)) then
#else
  if(docloud.AND.(doprecip.OR.dopredictNc)) then
#endif
     flag_micro3Dout = 1
  end if

  ! initialize factor for latent heat
  lfac(:) = 1. ! use one as default for number species
  lfac(iqv) = lcond
!bloss/qt  if(docloud) lfac(iqcl) = lcond
  if(doprecip) lfac(iqr) = lcond
  if(doicemicro) then
     lfac(iqci) = lsub
     lfac(iqs) = lsub
     if(dograupel) lfac(iqg) = lsub
  end if

  call graupel_init() ! call initialization routine within mphys module
#if (defined CRM && defined MODAL_AERO)
  call drop_activation_init
#endif

  if(nrestart.eq.0) then

! In SPCAM,  do not need this part. 
#ifndef CRM
 ! compute initial profiles of liquid water - M.K.
      call satadj_liquid(nzm,tabs0,q0,qc0,pres*100.)

     ! initialize microphysical quantities
     q0 = q0 + qc0
     do k = 1,nzm
        micro_field(:,:,k,iqv) = q0(k)
        cloudliq(:,:,k) = qc0(k)
        tabs(:,:,k) = tabs0(k)
     end do
     if(dopredictNc) then ! initialize concentration somehow...
       do k = 1,nzm
         if(q0(k).gt.0.) then
            micro_field(:,:,k,incl) = 0.5*ccnconst*1.e6
         end if
       end do
     end if
#endif  ! CRM

#ifdef CLUBB_CRM
     if(docloud.or.doclubb)  call micro_diagnose()   ! leave this line here
#else
     if(docloud) call micro_diagnose()   ! leave this here
#endif


  end if

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

use vars, only: fluxbq, fluxtq
#ifdef CLUBB_CRM
use params, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
#endif

fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
#ifdef CLUBB_CRM
if ( doclubb .and. (doclubb_sfc_fluxes.or.docam_sfc_fluxes) ) then
  fluxbmk(:,:,index_water_vapor) = 0.0 ! surface qv (latent heat) flux
else
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
end if
#else
fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
#endif
fluxtmk(:,:,index_water_vapor) = fluxtq(:,:) ! top of domain qv flux

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., that is  all the microphysics processes except for the spatial transport happen.

! IMPORTANT: You need to use the thermodynamic constants like specific heat, or
! specific heat of condensation, gas constant, etc, the same as in file params.f90
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg) 
! It should not be changed during all of your point microphysical processes!

subroutine micro_proc()

use params, only: fac_cond, fac_sub, rgas, rv
use grid, only: z, zi

#ifdef CRM
use vars, only: t, gamaz, precsfc, precssfc, precflux, qpfall, tlat, prec_xy, &
#else
use vars, only: t, gamaz, precsfc, precflux, qpfall, tlat, prec_xy, &
#endif /*CRM*/
     nstep, nstatis, icycle, total_water_prec

#ifdef ECPP
use ecppvars, only: qlsink, qlsink_bf, prain, precr, precsolid, rh, qcloud_bf
#endif

#ifdef CLUBB_CRM
use params, only: doclubb, docloud, dosmoke
use grid, only: nz
use error_code, only: clubb_at_least_debug_level
use fill_holes, only: fill_holes_driver
use clubbvars, only: wp2, cloud_frac, rho_ds_zt, rho_ds_zm ! are used, but not modified here
use vars, only: qcl ! Used here and updated in micro_diagnose
use vars, only: prespot ! exner^-1
use module_mp_GRAUPEL, only: &
  cloud_frac_thresh ! Threshold for using sgs cloud fraction to weight 
                    ! microphysical quantities [%]
use clubb_precision, only: core_rknd
use constants_clubb, only: T_freeze_K
use vars, only: CF3D
#ifdef CLDFRAC
use icecld_fraction_module, only: aist_single
use module_mp_GRAUPEL, only: POLYSVP
#endif
#endif


real, dimension(nzm) :: &
     tmpqcl, tmpqci, tmpqr, tmpqs, tmpqg, tmpqv, &
     tmpncl, tmpnci, tmpnr, tmpns, tmpng,  &
     tmpw, tmpwsub, tmppres, tmpdz, tmptabs, &
! hm 7/26/11, new output
     tmpaut,tmpacc,tmpevpc,tmpevpr,tmpmlt, &
     tmpsub,tmpdep,tmpcon, &
     tmtend1d, &
     mtendqcl, mtendqci, mtendqr, mtendqs, mtendqg, mtendqv, &
     mtendncl, mtendnci, mtendnr, mtendns, mtendng,  &
     stendqcl, stendqci, stendqr, stendqs, stendqg, stendqv, &
     stendncl, stendnci, stendnr, stendns, stendng,  &
     effg1d, effr1d, effs1d, effc1d, effi1d

#ifdef ECPP
real, dimension(nzm) :: C2PREC,QSINK_TMP, CSED,ISED,SSED,GSED,RSED,RH3D   ! used for cloud chemistry and wet deposition in ECPP
#endif

#ifdef CLUBB_CRM
real(kind=core_rknd), dimension(nz) :: &
     qv_clip, qcl_clip
real, dimension(nzm) :: cloud_frac_in, ice_cldfrac
#ifdef CLDFRAC
real, dimension(nzm) :: liquid_cldfrac, cldmax
real :: landfrac, snowh
real :: evs, qvs, rhw
#endif 
#endif /*CLUBB_CRM*/

real, dimension(nzm,nmicro_fields) :: stend1d, mtend1d
real :: tmpc, tmpr, tmpi, tmps, tmpg
integer :: i1, i2, j1, j2, i, j, k, m, n

real(8) :: tmp_total, tmptot

! call t_startf ('micro_proc')

#ifndef CRM
if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
      do i=1,nx
         precsfc(i,j)=0.    ! in SPCAM, done in crm.F90
      end do
   end do
   do k=1,nzm
      precflux(k) = 0.   ! in SPCAM, done in crm.F90
   end do
end if
#endif ! end CRM

if(dostatis) then ! initialize arrays for statistics
   mfrac(:,:) = 0.
   mtend(:,:) = 0.
   trtau(:,:) = 0.
!   qpfall(:)=0.     ! in SPCAM, done in crm.F90
   tlat(:) = 0.
   tmtend3d(:,:,:) = 0.
end if
stend(:,:) = 0.
mksed(:,:) = 0.

!!$if(doprecip) total_water_prec = total_water_prec + total_water()
 
do j = 1,ny
   do i = 1,nx

      ! zero out mixing ratios of microphysical species
      tmpqv(:) = 0.
      tmpqcl(:) = 0.
      tmpncl(:) = 0.
      tmpqr(:) = 0.
      tmpnr(:) = 0.
      tmpqci(:) = 0.
      tmpnci(:) = 0.
      tmpqs(:) = 0.
      tmpns(:) = 0.
      tmpqg(:) = 0.
      tmpng(:) = 0.

      ! get microphysical quantities in this grid column
      tmpqv(:) = micro_field(i,j,:,iqv) !bloss/qt: This is total water (qv+qcl)
!bloss/qt: compute below from saturation adjustment.
!bloss/qt      tmpqcl(:) = micro_field(i,j,:,iqcl)
      if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)
      if(doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
         tmpnr(:) = micro_field(i,j,:,inr)
      end if

      if(doicemicro) then
         tmpqci(:) = micro_field(i,j,:,iqci)
         tmpnci(:) = micro_field(i,j,:,inci)
         tmpqs(:) = micro_field(i,j,:,iqs)
         tmpns(:) = micro_field(i,j,:,ins)
         if(dograupel) then
            tmpqg(:) = micro_field(i,j,:,iqg)
            tmpng(:) = micro_field(i,j,:,ing)
         end if
      end if

      ! get absolute temperature in this column
      !bloss/qt: before saturation adjustment for liquid,
      !          this is Tcl = T - (L/Cp)*qcl (the cloud liquid water temperature)
      tmptabs(:) = t(i,j,:)  &           ! liquid water-ice static energy over Cp
           - gamaz(:) &                                   ! potential energy
           + fac_cond * (tmpqr(:)) &    ! bloss/qt: liquid latent energy due to rain only
           + fac_sub  * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy

      tmpdz = adz(:)*dz
!      tmpw = 0.5*(w(i,j,1:nzm) + w(i,j,2:nz))  ! MK: changed for stretched grids 
      tmpw = ((zi(2:nz)-z(1:nzm))*w(i,j,1:nzm)+ &
             (z(1:nzm)-zi(1:nzm))*w(i,j,2:nz))/(zi(2:nz)-zi(1:nzm))
#ifdef CLUBB_CRM
      ! Added by dschanen on 4 Nov 2008 to account for w_sgs 
      if ( doclubb .and. dosubgridw ) then
        ! Compute w_sgs.  Formula is consistent with that used with 
        ! TKE from MYJ pbl scheme in WRF (see module_mp_graupel.f90).
        tmpwsub = sqrt( LIN_INT( real( wp2(i,j,2:nz) ), real( wp2(i,j,1:nzm) ), &
                                  zi(2:nz), zi(1:nzm), z(1:nzm) ) )
      else
!        tmpwsub = 0.
! diagnose tmpwsub from tke.
! Notes: tke has to be already prognsotic or diagnostic.
        tmpwsub = sqrt(tke2(i,j,:)/3.)  ! diagnosed tmpwsub from tke
! diagnose tmpwsub from tk
!        tmpwsub = sqrt(2*3.141593)*tk(i,j,:)/(dz*adz(:))  ! from Ghan et al. (1997, JGR).
      end if
#else /* Old code */
!      tmpwsub = 0.
! diagnose tmpwsub from tke.
! Notes: tke has to be already prognsotic or diagnostic.
      tmpwsub = sqrt(tke2(i,j,:)/3.)  ! diagnosed tmpwsub from tke
! diagnose tmpwsub from tk
!      tmpwsub = sqrt(2*3.141593)*tk(i,j,:)/(dz*adz(:))  ! from Ghan et al. (1997, JGR).
#endif
      wvar(i,j,:) = tmpwsub(:)

      tmppres(:) = 100.*pres(1:nzm)

      !bloss/qt: saturation adjustment to compute cloud liquid water content.
      !          Note: tmpqv holds qv+qcl on input, qv on output.
      !                tmptabs hold T-(L/Cp)*qcl on input, T on output.
      !                tmpqcl hold qcl on output.
      !                tmppres is unchanged on output, should be in Pa.
#ifdef CLUBB_CRM
      ! In the CLUBB case, we want to call the microphysics on sub-saturated grid
      ! boxes and weight by cloud fraction, therefore we use the CLUBB value of 
      ! liquid water. -dschanen 23 Nov 2009
      if ( .not. ( docloud .or. dosmoke ) ) then
        if(.not.doclubb_tb) then
         tmpqcl  = cloudliq(i,j,:) ! Liquid updated by CLUBB just prior to this
         tmpqv   = tmpqv - tmpqcl ! Vapor
         tmptabs = tmptabs + fac_cond * tmpqcl ! Update temperature
        else
          call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
          cloudliq(i,j,:) = tmpqcl
        end if
      else
        call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
      end if
#else
      call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
#endif

#ifdef CLUBB_CRM
      if ( doclubb ) then
        cloud_frac_in(1:nzm) = cloud_frac(i,j,2:nz)

! Ice cloud fraction: 0 at 0 C, and 100% at -35C. 
        do k=1, nzm
!           ice_cldfrac(k) = -(tmptabs(k)-T_freeze_K)/35.0 
!           ice_cldfrac(k) = min(1.0, max(ice_cldfrac(k), 0.0))
!           if((tmpqci(k)).lt.1.0e-8) then 
!            ice_cldfrac(k) = 0.0   
!           end if
        
! diagnost ice cloud fraction using CAM5 formula (here option 5 is used, so landfrac and snowh are not used)
           landfrac = 0.0
           snowh = 0.0
           ice_cldfrac(k) = 0.0
           call aist_single( tmpqv(k), tmptabs(k), tmppres(k), tmpqci(k), landfrac, snowh, ice_cldfrac(k))

           if(doclubb_tb) then
               if(tmpqcl(k).gt.1.0e-14) then
                  liquid_cldfrac(k) = 1.0
               else
                  liquid_cldfrac(k) = 0.0
               end if
               ice_cldfrac(k) = max(ice_cldfrac(k), liquid_cldfrac(k))  ! ice cloud fraction is not less than liquid cloud fraction
               cloud_frac_in(k) = max(liquid_cldfrac(k), ice_cldfrac(K))
           end if
        end do
#ifdef CLDFRAC
        if(.not. doclubb_tb) then
           liquid_cldfrac(1:nzm) = cloud_frac(i, j, 2:nz)
        end if
#endif
      else  ! not clubb
#ifdef CLDFRAC
        do k=1, nzm
         if(tmpqcl(k).gt.1.0e-14) then
           liquid_cldfrac(k) = 1.0
         else
           liquid_cldfrac(k) = 0.0
         end if 

! diagnost ice cloud fraction using CAM5 formula (here option 5 is used, so landfrac and snowh are not used)
         landfrac = 0.0
         snowh = 0.0
         ice_cldfrac(k) = 0.0
         call aist_single( tmpqv(k), tmptabs(k), tmppres(k), tmpqci(k), landfrac, snowh, ice_cldfrac(k))
         ice_cldfrac(k) = max(ice_cldfrac(k), liquid_cldfrac(k))  ! ice cloud fraction is not less than liquid cloud fraction
         cloud_frac_in(k) = max(liquid_cldfrac(k), ice_cldfrac(K))
        end do
#else
        cloud_frac_in(1:nzm) = 0.0
#endif
      end if
#endif  /*CLUBB_CRM*/


#ifdef ECPP
! save cloud water before microphysics process for the calculation 
! of qlsink in ECPP
      qcloud_bf(i,j,:) =  tmpqcl(:)
#endif /*ECPP*/

      i1 = 1 ! dummy variables used by WRF convention in subroutine call
      i2 = 1
      j1 = 1
      j2 = 1

! hm 7/26/11, initialize new output
      tmpaut=0.
      tmpacc=0.
      tmpevpc=0.
      tmpevpr=0.
      tmpmlt=0.
      tmpsub=0.
      tmpdep=0.
      tmpcon=0.

      mtendqv = 0.
      mtendqcl = 0.
      mtendqr = 0.
      mtendqci = 0.
      mtendqs = 0.
      mtendqg = 0.
      mtendncl = 0.
      mtendnr = 0.
      mtendnci = 0.
      mtendns = 0.
      mtendng = 0.

      tmtend1d = 0.

      sfcpcp = 0.
      sfcicepcp = 0.

      effc1d(:) = 10. ! default liquid and ice effective radii
      effi1d(:) = 75.

#ifdef CLUBB_CRM
      if ( doclubb ) then
        if ( any( tmpqv < 0. ) ) then
          qv_clip(2:nz) = tmpqv(1:nzm)
          qv_clip(1) = 0.0_core_rknd
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(0,*) "M2005 has received a negative water vapor"
          end if
          call fill_holes_driver( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
          tmpqv = qv_clip(2:nz)
        end if
        if ( any( tmpqcl < 0. ) ) then
          qcl_clip(2:nz) = tmpqcl(1:nzm)
          qcl_clip(1) = 0.0_core_rknd
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(0,*) "M2005 has received a negative liquid water"
          end if
          call fill_holes_driver( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
          tmpqcl = qcl_clip(2:nz)
        end if
      end if ! doclubb

#ifdef CLDFRAC
      if(doclubb .and. (.not.doclubb_tb)) then
        do k=1, nzm   
          if(tmpqcl(k).lt.1.0e-14) then 
            liquid_cldfrac(k)=0.0
          end if
!
! set a RH limit for liquid cloud fraction 
! if RHw is below 0.3, no liquid water. 
!          evs = polysvp(tmptabs(k), 0) 
!          qvs = rgas/rv*evs/(tmppres(k)-evs)
!          rhw = tmpqv(k)/qvs
!          if(rhw.lt.0.3) then
!            liquid_cldfrac(k)= 0.0 
!            tmpqv(k) = tmpqv(k) + tmpqcl(k)
!            tmptabs(k) = tmptabs(k) - fac_cond * tmpqcl(k) ! Update temperature
!            tmpqcl(k) = 0.0
!          end if

          ice_cldfrac(k) = max(liquid_cldfrac(k), ice_cldfrac(k))
          cloud_frac_in(k) = max(liquid_cldfrac(k), ice_cldfrac(k))
        end do
      end if
#endif

      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/o sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho,tmpdz,tmpw,tmpwsub, &
! hm 7/26/11, new output
           tmpacc,tmpaut,tmpevpc,tmpevpr,tmpmlt, &
           tmpsub,tmpdep,tmpcon, &
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
          stendqr,stendqci,stendqs,stendqcl,cloud_frac_in  & ! cloud_frac added by dschanen UWM
#ifdef CLDFRAC
           ,liquid_cldfrac, ice_cldfrac,  cldmax  &  ! cloud fraction for liquid and ice condensate 
#endif
#ifdef ECPP
           ,C2PREC,QSINK_TMP,CSED,ISED,SSED,GSED,RSED,RH3D &        ! mhwang add, for ECPP
#endif
                                        )
 
      if ( doclubb ) then
        if ( any( tmpqv < 0. ) ) then
          qv_clip(2:nz) = tmpqv(1:nzm)
          qv_clip(1) = 0.0_core_rknd
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(0,*) "M2005 has produced a negative water vapor"
          end if
          call fill_holes_driver( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
          tmpqv = qv_clip(2:nz)
        end if
        if ( any( tmpqcl < 0. ) ) then
          qcl_clip(2:nz) = tmpqcl(1:nzm)
          qcl_clip(1) = 0.0_core_rknd
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(0,*) "M2005 has produced a negative liquid water"
          end if
          call fill_holes_driver( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
          tmpqcl = qcl_clip(2:nz)
        end if
      end if ! doclubb
#else
      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/o sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho,tmpdz,tmpw,tmpwsub, &
! hm 7/26/11, new output
           tmpacc,tmpaut,tmpevpc,tmpevpr,tmpmlt, &
           tmpsub,tmpdep,tmpcon, &
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
           stendqr,stendqci,stendqs,stendqcl           &
#ifdef ECPP
           ,C2PREC,QSINK_TMP,CSED,ISED,SSED,GSED,RSED,RH3D &        ! mhwang add, for ECPP
#endif
                                        )
#endif

#ifdef CRM
! hm 7/26/11, new output
      aut1(i,j,:) = tmpaut(:)      
      acc1(i,j,:) = tmpacc(:)      
      evpc1(i,j,:) = tmpevpc(:)      
      evpr1(i,j,:) = tmpevpr(:)      
      mlt1(i,j,:) = tmpmlt(:)      
      sub1(i,j,:) = tmpsub(:)      
      dep1(i,j,:) = tmpdep(:)      
      con1(i,j,:) = tmpcon(:)      

! hm 8/31/11, new output for gcm-grid and time-step avg
! rates are summed here over the icycle loop
! note: rates are multiplied by time step, and then
! divided by dt in crm.F90 to get mean rates
      aut1a(i,j,:) = aut1a(i,j,:) + aut1(i,j,:)*dtn 
      acc1a(i,j,:) = acc1a(i,j,:) + acc1(i,j,:)*dtn
      evpc1a(i,j,:) = evpc1a(i,j,:) + evpc1(i,j,:)*dtn
      evpr1a(i,j,:) = evpr1a(i,j,:) + evpr1(i,j,:)*dtn
      mlt1a(i,j,:) = mlt1a(i,j,:) + mlt1(i,j,:)*dtn
      sub1a(i,j,:) = sub1a(i,j,:) + sub1(i,j,:)*dtn
      dep1a(i,j,:) = dep1a(i,j,:) + dep1(i,j,:)*dtn
      con1a(i,j,:) = con1a(i,j,:) + con1(i,j,:)*dtn
#endif

     ! update microphysical quantities in this grid column
      if(doprecip) then
         total_water_prec = total_water_prec + sfcpcp

         ! take care of surface precipitation
         precsfc(i,j) = precsfc(i,j) + sfcpcp/dz
         prec_xy(i,j) = prec_xy(i,j) + sfcpcp/dtn/dz
#ifdef CRM
         precssfc(i,j) = precssfc(i,j) + sfcicepcp/dz    ! the corect unit of precssfc should be mm/dz +++mhwang
#endif
         ! update rain
         micro_field(i,j,:,iqr) = tmpqr(:)
         micro_field(i,j,:,inr) = tmpnr(:)
      else
         ! add rain to cloud
         tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
         tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

         ! zero out rain 
         tmpqr(:) = 0.
         tmpnr(:) = 0.

         ! add rain tendencies to cloud
         stendqcl(:) = stendqcl(:) + stendqr(:)
         mtendqcl(:) = mtendqcl(:) + mtendqr(:)
         mtendncl(:) = mtendncl(:) + mtendnr(:)

         ! zero out rain tendencies
         stendqr(:) = 0.
         mtendqr(:) = 0.
         mtendnr(:) = 0.
      end if

      !bloss/qt: update total water and cloud liquid.
      !          Note: update of total water moved to after if(doprecip),
      !                  since no precip moves rain --> cloud liq.
      micro_field(i,j,:,iqv) = tmpqv(:) + tmpqcl(:) !bloss/qt: total water
      cloudliq(i,j,:) = tmpqcl(:) !bloss/qt: auxilliary cloud liquid water variable
      if(dopredictNc) micro_field(i,j,:,incl) = tmpncl(:)

      reffc(i,j,:) = effc1d(:)

      if(doicemicro) then
         micro_field(i,j,:,iqci) = tmpqci(:)
         micro_field(i,j,:,inci) = tmpnci(:)
         micro_field(i,j,:,iqs) = tmpqs(:)
         micro_field(i,j,:,ins) = tmpns(:)
         if(dograupel) then
            micro_field(i,j,:,iqg) = tmpqg(:)
            micro_field(i,j,:,ing) = tmpng(:)
         end if
         reffi(i,j,:) = effi1d(:)  
      end if

      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================

      if(dostatis) then
!bloss/qt: total water microphysical tendency includes qv and qcl
         mtend(:,iqv) = mtend(:,iqv) + mtendqv + mtendqcl
!bloss/qt         mtend(:,iqcl) = mtend(:,iqcl) + mtendqcl
         if(dopredictNc) mtend(:,incl) = mtend(:,incl) + mtendncl
         if(doprecip) then
            mtend(:,iqr) = mtend(:,iqr) + mtendqr
            mtend(:,inr) = mtend(:,inr) + mtendnr
         end if

         if(doicemicro) then
            mtend(:,iqci) = mtend(:,iqci) + mtendqci
            mtend(:,inci) = mtend(:,inci) + mtendnci
            !bloss            stend(:,inci) = stend(:,inci) + stendnci

            mtend(:,iqs) = mtend(:,iqs) + mtendqs
            mtend(:,ins) = mtend(:,ins) + mtendns
            !bloss            stend(:,ins) = stend(:,ins) + stendns

            if(dograupel) then
               mtend(:,iqg) = mtend(:,iqg) + mtendqg
               mtend(:,ing) = mtend(:,ing) + mtendng
               !bloss            stend(:,ing) = stend(:,ing) + stendng
            end if
         end if

         do n = 1,nmicro_fields
            do k = 1,nzm
               if(micro_field(i,j,k,n).ge.1.e-6) mfrac(k,n) = mfrac(k,n)+1.
            end do
         end do

         ! approximate optical depth = 0.0018*lwp/effrad
         !  integrated up to level at which output
         tmpc = 0.
         tmpr = 0.
         tmpi = 0.
         tmps = 0.
         tmpg = 0.

         do k = 1,nzm
            tmpc = tmpc + 0.0018*rho(k)*dz*adz(k)*tmpqcl(k)/(1.e-20+1.e-6*effc1d(k))
            tmpr = tmpr + 0.0018*rho(k)*dz*adz(k)*tmpqr(k)/(1.e-20+1.e-6*effr1d(k))
            !bloss/qt: put cloud liquid optical depth in trtau(:,iqv)
            trtau(k,iqv) = trtau(k,iqv) + tmpc
            if(doprecip) trtau(k,iqr) = trtau(k,iqr) + tmpr

            if(doicemicro) then
               tmpi = tmpi + 0.0018*rho(k)*dz*adz(k)*tmpqci(k)/(1.e-20+1.e-6*effi1d(k))
               tmps = tmps + 0.0018*rho(k)*dz*adz(k)*tmpqs(k)/(1.e-20+1.e-6*effs1d(k))
               tmpg = tmpg + 0.0018*rho(k)*dz*adz(k)*tmpqg(k)/(1.e-20+1.e-6*effg1d(k))

               trtau(k,iqci) = trtau(k,iqci) + tmpi
               trtau(k,iqs) = trtau(k,iqs) + tmps
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
               if ( dograupel ) then
                 trtau(k,iqg) = trtau(k,iqg) + tmpg
               end if
#else
               trtau(k,iqg) = trtau(k,iqg) + tmpg
#endif /* CLUBB */
            end if
         end do

         tlat(1:nzm) = tlat(1:nzm) &
              - dtn*fac_cond*(stendqcl+stendqr) &
              - dtn*fac_sub*(stendqci+stendqs+stendqg)
         qpfall(1:nzm) = qpfall(1:nzm) + dtn*(stendqr+stendqs+stendqg)

#ifdef CRM
         qpsrc(1:nzm) = qpsrc(1:nzm) + dtn*(mtendqr+mtendqs+mtendqg)
         qpevp(1:nzm) = 0.0
#endif 

         !bloss: temperature tendency (sensible heating) due to phase changes
         tmtend3d(i,j,1:nzm) = tmtend1d(1:nzm)

      end if ! dostatis

      stend(:,iqv) = stend(:,iqv) + stendqcl !bloss/qt: iqcl --> iqv
      if(doprecip) then
         stend(:,iqr) = stend(:,iqr) + stendqr
      end if

      if(doicemicro) then
         stend(:,iqci) = stend(:,iqci) + stendqci
         stend(:,iqs) = stend(:,iqs) + stendqs
         if(dograupel) stend(:,iqg) = stend(:,iqg) + stendqg
      end if

#ifdef ECPP
      do k=1, nzm
        qlsink_bf(i,j,k)  = min(1.0/dt, QSINK_TMP(k))     ! /s
        rh(i,j,k)      = RH3D(k)       !0-1
        prain(i,j,k) = C2PREC(K)    ! kg/kg/s
        if(cloudliq(i,j,k).gt.1.0e-10) then
          qlsink(i,j,k) = min(1.0/dt, C2PREC(k)/cloudliq(i,j,k))                   
        else
          qlsink(i,j,k) = 0.0 
        end if 
      end do
      precr(i,j,:)=(RSED(:))    ! kg/m2/s
      precsolid(i,j,:)=(SSED(:)+GSED(:))  !kg/m2/s leave ISED out for the momenent, and we may want to
                                    ! test it effects in the future. +++mhwang
#endif /*ECPP*/
#ifdef CLUBB_CRM
if(doclubb) then
   CF3D(i, j, 1:nzm) = cloud_frac_in(1:nzm)
endif
#endif
   end do ! i = 1,nx
end do ! j = 1,ny

! back sedimentation flux out from sedimentation tendencies
tmpc = 0.
do k = 1,nzm
   m = nz-k
   tmpc = tmpc + stend(m,iqv)*rho(m)*dz*adz(m)  !bloss/qt: iqcl --> iqv
   mksed(m,iqv) = tmpc
end do
precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqv)*dtn/dz

if(doprecip) then
   tmpr = 0.
   do k = 1,nzm
      m = nz-k
      tmpr = tmpr + stend(m,iqr)*rho(m)*dz*adz(m)
      mksed(m,iqr) = tmpr
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqr)*dtn/dz
end if

if(doicemicro) then
   tmpi = 0.
   tmps = 0.
   tmpg = 0.
   do k = 1,nzm
      m = nz-k
      tmpi = tmpi + stend(m,iqci)*rho(m)*dz*adz(m)
      tmps = tmps + stend(m,iqs)*rho(m)*dz*adz(m)
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
      else
        tmpg = 0.
      end if
#else
      tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
#endif
      mksed(m,iqci) = tmpi
      mksed(m,iqs) = tmps
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        mksed(m,iqg) = tmpg
      end if
#else
      mksed(m,iqg) = tmpg
#endif
   end do
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
   if ( dograupel ) then
     precflux(1:nzm) = precflux(1:nzm) &
          - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
   else
     precflux(1:nzm) = precflux(1:nzm) &
          - (mksed(:,iqci) + mksed(:,iqs))*dtn/dz
   end if
#else
   precflux(1:nzm) = precflux(1:nzm) &
        - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
#endif
end if

!!$if(doprecip) total_water_prec = total_water_prec - total_water()

#ifdef CLUBB_CRM
if (docloud.or.doclubb)  call micro_diagnose()   ! leave this line here
#else
if (docloud)  call micro_diagnose()   ! leave this line here
#endif

! call t_stopf ('micro_proc')

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
!  This is the pace where the microphysics field that SAM actually cares about
!  are diagnosed.

subroutine micro_diagnose()

use vars
#ifdef CLUBB_CRM
use error_code, only: clubb_at_least_debug_level ! Procedure
use constants_clubb, only: fstderr, zero_threshold
implicit none
#endif

real omn, omp
integer i,j,k

! water vapor = total water - cloud liquid
qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqv) &
     - cloudliq(1:nx,1:ny,1:nzm)

#ifdef CLUBB_CRM
do i = 1, nx
  do j = 1, ny
    do k = 1, nzm
      ! Apply local hole-filling to vapor by converting liquid to vapor. Moist
      ! static energy should be conserved, so updating temperature is not
      ! needed here. -dschanen 31 August 2011
      if ( qv(i,j,k) < zero_threshold ) then
        cloudliq(i,j,k) = cloudliq(i,j,k) + qv(i,j,k)
        qv(i,j,k) = zero_threshold
        if ( cloudliq(i,j,k) < zero_threshold ) then
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Total water at", "i =", i, "j =", j, "k =", k, "is negative.", &
              "Applying non-conservative hard clipping."
          end if
          cloudliq(i,j,k) = zero_threshold
        end if ! cloud_liq < 0
      end if ! qv < 0
    end do ! 1.. nzm
  end do ! 1.. ny
end do ! 1.. nx
#endif /* CLUBB_CRM */
! cloud liquid water
qcl(1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm)

! rain water
if(doprecip) qpl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqr)

! cloud ice 
if(doicemicro) then
   qci(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqci)

   if(dograupel) then
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs) &
           + micro_field(1:nx,1:ny,1:nzm,iqg)
   else
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs)
   end if
end if

end subroutine micro_diagnose

#ifdef CLUBB_CRM
!---------------------------------------------------------------------
subroutine micro_update()

! Description:
! This subroutine essentially does what micro_proc does but does not
! call any microphysics subroutines.  We need to do this for the 
! single-moment bulk microphysics (SAM1MOM) so that CLUBB gets a
! properly updated value of ice fed in. 
!
! -dschanen UWM
!---------------------------------------------------------------------

  ! Update the dynamical core variables (e.g. qv, qcl) with the value in
  ! micro_field.  Diffusion, advection, and other processes are applied to
  ! micro_field but not the variables in vars.f90
  call micro_diagnose()

  return
end subroutine micro_update

!---------------------------------------------------------------------
subroutine micro_adjust( new_qv, new_qc )
! Description:
!   Adjust total water in SAM based on values from CLUBB.
! References:
!   None
!---------------------------------------------------------------------

  use vars, only: qci

  implicit none

  real, dimension(nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  ! Total water mixing ratio
  micro_field(1:nx,1:ny,1:nzm,iqv) = new_qv(1:nx,1:ny,1:nzm) &
                                   + new_qc(1:nx,1:ny,1:nzm)

  ! Cloud water mixing ratio
  cloudliq(1:nx,1:ny,1:nzm) = new_qc(1:nx,1:ny,1:nzm) 

  return
end subroutine micro_adjust

#endif /*CLUBB_CRM*/

!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need supply functions to compute terminal velocity for all of your 
! precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and water vapor (single values, not arrays). Var1 and var2 
! are some microphysics variables like water content and concentration.
! Don't change the number of arguments or their meaning!

!!$real function term_vel_qr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_qr
!!$
!!$real function term_vel_Nr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_Nr
!!$
!!$real function term_vel_qs(qs,ns,tabs,rho)
!!$! .......  
!!$end function term_vel_qs

! etc.

!----------------------------------------------------------------------
!!! compute sedimentation 
!
!  The perpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles
subroutine micro_precip_fall()

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for hour hypothetical case, there is no mixed hydrometeor, so omega is not actually used.

integer hydro_type
real omega(nx,ny,nzm) 

integer i,j,k

return ! do not need this routine -- sedimentation done in m2005micro.

!!$! Initialize arrays that accumulate surface precipitation flux
!!$
!!$ if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
!!$   do j=1,ny
!!$    do i=1,nx
!!$     precsfc(i,j)=0.
!!$    end do
!!$   end do
!!$   do k=1,nzm
!!$    precflux(k) = 0.
!!$   end do
!!$ end if
!!$
!!$ do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
!!$    qpfall(k)=0.
!!$    tlat(k) = 0.
!!$ end do
!!$   
!!$! Compute sedimentation of falling variables:
!!$
!!$ hydro_type=0
!!$ call precip_fall(qr, term_vel_qr, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Nr, term_vel_Nr, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qs, term_vel_qs, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ns, term_vel_Ns, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qg, term_vel_qg, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ng, term_vel_Ng, hydro_type, omega)
!!$


end subroutine micro_precip_fall

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
  implicit none
  integer :: k

  ! print out min/max values of all microphysical variables
  do k=1,nmicro_fields
     call fminmax_print(trim(mkname(k))//':', &
          micro_field(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  end do

end subroutine micro_print

!-----------------------------------------
subroutine satadj_liquid(nzm,tabs,qt,qc,pres)
  !bloss/qt: Utility routine based on cloud.f90 in 
  !  MICRO_SAM1MOM that was written by Marat Khairoutdinov.
  !  This routine performs a saturation adjustment for
  !  cloud liquid water only using a Newton method.
  !  While 20 iterations are allowed, most often this
  !  routine should exit in five iterations or less.
  !  Only a single calculation of the saturation vapor
  !  pressure is required in subsaturated air.

  use module_mp_GRAUPEL, only: polysvp
  use params, only: cp, lcond, rv, fac_cond
  implicit none

  integer, intent(in) :: nzm
  real, intent(inout), dimension(nzm) :: tabs ! absolute temperature, K
  real, intent(inout), dimension(nzm) :: qt  ! on input: qt; on output: qv
  real, intent(out), dimension(nzm) :: qc ! cloud liquid water, kg/kg
  real, intent(in), dimension(nzm) :: pres ! pressure, Pa

  real tabs1, dtabs, thresh, esat1, qsat1, fff, dfff
  integer k, niter

  integer, parameter :: maxiter = 20

  !bloss/qt: quick saturation adjustment to compute cloud liquid water content.
  do k = 1,nzm
    tabs1 = tabs(k) 
    esat1 = polysvp(tabs1,0)
    qsat1 = 0.622*esat1/ (pres(k) - esat1)
    qc(k) = 0. ! no cloud unless qt > qsat
    
    if (qt(k).gt.qsat1) then

      ! if unsaturated, nothing to do (i.e., qv=qt, T=Tl) --> just exit.
      ! if saturated, do saturation adjustment 
      !    (modeled after Marat's cloud.f90).

      ! generate initial guess based on above calculation of qsat
      dtabs = + fac_cond*MAX(0.,qt(k) - qsat1) &
           / ( 1. + lcond**2*qsat1/(cp*rv*tabs1**2) )
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(0.01, 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        esat1 = polysvp(tabs1,0)
        qsat1 = 0.622*esat1/ (pres(k) - esat1) ! saturation mixing ratio

        fff = tabs(k) - tabs1 + fac_cond*MAX(0.,qt(k) - qsat1)
        dfff = 1. + lcond**2*qsat1/(cp*rv*tabs1**2)
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      qc(k) = MAX( 0.,tabs1 - tabs(k) )/fac_cond ! cloud liquid mass mixing ratio
      qt(k) = qt(k) - qc(k) ! This now holds the water vapor mass mixing ratio.
      tabs(k) = tabs1 ! update temperature.
      
      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,nzm

end subroutine satadj_liquid

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi
#ifdef CLUBB_CRM
!-------------------------------------------------------------------------------
ELEMENTAL FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

! This function computes a linear interpolation of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height between those two height levels (rather 
! than a height outside of those two height levels) is computed.
!
! Here is a diagram:
!
!  ################################ Height high, know variable value
!
!
!
!  -------------------------------- Height to be interpolated to; linear interpolation
!
!
!
!
!
!  ################################ Height low, know variable value
!
!
! FORMULA:
!
! variable(@ Height interpolation) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height interpolation - Height low)  +  variable(@ Height low)

! Author: Brian Griffin, UW-Milwaukee
! Modifications: Dave Schanen added the elemental attribute 4 Nov 2008
! References: None

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_int

! Output Variable
REAL:: LIN_INT

LIN_INT = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT
#endif /*CLUBB_CRM*/
!------------------------------------------------------------------------------

end module microphysics



