module microphysics

! main interface to Morrison microphysics.
! original implementation by Peter Blossey, UW

use params, only: lcond, lsub, fac_cond, fac_sub, ggr, crm_rknd

use grid, only: nx,ny,nzm,nz, &  !grid dimensions; nzm = nz-1 # of scalar lvls
     dimx1_s,dimx2_s,dimy1_s,dimy2_s, & ! actual scalar-array dimensions in x,y
     dz, adz, dostatis, masterproc, &
     doSAMconditionals, dosatupdnconditionals

use vars, only: pres, rho, dt, dtn, w, t, tlatqi
use vars, only: tke2, tk2
use params, only: doprecip, docloud, doclubb
use task_util_mod
use utils

use module_mp_GRAUPEL, only: GRAUPEL_INIT, M2005MICRO_GRAUPEL, &
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! use graupel
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
#if (defined MODAL_AERO)
      domodal_aero,     &   ! use modal aerosol from the CAM
#endif
#ifdef CLUBB_CRM
      doclubb_tb,       &   ! use CLUBB as turbulence scheme only, but not cloud scheme,
                            ! so liquid water is diagnosed from saturation adjustment
      doclubb_gridmean,  &  ! feed grid-mean CLUBB values into Morrision microphysics
      doclubb_autoin,    &  ! use in-cloud values for autoconversion calculations
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

  use cam_abortutils,     only: endrun

implicit none

integer :: nmicro_fields ! total number of prognostic water vars

real(crm_rknd), allocatable, dimension(:,:,:,:,:) :: micro_field  ! holds mphys quantities

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing
integer :: index_water_vapor ! separate water vapor index used by SAM

real(crm_rknd), allocatable, dimension(:,:) :: lfac
integer, allocatable, dimension(:) :: flag_precip
integer, allocatable, dimension(:,:) :: flag_wmass, flag_number
integer, allocatable, dimension(:,:) :: flag_micro3Dout

integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

real(crm_rknd), allocatable, dimension(:,:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
real(crm_rknd), allocatable, dimension(:,:,:,:) :: reffc, reffi
real(crm_rknd), allocatable, dimension(:,:,:,:) :: cloudliq

! statistical arrays
real(crm_rknd), allocatable, dimension(:,:,:) :: mkwle ! resolved vertical flux
real(crm_rknd), allocatable, dimension(:,:,:) :: mkwsb ! SGS vertical flux
real(crm_rknd), allocatable, dimension(:,:,:) :: mksed ! sedimentation vertical flux
real(crm_rknd), allocatable, dimension(:,:,:) :: mkadv ! tendency due to vertical advection
real(crm_rknd), allocatable, dimension(:,:,:) :: mkdiff! tendency due to vertical diffusion
real(crm_rknd), allocatable, dimension(:,:,:) :: mklsadv ! tendency due to large-scale vertical advection
real(crm_rknd), allocatable, dimension(:,:,:) :: mfrac ! fraction of domain with microphysical quantity > 1.e-6
real(crm_rknd), allocatable, dimension(:,:,:) :: stend ! tendency due to sedimentation
real(crm_rknd), allocatable, dimension(:,:,:) :: mtend ! tendency due to microphysical processes (other than sedimentation)
real(crm_rknd), allocatable, dimension(:,:,:) :: mstor ! storage terms of microphysical variables
real(crm_rknd), allocatable, dimension(:,:,:) :: trtau    ! optical depths of various species

real(crm_rknd), allocatable, dimension(:,:) :: tmtend

real(crm_rknd) :: sfcpcp, sfcicepcp

! arrays with names/units for microphysical outputs in statistics.
character*3, allocatable, dimension(:) :: mkname
character*80, allocatable, dimension(:) :: mklongname
character*10, allocatable, dimension(:) :: mkunits
real(crm_rknd), allocatable, dimension(:) :: mkoutputscale
logical douse_reffc, douse_reffi

! You can also have some additional, diagnostic, arrays, for example, total
! nonprecipitating cloud water, etc:

!bloss: array which holds temperature tendency due to microphysics
real(crm_rknd), allocatable, dimension(:,:,:,:) :: tmtend3d

real(crm_rknd), allocatable, dimension(:,:) ::  qpevp   !sink of precipitating water due to evaporation (set to zero here)
real(crm_rknd), allocatable, dimension(:,:) ::  qpsrc   !source of precipitation microphysical processes (set to mtend)

real(crm_rknd), allocatable, dimension(:,:,:,:)  :: wvar  ! the vertical velocity variance from subgrid-scale motion,
                                              ! which is needed in droplet activation.
! hm 7/26/11 new output
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: aut1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: acc1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: evpc1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: evpr1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: mlt1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: sub1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: dep1  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: con1  !

real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: aut1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: acc1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: evpc1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: evpr1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: mlt1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: sub1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: dep1a  !
real(crm_rknd), public, allocatable, dimension(:,:,:,:)  :: con1a  !

!+++mhwangtest
! test water conservation
real(crm_rknd), public, allocatable, dimension(:,:,:) ::  sfcpcp2D  ! surface precipitation
!---mhwangtest

CONTAINS


subroutine allocate_micro(ncrms)
  implicit none
  integer, intent(in) :: ncrms
  ! allocate microphysical variables
  allocate(micro_field(ncrms,dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,nmicro_fields))
  allocate(fluxbmk(ncrms,nx,ny,nmicro_fields))
  allocate(fluxtmk(ncrms,nx,ny,nmicro_fields))
  allocate(reffc(nx,ny,nzm,ncrms))
  allocate(reffi(nx,ny,nzm,ncrms))
  allocate(mkwle(ncrms,nz,nmicro_fields))
  allocate(mkwsb(ncrms,nz,nmicro_fields))
  allocate(mkadv(ncrms,nz,nmicro_fields))
  allocate(mkdiff(ncrms,nz,nmicro_fields))
  allocate(mklsadv(nz,nmicro_fields,ncrms))
  allocate(stend(nzm,nmicro_fields,ncrms))
  allocate(mtend(nzm,nmicro_fields,ncrms))
  allocate(mfrac(nzm,nmicro_fields,ncrms))
  allocate(trtau(nzm,nmicro_fields,ncrms))
  allocate(mksed(nzm,nmicro_fields,ncrms))
  allocate(tmtend(nzm,ncrms))
  allocate(mstor(nzm,nmicro_fields,ncrms))
  allocate(cloudliq(ncrms,nx,ny,nzm))
  allocate(tmtend3d(nx,ny,nzm,ncrms))
  allocate(flag_micro3Dout(nmicro_fields,ncrms))
  allocate(flag_wmass(nmicro_fields,ncrms))
  allocate(flag_precip(nmicro_fields))
  allocate(flag_number(nmicro_fields,ncrms))
  allocate(lfac(nmicro_fields,ncrms))
  allocate(mkname(nmicro_fields))
  allocate(mklongname(nmicro_fields))
  allocate(mkunits(nmicro_fields))
  allocate(mkoutputscale(nmicro_fields))
  allocate(wvar(ncrms,nx,ny,nzm))
  allocate(qpevp(ncrms,nz))
  allocate(qpsrc(ncrms,nz))
  allocate(aut1(ncrms,nx,ny,nzm))
  allocate(acc1(ncrms,nx,ny,nzm))
  allocate(evpc1(ncrms,nx,ny,nzm))
  allocate(evpr1(ncrms,nx,ny,nzm))
  allocate(mlt1(ncrms,nx,ny,nzm))
  allocate(sub1(ncrms,nx,ny,nzm))
  allocate(dep1(ncrms,nx,ny,nzm))
  allocate(con1(ncrms,nx,ny,nzm))
  allocate(aut1a(ncrms,nx,ny,nzm))
  allocate(acc1a(ncrms,nx,ny,nzm))
  allocate(evpc1a(ncrms,nx,ny,nzm))
  allocate(evpr1a(ncrms,nx,ny,nzm))
  allocate(mlt1a(ncrms,nx,ny,nzm))
  allocate(sub1a(ncrms,nx,ny,nzm))
  allocate(dep1a(ncrms,nx,ny,nzm))
  allocate(con1a(ncrms,nx,ny,nzm))
  allocate(sfcpcp2D(nx,ny,ncrms))

  ! initialize these arrays
  micro_field = 0.
  cloudliq = 0. !bloss/qt: auxially cloud liquid water variable, analogous to qn in MICRO_SAM1MOM
  fluxbmk = 0.
  fluxtmk = 0.
  mkwle = 0.
  reffc = 0.
  reffi = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.
  mklsadv = 0.
  stend = 0.
  mtend = 0.
  mfrac = 0.
  trtau = 0.
  mksed = 0.
  tmtend = 0.
  tmtend3d = 0.
  mstor =0.
  wvar = 0.
  lfac = 0.
  ! initialize flag arrays to all mass, no number, no precip
  flag_wmass = 1
  flag_number = 0
  flag_precip = 0
  flag_micro3Dout = 0
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
  ! crjones notes: qpsrc must be initialized to zero before start of crm main time loop
  qpsrc = 0.
  qpevp = 0.
end subroutine allocate_micro


subroutine deallocate_micro()
  implicit none
  deallocate(micro_field)
  deallocate(fluxbmk)
  deallocate(fluxtmk)
  deallocate(reffc)
  deallocate(reffi)
  deallocate(mkwle)
  deallocate(mkwsb)
  deallocate(mkadv)
  deallocate(mkdiff)
  deallocate(mklsadv)
  deallocate(stend)
  deallocate(mtend)
  deallocate(mfrac)
  deallocate(trtau)
  deallocate(mksed)
  deallocate(tmtend)
  deallocate(mstor)
  deallocate(cloudliq)
  deallocate(tmtend3d)
  deallocate(flag_micro3Dout)
  deallocate(flag_wmass)
  deallocate(flag_precip)
  deallocate(flag_number)
  deallocate(lfac)
  deallocate(mkname)
  deallocate(mklongname)
  deallocate(mkunits)
  deallocate(mkoutputscale)
  deallocate (wvar)
  deallocate (qpevp)
  deallocate (qpsrc)
  deallocate (aut1)
  deallocate (acc1)
  deallocate (evpc1)
  deallocate (evpr1)
  deallocate (mlt1)
  deallocate (sub1)
  deallocate (dep1)
  deallocate (con1)
  deallocate (aut1a)
  deallocate (acc1a)
  deallocate (evpc1a)
  deallocate (evpr1a)
  deallocate (mlt1a)
  deallocate (sub1a)
  deallocate (dep1a)
  deallocate (con1a)
  deallocate (sfcpcp2D)
end subroutine deallocate_micro

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
     docloudedgeactivation  = .false. ! activate droplets at cloud base, and edges
     doclubb_tb = .false.
     doclubb_gridmean = .true.
     doclubb_autoin = .false.
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

  compute_reffc = douse_reffc
  compute_reffi = douse_reffi

end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the
!   beginning of each run, initial or restart:
subroutine micro_init(ncrms)

  use vars
#if (defined MODAL_AERO)
  use drop_activation, only: drop_activation_init
#endif

  implicit none
  integer, intent(in) :: ncrms

  real(crm_rknd), dimension(nzm) :: qc0, qi0

! Commented out by dschanen UWM 23 Nov 2009 to avoid a linking error
! real, external :: satadj_water
  integer :: k,icrm

  do icrm = 1 , ncrms

  ! initialize flag arrays
  if(dopredictNc) then
     ! Cloud droplet number concentration is a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass (:,icrm) = (/1,0,1,0,1,0,1,0,1,0/)
           flag_precip(:) = (/0,0,1,1,0,0,1,1,1,1/)
           flag_number(:,icrm) = (/0,1,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns
           flag_wmass (:,icrm) = (/1,0,1,0,1,0,1,0/)
           flag_precip(:) = (/0,0,1,1,0,0,1,1/)
           flag_number(:,icrm) = (/0,1,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, Nc, qr, Nr
           flag_wmass (:,icrm) = (/1,0,1,0/)
           flag_precip(:) = (/0,0,1,1/)
           flag_number(:,icrm) = (/0,1,0,1/)
        else
          !bloss/qt: qt, Nc
           flag_wmass (:,icrm) = (/1,0/)
           flag_precip(:) = (/0,0/)
           flag_number(:,icrm) = (/0,1/)
        end if
     end if
  else
     ! Cloud droplet number concentration is NOT a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass (:,icrm) = (/1,1,0,1,0,1,0,1,0/)
           flag_precip(:) = (/0,1,1,0,0,1,1,1,1/)
           flag_number(:,icrm) = (/0,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns
           flag_wmass (:,icrm) = (/1,1,0,1,0,1,0/)
           flag_precip(:) = (/0,1,1,0,0,1,1/)
           flag_number(:,icrm) = (/0,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, qr, Nr
           flag_wmass (:,icrm) = (/1,1,0/)
           flag_precip(:) = (/0,1,1/)
           flag_number(:,icrm) = (/0,0,1/)
        else
          !bloss/qt: only total water variable is needed for no-precip,
          !            fixed droplet number, warm cloud and no cloud simulations.
           flag_wmass (:,icrm) = (/1/)
           flag_precip(:) = (/0/)
           flag_number(:,icrm) = (/0/)
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
     flag_micro3Dout(:,icrm) = 1
  end if

  ! initialize factor for latent heat
  lfac(:,icrm) = 1. ! use one as default for number species
  lfac(iqv,icrm) = lcond
!bloss/qt  if(docloud) lfac(iqcl,icrm) = lcond
  if(doprecip) lfac(iqr,icrm) = lcond
  if(doicemicro) then
     lfac(iqci,icrm) = lsub
     lfac(iqs,icrm) = lsub
     if(dograupel) lfac(iqg,icrm) = lsub
  end if

  call graupel_init() ! call initialization routine within mphys module
#if (defined MODAL_AERO)
  call drop_activation_init
#endif

  if(nrestart.eq.0) then

#ifdef CLUBB_CRM
     if(docloud.or.doclubb)  call micro_diagnose(ncrms,icrm)   ! leave this line here
#else
     if(docloud) call micro_diagnose(ncrms,icrm)   ! leave this here
#endif


  end if

  enddo

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux(ncrms)
  use vars, only: fluxbq, fluxtq
#ifdef CLUBB_CRM
  use params, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
#endif
  implicit none
  integer, intent(in) :: ncrms
  integer :: icrm

  do icrm = 1 , ncrms
    fluxbmk(icrm,:,:,:) = 0. ! initialize all fluxes at surface to zero
    fluxtmk(icrm,:,:,:) = 0. ! initialize all fluxes at top of domain to zero
#ifdef CLUBB_CRM
    if ( doclubb .and. (doclubb_sfc_fluxes.or.docam_sfc_fluxes) ) then
      fluxbmk(icrm,:,:,index_water_vapor) = 0.0 ! surface qv(icrm,latent heat) flux
    else
      fluxbmk(icrm,:,:,index_water_vapor) = fluxbq(icrm,:,:) ! surface qv(icrm,latent heat) flux
    end if
#else
    fluxbmk(icrm,:,:,index_water_vapor) = fluxbq(icrm,:,:) ! surface qv(icrm,latent heat) flux
#endif
    fluxtmk(icrm,:,:,index_water_vapor) = fluxtq(icrm,:,:) ! top of domain qv flux
  enddo
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

subroutine micro_proc(ncrms)

use params, only: fac_cond, fac_sub, rgas
use grid, only: z, zi
use vars, only: t, gamaz, precsfc, precssfc, precflux, qpfall, tlat, prec_xy, &
                nstep, nstatis, icycle, total_water_prec
#ifdef ECPP
use ecppvars, only: qlsink, qlsink_bf, prain, precr, precsolid, rh, qcloud_bf
#endif

#ifdef CLUBB_CRM
use params, only: doclubb, docloud, dosmoke
use grid, only: nz
use error_code, only: clubb_at_least_debug_level
use fill_holes, only: fill_holes_driver
use clubbvars, only: wp2, cloud_frac, rho_ds_zt, rho_ds_zm, relvarg, accre_enhang ! are used, but not modified here
use vars, only: qcl ! Used here and updated in micro_diagnose
use vars, only: prespot ! exner^-1
use module_mp_GRAUPEL, only: &
  cloud_frac_thresh ! Threshold for using sgs cloud fraction to weight
                    ! microphysical quantities [%]
use clubb_precision, only: core_rknd
use constants_clubb, only: T_freeze_K
use vars, only: CF3D
#endif
implicit none
  integer, intent(in) :: ncrms


real(crm_rknd), dimension(nzm) :: &
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
real(crm_rknd), dimension(nzm) :: C2PREC,QSINK_TMP, CSED,ISED,SSED,GSED,RSED,RH3D   ! used for cloud chemistry and wet deposition in ECPP
#endif

#ifdef CLUBB_CRM
real(kind=core_rknd), dimension(nz) :: &
     qv_clip, qcl_clip
real(crm_rknd), dimension(nzm) :: cloud_frac_in, ice_cldfrac
real(crm_rknd), dimension(nzm) :: liq_cldfrac
real(crm_rknd), dimension(nzm) :: relvar        ! relative cloud water variance
real(crm_rknd), dimension(nzm) :: accre_enhan   ! optional accretion enhancement factor for MG
#endif /*CLUBB_CRM*/

real(crm_rknd), dimension(nzm,nmicro_fields) :: stend1d, mtend1d
real(crm_rknd) :: tmpc, tmpr, tmpi, tmps, tmpg
integer :: i1, i2, j1, j2, i, j, k, m, n, icrm

real(8) :: tmp_total, tmptot

do icrm = 1 , ncrms

! call t_startf ('micro_proc')

if(dostatis) then ! initialize arrays for statistics
   mfrac(:,:,icrm) = 0.
   mtend(:,:,icrm) = 0.
   trtau(:,:,icrm) = 0.
!   qpfall(icrm,:)=0.     ! in SPCAM, done in crm.F90
   tlat(icrm,:) = 0.
   tmtend3d(:,:,:,icrm) = 0.
end if
stend(:,:,icrm) = 0.
mksed(:,:,icrm) = 0.

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
      tmpqv(:) = micro_field(icrm,i,j,:,iqv) !bloss/qt: This is total water (qv+qcl)
!bloss/qt: compute below from saturation adjustment.
      if(dopredictNc) tmpncl(:) = micro_field(icrm,i,j,:,incl)
      if(doprecip) then
         tmpqr(:) = micro_field(icrm,i,j,:,iqr)
         tmpnr(:) = micro_field(icrm,i,j,:,inr)
      end if

      if(doicemicro) then
         tmpqci(:) = micro_field(icrm,i,j,:,iqci)
         tmpnci(:) = micro_field(icrm,i,j,:,inci)
         tmpqs(:) = micro_field(icrm,i,j,:,iqs)
         tmpns(:) = micro_field(icrm,i,j,:,ins)
         if(dograupel) then
            tmpqg(:) = micro_field(icrm,i,j,:,iqg)
            tmpng(:) = micro_field(icrm,i,j,:,ing)
         end if
      end if

      ! get absolute temperature in this column
      !bloss/qt: before saturation adjustment for liquid,
      !          this is Tcl = T - (L/Cp)*qcl(icrm,the cloud liquid water temperature)
      tmptabs(:) = t(icrm,i,j,:)  &           ! liquid water-ice static energy over Cp
           - gamaz(icrm,:) &                                   ! potential energy
           + fac_cond * (tmpqr(:)) &    ! bloss/qt: liquid latent energy due to rain only
           + fac_sub  * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy

      tmpdz = adz(icrm,:)*dz(icrm)
!      tmpw = 0.5*(w(icrm,i,j,1:nzm) + w(icrm,i,j,2:nz))  ! MK: changed for stretched grids
      tmpw = ((zi(icrm,2:nz)-z(icrm,1:nzm))*w(icrm,i,j,1:nzm)+ &
             (z(icrm,1:nzm)-zi(icrm,1:nzm))*w(icrm,i,j,2:nz))/(zi(icrm,2:nz)-zi(icrm,1:nzm))
#ifdef CLUBB_CRM
      ! Added by dschanen on 4 Nov 2008 to account for w_sgs
      if ( doclubb .and. dosubgridw ) then
        ! Compute w_sgs.  Formula is consistent with that used with
        ! TKE from MYJ pbl scheme in WRF (see module_mp_graupel.f90).
        tmpwsub = sqrt( LIN_INT( real( wp2(i,j,2:nz) ,crm_rknd), real( wp2(i,j,1:nzm) ,crm_rknd), &
                                  zi(icrm,2:nz), zi(icrm,1:nzm), z(icrm,1:nzm) ) )
      else
!        tmpwsub = 0.
! diagnose tmpwsub from tke.
! Notes: tke has to be already prognsotic or diagnostic.
        tmpwsub = sqrt(tke2(icrm,i,j,:)/3.)  ! diagnosed tmpwsub from tke
! diagnose tmpwsub from tk
!        tmpwsub = sqrt(2*3.141593)*tk(icrm,i,j,:)/(dz(icrm)*adz(icrm,:))  ! from Ghan et al. (1997, JGR).
      end if

      if ( doclubb ) then
        cloud_frac_in(1:nzm) = cloud_frac(i,j,2:nz)
        liq_cldfrac(1:nzm) = cloud_frac(i,j,2:nz)
      else
        cloud_frac_in(1:nzm) = 0.0
      end if

#else /* Old code */
!      tmpwsub = 0.
! diagnose tmpwsub from tke.
! Notes: tke has to be already prognsotic or diagnostic.
      tmpwsub = sqrt(tke2(icrm,i,j,:)/3.)  ! diagnosed tmpwsub from tke
! diagnose tmpwsub from tk
!      tmpwsub = sqrt(2*3.141593)*tk(icrm,i,j,:)/(dz(icrm)*adz(icrm,:))  ! from Ghan et al. (1997, JGR).
#endif
      wvar(icrm,i,j,:) = tmpwsub(:)

      tmppres(:) = 100.*pres(icrm,1:nzm)

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
         tmpqcl  = cloudliq(icrm,i,j,:) ! Liquid updated by CLUBB just prior to this
         tmpqv   = tmpqv - tmpqcl ! Vapor
         tmptabs = tmptabs + fac_cond * tmpqcl ! Update temperature
         if(doclubb_gridmean) then
           cloud_frac_in(1:nzm) = 0.0  ! to use grid mean for Morrison microphysics, just
                                       ! simply set cloud_frac_in to be zero.
           liq_cldfrac(1:nzm) = cloud_frac(i,j,2:nz)

           cf3d(icrm,i, j, 1:nzm) = cloud_frac(i, j, 2:nz)
           ice_cldfrac(:)= 0.0
           if(doicemicro) then
             do k=1, nzm
               if(tmpqci(k).gt.1.0e-8) then
                 ice_cldfrac(k) = 1.0
               end if
               if((tmpqcl(k) + tmpqci(k)).gt.1.0e-9) then
                 cf3d(icrm,i,j,k) = (cf3d(icrm,i,j,k) * tmpqcl(k) + ice_cldfrac(k) * tmpqci(k))  &
                           / (tmpqcl(k) + tmpqci(k))
               else
                 cf3d(icrm,i,j,k) = 0.0
               end if
               ice_cldfrac(k) = max(cf3d(icrm,i,j,k), liq_cldfrac(k))
             end do
           endif
         end if
        else
          call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
          cloudliq(icrm,i,j,:) = tmpqcl
          cloud_frac_in(1:nzm) = 0.0
        end if
      else
        call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
      end if
#else
      call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
#endif


#ifdef ECPP
! save cloud water before microphysics process for the calculation
! of qlsink in ECPP
      qcloud_bf(i,j,:,icrm) =  tmpqcl(:)
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

      sfcpcp2D(:,:,icrm) = 0.0  !+++mhwangtest

      effc1d(:) = 10. ! default liquid and ice effective radii
      effi1d(:) = 75.

#ifdef CLUBB_CRM
      relvar(:) = 8.
      accre_enhan(:) = 1.
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

   ! ------------------------------------------------- !
   ! Diagnose relative cloud water variance            !
   ! ------------------------------------------------- !
!        relvar(:) = 1.0  ! default
!        where (tmpqcl(:) /= 0. .and. qclvar(i,j, :) /= 0.) &
!          relvar(:) = min(8.0, max(0.35,tmpqcl(:)**2/qclvar(i,j,:)))
!          relvar(:) = min(8.0, max(0.35,tmpqcl(:)**2/qclvar(i,j,:)))

   ! ------------------------------------------------- !
   ! Optional Accretion enhancement factor             !
   ! ------------------------------------------------- !
!        accre_enhan(:) = 1.+0.65*(1.0/relvar(:))
        relvar(:) = relvarg(i,j,:)
        accre_enhan(:) = accre_enhang(i,j,:)
      end if ! doclubb

      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/o sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           ncrms,icrm,mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho(icrm,:),tmpdz,tmpw,tmpwsub, &
! hm 7/26/11, new output
           tmpacc,tmpaut,tmpevpc,tmpevpr,tmpmlt, &
           tmpsub,tmpdep,tmpcon, &
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
          stendqr,stendqci,stendqs,stendqcl,cloud_frac_in, liq_cldfrac, ice_cldfrac, relvar, accre_enhan  & ! cloud_frac added by dschanen UWM
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
           ncrms,icrm,mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho(icrm,:),tmpdz,tmpw,tmpwsub, &
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

! hm 7/26/11, new output
      aut1(icrm,i,j,:) = tmpaut(:)
      acc1(icrm,i,j,:) = tmpacc(:)
      evpc1(icrm,i,j,:) = tmpevpc(:)
      evpr1(icrm,i,j,:) = tmpevpr(:)
      mlt1(icrm,i,j,:) = tmpmlt(:)
      sub1(icrm,i,j,:) = tmpsub(:)
      dep1(icrm,i,j,:) = tmpdep(:)
      con1(icrm,i,j,:) = tmpcon(:)

! hm 8/31/11, new output for gcm-grid and time-step avg
! rates are summed here over the icycle loop
! note: rates are multiplied by time step, and then
! divided by dt in crm.F90 to get mean rates
      aut1a(icrm,i,j,:) = aut1a(icrm,i,j,:) + aut1(icrm,i,j,:)*dtn
      acc1a(icrm,i,j,:) = acc1a(icrm,i,j,:) + acc1(icrm,i,j,:)*dtn
      evpc1a(icrm,i,j,:) = evpc1a(icrm,i,j,:) + evpc1(icrm,i,j,:)*dtn
      evpr1a(icrm,i,j,:) = evpr1a(icrm,i,j,:) + evpr1(icrm,i,j,:)*dtn
      mlt1a(icrm,i,j,:) = mlt1a(icrm,i,j,:) + mlt1(icrm,i,j,:)*dtn
      sub1a(icrm,i,j,:) = sub1a(icrm,i,j,:) + sub1(icrm,i,j,:)*dtn
      dep1a(icrm,i,j,:) = dep1a(icrm,i,j,:) + dep1(icrm,i,j,:)*dtn
      con1a(icrm,i,j,:) = con1a(icrm,i,j,:) + con1(icrm,i,j,:)*dtn

     ! update microphysical quantities in this grid column
      if(doprecip) then
         total_water_prec(icrm) = total_water_prec(icrm) + sfcpcp

         ! take care of surface precipitation
         precsfc(icrm,i,j) = precsfc(icrm,i,j) + sfcpcp/dz(icrm)
         prec_xy(icrm,i,j) = prec_xy(icrm,i,j) + sfcpcp/dtn/dz(icrm)
!+++mhwang
         sfcpcp2D(i,j,icrm) = sfcpcp/dtn/dz(icrm)
!---mhwang
         precssfc(icrm,i,j) = precssfc(icrm,i,j) + sfcicepcp/dz(icrm)    ! the corect unit of precssfc should be mm/dz +++mhwang
         ! update rain
         micro_field(icrm,i,j,:,iqr) = tmpqr(:)
         micro_field(icrm,i,j,:,inr) = tmpnr(:)
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
      micro_field(icrm,i,j,:,iqv) = tmpqv(:) + tmpqcl(:) !bloss/qt: total water
      cloudliq(icrm,i,j,:) = tmpqcl(:) !bloss/qt: auxilliary cloud liquid water variable
      if(dopredictNc) micro_field(icrm,i,j,:,incl) = tmpncl(:)

      reffc(i,j,:,icrm) = effc1d(:)

      if(doicemicro) then
         micro_field(icrm,i,j,:,iqci) = tmpqci(:)
         micro_field(icrm,i,j,:,inci) = tmpnci(:)
         micro_field(icrm,i,j,:,iqs) = tmpqs(:)
         micro_field(icrm,i,j,:,ins) = tmpns(:)
         if(dograupel) then
            micro_field(icrm,i,j,:,iqg) = tmpqg(:)
            micro_field(icrm,i,j,:,ing) = tmpng(:)
         end if
         reffi(i,j,:,icrm) = effi1d(:)
      end if

      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(icrm,i,j,:) = t(icrm,i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================

      if(dostatis) then
!bloss/qt: total water microphysical tendency includes qv and qcl
         mtend(:,iqv,icrm) = mtend(:,iqv,icrm) + mtendqv + mtendqcl
!bloss/qt         mtend(:,iqcl,icrm) = mtend(:,iqcl,icrm) + mtendqcl
         if(dopredictNc) mtend(:,incl,icrm) = mtend(:,incl,icrm) + mtendncl
         if(doprecip) then
            mtend(:,iqr,icrm) = mtend(:,iqr,icrm) + mtendqr
            mtend(:,inr,icrm) = mtend(:,inr,icrm) + mtendnr
         end if

         if(doicemicro) then
            mtend(:,iqci,icrm) = mtend(:,iqci,icrm) + mtendqci
            mtend(:,inci,icrm) = mtend(:,inci,icrm) + mtendnci
            !bloss            stend(:,inci,icrm) = stend(:,inci,icrm) + stendnci

            mtend(:,iqs,icrm) = mtend(:,iqs,icrm) + mtendqs
            mtend(:,ins,icrm) = mtend(:,ins,icrm) + mtendns
            !bloss            stend(:,ins,icrm) = stend(:,ins,icrm) + stendns

            if(dograupel) then
               mtend(:,iqg,icrm) = mtend(:,iqg,icrm) + mtendqg
               mtend(:,ing,icrm) = mtend(:,ing,icrm) + mtendng
               !bloss            stend(:,ing,icrm) = stend(:,ing,icrm) + stendng
            end if
         end if

         do n = 1,nmicro_fields
            do k = 1,nzm
               if(micro_field(icrm,i,j,k,n).ge.1.e-6) mfrac(k,n,icrm) = mfrac(k,n,icrm)+1.
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
            tmpc = tmpc + 0.0018*rho(icrm,k)*dz(icrm)*adz(icrm,k)*tmpqcl(k)/(1.e-20+1.e-6*effc1d(k))
            tmpr = tmpr + 0.0018*rho(icrm,k)*dz(icrm)*adz(icrm,k)*tmpqr(k)/(1.e-20+1.e-6*effr1d(k))
            !bloss/qt: put cloud liquid optical depth in trtau(:,iqv,icrm)
            trtau(k,iqv,icrm) = trtau(k,iqv,icrm) + tmpc
            if(doprecip) trtau(k,iqr,icrm) = trtau(k,iqr,icrm) + tmpr

            if(doicemicro) then
               tmpi = tmpi + 0.0018*rho(icrm,k)*dz(icrm)*adz(icrm,k)*tmpqci(k)/(1.e-20+1.e-6*effi1d(k))
               tmps = tmps + 0.0018*rho(icrm,k)*dz(icrm)*adz(icrm,k)*tmpqs(k)/(1.e-20+1.e-6*effs1d(k))
               tmpg = tmpg + 0.0018*rho(icrm,k)*dz(icrm)*adz(icrm,k)*tmpqg(k)/(1.e-20+1.e-6*effg1d(k))

               trtau(k,iqci,icrm) = trtau(k,iqci,icrm) + tmpi
               trtau(k,iqs,icrm) = trtau(k,iqs,icrm) + tmps
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
               if ( dograupel ) then
                 trtau(k,iqg,icrm) = trtau(k,iqg,icrm) + tmpg
               end if
#else
               trtau(k,iqg,icrm) = trtau(k,iqg,icrm) + tmpg
#endif /* CLUBB */
            end if
         end do

         tlat(icrm,1:nzm) = tlat(icrm,1:nzm) &
              - dtn*fac_cond*(stendqcl+stendqr) &
              - dtn*fac_sub*(stendqci+stendqs+stendqg)
         qpfall(icrm,1:nzm) = qpfall(icrm,1:nzm) + dtn*(stendqr+stendqs+stendqg)

         qpsrc(icrm,1:nzm) = qpsrc(icrm,1:nzm) + dtn*(mtendqr+mtendqs+mtendqg)
         qpevp(icrm,1:nzm) = 0.0

         !bloss: temperature tendency (sensible heating) due to phase changes
         tmtend3d(i,j,1:nzm,icrm) = tmtend1d(1:nzm)

      end if ! dostatis

      stend(:,iqv,icrm) = stend(:,iqv,icrm) + stendqcl !bloss/qt: iqcl --> iqv
      if(doprecip) then
         stend(:,iqr,icrm) = stend(:,iqr,icrm) + stendqr
      end if

      if(doicemicro) then
         stend(:,iqci,icrm) = stend(:,iqci,icrm) + stendqci
         stend(:,iqs,icrm) = stend(:,iqs,icrm) + stendqs
         if(dograupel) stend(:,iqg,icrm) = stend(:,iqg,icrm) + stendqg
      end if

#ifdef ECPP
      do k=1, nzm
        qlsink_bf(i,j,k,icrm)  = min(1.0/dt, QSINK_TMP(k))     ! /s
        rh(i,j,k,icrm)      = RH3D(k)       !0-1
        prain(i,j,k,icrm) = C2PREC(K)    ! kg/kg/s
        if(cloudliq(icrm,i,j,k).gt.1.0e-10) then
          qlsink(i,j,k,icrm) = min(1.0/dt, C2PREC(k)/cloudliq(icrm,i,j,k))
        else
          qlsink(i,j,k,icrm) = 0.0
        end if
      end do
      precr(i,j,:,icrm)=(RSED(:))    ! kg/m2/s
      precsolid(i,j,:,icrm)=(SSED(:)+GSED(:))  !kg/m2/s leave ISED out for the momenent, and we may want to
                                    ! test it effects in the future. +++mhwang
#endif /*ECPP*/

   end do ! i = 1,nx
end do ! j = 1,ny

! back sedimentation flux out from sedimentation tendencies
tmpc = 0.
do k = 1,nzm
   m = nz-k
   tmpc = tmpc + stend(m,iqv,icrm)*rho(icrm,m)*dz(icrm)*adz(icrm,m)  !bloss/qt: iqcl --> iqv
   mksed(m,iqv,icrm) = tmpc
end do
precflux(icrm,1:nzm) = precflux(icrm,1:nzm) - mksed(:,iqv,icrm)*dtn/dz(icrm)

if(doprecip) then
   tmpr = 0.
   do k = 1,nzm
      m = nz-k
      tmpr = tmpr + stend(m,iqr,icrm)*rho(icrm,m)*dz(icrm)*adz(icrm,m)
      mksed(m,iqr,icrm) = tmpr
   end do
   precflux(icrm,1:nzm) = precflux(icrm,1:nzm) - mksed(:,iqr,icrm)*dtn/dz(icrm)
end if

if(doicemicro) then
   tmpi = 0.
   tmps = 0.
   tmpg = 0.
   do k = 1,nzm
      m = nz-k
      tmpi = tmpi + stend(m,iqci,icrm)*rho(icrm,m)*dz(icrm)*adz(icrm,m)
      tmps = tmps + stend(m,iqs,icrm)*rho(icrm,m)*dz(icrm)*adz(icrm,m)
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        tmpg = tmpg + stend(m,iqg,icrm)*rho(icrm,m)*dz(icrm)*adz(icrm,m)
      else
        tmpg = 0.
      end if
#else
      tmpg = tmpg + stend(m,iqg,icrm)*rho(icrm,m)*dz(icrm)*adz(icrm,m)
#endif
      mksed(m,iqci,icrm) = tmpi
      mksed(m,iqs,icrm) = tmps
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        mksed(m,iqg,icrm) = tmpg
      end if
#else
      mksed(m,iqg,icrm) = tmpg
#endif
   end do
#ifdef CLUBB_CRM /* Bug fix -dschanen 9 Mar 2012 */
   if ( dograupel ) then
     precflux(icrm,1:nzm) = precflux(icrm,1:nzm) &
          - (mksed(:,iqci,icrm) + mksed(:,iqs,icrm) + mksed(:,iqg,icrm))*dtn/dz(icrm)
   else
     precflux(icrm,1:nzm) = precflux(icrm,1:nzm) &
          - (mksed(:,iqci,icrm) + mksed(:,iqs,icrm))*dtn/dz(icrm)
   end if
#else
   precflux(icrm,1:nzm) = precflux(icrm,1:nzm) &
        - (mksed(:,iqci,icrm) + mksed(:,iqs,icrm) + mksed(:,iqg,icrm))*dtn/dz(icrm)
#endif
end if

!!$if(doprecip) total_water_prec = total_water_prec - total_water()

#ifdef CLUBB_CRM
if (docloud.or.doclubb)  call micro_diagnose(ncrms,icrm)   ! leave this line here
if(doclubb) then
   cf3d(icrm,1:nx, 1:ny, 1:nzm) = cloud_frac(1:nx, 1:ny, 2:nzm+1)
   if(doicemicro) then
     do i=1, nx
       do j=1, ny
          ice_cldfrac(:) = 0.0
          do k=1, nzm
! Ice cloud fraction: 0 at 0 C, and 100% at -35C.
!           ice_cldfrac(k) = -(tmptabs(k)-T_freeze_K)/35.0
!           ice_cldfrac(k) = min(1.0, max(ice_cldfrac(k), 0.0))
           if(micro_field(icrm,i,j,k,iqci) .gt. 1.0e-8) then
             ice_cldfrac(k) = 1.0
           end if
           if(cloudliq(icrm,i,j,k) + micro_field(icrm,i,j,k,iqci) .gt.1.0e-9) then
             cf3d(icrm,i,j,k) = (cf3d(icrm,i,j,k)* cloudliq(icrm,i,j,k) + ice_cldfrac(k) * micro_field(icrm,i,j,k,iqci))  &
                           / (cloudliq(icrm,i,j,k) + micro_field(icrm,i,j,k,iqci))
           else
             cf3d(icrm,i,j,k) = 0.0
           end if
          end do
        end do
     end do
   endif
endif
#else
if (docloud)  call micro_diagnose(ncrms,icrm)   ! leave this line here
#endif

! call t_stopf ('micro_proc')

enddo !icrm

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
!  This is the pace where the microphysics field that SAM actually cares about
!  are diagnosed.

subroutine micro_diagnose(ncrms,icrm)

use vars
#ifdef CLUBB_CRM
use error_code, only: clubb_at_least_debug_level ! Procedure
use constants_clubb, only: fstderr, zero_threshold
#endif
implicit none
integer, intent(in) :: ncrms,icrm

real(crm_rknd) omn, omp
integer i,j,k

! water vapor = total water - cloud liquid
qv(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqv) &
     - cloudliq(icrm,1:nx,1:ny,1:nzm)

#ifdef CLUBB_CRM
do i = 1, nx
  do j = 1, ny
    do k = 1, nzm
      ! Apply local hole-filling to vapor by converting liquid to vapor. Moist
      ! static energy should be conserved, so updating temperature is not
      ! needed here. -dschanen 31 August 2011
      if ( qv(icrm,i,j,k) < zero_threshold ) then
        cloudliq(icrm,i,j,k) = cloudliq(icrm,i,j,k) + qv(icrm,i,j,k)
        qv(icrm,i,j,k) = zero_threshold
        if ( cloudliq(icrm,i,j,k) < zero_threshold ) then
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Total water at", "i =", i, "j =", j, "k =", k, "is negative.", &
              "Applying non-conservative hard clipping."
          end if
          cloudliq(icrm,i,j,k) = zero_threshold
        end if ! cloud_liq < 0
      end if ! qv < 0
    end do ! 1.. nzm
  end do ! 1.. ny
end do ! 1.. nx
#endif /* CLUBB_CRM */
! cloud liquid water
qcl(icrm,1:nx,1:ny,1:nzm) = cloudliq(icrm,1:nx,1:ny,1:nzm)

! rain water
if(doprecip) qpl(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqr)

! cloud ice
if(doicemicro) then
   qci(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqci)

   if(dograupel) then
      qpi(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqs) &
           + micro_field(icrm,1:nx,1:ny,1:nzm,iqg)
   else
      qpi(icrm,1:nx,1:ny,1:nzm) = micro_field(icrm,1:nx,1:ny,1:nzm,iqs)
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
  call micro_diagnose(ncrms,icrm)

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

  real(crm_rknd), dimension(nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  ! Total water mixing ratio
  micro_field(icrm,1:nx,1:ny,1:nzm,iqv) = new_qv(1:nx,1:ny,1:nzm) &
                                   + new_qc(1:nx,1:ny,1:nzm)

  ! Cloud water mixing ratio
  cloudliq(icrm,1:nx,1:ny,1:nzm) = new_qc(1:nx,1:ny,1:nzm)

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

!----------------------------------------------------------------------
!!! compute sedimentation
subroutine micro_precip_fall(ncrms)
  implicit none
  integer, intent(in) :: ncrms
  return ! do not need this routine -- sedimentation done in m2005micro.
end subroutine micro_precip_fall

!----------------------------------------------------------------------
! called when stepout() called

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
  real(crm_rknd), intent(inout), dimension(nzm) :: tabs ! absolute temperature, K
  real(crm_rknd), intent(inout), dimension(nzm) :: qt  ! on input: qt; on output: qv
  real(crm_rknd), intent(out), dimension(nzm) :: qc ! cloud liquid water, kg/kg
  real(crm_rknd), intent(in), dimension(nzm) :: pres ! pressure, Pa

  real(crm_rknd) tabs1, dtabs, thresh, esat1, qsat1, fff, dfff
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
      dtabs = + fac_cond*MAX(real(0.,crm_rknd),qt(k) - qsat1) &
           / ( 1. + lcond**2*qsat1/(cp*rv*tabs1**2) )
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(real(0.01,crm_rknd), 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        esat1 = polysvp(tabs1,0)
        qsat1 = 0.622*esat1/ (pres(k) - esat1) ! saturation mixing ratio

        fff = tabs(k) - tabs1 + fac_cond*MAX(real(0.,crm_rknd),qt(k) - qsat1)
        dfff = 1. + lcond**2*qsat1/(cp*rv*tabs1**2)
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      qc(k) = MAX( real(0.,crm_rknd), tabs1 - tabs(k) )/fac_cond ! cloud liquid mass mixing ratio
      qt(k) = qt(k) - qc(k) ! This now holds the water vapor mass mixing ratio.
      tabs(k) = tabs1 ! update temperature.

      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,nzm

end subroutine satadj_liquid

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water(ncrms,icrm)
  use vars, only : nstep,adz,dz,rho
  implicit none
  integer, intent(in) :: ncrms,icrm
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m,icrm).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(icrm,i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(icrm,k)*dz(icrm)*rho(icrm,k)
    end do
   end if
  end do

end function total_water
#ifdef CLUBB_CRM
!-------------------------------------------------------------------------------
ELEMENTAL REAL(crm_rknd) FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

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
REAL(crm_rknd), INTENT(IN):: var_high
REAL(crm_rknd), INTENT(IN):: var_low
REAL(crm_rknd), INTENT(IN):: height_high
REAL(crm_rknd), INTENT(IN):: height_low
REAL(crm_rknd), INTENT(IN):: height_int

! Output Variable
REAL(crm_rknd) :: LIN_INT

LIN_INT = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT
#endif /*CLUBB_CRM*/
!------------------------------------------------------------------------------

end module microphysics
