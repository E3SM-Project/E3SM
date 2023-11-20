module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
! Jul 2011  F Vitt -- made avaliable to be used in a prescribed modal aerosol mode (no prognostic MAM)
! Mar 2012  F Vitt -- made changes for to prevent abort when 7-mode aeroslol model is used
!                     some of the needed consituents do not exist in 7-mode so bin_fluxes will be false
! May 2014  F Vitt -- included contributions from MAM4 aerosols and added soa_a2 to the ocphiwet fluxes
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch,       only: cam_out_t     
use constituents,     only: pcnst, cnst_get_ind
use ppgrid,           only: pcols
use cam_abortutils,       only: endrun
use modal_aero_data,  only: nso4, nsoa, npoa, nbc, so4SpecName, soaSpecName, poaSpecName, bcSpecName

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data
integer :: idx_bc1(1:nbc)   = -1
integer :: idx_pom1(1:npoa) = -1
integer :: idx_soa1(1:nsoa) = -1
integer :: idx_soa2(1:nsoa) = -1
integer :: idx_dst1 = -1
integer :: idx_dst3 = -1
integer :: idx_ncl3 = -1
integer :: idx_so43(1:nso4) = -1
integer :: idx_bc4(1:nbc)  = -1
integer :: idx_pom4(1:npoa) = -1

!mgf++ MAM7
integer :: idx_bc3(1:nbc)  = -1
integer :: idx_pom3(1:npoa) = -1
integer :: idx_dst5 = -1
integer :: idx_dst7 = -1
!mgf--

logical :: bin_fluxes = .false.

logical :: initialized = .false.

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init()

! set aerosol indices for re-mapping surface deposition fluxes:
! *_a1 = accumulation mode
! *_a2 = aitken mode
! *_a3 = coarse mode

   ! can be initialized with user specified indices
   ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified

   integer :: i

   ! if already initialized abort the run
   if (initialized) then
     call endrun('modal_aero_deposition_init is already initialized')
   endif

   do i = 1, nbc
      call cnst_get_ind(trim(bcSpecName(i))//'_a1',  idx_bc1(i))
   enddo

   do i = 1, npoa
      call cnst_get_ind(trim(poaSpecName(i))//'_a1',  idx_pom1(i))
   enddo

   do i = 1, nsoa
      call cnst_get_ind(trim(soaSpecName(i))//'_a1',  idx_soa1(i))
      call cnst_get_ind(trim(soaSpecName(i))//'_a2',  idx_soa2(i))
   enddo

   do i = 1, nso4
      call cnst_get_ind(trim(so4SpecName(i))//'_a3',  idx_so43(i))
   enddo

   call cnst_get_ind('dst_a1', idx_dst1)
   call cnst_get_ind('dst_a3', idx_dst3)
   call cnst_get_ind('ncl_a3', idx_ncl3)

#ifdef MODAL_AER
#if( (defined MODAL_AERO_4MODE) || (defined MODAL_AERO_4MODE_MOM) )
   do i = 1, nbc
      call cnst_get_ind(trim(bcSpecName(i))//'_a4',  idx_bc4(i))
   enddo

   do i = 1, npoa
      call cnst_get_ind(trim(poaSpecName(i))//'_a4',  idx_pom4(i))
   enddo
#endif
#endif

#ifdef MODAL_AER
!mgf++ incorporate MAM7 fluxes
   
   ! unsure of the consequences of this being true by default, but it
   ! is needed for MAM7 prognostic fluxes to be active
   bin_fluxes = .true.  

#if( (defined MODAL_AERO_7MODE) || (defined MODAL_AERO_9MODE) )
      ! assign additional indices for MAM7 species:
      call cnst_get_ind('bc_a3',  idx_bc3)
      call cnst_get_ind('pom_a3', idx_pom3)
      call cnst_get_ind('dst_a5', idx_dst5)
      call cnst_get_ind('dst_a7', idx_dst7)
#endif
!mgf--

#else
!  for 7 mode bin_fluxes will be false
   bin_fluxes = idx_dst1>0 .and. idx_dst3>0 .and.idx_ncl3>0 .and. idx_so43(1)>0
#endif

#ifdef RAIN_EVAP_TO_COARSE_AERO
      ! assign additional indices for resuspended BC and POM to coarse mode:

   do i = 1, nbc
      call cnst_get_ind(trim(bcSpecName(i))//'_a3',  idx_bc3(i))
   enddo

   do i = 1, npoa
      call cnst_get_ind(trim(poaSpecName(i))//'_a3',  idx_pom3(i))
   enddo
#endif

   initialized = .true.

end subroutine modal_aero_deposition_init

!==============================================================================
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec
   integer :: ncol                      ! number of columns
   real(r8) :: bcphiwet_sum, bcphidry_sum, ocphiwet_sum, ocphidry_sum
   !----------------------------------------------------------------------------

   ncol = cam_out%ncol

#ifdef MODAL_AER
   do i = 1, ncol

!mgf++
#ifdef MODAL_AERO_3MODE
      ! MAM3
 
!      ! in SNICAR+MAM, bcphiwet represents BC mixed internally within hydrometeors
!      ! bcphidry represents BC mixed externally to hydrometeors
!      ! ocphiwet represents OC mixed internally within hydrometeors
!      ! ocphidry represents OC mixed externally to hydrometeors

      bcphiwet_sum = 0.0_r8
      do ispec = 1, nbc
         bcphiwet_sum = bcphiwet_sum + aerdepwetcw(i,idx_bc1(ispec))
      enddo
      cam_out%bcphiwet(i) = -(bcphiwet_sum)

      bcphidry_sum = 0.0_r8
      do ispec = 1, nbc
         bcphidry_sum = bcphidry_sum + aerdepwetis(i,idx_bc1(ispec))
      enddo
      cam_out%bcphidry(i) = -(bcphidry_sum)

      ocphiwet_sum = 0.0_r8
      do ispec = 1, npoa
         ocphiwet_sum = ocphiwet_sum + aerdepwetcw(i,idx_pom1(ispec))
      enddo
      do ispec = 1, nsoa
         ocphiwet_sum = ocphiwet_sum + aerdepwetcw(i,idx_soa1(ispec)) + aerdepwetcw(i,idx_soa2(ispec))
      enddo
      cam_out%ocphiwet(i)= -(ocphiwet_sum)

      ocphidry_sum = 0.0_r8
      do ispec = 1, npoa
         ocphidry_sum = ocphidry_sum + aerdepwetis(i,idx_pom1(ispec))
      enddo
      do ispec = 1, nsoa
         ocphidry_sum = ocphidry_sum + aerdepwetis(i,idx_soa1(ispec)) + aerdepwetis(i,idx_soa2(ispec))
      enddo
      cam_out%ocphidry(i)= -(ocphidry_sum)

#ifdef RAIN_EVAP_TO_COARSE_AERO
      do ispec = 1, nbc
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -aerdepwetcw(i,idx_bc3(ispec))
         cam_out%bcphidry(i) = cam_out%bcphidry(i) -aerdepwetis(i,idx_bc3(ispec))
      enddo

      do ispec = 1, npoa
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -aerdepwetcw(i,idx_pom3(ispec))
         cam_out%ocphidry(i) = cam_out%ocphidry(i) -aerdepwetis(i,idx_pom3(ispec))
      enddo
#endif

      
      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).  In MAM3,
      ! accumulation-mode dust is assumed to resemble fine dust
      cam_out%dstwet1(i) = -(0.625_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.015_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

      cam_out%dstwet2(i) = -(0.345_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.252_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

      cam_out%dstwet3(i) = -(0.029_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.444_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

      cam_out%dstwet4(i) = -(0.001_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.289_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

#endif

#if( (defined MODAL_AERO_4MODE) || (defined MODAL_AERO_4MODE_MOM) )
      ! MAM4
!      ! in SNICAR+MAM, bcphiwet represents BC mixed internally within
!      ! hydrometeors
!      ! bcphidry represents BC mixed externally to hydrometeors
!      ! ocphiwet represents OC mixed internally within hydrometeors
!      ! ocphidry represents OC mixed externally to hydrometeors

      bcphiwet_sum = 0.0_r8
      do ispec = 1, nbc
         bcphiwet_sum = bcphiwet_sum + aerdepwetcw(i,idx_bc1(ispec)) + aerdepwetcw(i,idx_bc4(ispec))
      enddo
      cam_out%bcphiwet(i) = -(bcphiwet_sum)

      bcphidry_sum = 0.0_r8
      do ispec = 1, nbc
         bcphidry_sum = bcphidry_sum + aerdepwetis(i,idx_bc1(ispec)) + aerdepwetis(i,idx_bc4(ispec))
      enddo
      cam_out%bcphidry(i) = -(bcphidry_sum)

      ocphiwet_sum = 0.0_r8
      do ispec = 1, npoa
         ocphiwet_sum = ocphiwet_sum + aerdepwetcw(i,idx_pom1(ispec)) + aerdepwetcw(i,idx_pom4(ispec))
      enddo
      do ispec = 1, nsoa
         ocphiwet_sum = ocphiwet_sum + aerdepwetcw(i,idx_soa1(ispec)) + aerdepwetcw(i,idx_soa2(ispec))
      enddo
      cam_out%ocphiwet(i)= -(ocphiwet_sum)

      ocphidry_sum = 0.0_r8
      do ispec = 1, npoa
         ocphidry_sum = ocphidry_sum + aerdepwetis(i,idx_pom1(ispec)) + aerdepwetis(i,idx_pom4(ispec))
      enddo
      do ispec = 1, nsoa
         ocphidry_sum = ocphidry_sum + aerdepwetis(i,idx_soa1(ispec)) + aerdepwetis(i,idx_soa2(ispec))
      enddo
      cam_out%ocphidry(i)= -(ocphidry_sum)

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC

      do ispec = 1, nbc
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -aerdepwetcw(i,idx_bc3(ispec))
         cam_out%bcphidry(i) = cam_out%bcphidry(i) -aerdepwetis(i,idx_bc3(ispec))
      enddo

      do ispec = 1, npoa
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -aerdepwetcw(i,idx_pom3(ispec))
         cam_out%ocphidry(i) = cam_out%ocphidry(i) -aerdepwetis(i,idx_pom3(ispec))
      enddo
#endif

      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).  In MAM3,
      ! accumulation-mode dust is assumed to resemble fine dust
      cam_out%dstwet1(i) = -(0.625_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.015_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

      cam_out%dstwet2(i) = -(0.345_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.252_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

      cam_out%dstwet3(i) = -(0.029_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.444_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

      cam_out%dstwet4(i) = -(0.001_r8*(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))+ &
                             0.289_r8*(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3)))

#endif

#if( (defined MODAL_AERO_7MODE) || (defined MODAL_AERO_9MODE) )
      ! MAM7
      
      ! in SNICAR+MAM, bcphiwet represents BC mixed internally within hydrometeors
      cam_out%bcphiwet(i) = -(aerdepwetcw(i,idx_bc1)+aerdepwetcw(i,idx_bc3))

      ! bcphidry represents BC mixed externally to hydrometeors
      cam_out%bcphidry(i) = -(aerdepwetis(i,idx_bc1)+aerdepwetis(i,idx_bc3))

      ! ocphiwet represents OC mixed internally within hydrometeors
      cam_out%ocphiwet(i) = -(aerdepwetcw(i,idx_pom1)+aerdepwetcw(i,idx_pom3)+ &
                              aerdepwetcw(i,idx_soa1)+aerdepwetcw(i,idx_soa2))

      ! ocphidry represents OC mixed externally to hydrometeors
      cam_out%ocphidry(i) = -(aerdepwetis(i,idx_pom1)+aerdepwetis(i,idx_pom3)+ &
                              aerdepwetis(i,idx_soa1)+aerdepwetis(i,idx_soa2))


      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).
      cam_out%dstwet1(i) = -(0.625_r8*(aerdepwetis(i,idx_dst5)+aerdepwetcw(i,idx_dst5))+ &
                             0.015_r8*(aerdepwetis(i,idx_dst7)+aerdepwetcw(i,idx_dst7)))

      cam_out%dstwet2(i) = -(0.345_r8*(aerdepwetis(i,idx_dst5)+aerdepwetcw(i,idx_dst5))+ &
                             0.252_r8*(aerdepwetis(i,idx_dst7)+aerdepwetcw(i,idx_dst7)))

      cam_out%dstwet3(i) = -(0.029_r8*(aerdepwetis(i,idx_dst5)+aerdepwetcw(i,idx_dst5))+ &
                             0.444_r8*(aerdepwetis(i,idx_dst7)+aerdepwetcw(i,idx_dst7)))

      cam_out%dstwet4(i) = -(0.001_r8*(aerdepwetis(i,idx_dst5)+aerdepwetcw(i,idx_dst5))+ &
                             0.289_r8*(aerdepwetis(i,idx_dst7)+aerdepwetcw(i,idx_dst7)))

#endif
      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet2(i)  .lt. 0._r8) cam_out%dstwet2(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
      if (cam_out%dstwet4(i)  .lt. 0._r8) cam_out%dstwet4(i)  = 0._r8
   enddo

!mgf--

#else

   cam_out%bcphiwet(:) = 0._r8
   cam_out%ocphiwet(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol
      ! black carbon fluxes
      bcphiwet_sum = 0.0_r8
      do ispec = 1, nbc
         if (idx_bc1(1)>0) then
            bcphiwet_sum = bcphiwet_sum +aerdepwetis(i,idx_bc1(ispec))+aerdepwetcw(i,idx_bc1(ispec))
         endif
         if (idx_bc4(1)>0) then
            bcphiwet_sum = bcphiwet_sum +aerdepwetis(i,idx_bc4(ispec))+aerdepwetcw(i,idx_bc4(ispec))
         endif
      enddo
      cam_out%bcphiwet(i) = -(bcphiwet_sum)

      ! organic carbon fluxes
      ocphiwet_sum = 0.0_r8
      do ispec = 1, npoa
         if (idx_pom1(1)>0) then
            ocphiwet_sum = ocphiwet_sum +aerdepwetis(i,idx_pom1(ispec))+aerdepwetcw(i,idx_pom1(ispec))
         endif
         if (idx_pom4(1)>0) then
            ocphiwet_sum = ocphiwet_sum +aerdepwetis(i,idx_pom4(ispec))+aerdepwetcw(i,idx_pom4(ispec))
         endif
      enddo
      do ispec = 1, nsoa
         if (idx_soa1(1)>0) then
            ocphiwet_sum = ocphiwet_sum +aerdepwetis(i,idx_soa1(ispec))+aerdepwetcw(i,idx_soa1(ispec))
         endif
         if (idx_soa2(1)>0) then
            ocphiwet_sum = ocphiwet_sum +aerdepwetis(i,idx_soa2(ispec))+aerdepwetcw(i,idx_soa2(ispec))
         endif
      enddo
      cam_out%ocphiwet(i)= -(ocphiwet_sum)

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC

      do ispec = 1, nbc
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -(aerdepwetis(i,idx_bc3(ispec))+aerdepwetcw(i,idx_bc3(ispec)))
      enddo
      do ispec = 1, npoa
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_pom3(ispec))+aerdepwetcw(i,idx_pom3(ispec)))
      enddo
#endif

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstwet1(i) = -(aerdepwetis(i,idx_dst1)+aerdepwetcw(i,idx_dst1))
      
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstwet2(i) = 0._r8
      cam_out%dstwet3(i) = -(aerdepwetis(i,idx_dst3)+aerdepwetcw(i,idx_dst3))
      cam_out%dstwet4(i) = 0._r8

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
   enddo
#endif

end subroutine set_srf_wetdep

!==============================================================================

subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)

! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec
   integer :: ncol                      ! number of columns
   real(r8):: bcphidry_sum, bcphodry_sum, ocphidry_sum, ocphodry_sum
   !----------------------------------------------------------------------------

   ncol = cam_out%ncol

#ifdef MODAL_AER

   do i = 1, ncol
      
!mgf++
#ifdef MODAL_AERO_3MODE
      ! MAM3
      ! in SNICAR+MAM, bcphodry represents BC mixed external to hydrometeors
      bcphodry_sum = 0.0_r8
      do ispec = 1, nbc
         bcphodry_sum = bcphodry_sum + aerdepdryis(i,idx_bc1(ispec))+aerdepdrycw(i,idx_bc1(ispec))
      enddo
      cam_out%bcphodry(i) = bcphodry_sum
      
      ! ocphodry represents OC mixed external to hydrometeors
      ocphodry_sum = 0.0_r8
      do ispec = 1, npoa
         ocphodry_sum = ocphodry_sum+aerdepdryis(i,idx_pom1(ispec))+aerdepdrycw(i,idx_pom1(ispec))
      enddo
      do ispec = 1, nsoa
         ocphodry_sum = ocphodry_sum+aerdepdryis(i,idx_soa1(ispec))+aerdepdrycw(i,idx_soa1(ispec))+ &
                                     aerdepdryis(i,idx_soa2(ispec))+aerdepdrycw(i,idx_soa2(ispec))
      enddo
      cam_out%ocphodry(i) = ocphodry_sum

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC
      do ispec = 1, nbc
         cam_out%bcphodry(i) = cam_out%bcphodry(i) +(aerdepdryis(i,idx_bc3(ispec))+aerdepdrycw(i,idx_bc3(ispec)))
      enddo
      do ispec = 1, npoa
         cam_out%ocphodry(i) = cam_out%ocphodry(i) +(aerdepdryis(i,idx_pom3(ispec))+aerdepdrycw(i,idx_pom3(ispec)))
      enddo
#endif


      ! NOTE: drycw fluxes shown above ideally would be included as
      ! within-hydrometeor species, but this would require passing
      ! additional species through the coupler.  drycw fluxes are
      ! extremely small in the global-mean, so this will make little
      ! difference.


      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).  In MAM3,
      ! accumulation-mode dust is assumed to resemble fine dust
      cam_out%dstdry1(i) = (0.625_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.015_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

      cam_out%dstdry2(i) = (0.345_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.252_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

      cam_out%dstdry3(i) = (0.029_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.444_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

      cam_out%dstdry4(i) = (0.001_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.289_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

#endif

#if( (defined MODAL_AERO_4MODE) || (defined MODAL_AERO_4MODE_MOM) )
      ! MAM4

      ! in SNICAR+MAM, bcphodry represents BC mixed external to hydrometeors
      bcphodry_sum = 0.0_r8
      do ispec = 1, nbc
         bcphodry_sum = bcphodry_sum +aerdepdryis(i,idx_bc1(ispec))+aerdepdrycw(i,idx_bc1(ispec))+ &
                                      aerdepdryis(i,idx_bc4(ispec))+aerdepdrycw(i,idx_bc4(ispec))
      enddo
      cam_out%bcphodry(i) = bcphodry_sum

      ! ocphodry represents OC mixed external to hydrometeors
      ocphodry_sum = 0.0_r8
      do ispec = 1, npoa
         ocphodry_sum = ocphodry_sum+aerdepdryis(i,idx_pom1(ispec))+aerdepdrycw(i,idx_pom1(ispec))+ &
                                     aerdepdryis(i,idx_pom4(ispec))+aerdepdrycw(i,idx_pom4(ispec))
      enddo
      do ispec = 1, nsoa
         ocphodry_sum = ocphodry_sum+aerdepdryis(i,idx_soa1(ispec))+aerdepdrycw(i,idx_soa1(ispec))+ &
                                     aerdepdryis(i,idx_soa2(ispec))+aerdepdrycw(i,idx_soa2(ispec))
      enddo
      cam_out%ocphodry(i) = ocphodry_sum

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC

      do ispec = 1, nbc
         cam_out%bcphodry(i) = cam_out%bcphodry(i) +(aerdepdryis(i,idx_bc3(ispec))+aerdepdrycw(i,idx_bc3(ispec)))
      enddo
      do ispec = 1, npoa
         cam_out%ocphodry(i) = cam_out%ocphodry(i) +(aerdepdryis(i,idx_pom3(ispec))+aerdepdrycw(i,idx_pom3(ispec)))
      enddo
#endif

      ! NOTE: drycw fluxes shown above ideally would be included as
      ! within-hydrometeor species, but this would require passing
      ! additional species through the coupler.  drycw fluxes are
      ! extremely small in the global-mean, so this will make little
      ! difference.


      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).  In MAM3,
      ! accumulation-mode dust is assumed to resemble fine dust
      cam_out%dstdry1(i) = (0.625_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.015_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

      cam_out%dstdry2(i) = (0.345_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.252_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

      cam_out%dstdry3(i) = (0.029_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.444_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

      cam_out%dstdry4(i) = (0.001_r8*(aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1))+ &
                            0.289_r8*(aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)))

#endif


#if( (defined MODAL_AERO_7MODE) || (defined MODAL_AERO_9MODE) )
      ! MAM7

      ! in SNICAR+MAM, bcphodry represents BC mixed external to hydrometeors
      cam_out%bcphodry(i) = aerdepdryis(i,idx_bc1)+aerdepdryis(i,idx_bc3)+ &
                            aerdepdrycw(i,idx_bc1)+aerdepdrycw(i,idx_bc3)
      
      ! ocphodry represents OC mixed external to hydrometeors
      cam_out%ocphodry(i) = aerdepdryis(i,idx_pom1)+aerdepdryis(i,idx_pom3)+aerdepdryis(i,idx_soa1)+aerdepdryis(i,idx_soa2)+ &
                            aerdepdrycw(i,idx_pom1)+aerdepdrycw(i,idx_pom3)+aerdepdrycw(i,idx_soa1)+aerdepdrycw(i,idx_soa2)

      ! NOTE: drycw fluxes ideally would be included as
      ! within-hydrometeor species, but this would require passing
      ! additional species through the coupler.  drycw fluxes are
      ! extremely small in the global-mean, so this will make little
      ! difference.


      ! Four dust bins in SNICAR represent dust with dry diameters of
      ! 0.1-1.0um, 1.0-2.5um, 2.5-5.0um, 5.0-10um, respectively.  Dust
      ! mass is partitioned into these bins based on global-mean size
      ! distributions of MAM7 fine dust and coarse dust shown in Table
      ! 1 of Liu et al (2012, doi:10.5194/gmd-5-709-2012).
      cam_out%dstdry1(i) = (0.625_r8*(aerdepdryis(i,idx_dst5)+aerdepdrycw(i,idx_dst5))+ &
                            0.015_r8*(aerdepdryis(i,idx_dst7)+aerdepdrycw(i,idx_dst7)))

      cam_out%dstdry2(i) = (0.345_r8*(aerdepdryis(i,idx_dst5)+aerdepdrycw(i,idx_dst5))+ &
                            0.252_r8*(aerdepdryis(i,idx_dst7)+aerdepdrycw(i,idx_dst7)))

      cam_out%dstdry3(i) = (0.029_r8*(aerdepdryis(i,idx_dst5)+aerdepdrycw(i,idx_dst5))+ &
                            0.444_r8*(aerdepdryis(i,idx_dst7)+aerdepdrycw(i,idx_dst7)))

      cam_out%dstdry4(i) = (0.001_r8*(aerdepdryis(i,idx_dst5)+aerdepdrycw(i,idx_dst5))+ &
                            0.289_r8*(aerdepdryis(i,idx_dst7)+aerdepdrycw(i,idx_dst7)))

#endif
      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%bcphodry(i) .lt. 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%ocphodry(i) .lt. 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  .lt. 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry2(i)  .lt. 0._r8) cam_out%dstdry2(i)  = 0._r8
      if (cam_out%dstdry3(i)  .lt. 0._r8) cam_out%dstdry3(i)  = 0._r8
      if (cam_out%dstdry4(i)  .lt. 0._r8) cam_out%dstdry4(i)  = 0._r8
   enddo

#else

   cam_out%bcphidry(:) = 0._r8
   cam_out%bcphodry(:) = 0._r8
   cam_out%ocphidry(:) = 0._r8
   cam_out%ocphodry(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol
      ! black carbon fluxes
      bcphidry_sum = 0.0_r8
      if (idx_bc1(1)>0) then
         do ispec = 1, nbc
            bcphidry_sum = bcphidry_sum + aerdepdryis(i,idx_bc1(ispec))+aerdepdrycw(i,idx_bc1(ispec))
         enddo
         cam_out%bcphidry(i) = bcphidry_sum
      endif

      bcphodry_sum = 0.0_r8
      if (idx_bc4(1)>0) then
         do ispec = 1, nbc
            bcphodry_sum = bcphodry_sum + aerdepdryis(i,idx_bc4(ispec))+aerdepdrycw(i,idx_bc4(ispec))
         enddo
         cam_out%bcphodry(i) = bcphodry_sum
      endif

      ! organic carbon fluxes
      ocphidry_sum = 0.0_r8
      if (idx_pom1(1)>0) then
         do ispec = 1, npoa
            ocphidry_sum = ocphidry_sum + aerdepdryis(i,idx_pom1(ispec))+aerdepdrycw(i,idx_pom1(ispec))
         enddo
      endif
      if (idx_soa1(1)>0) then
         do ispec = 1, nsoa
            ocphidry_sum = ocphidry_sum + aerdepdryis(i,idx_soa1(ispec)) + aerdepdrycw(i,idx_soa1(ispec))
         enddo
      endif
      cam_out%ocphidry(i) = ocphidry_sum

      ocphodry_sum = 0.0_r8
      if (idx_pom4(1)>0) then
         do ispec = 1, npoa
            ocphodry_sum = ocphodry_sum + aerdepdryis(i,idx_pom4(ispec))+aerdepdrycw(i,idx_pom4(ispec))
         enddo
      endif
      if (idx_soa2(1)>0) then
         do ispec = 1, nsoa
            ocphodry_sum = ocphodry_sum + aerdepdryis(i,idx_soa2(ispec))+aerdepdrycw(i,idx_soa2(ispec))
         enddo
      endif
      cam_out%ocphodry(i) = ocphodry_sum
#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC to xxphidry
      do ispec = 1, nbc
         cam_out%bcphidry(i) = cam_out%bcphidry(i) +(aerdepdryis(i,idx_bc3(ispec))+aerdepdrycw(i,idx_bc3(ispec)))
      enddo
      do ispec = 1, npoa
         cam_out%ocphidry(i) = cam_out%ocphidry(i) +(aerdepdryis(i,idx_pom3(ispec))+aerdepdrycw(i,idx_pom3(ispec)))
      enddo
#endif

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstdry1(i) = aerdepdryis(i,idx_dst1)+aerdepdrycw(i,idx_dst1)
      
      ! Two options for partitioning deposition into bins 2-4:
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstdry2(i) = 0._r8
      cam_out%dstdry3(i) = aerdepdryis(i,idx_dst3)+aerdepdrycw(i,idx_dst3)
      cam_out%dstdry4(i) = 0._r8

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) .lt. 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%bcphodry(i) .lt. 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphidry(i) .lt. 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%ocphodry(i) .lt. 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  .lt. 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry3(i)  .lt. 0._r8) cam_out%dstdry3(i)  = 0._r8
   enddo
#endif

end subroutine set_srf_drydep


!==============================================================================

end module modal_aero_deposition
