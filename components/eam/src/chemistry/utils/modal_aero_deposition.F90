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

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data
integer :: idx_bc1  = -1
integer :: idx_pom1 = -1
integer :: idx_soa1 = -1
integer :: idx_soa2 = -1
integer :: idx_dst1 = -1
integer :: idx_dst3 = -1
integer :: idx_ncl3 = -1
integer :: idx_so43 = -1
integer :: idx_bc4  = -1
integer :: idx_pom4 = -1

!mgf++ MAM7
integer :: idx_bc3  = -1
integer :: idx_pom3 = -1
integer :: idx_dst5 = -1
integer :: idx_dst7 = -1
!mgf--

logical :: bin_fluxes = .false.

logical :: initialized = .false.

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init(bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx, &
                            dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx,bc4_ndx,pom4_ndx)

! set aerosol indices for re-mapping surface deposition fluxes:
! *_a1 = accumulation mode
! *_a2 = aitken mode
! *_a3 = coarse mode

   ! can be initialized with user specified indices
   ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified

   integer, optional, intent(in) :: bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx,dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx
   integer, optional, intent(in) :: bc4_ndx,pom4_ndx

   ! if already initialized abort the run
   if (initialized) then
     call endrun('modal_aero_deposition_init is already initialized')
   endif

   if (present(bc1_ndx)) then
      idx_bc1  = bc1_ndx
   else
      call cnst_get_ind('bc_a1',  idx_bc1)
   endif
   if (present(pom1_ndx)) then
      idx_pom1 = pom1_ndx
   else
      call cnst_get_ind('pom_a1', idx_pom1)
   endif
   if (present(soa1_ndx)) then
      idx_soa1 = soa1_ndx
   else
      call cnst_get_ind('soa_a1', idx_soa1)
   endif
   if (present(soa2_ndx)) then
      idx_soa2 = soa2_ndx
   else
      call cnst_get_ind('soa_a2', idx_soa2)
   endif
   if (present(dst1_ndx)) then
      idx_dst1 = dst1_ndx
   else
      call cnst_get_ind('dst_a1', idx_dst1,abort=.false.)
   endif
   if (present(dst3_ndx)) then
      idx_dst3 = dst3_ndx
   else
      call cnst_get_ind('dst_a3', idx_dst3,abort=.false.)
   endif
   if (present(ncl3_ndx)) then
      idx_ncl3 = ncl3_ndx
   else
      call cnst_get_ind('ncl_a3', idx_ncl3,abort=.false.)
   endif
   if (present(so43_ndx)) then
      idx_so43 = so43_ndx
   else
      call cnst_get_ind('so4_a3', idx_so43,abort=.false.)
   endif
   if (present(bc4_ndx)) then
      idx_bc4 = bc4_ndx
   else
      call cnst_get_ind('bc_a4', idx_bc4,abort=.false.)
   endif
   if (present(pom4_ndx)) then
      idx_pom4 = pom4_ndx
   else
      call cnst_get_ind('pom_a4', idx_pom4,abort=.false.)   
   endif

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
   bin_fluxes = idx_dst1>0 .and. idx_dst3>0 .and.idx_ncl3>0 .and. idx_so43>0
#endif

#ifdef RAIN_EVAP_TO_COARSE_AERO
      ! assign additional indices for resuspended BC and POM to coarse mode:
      call cnst_get_ind('bc_a3',  idx_bc3)
      call cnst_get_ind('pom_a3', idx_pom3)
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
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   if (.not.bin_fluxes) return

   ncol = cam_out%ncol

#ifdef MODAL_AER
   do i = 1, ncol

!mgf++
#ifdef MODAL_AERO_3MODE
      ! MAM3
    
      ! in SNICAR+MAM, bcphiwet represents BC mixed internally within hydrometeors
      cam_out%bcphiwet(i) = -(aerdepwetcw(i,idx_bc1))

      ! bcphidry represents BC mixed externally to hydrometeors
      cam_out%bcphidry(i) = -(aerdepwetis(i,idx_bc1))

      ! ocphiwet represents OC mixed internally within hydrometeors
      cam_out%ocphiwet(i) = -(aerdepwetcw(i,idx_pom1)+aerdepwetcw(i,idx_soa1)+aerdepwetcw(i,idx_soa2))

      ! ocphidry represents OC mixed externally to hydrometeors
      cam_out%ocphidry(i) = -(aerdepwetis(i,idx_pom1)+aerdepwetis(i,idx_soa1)+aerdepwetis(i,idx_soa2))

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -aerdepwetcw(i,idx_bc3)
         cam_out%bcphidry(i) = cam_out%bcphidry(i) -aerdepwetis(i,idx_bc3)
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -aerdepwetcw(i,idx_pom3)
         cam_out%ocphidry(i) = cam_out%ocphidry(i) -aerdepwetis(i,idx_pom3)
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

      ! in SNICAR+MAM, bcphiwet represents BC mixed internally within
      ! hydrometeors
      cam_out%bcphiwet(i) = -(aerdepwetcw(i,idx_bc1)+aerdepwetcw(i,idx_bc4))

      ! bcphidry represents BC mixed externally to hydrometeors
      cam_out%bcphidry(i) = -(aerdepwetis(i,idx_bc1)+aerdepwetis(i,idx_bc4))

      ! ocphiwet represents OC mixed internally within hydrometeors
      cam_out%ocphiwet(i) = -(aerdepwetcw(i,idx_pom1)+aerdepwetcw(i,idx_pom4)+ &
                              aerdepwetcw(i,idx_soa1)+aerdepwetcw(i,idx_soa2))

      ! ocphidry represents OC mixed externally to hydrometeors
      cam_out%ocphidry(i) = -(aerdepwetis(i,idx_pom1)+aerdepwetis(i,idx_pom4)+ &
                              aerdepwetis(i,idx_soa1)+aerdepwetis(i,idx_soa2))

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -aerdepwetcw(i,idx_bc3)
         cam_out%bcphidry(i) = cam_out%bcphidry(i) -aerdepwetis(i,idx_bc3)
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -aerdepwetcw(i,idx_pom3)
         cam_out%ocphidry(i) = cam_out%ocphidry(i) -aerdepwetis(i,idx_pom3)
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
      if (idx_bc1>0) &
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -(aerdepwetis(i,idx_bc1)+aerdepwetcw(i,idx_bc1))
      if (idx_bc4>0) &
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -(aerdepwetis(i,idx_bc4)+aerdepwetcw(i,idx_bc4))

      ! organic carbon fluxes
      if (idx_soa1>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_soa1)+aerdepwetcw(i,idx_soa1))
      if (idx_soa2>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_soa2)+aerdepwetcw(i,idx_soa2))
      if (idx_pom1>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_pom1)+aerdepwetcw(i,idx_pom1))
      if (idx_pom4>0) &
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_pom4)+aerdepwetcw(i,idx_pom4))
#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) -(aerdepwetis(i,idx_bc3)+aerdepwetcw(i,idx_bc3))
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) -(aerdepwetis(i,idx_pom3)+aerdepwetcw(i,idx_pom3))
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
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------
   if (.not.bin_fluxes) return

   ncol = cam_out%ncol

#ifdef MODAL_AER

   do i = 1, ncol
      
!mgf++
#ifdef MODAL_AERO_3MODE
      ! MAM3
     
      ! in SNICAR+MAM, bcphodry represents BC mixed external to hydrometeors
      cam_out%bcphodry(i) = aerdepdryis(i,idx_bc1)+aerdepdrycw(i,idx_bc1)
      
      ! ocphodry represents OC mixed external to hydrometeors
      cam_out%ocphodry(i) = aerdepdryis(i,idx_pom1)+aerdepdryis(i,idx_soa1)+aerdepdryis(i,idx_soa2)+ &
                            aerdepdrycw(i,idx_pom1)+aerdepdrycw(i,idx_soa1)+aerdepdrycw(i,idx_soa2)
#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC
         cam_out%bcphodry(i) = cam_out%bcphodry(i) + (aerdepdryis(i,idx_bc3)+aerdepdrycw(i,idx_bc3))
         cam_out%ocphodry(i) = cam_out%ocphodry(i) + (aerdepdryis(i,idx_pom3)+aerdepdrycw(i,idx_pom3))
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
      cam_out%bcphodry(i) = aerdepdryis(i,idx_bc1)+aerdepdryis(i,idx_bc4)+ &
                            aerdepdrycw(i,idx_bc1)+aerdepdrycw(i,idx_bc4)

      ! ocphodry represents OC mixed external to hydrometeors
      cam_out%ocphodry(i) = aerdepdryis(i,idx_pom1)+aerdepdryis(i,idx_pom4)+aerdepdryis(i,idx_soa1)+aerdepdryis(i,idx_soa2)+ &
                            aerdepdrycw(i,idx_pom1)+aerdepdrycw(i,idx_pom4)+aerdepdrycw(i,idx_soa1)+aerdepdrycw(i,idx_soa2)

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC
         cam_out%bcphodry(i) = cam_out%bcphodry(i) + (aerdepdryis(i,idx_bc3)+aerdepdrycw(i,idx_bc3))
         cam_out%ocphodry(i) = cam_out%ocphodry(i) + (aerdepdryis(i,idx_pom3)+aerdepdrycw(i,idx_pom3))
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
      if (idx_bc1>0) &
           cam_out%bcphidry(i) = cam_out%bcphidry(i) + aerdepdryis(i,idx_bc1)+aerdepdrycw(i,idx_bc1)
      if (idx_bc4>0) &
           cam_out%bcphodry(i) = cam_out%bcphodry(i) + aerdepdryis(i,idx_bc4)+aerdepdrycw(i,idx_bc4)

      ! organic carbon fluxes
      if (idx_pom1>0) &
           cam_out%ocphidry(i) = cam_out%ocphidry(i) + aerdepdryis(i,idx_pom1)+aerdepdrycw(i,idx_pom1)
      if( idx_pom4>0) &
           cam_out%ocphodry(i) = cam_out%ocphodry(i) + aerdepdryis(i,idx_pom4)+aerdepdrycw(i,idx_pom4)
      if (idx_soa1>0) &
           cam_out%ocphidry(i) = cam_out%ocphidry(i) + aerdepdryis(i,idx_soa1)+aerdepdrycw(i,idx_soa1)
      if (idx_soa2>0) &
           cam_out%ocphodry(i) = cam_out%ocphodry(i) + aerdepdryis(i,idx_soa2)+aerdepdrycw(i,idx_soa2)

#ifdef RAIN_EVAP_TO_COARSE_AERO
       ! add resuspended coarse-mode BC and OC to xxphidry
         cam_out%bcphidry(i) = cam_out%bcphidry(i) + (aerdepdryis(i,idx_bc3)+aerdepdrycw(i,idx_bc3))
         cam_out%ocphidry(i) = cam_out%ocphidry(i) + (aerdepdryis(i,idx_pom3)+aerdepdrycw(i,idx_pom3))
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
