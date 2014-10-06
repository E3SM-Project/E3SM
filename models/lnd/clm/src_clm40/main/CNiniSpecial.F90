!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNiniSpecial
!
! !INTERFACE:
subroutine CNiniSpecial ()
!
! !DESCRIPTION:
! One-time initialization of CN variables for special landunits
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pftvarcon   , only: noveg
   use decompMod   , only: get_proc_bounds
   use clm_varcon  , only: spval
   use clm_varctl  , only: iulog, use_c13, use_cn
   use clmtype
   use CNSetValueMod
!
! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:
! subroutine iniTimeConst in file iniTimeConst.F90
!
! !REVISION HISTORY:
! 11/13/03: Created by Peter Thornton
!
!
! local pointers to implicit in arguments
!
  integer , pointer :: clandunit(:)    ! landunit index of corresponding column
  integer , pointer :: plandunit(:)    ! landunit index of corresponding pft
  logical , pointer :: ifspecial(:)    ! BOOL: true=>landunit is wetland,ice,lake, or urban
!
! local pointers to implicit out arguments
!
! !LOCAL VARIABLES:
!EOP
   integer :: fc,fp,l,c,p  ! indices
   integer :: begp, endp   ! per-clump/proc beginning and ending pft indices
   integer :: begc, endc   ! per-clump/proc beginning and ending column indices
   integer :: begl, endl   ! per-clump/proc beginning and ending landunit indices
   integer :: begg, endg   ! per-clump/proc gridcell ending gridcell indices
   integer :: num_specialc ! number of good values in specialc filter
   integer :: num_specialp ! number of good values in specialp filter
   integer, allocatable :: specialc(:) ! special landunit filter - columns
   integer, allocatable :: specialp(:) ! special landunit filter - pfts
!-----------------------------------------------------------------------
   ! assign local pointers at the landunit level
   ifspecial => lun%ifspecial

   ! assign local pointers at the column level
   clandunit => col%landunit

   ! assign local pointers at the pft level
   plandunit => pft%landunit

   ! Determine subgrid bounds on this processor
   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

   ! allocate special landunit filters
   allocate(specialc(endc-begc+1))
   allocate(specialp(endp-begp+1))

   ! fill special landunit filters
   num_specialc = 0
   do c = begc, endc
      l = clandunit(c)
      if (ifspecial(l)) then
         num_specialc = num_specialc + 1
         specialc(num_specialc) = c
      end if
   end do

   num_specialp = 0
   do p = begp, endp
      l = plandunit(p)
      if (ifspecial(l)) then
         num_specialp = num_specialp + 1
         specialp(num_specialp) = p
      end if
   end do

   ! initialize column-level fields
   call CNSetCps(num_specialc, specialc, spval, cps)
   call CNSetCcs(num_specialc, specialc, 0._r8, ccs)
   call CNSetCns(num_specialc, specialc, 0._r8, cns)
   call CNSetCcf(num_specialc, specialc, 0._r8, ccf)
   call CNSetCnf(num_specialc, specialc, 0._r8, cnf)
   if (use_c13) then
      ! 4/14/05: PET
      ! adding isotope code
      call CNSetCcs(num_specialc, specialc, 0._r8, cc13s)
      call CNSetCcf(num_specialc, specialc, 0._r8, cc13f)
   endif

   ! initialize column-average pft fields
   call CNSetPps(num_specialc, specialc, spval, pps_a)
   call CNSetPcs(num_specialc, specialc, 0._r8, pcs_a)
   call CNSetPns(num_specialc, specialc, 0._r8, pns_a)
   call CNSetPcf(num_specialc, specialc, 0._r8, pcf_a)
   call CNSetPnf(num_specialc, specialc, 0._r8, pnf_a)

   ! initialize pft-level fields
   call CNSetPepv(num_specialp, specialp, spval, pepv)
   call CNSetPps(num_specialp, specialp, spval, pps)
   call CNSetPcs(num_specialp, specialp, 0._r8, pcs)
   call CNSetPns(num_specialp, specialp, 0._r8, pns)
   call CNSetPcf(num_specialp, specialp, 0._r8, pcf)
   call CNSetPnf(num_specialp, specialp, 0._r8, pnf)
   if (use_c13) then
      ! 4/14/05: PET
      ! adding isotope code
      call CNSetPcs(num_specialp, specialp, 0._r8, pc13s)
      call CNSetPcf(num_specialp, specialp, 0._r8, pc13f)
   endif

   ! now loop through special filters and explicitly set the variables that
   ! have to be in place for SurfaceAlbedo and biogeophysics
   ! also set pcf%psnsun and pcf%psnsha to 0 (not included in CNSetPcf())

   do fp = 1,num_specialp
      p = specialp(fp)
      pps%tlai(p) = 0._r8
      pps%tsai(p) = 0._r8
      pps%elai(p) = 0._r8
      pps%esai(p) = 0._r8
      pps%htop(p) = 0._r8
      pps%hbot(p) = 0._r8
      pps%fwet(p) = 0._r8
      pps%fdry(p) = 0._r8
      pps%frac_veg_nosno_alb(p) = 0._r8
      pps%frac_veg_nosno(p) = 0._r8
      pcf%psnsun(p) = 0._r8
      pcf%psnsha(p) = 0._r8
      if (use_c13) then
         ! 4/14/05: PET
         ! Adding isotope code
         pc13f%psnsun(p) = 0._r8
         pc13f%psnsha(p) = 0._r8
      endif
      
   end do

   do fc = 1,num_specialc
      c = specialc(fc)
      pcf_a%psnsun(c) = 0._r8
      pcf_a%psnsha(c) = 0._r8
      if (use_c13) then
         ! 8/17/05: PET
         ! Adding isotope code
         pcf_a%psnsun(c) = 0._r8
         pcf_a%psnsha(c) = 0._r8
      endif
      
	  ! adding dynpft code
	  ccs%seedc(c) = 0._r8
	  ccs%prod10c(c) = 0._r8	  
	  ccs%prod100c(c) = 0._r8	  
	  ccs%totprodc(c) = 0._r8	  
          if (use_c13) then
             cc13s%seedc(c) = 0._r8
             cc13s%prod10c(c) = 0._r8	  
             cc13s%prod100c(c) = 0._r8	  
             cc13s%totprodc(c) = 0._r8	  
          endif
	  cns%seedn(c) = 0._r8
	  cns%prod10n(c) = 0._r8	  
	  cns%prod100n(c) = 0._r8	  
	  cns%totprodn(c) = 0._r8	  
	  ccf%dwt_seedc_to_leaf(c) = 0._r8
	  ccf%dwt_seedc_to_deadstem(c) = 0._r8
	  ccf%dwt_conv_cflux(c) = 0._r8
	  ccf%dwt_prod10c_gain(c) = 0._r8
	  ccf%prod10c_loss(c) = 0._r8
	  ccf%dwt_prod100c_gain(c) = 0._r8
	  ccf%prod100c_loss(c) = 0._r8
	  ccf%dwt_frootc_to_litr1c(c) = 0._r8
	  ccf%dwt_frootc_to_litr2c(c) = 0._r8
	  ccf%dwt_frootc_to_litr3c(c) = 0._r8
	  ccf%dwt_livecrootc_to_cwdc(c) = 0._r8
	  ccf%dwt_deadcrootc_to_cwdc(c) = 0._r8
	  ccf%dwt_closs(c) = 0._r8
	  ccf%landuseflux(c) = 0._r8
	  ccf%landuptake(c) = 0._r8
          if (use_c13) then
             cc13f%dwt_seedc_to_leaf(c) = 0._r8
             cc13f%dwt_seedc_to_deadstem(c) = 0._r8
             cc13f%dwt_conv_cflux(c) = 0._r8
             cc13f%dwt_prod10c_gain(c) = 0._r8
             cc13f%prod10c_loss(c) = 0._r8
             cc13f%dwt_prod100c_gain(c) = 0._r8
             cc13f%prod100c_loss(c) = 0._r8
             cc13f%dwt_frootc_to_litr1c(c) = 0._r8
             cc13f%dwt_frootc_to_litr2c(c) = 0._r8
             cc13f%dwt_frootc_to_litr3c(c) = 0._r8
             cc13f%dwt_livecrootc_to_cwdc(c) = 0._r8
             cc13f%dwt_deadcrootc_to_cwdc(c) = 0._r8
             cc13f%dwt_closs(c) = 0._r8
          endif
	  cnf%dwt_seedn_to_leaf(c) = 0._r8
	  cnf%dwt_seedn_to_deadstem(c) = 0._r8
	  cnf%dwt_conv_nflux(c) = 0._r8
	  cnf%dwt_prod10n_gain(c) = 0._r8
	  cnf%prod10n_loss(c) = 0._r8
	  cnf%dwt_prod100n_gain(c) = 0._r8
	  cnf%prod100n_loss(c) = 0._r8
	  cnf%dwt_frootn_to_litr1n(c) = 0._r8
	  cnf%dwt_frootn_to_litr2n(c) = 0._r8
	  cnf%dwt_frootn_to_litr3n(c) = 0._r8
	  cnf%dwt_livecrootn_to_cwdn(c) = 0._r8
	  cnf%dwt_deadcrootn_to_cwdn(c) = 0._r8
	  cnf%dwt_nloss(c) = 0._r8
      
   end do

   ! deallocate special landunit filters
   deallocate(specialc)
   deallocate(specialp)

end subroutine CNiniSpecial
