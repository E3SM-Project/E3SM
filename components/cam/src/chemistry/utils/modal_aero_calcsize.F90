module modal_aero_calcsize

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use physconst,        only: pi, rhoh2o, gravit

use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state, physics_ptend
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use phys_control,     only: phys_getopts
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use shr_log_mod ,     only: errMsg => shr_log_errMsg
use cam_history,      only: addfld, horiz_only, add_default, fieldname_len, outfld
use constituents,     only: pcnst, cnst_name, cnst_get_ind

use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use time_manager,  only: is_first_step !DELETE IT BALLI


#ifdef MODAL_AERO

! these are the variables needed for the diagnostic calculation of dry radius
use modal_aero_data, only: ntot_amode, nspec_amode, &
                           numptr_amode, &
                           alnsg_amode, &
                           voltonumbhi_amode, voltonumblo_amode, &
                           dgnum_amode, dgnumhi_amode, dgnumlo_amode


! these variables are needed for the prognostic calculations to exchange mass
! between modes
use modal_aero_data,  only: numptrcw_amode, mprognum_amode, qqcw_get_field, lmassptrcw_amode, &
           lmassptr_amode, modeptr_accum, modeptr_aitken, ntot_aspectype, &
           lspectype_amode, specmw_amode, specdens_amode, voltonumb_amode, &
           cnst_name_cw

#endif


implicit none
private
save

public modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
public :: modal_aero_calcsize_reg

logical :: do_adjust_default
logical :: do_aitacc_transfer_default

!Mimic enumerators for aerosol types
integer, parameter:: inter_aero   = 1
integer, parameter:: cld_brn_aero = 2

integer :: dgnum_idx = -1

#ifdef MODAL_AERO
integer, parameter :: maxspec_csizxf = ntot_aspectype
#else
! TODO: this is a kludge.  This value should probably be assigned
! elsewhere for the non-modal case.  S.M. Burrows.
integer, parameter :: maxspec_csizxf = 8
#endif

real(r8), parameter :: third = 1.0_r8/3.0_r8

integer :: npair_csizxf = -123456789
integer :: modefrm_csizxf
integer :: modetoo_csizxf
integer :: nspecfrm_csizxf
integer :: lspecfrmc_csizxf(maxspec_csizxf)
integer :: lspecfrma_csizxf(maxspec_csizxf)
integer :: lspectooc_csizxf(maxspec_csizxf)
integer :: lspectooa_csizxf(maxspec_csizxf)

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_calcsize_reg()
  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

  integer :: nmodes

  call rad_cnst_get_info(0, nmodes=nmodes)

  call pbuf_add_field('DGNUM', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnum_idx)

end subroutine modal_aero_calcsize_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_init( pbuf2d, species_class)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   !    set do_adjust_default and do_aitacc_transfer_default flags
   !    create history fields for column tendencies associated with
   !       modal_aero_calcsize
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer, intent(in) :: species_class(:)

   ! local
   integer  :: ipair, iq, iqfrm, iqtoo
   integer  :: aer_type
   integer  :: lsfrm, lstoo, lsfrma, lsfrmc, lstooa, lstooc, lunout
   integer  :: mfrm, mtoo
   integer  :: n, nacc, nait, nspec
   integer  :: nchfrma, nchfrmc, nchfrmskip, nchtooa, nchtooc, nchtooskip
   logical  :: history_aerosol
   logical  :: history_verbose

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
   !-----------------------------------------------------------------------

   call phys_getopts(history_aerosol_out=history_aerosol, &
                     history_verbose_out=history_verbose)

   ! init entities required for both prescribed and prognostic modes

   if (is_first_step()) then ! bsingh- why do we need this conditional in init routine?
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnum_idx, 0.0_r8)
   endif

   npair_csizxf = 0
   modefrm_csizxf = 0
   modetoo_csizxf = 0

#ifndef MODAL_AERO
   do_adjust_default          = .false.
   do_aitacc_transfer_default = .false.

#else
   !  do_adjust_default allows adjustment to be turned on/off
   do_adjust_default = .true.

   !  do_aitacc_transfer_default allows aitken <--> accum mode transfer to be turned on/off
   !  *** it can only be true when aitken & accum modes are both present
   !      and have prognosed number and diagnosed surface/sigmag
   nait = modeptr_aitken
   nacc = modeptr_accum
   do_aitacc_transfer_default = .false.
   if ((modeptr_aitken > 0) .and.   &
      (modeptr_accum  > 0) .and.   &
      (modeptr_aitken /= modeptr_accum)) then
      do_aitacc_transfer_default = .true.
      if (mprognum_amode(nait) <= 0) do_aitacc_transfer_default = .false.
      if (mprognum_amode(nacc) <= 0) do_aitacc_transfer_default = .false.
   end if

   if ( .not. do_adjust_default ) return

!
!  define history fields for number-adjust source-sink for all modes
!

do_aitacc_transfer_if_block1: &
      if ( do_aitacc_transfer_default ) then
!
!   compute pointers for aitken <--> accum mode transfer
!	(a2 <--> a1 transfer)
!   transfers include number_a, number_c, mass_a, mass_c
!
      npair_csizxf = 1
      modefrm_csizxf = nait
      modetoo_csizxf = nacc

!
!   define species involved in each transfer pairing
!
aa_ipair: do ipair = 1, npair_csizxf

      mfrm = modefrm_csizxf
      mtoo = modetoo_csizxf
      if (mfrm < 10) then
          nchfrmskip = 1
      else if (mfrm < 100) then
          nchfrmskip = 2
      else
          nchfrmskip = 3
      end if
      if (mtoo < 10) then
          nchtooskip = 1
      else if (mtoo < 100) then
          nchtooskip = 2
      else
          nchtooskip = 3
      end if
      nspec = 0

aa_iqfrm: do iqfrm = -1, nspec_amode(mfrm)

         if (iqfrm == -1) then
            lsfrma = numptr_amode(mfrm)
            lstooa = numptr_amode(mtoo)
            lsfrmc = numptrcw_amode(mfrm)
            lstooc = numptrcw_amode(mtoo)
         else if (iqfrm == 0) then
!   bypass transfer of aerosol water due to calcsize transfer
            cycle aa_iqfrm
         else
            lsfrma = lmassptr_amode(iqfrm,mfrm)
            lsfrmc = lmassptrcw_amode(iqfrm,mfrm)
            lstooa = 0
            lstooc = 0
         end if

         if ((lsfrma < 1) .or. (lsfrma > pcnst)) then
            write(iulog,9100) mfrm, iqfrm, lsfrma
            call endrun( 'modal_aero_calcsize_init error aa' )
         end if
         if ((lsfrmc < 1) .or. (lsfrmc > pcnst)) then
            write(iulog,9102) mfrm, iqfrm, lsfrmc
            call endrun( 'modal_aero_calcsize_init error bb' )
         end if

         if (iqfrm > 0) then
            nchfrma = len( trim( cnst_name(lsfrma) ) ) - nchfrmskip

! find "too" species having same cnst_name as the "frm" species
! (except for last 1/2/3 characters which are the mode index)
            do iqtoo = 1, nspec_amode(mtoo)
               lstooa = lmassptr_amode(iqtoo,mtoo)
               nchtooa = len( trim( cnst_name(lstooa) ) ) - nchtooskip
               if (cnst_name(lsfrma)(1:nchfrma) == cnst_name(lstooa)(1:nchtooa)) then
               ! interstitial names match, so check cloudborne names too
                   nchfrmc = len( trim( cnst_name_cw(lsfrmc) ) ) - nchfrmskip
                   lstooc = lmassptrcw_amode(iqtoo,mtoo)
                   nchtooc = len( trim( cnst_name_cw(lstooc) ) ) - nchtooskip
                   if (cnst_name_cw(lsfrmc)(1:nchfrmc) /= &
                       cnst_name_cw(lstooc)(1:nchtooc)) lstooc = 0
                   exit
               else
                   lstooa = 0
               end if
            end do
         end if ! (iqfrm > 0)

         if ((lstooc < 1) .or. (lstooc > pcnst)) lstooc = 0
         if ((lstooa < 1) .or. (lstooa > pcnst)) lstooa = 0
         if (lstooa == 0) then
            write(iulog,9104) mfrm, iqfrm, lsfrma, iqtoo, lstooa
            call endrun( 'modal_aero_calcsize_init error cc' )
         end if
         if ((lstooc == 0) .and. (iqfrm /= 0)) then
            write(iulog,9104) mfrm, iqfrm, lsfrmc, iqtoo, lstooc
            call endrun( 'modal_aero_calcsize_init error dd' )
         end if

         nspec = nspec + 1
         lspecfrma_csizxf(nspec) = lsfrma
         lspectooa_csizxf(nspec) = lstooa
         lspecfrmc_csizxf(nspec) = lsfrmc
         lspectooc_csizxf(nspec) = lstooc
      end do aa_iqfrm

      nspecfrm_csizxf = nspec
      end do aa_ipair

9100  format( / '*** subr. modal_aero_calcsize_init' /   &
         'lspecfrma out of range' /   &
         'modefrm, ispecfrm, lspecfrma =', 3i6 / )
9102  format( / '*** subr. modal_aero_calcsize_init' /   &
         'lspecfrmc out of range' /   &
         'modefrm, ispecfrm, lspecfrmc =', 3i6 / )
9104  format( / '*** subr. modal_aero_calcsize_init' /   &
         'lspectooa out of range' /   &
         'modefrm, ispecfrm, lspecfrma, ispectoo, lspectooa =', 5i6 / )
9106  format( / '*** subr. modal_aero_calcsize_init' /   &
         'lspectooc out of range' /   &
         'modefrm, ispecfrm, lspecfrmc, ispectoo, lspectooc =', 5i6 / )

!
!   output results
!
      if ( masterproc ) then

      write(iulog,9310) do_adjust_default, do_aitacc_transfer_default

      do ipair = 1, npair_csizxf
      mfrm = modefrm_csizxf
      mtoo = modetoo_csizxf
      write(iulog,9320) ipair, mfrm, mtoo

      do iq = 1, nspecfrm_csizxf
         lsfrma = lspecfrma_csizxf(iq)
         lstooa = lspectooa_csizxf(iq)
         lsfrmc = lspecfrmc_csizxf(iq)
         lstooc = lspectooc_csizxf(iq)
         if (lstooa .gt. 0) then
            write(iulog,9330) lsfrma, cnst_name(lsfrma),   &
                               lstooa, cnst_name(lstooa)
         else
            write(iulog,9340) lsfrma, cnst_name(lsfrma)
         end if
         if (lstooc .gt. 0) then
            write(iulog,9330) lsfrmc, cnst_name_cw(lsfrmc),   &
                               lstooc, cnst_name_cw(lstooc)
         else if (lsfrmc .gt. 0) then
            write(iulog,9340) lsfrmc, cnst_name_cw(lsfrmc)
         else
            write(iulog,9350)
         end if
      end do ! iq

      end do ! ipair
      write(iulog,*)

      end if ! ( masterproc )


      else ! do_aitacc_transfer_if_block1

      npair_csizxf = 0
      if ( masterproc ) then
      write(iulog,9310) do_adjust_default, do_aitacc_transfer_default
      write(iulog,9320) 0, 0, 0
      end if

      end if do_aitacc_transfer_if_block1

9310  format( / 'subr. modal_aero_calcsize_init' / &
         'do_adjust_default, do_aitacc_transfer_default = ', 2l10 )
9320  format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330  format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340  format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350  format( 5x, 'no corresponding activated species' )



!  define history fields for number-adjust source-sink for all modes
do_adjust_if_block2: &
      if ( do_adjust_default ) then

   do n = 1, ntot_amode
      if (mprognum_amode(n) <= 0) cycle

      do aer_type = 1, 2
         if (aer_type == inter_aero) then
            tmpnamea = cnst_name(numptr_amode(n))
         else
            tmpnamea = cnst_name_cw(numptrcw_amode(n))
         end if
         unit = '#/m2/s'
         fieldname = trim(tmpnamea) // '_sfcsiz1'
         long_name = trim(tmpnamea) // ' calcsize number-adjust column source'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol .and. history_verbose) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(iulog,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnamea) // '_sfcsiz2'
         long_name = trim(tmpnamea) // ' calcsize number-adjust column sink'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if (history_aerosol .and. history_verbose) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(iulog,'(2a)') 'calcsize addfld - ', fieldname
      end do   ! aer_type = ...

   end do   ! n = ...


!  define history fields for aitken-accum transfer
do_aitacc_transfer_if_block2: &
      if ( do_aitacc_transfer_default ) then

! check that calcsize transfer ipair=1 is aitken-->accum
      ipair = 1
      if ((modefrm_csizxf .ne. nait) .or.   &
          (modetoo_csizxf .ne. nacc)) then
         write( iulog, '(//2a//)' )   &
            '*** modal_aero_calcaersize_init error -- ',   &
            'modefrm/too_csizxf(1) are wrong'
         call endrun( 'modal_aero_calcaersize_init error' )
      end if

      do iq = 1, nspecfrm_csizxf

! aer_type=1 does interstitial ("_a"); aer_type=2 does activated ("_c");
         do aer_type = 1, 2

! the lspecfrma_csizxf (and lspecfrmc_csizxf) are aitken species
! the lspectooa_csizxf (and lspectooc_csizxf) are accum  species
            if (aer_type .eq. inter_aero) then
               lsfrm = lspecfrma_csizxf(iq)
               lstoo = lspectooa_csizxf(iq)
            else
               lsfrm = lspecfrmc_csizxf(iq)
               lstoo = lspectooc_csizxf(iq)
            end if
            if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

            if (aer_type .eq. inter_aero) then
               tmpnamea = cnst_name(lsfrm)
               tmpnameb = cnst_name(lstoo)
            else
               tmpnamea = cnst_name_cw(lsfrm)
               tmpnameb = cnst_name_cw(lstoo)
            end if

            unit = 'kg/m2/s'
            if ((tmpnamea(1:3) == 'num') .or. &
               (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
            fieldname = trim(tmpnamea) // '_sfcsiz3'
            long_name = trim(tmpnamea) // ' calcsize aitken-to-accum adjust column tendency'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if (history_aerosol .and. history_verbose) then
               call add_default(fieldname, 1, ' ')
            end if
            if ( masterproc ) write(iulog,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnameb) // '_sfcsiz3'
            long_name = trim(tmpnameb) // ' calcsize aitken-to-accum adjust column tendency'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if (history_aerosol .and. history_verbose) then
               call add_default(fieldname, 1, ' ')
            end if
            if ( masterproc ) write(iulog,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnamea) // '_sfcsiz4'
            long_name = trim(tmpnamea) // ' calcsize accum-to-aitken adjust column tendency'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if (history_aerosol .and. history_verbose) then
               call add_default(fieldname, 1, ' ')
            end if
            if ( masterproc ) write(iulog,'(2a)') 'calcsize addfld - ', fieldname

            fieldname = trim(tmpnameb) // '_sfcsiz4'
            long_name = trim(tmpnameb) // ' calcsize accum-to-aitken adjust column tendency'
            call addfld( fieldname, horiz_only, 'A', unit, long_name )
            if (history_aerosol .and. history_verbose) then
               call add_default(fieldname, 1, ' ')
            end if
            if ( masterproc ) write(iulog,'(2a)') 'calcsize addfld - ', fieldname

         end do   ! aer_type = ...
      end do   ! iq = ...

      end if do_aitacc_transfer_if_block2

      end if do_adjust_if_block2


      if ( masterproc ) then
         write(iulog,'(/a)') 'l, species_class, name'
         do n = 1, pcnst
            write(iulog,'(2i4,2x,a)') n, species_class(n), cnst_name(n)
         end do
      end if
   if ( masterproc ) write(iulog,'(a)') 'modal_aero_calcsize_init ALL DONE'

#endif

end subroutine modal_aero_calcsize_init

!===============================================================================

subroutine modal_aero_calcsize_sub(state, pbuf, ptend, deltat, do_adjust_in, &
   do_aitacc_transfer_in, list_idx_in, update_mmr_in, dgnumdry_m,cp_buf,cp_num_buf)

  implicit none
   !-----------------------------------------------------------------------
   !
   ! Calculates aerosol size distribution parameters
   !    mprognum_amode >  0
   !       calculate Dgnum from mass, number, and fixed sigmag
   !    mprognum_amode <= 0
   !       calculate number from mass, fixed Dgnum, and fixed sigmag
   !
   ! Also (optionally) adjusts prognostic number to
   !    be within bounds determined by mass, Dgnum bounds, and sigma bounds
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), target, intent(in)    :: state       ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer
   type(physics_ptend), target, optional, intent(inout) :: ptend       ! indivdual parameterization tendencies
   real(r8),                    optional, intent(in)    :: deltat, cp_buf(pcols,pver,7,4),cp_num_buf(pcols,pver,4)      ! model time-step size (s)


   logical,  optional, intent(in) :: do_adjust_in
   logical,  optional, intent(in) :: do_aitacc_transfer_in
   logical,  optional, intent(in) :: update_mmr_in
   integer,  optional, intent(in) :: list_idx_in    ! diagnostic list index
   real(r8), optional, pointer    :: dgnumdry_m(:,:,:) ! interstital aerosol dry number mode radius (m)

#ifdef MODAL_AERO

   ! local

   logical :: do_adjust
   logical :: do_aitacc_transfer
   logical :: is_list_present, update_mmr

   integer  :: lchnk                ! chunk identifier
   integer  :: ncol                 ! number of columns
   integer  :: list_idx_local       !list idx local to this subroutine

   real(r8), parameter :: close_to_one = 1.0_r8 + 1.0e-15_r8
   real(r8), pointer :: pdel(:,:)   ! pressure thickness of levels
   real(r8), pointer :: state_q(:,:,:)    ! Tracer MR array

   logical,  pointer :: dotend(:)   ! flag for doing tendency
   real(r8), pointer :: dqdt(:,:,:) ! TMR tendency array

   real(r8), pointer :: dgncur_a(:,:,:)

   integer  :: iq, nspec, imode, klev, icol, idx
   integer  :: num_idx
   integer  :: num_mode_idx, num_cldbrn_mode_idx, lsfrm, lstoo
   integer  :: stat

   logical  :: dotendqqcw(pcnst)

   character(len=fieldname_len)   :: tmpnamea, tmpnameb, name
   character(len=fieldname_len+3) :: fieldname
   real(r8) :: deltatinv                     ! 1/deltat
   real(r8) :: dgncur_c(pcols,pver,ntot_amode)
   real(r8) :: dqqcwdt(pcols,pver,pcnst)     ! cloudborne TMR tendency array
   real(r8) :: drv_a_accsv(pcols,pver), drv_c_accsv(pcols,pver)
   real(r8) :: drv_a_aitsv(pcols,pver), drv_c_aitsv(pcols,pver)
   real(r8) :: drv_a_sv(pcols,pver,ntot_amode), drv_c_sv(pcols,pver,ntot_amode)
   real(r8) :: dryvol_a(pcols,pver)          ! interstital aerosol dry
   ! volume (cm^3/mol_air)
   real(r8) :: dryvol_c(pcols,pver)          ! activated aerosol dry volume
   ! to size bounds
   real(r8) :: fracadj                       ! deltat/tadj
   real(r8) :: num_a_accsv(pcols,pver), num_c_accsv(pcols,pver)
   real(r8) :: num_a_aitsv(pcols,pver), num_c_aitsv(pcols,pver)
   real(r8) :: num_a_sv(pcols,pver,ntot_amode), num_c_sv(pcols,pver,ntot_amode)
   real(r8) :: tadj                          ! adjustment time scale
   real(r8) :: tadjinv                       ! 1/tadj
   real(r8) :: v2ncur_a(pcols,pver,ntot_amode)
   real(r8) :: v2ncur_c(pcols,pver,ntot_amode)

   !---------------------------------------------------------------------
   ! "qsrflx" array:
   !----------------
   !process-specific column tracer tendencies
   ! 3rd index --
   !    1="standard" number adjust gain;
   !    2="standard" number adjust loss;
   !    3=aitken-->accum transfer;
   !    4=accum-->aitken transferk
   ! 4th index --
   !    1="a" species (interstitial); 2="c" species(cloud borne)
   !-----------------------------------------------------------------------
   integer, parameter :: nsrflx = 4   ! last dimension of qsrflx
   real(r8) :: qsrflx(pcols,pcnst,nsrflx,2)

   !BSINGH: THINGS TO DO:
   !2. In compute_dry_volume, find if _c spec for radiation diags is same as _a spec
   !3. In compute_dry_volume, see if two nested do loops for _a and _c can be combined
   !4. state_q should not be used as some specs may be missing in the diganostic calls
   !5. fldcw, qqcw_get_field should be called with the right index as some species may be missing from radiation diagnostics
   !4. deltat should NOT be optional
   !5. remove unused vars
   !6. update mmr should not happen for diagnostics, add a check!!!
   !7. add comments and author, date , called by and calls
   !9. why call master proc......is it tp print something during first time step? (aitken_accum_exchange)
   !10. cloud borne processes and interstital can be combined and refactored into one subroutine.(aitken_accum_exchange)
   !11. init and register routine should be re-worked/streamlined to remove redundant logics
   !12. improve update_mmr logic in refactored subroutines


   !-----------------------------------------------------------------------------------
   !Extract info about optional variables and initialize local variables accordingly
   !------------------------------------------------------------------------------------
   !
   !The default behavior is to update species mass mixing ratios (interstitial and
   !cloud borne), but there are other calls (such as radiation diagnostics and otherwise)
   !which might just require updated aerosol sizes without any update to species mmr
   !-------------------------------------------------------------------------------------

   update_mmr = .true. !update mmr for both interstitial and cloud borne aerosols
   if(present(update_mmr_in)) update_mmr = update_mmr_in

   !For adjusting aerosol sizes
   do_adjust = do_adjust_default
   if (present(do_adjust_in)) do_adjust = do_adjust_in


   !For transerfering aerosols between modes (Aitkem<-->Accumulation) based on their new size
   do_aitacc_transfer = do_aitacc_transfer_default
   if (present(do_aitacc_transfer_in))  do_aitacc_transfer = do_aitacc_transfer_in

   !Set list_idx_local and other sanity checks
   list_idx_local  = 0
   if(present(list_idx_in)) then
      list_idx_local  = list_idx_in
      if (.not. present(dgnumdry_m)) &
           call endrun('list_idx_in is present but dgnumdry_m pointer is missing'//errmsg(__FILE__,__LINE__))
      dgncur_a => dgnumdry_m(:,:,:)
   else
      call pbuf_get_field(pbuf, dgnum_idx, dgncur_a)
   endif

   !Set tendency variable based on the presence of ptend
   if(update_mmr) then
      if(.not. present(ptend)) call endrun('ptend should be present if update_mmr is true'//errmsg(__FILE__,__LINE__))
      dotend => ptend%lq
      dqdt   => ptend%q
      dotendqqcw(:)   = .false.
      dqqcwdt(:,:,:)  = 0.0_r8
      qsrflx(:,:,:,:) = 0.0_r8
   else
      dqqcwdt(:,:,:)  = huge(dqqcwdt)
      qsrflx(:,:,:,:) = huge(qsrflx)
   endif

   pdel     => state%pdel !balli add if for only update_mmr
   state_q  => state%q    !BSINGH - it is okay to use it for num mmr but not for specie mmr (as diagnostic call may miss some species)

   !----------------------------------------------------------------------------
   ! tadj = adjustment time scale for number, surface when they are prognosed
   !----------------------------------------------------------------------------
   tadj    = max( 86400.0_r8, deltat )
   tadjinv = 1.0_r8/(tadj*close_to_one)
   fracadj = max( 0.0_r8, min( 1.0_r8, deltat*tadjinv ) )

   !grid parameters
   ncol  = state%ncol  !# of columns
   lchnk = state%lchnk !chunk # !BALLI- add if for only for cld borne

   !inverse of time step
   deltatinv = 1.0_r8/(deltat*close_to_one)


   !Now compute dry diameter for both interstitial and cloud borne aerosols
   do imode = 1, ntot_amode

      !----------------------------------------------------------------------
      !Initialize all parameters to the default values for the mode
      !----------------------------------------------------------------------
      !interstitial
      call set_initial_sz_and_volumes(top_lev, pver, ncol, imode,    & !input
           dgncur_a, v2ncur_a, dryvol_a)                               !output

      !cloud borne
      call set_initial_sz_and_volumes(top_lev, pver, ncol, imode,    & !input
           dgncur_c, v2ncur_c, dryvol_c)                               !output

      !----------------------------------------------------------------------
      !Find # of species in this mode
      !----------------------------------------------------------------------
      !(for radiation diagnostics, # of species in a mode can be differnet from the default)
      call rad_cnst_get_info(list_idx_local, imode, nspec=nspec)

      !----------------------------------------------------------------------
      !Compute dry volume mixrats (aerosol diameter)
      !Current default: number mmr is prognosed
      !       Algorithm:calculate aerosol diameter from mass, number, and fixed sigmag
      !
      !sigmag ("sigma g") is "geometric standard deviation for aerosol mode"
      !
      !Volume = sum_over_components{ component_mass mixrat / density }
      !----------------------------------------------------------------------
      if(present(cp_buf)) then!<<TO BE REMOVED
         call compute_dry_volume(top_lev, pver, ncol, imode, nspec, state, pbuf, dryvol_a, dryvol_c, list_idx_local, cp_buf)
         !^^TO BE REMOVED
      else!<<TO BE REMOVED
         call compute_dry_volume(top_lev, pver, ncol, imode, nspec, state, pbuf, dryvol_a, dryvol_c, list_idx_local)
      endif

      ! do size adjustment based on computed dry diameter values
      call size_adjustment(top_lev, pver, ncol, lchnk, imode, dryvol_a, state_q, dryvol_c, pdel, & !input
           do_adjust, update_mmr, do_aitacc_transfer, deltatinv, fracadj, pbuf,                  & !input
           dgncur_a, dgncur_c, v2ncur_a, v2ncur_c,                                               & !output
           drv_a_accsv, drv_c_accsv, drv_a_aitsv, drv_c_aitsv, drv_a_sv, drv_c_sv,               & !output
           num_a_accsv, num_c_accsv, num_a_aitsv, num_c_aitsv, num_a_sv, num_c_sv,               & !output
           dotend, dotendqqcw, dqdt, dqqcwdt, qsrflx, cp_num_buf)                                  !output

   end do  ! do imode = 1, ntot_amode


   !
   !
   ! the following section (from here to label 49000)
   !    does aitken <--> accum mode transfer
   !
   ! when the aitken mode mean size is too big, the largest
   !    aitken particles are transferred into the accum mode
   !    to reduce the aitken mode mean size
   ! when the accum mode mean size is too small, the smallest
   !    accum particles are transferred into the aitken mode
   !    to increase the accum mode mean size
   !
   !
   if ( do_aitacc_transfer ) then
      if(present(cp_buf)) then!<<TO BE REMOVED
         call aitken_accum_exchange(ncol, lchnk, list_idx_local, tadjinv, &
              deltat, pdel, state_q, state, &
              pbuf, &
              drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv, &
              drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,  &
              dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
              dqdt, dqqcwdt, qsrflx, cp_buf)
      else!<<TO BE REMOVED
         call aitken_accum_exchange(ncol, lchnk, list_idx_local, tadjinv, &
              deltat, pdel, state_q, state, &
              pbuf, &
              drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv, &
              drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,  &
              dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
              dqdt, dqqcwdt, qsrflx)
      endif
   end if  !  do_aitacc_transfer

   lsfrm = -huge(lsfrm) !initialize


   !----------------------------------------------------------------------
   ! apply tendencies to cloud-borne species MRs
   ! Only if "update_mmr" is TRUE
   !----------------------------------------------------------------------
   if(update_mmr) then
      call update_cld_brn_mmr(top_lev, pver, ncol, lchnk, pcnst, deltat, pbuf, dotendqqcw, dqqcwdt)

      !----------------------------------------------------------------------
      ! do outfld calls
      !----------------------------------------------------------------------

      ! history fields for number-adjust source-sink for all modes
      if ( .not. do_adjust ) return ! return if adjustment not required

      do imode = 1, ntot_amode
         if (mprognum_amode(imode) <= 0) cycle

         !interstitial
         idx  = numptr_amode(imode)
         name = cnst_name(idx)
         call output_flds(name, idx, lchnk, inter_aero, qsrflx)

         !cloud borne
         idx = numptrcw_amode(imode)
         name = cnst_name_cw(idx)
         call output_flds(name, idx, lchnk, cld_brn_aero, qsrflx)
      end do   ! imode

      ! history fields for aitken-accum transfer
      if ( .not. do_aitacc_transfer ) return ! return if transfer of species not required


      do iq = 1, nspecfrm_csizxf

         ! aer_type=1 does interstitial ("_a"); aer_type=2 does activated ("_c");

         lsfrm = lspecfrma_csizxf(iq)
         lstoo = lspectooa_csizxf(iq)

         if ((lsfrm > 0) .and. (lstoo > 0)) then
            tmpnamea = cnst_name(lsfrm)
            tmpnameb = cnst_name(lstoo)
            fieldname = trim(tmpnamea) // '_sfcsiz3'
            call outfld( fieldname, qsrflx(:,lsfrm,3,inter_aero), pcols, lchnk)

            fieldname = trim(tmpnameb) // '_sfcsiz3'
            call outfld( fieldname, qsrflx(:,lstoo,3,inter_aero), pcols, lchnk)

            fieldname = trim(tmpnamea) // '_sfcsiz4'
            call outfld( fieldname, qsrflx(:,lsfrm,4,inter_aero), pcols, lchnk)

            fieldname = trim(tmpnameb) // '_sfcsiz4'
            call outfld( fieldname, qsrflx(:,lstoo,4,inter_aero), pcols, lchnk)
         endif

         lsfrm = lspecfrmc_csizxf(iq)
         lstoo = lspectooc_csizxf(iq)
         if ((lsfrm > 0) .and. (lstoo > 0)) then
            tmpnamea = cnst_name_cw(lsfrm)
            tmpnameb = cnst_name_cw(lstoo)

            fieldname = trim(tmpnamea) // '_sfcsiz3'
            call outfld( fieldname, qsrflx(:,lsfrm,3,cld_brn_aero), pcols, lchnk)

            fieldname = trim(tmpnameb) // '_sfcsiz3'
            call outfld( fieldname, qsrflx(:,lstoo,3,cld_brn_aero), pcols, lchnk)

            fieldname = trim(tmpnamea) // '_sfcsiz4'
            call outfld( fieldname, qsrflx(:,lsfrm,4,cld_brn_aero), pcols, lchnk)

            fieldname = trim(tmpnameb) // '_sfcsiz4'
            call outfld( fieldname, qsrflx(:,lstoo,4,cld_brn_aero), pcols, lchnk)
         endif
      end do   ! iq = ...

   endif!if(update_mmr)
#endif

return
end subroutine modal_aero_calcsize_sub

!---------------------------------------------------------------------------------------------

subroutine set_initial_sz_and_volumes(top_lev, pver, ncol, imode, & !input
     dgncur, v2ncur, dryvol )                                       !output

  !-----------------------------------------------------------------------------
  !Purpose: Set initial defaults for the dry diameter, volum to num
  ! and dry volume
  !
  !Called by: modal_aero_calcsize_sub
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------
  implicit none
  !inputs
  integer, intent(in) :: top_lev, pver !for model level loop
  integer, intent(in) :: ncol          !# of columns
  integer, intent(in) :: imode         !mode index

  !outputs
  real(r8), intent(out) :: dgncur(:,:,:) !diameter
  real(r8), intent(out) :: v2ncur(:,:,:) !volume to number
  real(r8), intent(out) :: dryvol(:,:)   !dry volume

  !local variables
  integer :: icol, klev

  do klev = top_lev, pver
     do icol = 1, ncol
        dgncur(icol,klev,imode) = dgnum_amode(imode)     !diameter
        v2ncur(icol,klev,imode) = voltonumb_amode(imode) !volume to number
        dryvol(icol,klev)       = 0.0_r8                 !initialize dry vol
     end do
  end do

  return

end subroutine set_initial_sz_and_volumes

!---------------------------------------------------------------------------------------------

subroutine compute_dry_volume(top_lev, pver, ncol, imode, nspec, state, pbuf, &
     dryvol_a, dryvol_c, list_idx_in, cp_buf)

  !-----------------------------------------------------------------------------
  !Purpose: Compute initial dry volume based on mmr and specie density
  ! volume = mmr/density
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : rad_cnst_get_aer_mmr, rad_cnst_get_aer_props
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: top_lev, pver !for model level loop
  integer,  intent(in) :: ncol          !# of columns
  integer,  intent(in) :: imode         !mode index
  integer,  intent(in) :: nspec         !# of species in mode "imode"
  integer,  intent(in) :: list_idx_in
  type(physics_state), target, intent(in)    :: state       ! Physics state variables
  type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer
  real(r8), intent(in), optional    :: cp_buf(:,:,:,:)

  !in-outs
  real(r8), intent(inout) :: dryvol_a(:,:)                    ! interstital aerosol dry
  real(r8), intent(inout) :: dryvol_c(:,:)                    ! interstital aerosol dry

  !local vars
  integer  :: ispec, icol, klev, lchnk, idx_cw
  real(r8) :: specdens  !specie density
  real(r8) :: dummwdens !density inverse
  real(r8), pointer :: specmmr(:,:)  !specie mmr (interstitial)
  real(r8), pointer :: fldcw(:,:)    !specie mmr (cloud borne)

  character(len=32) :: spec_name

  lchnk = state%lchnk !get chunk info for retrieving cloud borne aerosols mmr

  do ispec = 1, nspec
     ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air

     !Get mmr and density for a specie
     call rad_cnst_get_aer_mmr(list_idx_in, imode, ispec, 'a', state, pbuf, specmmr) !get mmr
     call rad_cnst_get_aer_props(list_idx_in, imode, ispec, density_aer=specdens)    !get density
     dummwdens = 1.0_r8 / specdens !inverse of density

     !compute dry volume as a function of space (i,k)
     do klev = top_lev, pver
        do icol = 1, ncol
           dryvol_a(icol,klev) = dryvol_a(icol,klev) + max(0.0_r8,specmmr(icol,klev))*dummwdens
        end do
     end do

     !For radiation diagnostic calls, some species might be missing. We accomodate that,
     !we first get the specie name and then find its index from cnst_name array
     !This index is exactly the same as pbuf index for that specie
     call rad_cnst_get_info(list_idx_in, imode, ispec,spec_name=spec_name)
     call cnst_get_ind(trim(adjustl(spec_name)), idx_cw)

     !fldcw => qqcw_get_field(pbuf,lmassptrcw_amode(ispec,imode),lchnk)
     fldcw => qqcw_get_field(pbuf,idx_cw,lchnk)
     !######################## TO BE REMOVED
     if(present(cp_buf)) then
        do klev=top_lev,pver
           do icol=1,ncol
              dryvol_c(icol,klev) = dryvol_c(icol,klev)    &
                   + max(0.0_r8,cp_buf(icol,klev,ispec,imode))*dummwdens
           end do
        end do
     else
     !######################## TO BE REMOVED ^^^^^^^^^^^^^^^
        do klev = top_lev, pver
           do icol = 1, ncol
              dryvol_c(icol,klev) = dryvol_c(icol,klev)    &
                   + max(0.0_r8,fldcw(icol,klev))*dummwdens
           end do
        end do
     endif
  end do

  return
end subroutine compute_dry_volume

!---------------------------------------------------------------------------------------------

subroutine size_adjustment(top_lev, pver, ncol, lchnk, imode, dryvol_a, state_q, dryvol_c, pdel, &
     do_adjust, update_mmr, do_aitacc_transfer, deltatinv, fracadj, pbuf, &
     dgncur_a, dgncur_c, v2ncur_a, v2ncur_c, &
     drv_a_accsv, drv_c_accsv, drv_a_aitsv, drv_c_aitsv, drv_a_sv, drv_c_sv, &
     num_a_accsv, num_c_accsv, num_a_aitsv, num_c_aitsv, num_a_sv, num_c_sv, &
     dotend, dotendqqcw, dqdt, dqqcwdt, qsrflx, cp_num_buf)

  !-----------------------------------------------------------------------------
  !Purpose: Do the aerosol size adjustment if needed
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : compute_dgn_vol_limits, adjust_num_sizes, update_dgn_voltonum
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: top_lev, pver !for model level loop
  integer,  intent(in) :: ncol          !# of columns
  integer,  intent(in) :: lchnk
  integer,  intent(in) :: imode         !mode index
  real(r8), intent(in) :: deltatinv
  real(r8), intent(in) :: fracadj
  real(r8), intent(in) :: pdel(:,:)
  real(r8), intent(in) :: dryvol_a(:,:), dryvol_c(:,:)
  real(r8), intent(in) :: state_q(:,:,:)

  logical, intent(in) :: do_adjust, update_mmr, do_aitacc_transfer
  type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer

  !TO BE REMOVED
  real(r8), intent(in), optional :: cp_num_buf(:,:,:)

  !outputs
  real(r8), intent(inout) :: dgncur_a(:,:,:) !BALLI- comments
  real(r8), intent(inout) :: dgncur_c(:,:,:)
  real(r8), intent(inout) :: v2ncur_a(:,:,:)
  real(r8), intent(inout) :: v2ncur_c(:,:,:)
  real(r8), intent(inout) :: drv_a_accsv(pcols,pver), drv_c_accsv(pcols,pver)
  real(r8), intent(inout) :: drv_a_aitsv(pcols,pver), drv_c_aitsv(pcols,pver)
  real(r8), intent(inout) :: drv_a_sv(pcols,pver,ntot_amode), drv_c_sv(pcols,pver,ntot_amode)
  real(r8), intent(inout) :: num_a_accsv(pcols,pver), num_c_accsv(pcols,pver)
  real(r8), intent(inout) :: num_a_aitsv(pcols,pver), num_c_aitsv(pcols,pver)
  real(r8), intent(inout) :: num_a_sv(pcols,pver,ntot_amode), num_c_sv(pcols,pver,ntot_amode)

  logical,  intent(inout), optional :: dotend(:), dotendqqcw(:)
  real(r8), intent(inout), optional :: dqdt(:,:,:), dqqcwdt(:,:,:), qsrflx(:,:,:,:)

  !local
  integer :: klev, icol
  integer :: num_mode_idx, num_cldbrn_mode_idx, nait, nacc
  real(r8) :: v2nyy, v2nxx                  ! voltonumblo/hi of current mode
  real(r8) :: v2nyyrl, v2nxxrl              ! relaxed voltonumblo/hi
  real(r8) :: dgnyy, dgnxx                  ! dgnumlo/hi of current mode
  real(r8) :: drv_a, drv_c                  ! dry volume (cm3/mol_air)
  real(r8) :: num_a0, num_c0                ! initial number (#/mol_air)
  real(r8) :: num_a, num_c                  ! working number (#/mol_air)
  real(r8) :: cmn_factor
  real(r8), pointer :: fldcw(:,:)           !specie mmr (cloud borne)

  !find state q array number mode indices for interstitial and cloud borne aerosols
  num_mode_idx        = numptr_amode(imode)
  num_cldbrn_mode_idx = numptrcw_amode(imode)

  nait = modeptr_aitken !aitken mode number
  nacc = modeptr_accum  !accumulation mode number

  !set hi/lo limits (min, max and their relaxed counterparts) for volumes and diameters
  call compute_dgn_vol_limits(imode, nait, nacc, do_aitacc_transfer, & !input
       v2nxx, v2nyy, v2nxxrl, v2nyyrl, dgnxx, dgnyy) !output

  !set tendency logicals to true if we need to update mmr
  if (update_mmr) then
     dotend(num_mode_idx)            = .true.
     dotendqqcw(num_cldbrn_mode_idx) = .true.
  end if

  !pointer to cloud borne number mmr for mode imode
  fldcw => qqcw_get_field(pbuf,num_cldbrn_mode_idx,lchnk,.true.)

  do  klev = top_lev, pver
     do  icol = 1, ncol

        drv_a  = dryvol_a(icol,klev)
        num_a0 = state_q(icol,klev,num_mode_idx)
        num_a  = max( 0.0_r8, num_a0 )
        drv_c  = dryvol_c(icol,klev)
        if(present(cp_num_buf)) then
           num_c0 = cp_num_buf(icol,klev,imode)
        else
           num_c0 = fldcw(icol,klev)
        endif
        num_c = max( 0.0_r8, num_c0 )

        if (do_adjust) then

           !-----------------------------------------------------------------
           ! Do number adjustment for interstitial and activated particles
           !-----------------------------------------------------------------
           !Adjustments that:
           !(1) make numbers non-negative or
           !(2) make numbers zero when volume is zero
           !are applied over time-scale deltat
           !Adjustments that bring numbers to within specified bounds are
           !applied over time-scale tadj
           !-----------------------------------------------------------------
           if(update_mmr) then
              call adjust_num_sizes(icol, klev, num_mode_idx, num_cldbrn_mode_idx, &                   !input
                   drv_a, num_a0, drv_c, num_c0, deltatinv, v2nxx, v2nxxrl, v2nyy, v2nyyrl, fracadj, & !input
                   num_a, num_c, dqdt, dqqcwdt)                                                        !output

           else
              call adjust_num_sizes(icol, klev, num_mode_idx, num_cldbrn_mode_idx, &                   !input
                   drv_a, num_a0, drv_c, num_c0, deltatinv, v2nxx, v2nxxrl, v2nyy, v2nyyrl, fracadj, & !input
                   num_a, num_c)                                                                       !output

           endif
        endif !do_adjust

        !Compute a common factor for size computations
        cmn_factor = exp(4.5_r8*alnsg_amode(imode)**2)*pi/6.0_r8

        !
        ! now compute current dgn and v2n
        !
        if(update_mmr) then
           call update_dgn_voltonum(icol, klev, imode, num_mode_idx, inter_aero, gravit, cmn_factor, drv_a, num_a, &
                v2nxx, v2nyy, dgnxx, dgnyy, pdel, &
                dgncur_a, v2ncur_a, qsrflx, dqdt)

           !for cloud borne aerosols
           call update_dgn_voltonum(icol, klev, imode, num_cldbrn_mode_idx, cld_brn_aero, gravit, cmn_factor, drv_c, num_c, &
                v2nxx, v2nyy, dgnumhi_amode(imode), dgnumlo_amode(imode), pdel, &
                dgncur_c, v2ncur_c, qsrflx, dqqcwdt)
        else
           call update_dgn_voltonum(icol, klev, imode, num_mode_idx, inter_aero, gravit, cmn_factor, drv_a, num_a, &
                v2nxx, v2nyy, dgnxx, dgnyy, pdel, &
                dgncur_a, v2ncur_a)

           !for cloud borne aerosols
           call update_dgn_voltonum(icol, klev, imode, num_cldbrn_mode_idx, cld_brn_aero, gravit, cmn_factor, drv_c, num_c, &
                v2nxx, v2nyy, dgnumhi_amode(imode), dgnumlo_amode(imode), pdel, &
                dgncur_c, v2ncur_c)
        endif

        ! save number and dryvol for aitken <--> accum transfer
        if ( do_aitacc_transfer ) then
           if (imode == nait) then
              drv_a_aitsv(icol,klev) = drv_a
              num_a_aitsv(icol,klev) = num_a
              drv_c_aitsv(icol,klev) = drv_c
              num_c_aitsv(icol,klev) = num_c
           else if (imode == nacc) then
              drv_a_accsv(icol,klev) = drv_a
              num_a_accsv(icol,klev) = num_a
              drv_c_accsv(icol,klev) = drv_c
              num_c_accsv(icol,klev) = num_c
           end if
        end if
        drv_a_sv(icol,klev,imode) = drv_a
        num_a_sv(icol,klev,imode) = num_a
        drv_c_sv(icol,klev,imode) = drv_c
        num_c_sv(icol,klev,imode) = num_c

     end do
  end do

  return

end subroutine size_adjustment

!---------------------------------------------------------------------------------------------

subroutine compute_dgn_vol_limits(imode, nait, nacc, do_aitacc_transfer, & !input
           v2nxx, v2nyy, v2nxxrl, v2nyyrl, dgnxx, dgnyy) !output

  !-----------------------------------------------------------------------------
  !Purpose: Compute hi and lo limits for diameter and volume based on the
  ! relaxation factor
  !
  !Called by: size_adjustment
  !Calls    : None
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

implicit none

!inputs
integer,  intent(in) :: imode
integer,  intent(in) :: nait, nacc
logical,  intent(in) :: do_aitacc_transfer

!outputs
real(r8), intent(out) :: v2nxx
real(r8), intent(out) :: v2nyy
real(r8), intent(out) :: v2nxxrl
real(r8), intent(out) :: v2nyyrl
real(r8), intent(out) :: dgnxx
real(r8), intent(out) :: dgnyy

!local
real(r8), parameter :: relax_factor = 27.0_r8 !relax_factor=3**3=27,
                                              !i.e. dgnumlo_relaxed = dgnumlo/3 and dgnumhi_relaxed = dgnumhi*3


!v2nxx = voltonumbhi_amode is proportional to dgnumhi**(-3),
!        and produces the minimum allowed number for a given volume
v2nxx   = voltonumbhi_amode(imode)

!v2nyy = voltonumblo_amode is proportional to dgnumlo**(-3),
!        and produces the maximum allowed number for a given volume
v2nyy   = voltonumblo_amode(imode)

!v2nxxrl and v2nyyrl are their "relaxed" equivalents.
v2nxxrl = v2nxx/relax_factor
v2nyyrl = v2nyy*relax_factor
dgnxx   = dgnumhi_amode(imode)
dgnyy   = dgnumlo_amode(imode)

if ( do_aitacc_transfer ) then
   !for n=nait, divide v2nxx by 1.0e6 to effectively turn off the
   !         adjustment when number is too small (size is too big)
   if (imode == nait) v2nxx = v2nxx/1.0e6_r8

   !for n=nacc, multiply v2nyy by 1.0e6 to effectively turn off the
   !         adjustment when number is too big (size is too small)
   if (imode == nacc) v2nyy = v2nyy*1.0e6_r8

   !Also change the v2nyyrl/v2nxxrl so that
   !the interstitial<-->activated adjustment is turned off
   v2nxxrl = v2nxx/relax_factor
   v2nyyrl = v2nyy*relax_factor
end if

return
end subroutine compute_dgn_vol_limits

!---------------------------------------------------------------------------------------------

subroutine adjust_num_sizes(icol, klev, num_mode_idx, num_cldbrn_mode_idx, &
     drv_a, num_a0, drv_c, num_c0, deltatinv, v2nxx, v2nxxrl, v2nyy, v2nyyrl, fracadj, &
     num_a, num_c, dqdt, dqqcwdt)

  !-----------------------------------------------------------------------------
  !Purpose: Adjust num sizes if needed
  !
  !Called by: size_adjustment
  !Calls    : None
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: icol, klev
  real(r8), intent(in) :: drv_a, num_a0, drv_c, num_c0
  real(r8), intent(in) :: deltatinv, v2nxx, v2nxxrl, v2nyy, v2nyyrl, fracadj
  integer,  intent(in) :: num_mode_idx, num_cldbrn_mode_idx

  !outputs
  real(r8), intent(inout) :: num_a, num_c
  real(r8), intent(inout), optional :: dqdt(:,:,:), dqqcwdt(:,:,:)

  !local
  logical  :: update_mmr
  real(r8) :: num_a1, num_c1, numbnd
  real(r8) :: num_a2, num_c2, delnum_a2, delnum_c2
  real(r8) :: delnum_a3, delnum_c3
  real(r8) :: num_t2, drv_t, delnum_t3

  update_mmr = .false.
  if(present(dqdt) .and. present(dqqcwdt)) update_mmr = .true.


  if ((drv_a <= 0.0_r8) .and. (drv_c <= 0.0_r8)) then
     ! both interstitial and activated volumes are zero
     ! adjust both numbers to zero
     num_a = 0.0_r8
     num_c = 0.0_r8
     if(update_mmr) then
        dqdt(icol,klev,num_mode_idx)           = -num_a0*deltatinv
        dqqcwdt(icol,klev,num_cldbrn_mode_idx) = -num_c0*deltatinv
     endif
  else if (drv_c <= 0.0_r8) then
     ! activated volume is zero, so interstitial number/volume == total/combined
     ! apply step 1 and 3, but skip the relaxed adjustment (step 2, see below)
     num_c = 0.0_r8

     num_a1 = num_a
     numbnd = max( drv_a*v2nxx, min( drv_a*v2nyy, num_a1 ) )
     num_a  = num_a1 + (numbnd - num_a1)*fracadj
     if(update_mmr) then
        dqdt(icol,klev,num_mode_idx)           = (num_a - num_a0)*deltatinv
        dqqcwdt(icol,klev,num_cldbrn_mode_idx) = -num_c0*deltatinv
     endif
  else if (drv_a <= 0.0_r8) then
     ! interstitial volume is zero, treat similar to above
     num_a = 0.0_r8
     num_c1 = num_c
     numbnd = max( drv_c*v2nxx, min( drv_c*v2nyy, num_c1 ) )
     num_c  = num_c1 + (numbnd - num_c1)*fracadj
     if(update_mmr) then
        dqdt(icol,klev,num_mode_idx)           = -num_a0*deltatinv
        dqqcwdt(icol,klev,num_cldbrn_mode_idx) = (num_c - num_c0)*deltatinv
     endif
  else
     ! both volumes are positive
     ! apply 3 adjustment steps
     ! step1:  num_a,c0 --> num_a,c1 forces non-negative values
     num_a1 = num_a
     num_c1 = num_c
     ! step2:  num_a,c1 --> num_a,c2 applies relaxed bounds to the interstitial
     !    and activated number (individually)
     !    if only a or c changes, adjust the other in the opposite direction
     !    as much as possible to conserve a+c
     numbnd = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl, num_a1 ) )
     delnum_a2 = (numbnd - num_a1)*fracadj
     num_a2 = num_a1 + delnum_a2
     numbnd = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl, num_c1 ) )
     delnum_c2 = (numbnd - num_c1)*fracadj
     num_c2 = num_c1 + delnum_c2
     if ((delnum_a2 == 0.0_r8) .and. (delnum_c2 /= 0.0_r8)) then
        num_a2 = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl,   &
             num_a1-delnum_c2 ) )
     else if ((delnum_a2 /= 0.0_r8) .and. (delnum_c2 == 0.0_r8)) then
        num_c2 = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl,   &
             num_c1-delnum_a2 ) )
     end if
     ! step3:  num_a,c2 --> num_a,c3 applies stricter bounds to the
     !    combined/total number
     drv_t = drv_a + drv_c
     num_t2 = num_a2 + num_c2
     delnum_a3 = 0.0_r8
     delnum_c3 = 0.0_r8
     if (num_t2 < drv_t*v2nxx) then
        delnum_t3 = (drv_t*v2nxx - num_t2)*fracadj
        ! if you are here then (num_a2 < drv_a*v2nxx) and/or
        !                      (num_c2 < drv_c*v2nxx) must be true
        if ((num_a2 < drv_a*v2nxx) .and. (num_c2 < drv_c*v2nxx)) then
           delnum_a3 = delnum_t3*(num_a2/num_t2)
           delnum_c3 = delnum_t3*(num_c2/num_t2)
        else if (num_c2 < drv_c*v2nxx) then
           delnum_c3 = delnum_t3
        else if (num_a2 < drv_a*v2nxx) then
           delnum_a3 = delnum_t3
        end if
     else if (num_t2 > drv_t*v2nyy) then
        delnum_t3 = (drv_t*v2nyy - num_t2)*fracadj
        ! if you are here then (num_a2 > drv_a*v2nyy) and/or
        !                      (num_c2 > drv_c*v2nyy) must be true
        if ((num_a2 > drv_a*v2nyy) .and. (num_c2 > drv_c*v2nyy)) then
           delnum_a3 = delnum_t3*(num_a2/num_t2)
           delnum_c3 = delnum_t3*(num_c2/num_t2)
        else if (num_c2 > drv_c*v2nyy) then
           delnum_c3 = delnum_t3
        else if (num_a2 > drv_a*v2nyy) then
           delnum_a3 = delnum_t3
        end if
     end if
     num_a = num_a2 + delnum_a3

     num_c = num_c2 + delnum_c3
     if(update_mmr) then
        dqdt(icol,klev,num_mode_idx)           = (num_a - num_a0)*deltatinv
        dqqcwdt(icol,klev,num_cldbrn_mode_idx) = (num_c - num_c0)*deltatinv
     endif
  end if

  return
end subroutine adjust_num_sizes

!---------------------------------------------------------------------------------------------

subroutine update_dgn_voltonum(icol, klev, imode, num_idx, aer_type, gravit, cmn_factor, drv, num, &
     v2nxx, v2nyy, dgnxx, dgnyy, pdel, &
     dgncur, v2ncur, qsrflx, dqdt)

  !-----------------------------------------------------------------------------
  !Purpose: updates diameter and colume to num based on limits
  !
  !Called by: size_adjustment
  !Calls    : None
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: icol, klev, imode, num_idx, aer_type
  real(r8), intent(in) :: gravit, cmn_factor
  real(r8), intent(in) :: drv, num
  real(r8), intent(in) :: v2nxx, v2nyy, dgnxx, dgnyy
  real(r8), intent(in) :: pdel(:,:)
  real(r8), intent(in), optional :: dqdt(:,:,:)
  !output
  real(r8), intent(inout) :: dgncur(:,:,:), v2ncur(:,:,:)
  real(r8), intent(inout), optional :: qsrflx(:,:,:,:)

  !local
  logical  :: update_mmr
  real(r8) :: pdel_fac

  update_mmr = .false.
  if(present(qsrflx) .and. present(dqdt)) update_mmr = .true.

  if (drv > 0.0_r8) then
     if (num <= drv*v2nxx) then
        dgncur(icol,klev,imode) = dgnxx
        v2ncur(icol,klev,imode) = v2nxx
     else if (num >= drv*v2nyy) then
        dgncur(icol,klev,imode) = dgnyy
        v2ncur(icol,klev,imode) = v2nyy
     else
        dgncur(icol,klev,imode) = (drv/(cmn_factor*num))**third
        v2ncur(icol,klev,imode) = num/drv
     end if
  end if
  if (update_mmr) then
     pdel_fac = pdel(icol,klev)/gravit   ! = rho*dz
     qsrflx(icol,num_idx,1,aer_type) = qsrflx(icol,num_idx,1,aer_type) + max(0.0_r8,dqdt(icol,klev,num_idx))*pdel_fac
     qsrflx(icol,num_idx,2,aer_type) = qsrflx(icol,num_idx,2,aer_type) + min(0.0_r8,dqdt(icol,klev,num_idx))*pdel_fac
  endif

  return
end subroutine update_dgn_voltonum

!---------------------------------------------------------------------------------------------

subroutine  aitken_accum_exchange(ncol, lchnk, list_idx, tadjinv, &
     deltat, pdel, state_q, state, &
     pbuf, &
     drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv,     &
     drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,      &
     dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
     dqdt, dqqcwdt, qsrflx, cp_buf )

  !-----------------------------------------------------------------------------
  !Purpose: Exchange aerosols between aitken and accumulation modes based on new
  ! sizes
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : endrun
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: ncol, lchnk, list_idx
  real(r8), intent(in) :: tadjinv
  real(r8), intent(in) :: deltat
  real(r8), intent(in) :: pdel(:,:)
  real(r8), intent(in) :: state_q(:,:,:)
  real(r8), intent(in) :: drv_a_aitsv(:,:), num_a_aitsv(:,:)
  real(r8), intent(in) :: drv_a_accsv(:,:), num_a_accsv(:,:)
  real(r8), intent(in) :: drv_c_aitsv(:,:), num_c_aitsv(:,:)
  real(r8), intent(in) :: drv_c_accsv(:,:), num_c_accsv(:,:)
  type(physics_state), intent(in)    :: state       ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer
  real(r8), intent(in), optional    :: cp_buf(:,:,:,:)! TO BE REMOVED

  !outputs
  real(r8), intent(inout) :: dgncur_a(:,:,:) !BALLI- comments
  real(r8), intent(inout) :: dgncur_c(:,:,:)
  real(r8), intent(inout) :: v2ncur_a(:,:,:)
  real(r8), intent(inout) :: v2ncur_c(:,:,:)
  logical,  intent(inout), optional :: dotend(:), dotendqqcw(:)
  real(r8), intent(inout), optional :: dqdt(:,:,:), dqqcwdt(:,:,:), qsrflx(:,:,:,:)

  !local
  integer  :: imode, jmode, aer_type, jsrflx, icol, klev, iq
  integer  :: lsfrm, lstoo, ispec, idx, n_use, m_use, nb, mb
  integer  :: ixfer_ait2acc, ixfer_acc2ait
  integer  :: nait, nacc, idx_cw
  integer, save  :: idiagaa = 1
  real(r8) :: num_a, drv_a, num_c, drv_c
  real(r8) :: num_a_acc, num_c_acc
  real(r8) :: drv_a_acc, drv_c_acc
  real(r8) :: v2nzz, pdel_fac, num_t, drv_t
  real(r8) :: drv_a_noxf, drv_c_noxf, drv_t_noxf, dummwdens
  real(r8) :: num_t0, num_t_noxf
  real(r8) :: duma, cmn_factor       ! work variables
  real(r8) :: xfercoef
  real(r8) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
  real(r8) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
  real(r8) :: xferfrac_num_acc2ait, xferfrac_vol_acc2ait
  real(r8) :: xferfrac_num_ait2acc, xferfrac_vol_ait2acc
  real(r8) :: xfertend, xfertend_num(2,2)
  real(r8), pointer :: specmmr(:,:)  !specie mmr (interstitial)
  real(r8), pointer :: fldcw(:,:)    !specie mmr (cloud borne)
  logical  :: update_mmr, noxf_acc2ait(ntot_aspectype)

  character(len=32) :: spec_name

  update_mmr = .false.
  if(present(dotend) .and. present(dotendqqcw)) update_mmr = .true.

  if (npair_csizxf .le. 0)call endrun('npair_csizxf <= 0'//errmsg(__FILE__,__LINE__))

  ! check that calcsize transfer aitken-->accum
  nait = modeptr_aitken !aitken mode number
  nacc = modeptr_accum  !accumulation mode number

  if ((modefrm_csizxf .ne. nait) .or.   &
       (modetoo_csizxf .ne. nacc)) call endrun('modefrm/too_csizxf(1) are wrong'//errmsg(__FILE__,__LINE__))

  ! set dotend() for species that will be transferred
  do iq = 1, nspecfrm_csizxf
     lsfrm = lspecfrma_csizxf(iq)
     lstoo = lspectooa_csizxf(iq)
     if ((lsfrm > 0) .and. (lstoo > 0)) then
        if(update_mmr)dotend(lsfrm) = .true.
        if(update_mmr)dotend(lstoo) = .true.
     end if
     lsfrm = lspecfrmc_csizxf(iq)
     lstoo = lspectooc_csizxf(iq)
     if ((lsfrm > 0) .and. (lstoo > 0)) then
        if(update_mmr)dotendqqcw(lsfrm) = .true.
        if(update_mmr)dotendqqcw(lstoo) = .true.
     end if
  end do

  ! identify accum species cannot be transferred to aitken mode
  noxf_acc2ait(:) = .true.
  do ispec = 1, nspec_amode(nacc)
     idx = lmassptr_amode(ispec,nacc)
     do iq = 1, nspecfrm_csizxf
        if (lspectooa_csizxf(iq) == idx) then
           noxf_acc2ait(ispec) = .false.
        end if
     end do
  end do

  ! v2nzz is voltonumb at the "geometrically-defined" mid-point
  ! between the aitken and accum modes
  v2nzz = sqrt(voltonumb_amode(nait)*voltonumb_amode(nacc))

  ! loop over columns and levels
  do  klev = top_lev, pver
     do  icol = 1, ncol

        pdel_fac = pdel(icol,klev)/gravit   ! = rho*dz
        xfertend_num(:,:) = 0.0_r8

        ! compute aitken --> accum transfer rates
        ixfer_ait2acc = 0
        xfercoef_num_ait2acc = 0.0_r8
        xfercoef_vol_ait2acc = 0.0_r8

        drv_t = drv_a_aitsv(icol,klev) + drv_c_aitsv(icol,klev)
        num_t = num_a_aitsv(icol,klev) + num_c_aitsv(icol,klev)
        if (drv_t > 0.0_r8) then
           if (num_t < drv_t*v2nzz) then
              ixfer_ait2acc = 1
              if (num_t < drv_t*voltonumb_amode(nacc)) then
                 xferfrac_num_ait2acc = 1.0_r8
                 xferfrac_vol_ait2acc = 1.0_r8
              else
                 xferfrac_vol_ait2acc = ((num_t/drv_t) - v2nzz)/   &
                      (voltonumb_amode(nacc) - v2nzz)
                 xferfrac_num_ait2acc = xferfrac_vol_ait2acc*   &
                      (drv_t*voltonumb_amode(nacc)/num_t)
                 if ((xferfrac_num_ait2acc <= 0.0_r8) .or.   &
                      (xferfrac_vol_ait2acc <= 0.0_r8)) then
                    xferfrac_num_ait2acc = 0.0_r8
                    xferfrac_vol_ait2acc = 0.0_r8
                 else if ((xferfrac_num_ait2acc >= 1.0_r8) .or.   &
                      (xferfrac_vol_ait2acc >= 1.0_r8)) then
                    xferfrac_num_ait2acc = 1.0_r8
                    xferfrac_vol_ait2acc = 1.0_r8
                 end if
              end if
              xfercoef_num_ait2acc = xferfrac_num_ait2acc*tadjinv
              xfercoef_vol_ait2acc = xferfrac_vol_ait2acc*tadjinv
              xfertend_num(1,1) = num_a_aitsv(icol,klev)*xfercoef_num_ait2acc
              xfertend_num(1,2) = num_c_aitsv(icol,klev)*xfercoef_num_ait2acc
           end if
        end if

        ! compute accum --> aitken transfer rates
        ! accum may have some species (seasalt, dust, poa, lll) that are
        !    not in aitken mode
        ! so first divide the accum drv & num into not-transferred (noxf) species
        !    and transferred species, and use the transferred-species
        !    portion in what follows
        ixfer_acc2ait = 0
        xfercoef_num_acc2ait = 0.0_r8
        xfercoef_vol_acc2ait = 0.0_r8

        drv_t = drv_a_accsv(icol,klev) + drv_c_accsv(icol,klev)
        num_t = num_a_accsv(icol,klev) + num_c_accsv(icol,klev)
        drv_a_noxf = 0.0_r8
        drv_c_noxf = 0.0_r8
        if (drv_t > 0.0_r8) then
           if (num_t > drv_t*v2nzz) then
              do ispec = 1, nspec_amode(nacc)

                 if ( noxf_acc2ait(ispec) ) then
                    ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
                    dummwdens = 1.0_r8 / specdens_amode(lspectype_amode(ispec,nacc))
                    idx = lmassptr_amode(ispec,nacc)
                    call rad_cnst_get_aer_mmr(list_idx, nacc, ispec, 'a', state, pbuf, specmmr) !get mmr
                    drv_a_noxf = drv_a_noxf    &
                         + max(0.0_r8,specmmr(icol,klev))*dummwdens
                    !if(masterproc)then
                    !   call rad_cnst_get_info(list_idx, nacc, ispec,spec_name=spec_name)
                    !   call cnst_get_ind(trim(adjustl(spec_name)), idx_cw)
                    !   if(lmassptrcw_amode(ispec,nacc) == idx_cw) then
                    !      write(iulog,*)'BALLI, it okay:', list_idx, lmassptrcw_amode(ispec,nacc), idx_cw, ispec,nacc, spec_name
                    !   else
                    !      write(iulog,*)'BALLI, it NOTokay:', list_idx, lmassptrcw_amode(ispec,nacc), idx_cw, ispec,nacc, spec_name
                    !   endif
                    !endif
                    fldcw => qqcw_get_field(pbuf,lmassptrcw_amode(ispec,nacc),lchnk)
                    if(present(cp_buf)) then !<< remove cp_buf stuff
                       drv_c_noxf = drv_c_noxf    &
                            + max(0.0_r8,cp_buf(icol,klev,ispec,nacc))*dummwdens
                    else
                       drv_c_noxf = drv_c_noxf    &
                            + max(0.0_r8,fldcw(icol,klev))*dummwdens
                    endif
                 end if
              end do
              drv_t_noxf = drv_a_noxf + drv_c_noxf
              num_t_noxf = drv_t_noxf*voltonumblo_amode(nacc)
              num_t0 = num_t
              num_t = max( 0.0_r8, num_t - num_t_noxf )
              drv_t = max( 0.0_r8, drv_t - drv_t_noxf )
           end if
        end if

        if (drv_t > 0.0_r8) then
           if (num_t > drv_t*v2nzz) then
              ixfer_acc2ait = 1
              if (num_t > drv_t*voltonumb_amode(nait)) then
                 xferfrac_num_acc2ait = 1.0_r8
                 xferfrac_vol_acc2ait = 1.0_r8
              else
                 xferfrac_vol_acc2ait = ((num_t/drv_t) - v2nzz)/   &
                      (voltonumb_amode(nait) - v2nzz)
                 xferfrac_num_acc2ait = xferfrac_vol_acc2ait*   &
                      (drv_t*voltonumb_amode(nait)/num_t)
                 if ((xferfrac_num_acc2ait <= 0.0_r8) .or.   &
                      (xferfrac_vol_acc2ait <= 0.0_r8)) then
                    xferfrac_num_acc2ait = 0.0_r8
                    xferfrac_vol_acc2ait = 0.0_r8
                 else if ((xferfrac_num_acc2ait >= 1.0_r8) .or.   &
                      (xferfrac_vol_acc2ait >= 1.0_r8)) then
                    xferfrac_num_acc2ait = 1.0_r8
                    xferfrac_vol_acc2ait = 1.0_r8
                 end if
              end if
              duma = 1.0e-37_r8
              xferfrac_num_acc2ait = xferfrac_num_acc2ait*   &
                   num_t/max( duma, num_t0 )
              xfercoef_num_acc2ait = xferfrac_num_acc2ait*tadjinv
              xfercoef_vol_acc2ait = xferfrac_vol_acc2ait*tadjinv
              xfertend_num(2,1) = num_a_accsv(icol,klev)*xfercoef_num_acc2ait
              xfertend_num(2,2) = num_c_accsv(icol,klev)*xfercoef_num_acc2ait
           end if
        end if

        ! jump to end-of-loop if no transfer is needed at current icol,klev
        if (ixfer_ait2acc+ixfer_acc2ait > 0) then
           !
           ! compute new dgncur & v2ncur for aitken & accum modes
           !
           ! currently inactive
           do imode = nait, nacc, (nacc-nait)
              if (imode .eq. nait) then
                 duma = (xfertend_num(1,1) - xfertend_num(2,1))*deltat
                 num_a     = max( 0.0_r8, num_a_aitsv(icol,klev) - duma )
                 num_a_acc = max( 0.0_r8, num_a_accsv(icol,klev) + duma )
                 duma = (drv_a_aitsv(icol,klev)*xfercoef_vol_ait2acc -   &
                      (drv_a_accsv(icol,klev)-drv_a_noxf)*xfercoef_vol_acc2ait)*deltat
                 drv_a     = max( 0.0_r8, drv_a_aitsv(icol,klev) - duma )
                 drv_a_acc = max( 0.0_r8, drv_a_accsv(icol,klev) + duma )
                 duma = (xfertend_num(1,2) - xfertend_num(2,2))*deltat
                 num_c     = max( 0.0_r8, num_c_aitsv(icol,klev) - duma )
                 num_c_acc = max( 0.0_r8, num_c_accsv(icol,klev) + duma )
                 duma = (drv_c_aitsv(icol,klev)*xfercoef_vol_ait2acc -   &
                      (drv_c_accsv(icol,klev)-drv_c_noxf)*xfercoef_vol_acc2ait)*deltat
                 drv_c     = max( 0.0_r8, drv_c_aitsv(icol,klev) - duma )
                 drv_c_acc = max( 0.0_r8, drv_c_accsv(icol,klev) + duma )
              else
                 num_a = num_a_acc
                 drv_a = drv_a_acc
                 num_c = num_c_acc
                 drv_c = drv_c_acc
              end if

              cmn_factor = exp(4.5_r8*alnsg_amode(4)**2)*pi/6.0_r8 !THIS SHOULD BE FIXED : replac 4 with imode

              if (drv_a > 0.0_r8) then
                 if (num_a <= drv_a*voltonumbhi_amode(imode)) then
                    dgncur_a(icol,klev,imode) = dgnumhi_amode(imode)
                    v2ncur_a(icol,klev,imode) = voltonumbhi_amode(imode)
                 else if (num_a >= drv_a*voltonumblo_amode(imode)) then
                    dgncur_a(icol,klev,imode) = dgnumlo_amode(imode)
                    v2ncur_a(icol,klev,imode) = voltonumblo_amode(imode)
                 else
                    dgncur_a(icol,klev,imode) = (drv_a/(cmn_factor*num_a))**third
                    v2ncur_a(icol,klev,imode) = num_a/drv_a
                 end if
              else
                 dgncur_a(icol,klev,imode) = dgnum_amode(imode)
                 v2ncur_a(icol,klev,imode) = voltonumb_amode(imode)
              end if

              if (drv_c > 0.0_r8) then
                 if (num_c <= drv_c*voltonumbhi_amode(imode)) then
                    dgncur_c(icol,klev,imode) = dgnumhi_amode(imode)
                    v2ncur_c(icol,klev,imode) = voltonumbhi_amode(imode)
                 else if (num_c >= drv_c*voltonumblo_amode(imode)) then
                    dgncur_c(icol,klev,imode) = dgnumlo_amode(imode)
                    v2ncur_c(icol,klev,imode) = voltonumblo_amode(imode)
                 else
                    dgncur_c(icol,klev,imode) = (drv_c/(cmn_factor*num_c))**third
                    v2ncur_c(icol,klev,imode) = num_c/drv_c
                 end if
              else
                 dgncur_c(icol,klev,imode) = dgnum_amode(imode)
                 v2ncur_c(icol,klev,imode) = voltonumb_amode(imode)
              end if

           end do


           !
           ! compute tendency amounts for aitken <--> accum transfer
           !

           if ( masterproc ) then
              if (idiagaa > 0) then
                 do jmode = 1, 2
                    do iq = 1, nspecfrm_csizxf
                       do aer_type = 1, 2
                          if (jmode .eq. 1) then
                             if (aer_type .eq. inter_aero) then
                                lsfrm = lspecfrma_csizxf(iq)
                                lstoo = lspectooa_csizxf(iq)
                             else
                                lsfrm = lspecfrmc_csizxf(iq)
                                lstoo = lspectooc_csizxf(iq)
                             end if
                          else
                             if (aer_type .eq. inter_aero) then
                                lsfrm = lspectooa_csizxf(iq)
                                lstoo = lspecfrma_csizxf(iq)
                             else
                                lsfrm = lspectooc_csizxf(iq)
                                lstoo = lspecfrmc_csizxf(iq)
                             end if
                          end if
                          write( iulog, '(a,3i3,2i4)' ) 'calcsize jmode,iq,aer_type, lsfrm,lstoo',   &
                               jmode,iq,aer_type, lsfrm,lstoo
                       end do
                    end do
                 end do
              end if
           end if
           idiagaa = -1


           ! jmode=1 does aitken-->accum; jmode=2 does accum-->aitken
           do  jmode = 1, 2

              if ((jmode .eq. 1 .and. ixfer_ait2acc > 0) .or. &
                   (jmode .eq. 2 .and. ixfer_acc2ait > 0)) then

                 jsrflx = jmode+2
                 if (jmode .eq. 1) then
                    xfercoef = xfercoef_vol_ait2acc
                 else
                    xfercoef = xfercoef_vol_acc2ait
                 end if

                 do  iq = 1, nspecfrm_csizxf

                    ! aer_type=1 does interstitial ("_a"); aer_type=2 does activated ("_c");
                    do  aer_type = 1, 2

                       ! the lspecfrma_csizxf (and lspecfrmc_csizxf) are aitken species
                       ! the lspectooa_csizxf (and lspectooc_csizxf) are accum  species
                       ! for jmode=1, want lsfrm=aitken species, lstoo=accum  species
                       ! for jmode=2, want lsfrm=accum  species,  lstoo=aitken species
                       if (jmode .eq. 1) then
                          if (aer_type .eq. inter_aero) then
                             lsfrm = lspecfrma_csizxf(iq)
                             lstoo = lspectooa_csizxf(iq)
                          else
                             lsfrm = lspecfrmc_csizxf(iq)
                             lstoo = lspectooc_csizxf(iq)
                          end if
                       else
                          if (aer_type .eq. inter_aero) then
                             lsfrm = lspectooa_csizxf(iq)
                             lstoo = lspecfrma_csizxf(iq)
                          else
                             lsfrm = lspectooc_csizxf(iq)
                             lstoo = lspecfrmc_csizxf(iq)
                          end if
                       end if

                       if ((lsfrm > 0) .and. (lstoo > 0)) then
                          if (aer_type .eq. inter_aero) then
                             if (iq .eq. 1) then
                                xfertend = xfertend_num(jmode,aer_type)
                             else
                                xfertend = max(0.0_r8,state_q(icol,klev,lsfrm))*xfercoef
                             end if
                             if(update_mmr)dqdt(icol,klev,lsfrm) = dqdt(icol,klev,lsfrm) - xfertend
                             if(update_mmr)dqdt(icol,klev,lstoo) = dqdt(icol,klev,lstoo) + xfertend
                          else
                             if (iq .eq. 1) then
                                xfertend = xfertend_num(jmode,aer_type)
                             else
                                fldcw => qqcw_get_field(pbuf,lsfrm,lchnk)
                                if(present(cp_buf)) then
                                   n_use = huge(1)
                                   m_use = huge(1)
                                   do nb = 1,4
                                      do mb = 1,nspec_amode(nb)
                                         if(lsfrm == lmassptr_amode(mb,nb)) then
                                            n_use = nb
                                            m_use = mb
                                         endif
                                      enddo
                                   enddo
                                   xfertend = max(0.0_r8,cp_buf(icol,klev,m_use,n_use))*xfercoef
                                else
                                   xfertend = max(0.0_r8,fldcw(icol,klev))*xfercoef
                                endif
                             end if
                             dqqcwdt(icol,klev,lsfrm) = dqqcwdt(icol,klev,lsfrm) - xfertend
                             dqqcwdt(icol,klev,lstoo) = dqqcwdt(icol,klev,lstoo) + xfertend
                          end if
                          qsrflx(icol,lsfrm,jsrflx,aer_type) = qsrflx(icol,lsfrm,jsrflx,aer_type) - xfertend*pdel_fac
                          qsrflx(icol,lstoo,jsrflx,aer_type) = qsrflx(icol,lstoo,jsrflx,aer_type) + xfertend*pdel_fac
                       end if

                    end do
                 end do
              end if
           end do

        end if
     end do
  end do

  return
end subroutine aitken_accum_exchange

!---------------------------------------------------------------------------------------------

subroutine update_cld_brn_mmr(top_lev, pver, ncol, lchnk, pcnst, deltat, pbuf, dotendqqcw, dqqcwdt)

  !-----------------------------------------------------------------------------
  !Purpose: updates mmr of cloud borne aerosols
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : None
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------

  implicit none

  !input
  integer,  intent(in) :: top_lev, pver !for model level loop
  integer,  intent(in) :: ncol          !# of columns
  integer,  intent(in) :: lchnk
  integer,  intent(in) :: pcnst         !# of constituents
  real(r8), intent(in) :: deltat        !time step

  logical,  intent(in) :: dotendqqcw(:)

  type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer

  !output
  real(r8), intent(inout) :: dqqcwdt(:,:,:)

  !local
  integer :: icnst, lc, klev, icol
  real(r8), pointer :: fldcw(:,:)    !specie mmr (cloud borne)

  do icnst = 1, pcnst
     lc = icnst
     if ( lc>0 .and. dotendqqcw(lc) ) then
        fldcw=> qqcw_get_field(pbuf,icnst,lchnk)
        do klev = top_lev, pver
           do icol = 1, ncol
              fldcw(icol,klev) = max( 0.0_r8, (fldcw(icol,klev) + dqqcwdt(icol,klev,lc)*deltat) )
           end do
        end do
     end if
  end do

  return
end subroutine update_cld_brn_mmr

!---------------------------------------------------------------------------------------------

subroutine output_flds(name, idx, lchnk, aer_type, qsrflx)

  !-----------------------------------------------------------------------------
  !Purpose: Output model fields in the history file
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : outfld
  !
  !Author: Balwinder Singh based on the code by Richard Easter (RCE)
  !-----------------------------------------------------------------------------
  implicit none

  !inputs
  character(len=*), intent(in) :: name

  integer,  intent(in) :: idx, lchnk, aer_type
  real(r8), intent(in) :: qsrflx(:,:,:,:)

  !local
  character(len=fieldname_len) :: fieldname

  fieldname = trim(name) // '_sfcsiz1'
  call outfld( fieldname, qsrflx(:,idx,1,aer_type), pcols, lchnk)

  fieldname = trim(name) // '_sfcsiz2'
  call outfld( fieldname, qsrflx(:,idx,2,aer_type), pcols, lchnk)

  return
end subroutine output_flds

!---------------------------------------------------------------------------------------------

subroutine modal_aero_calcsize_diag(state, pbuf, list_idx_in, dgnum_m)

   !-----------------------------------------------------------------------
   !
   ! Calculate aerosol size distribution parameters
   !
   ! ***N.B.*** DGNUM for the modes in the climate list are put directly into
   !            the physics buffer.  For diagnostic list calculations use the
   !            optional list_idx and dgnum args.
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), intent(in), target :: state   ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)      ! physics buffer

   integer,  optional, intent(in)   :: list_idx_in    ! diagnostic list index
   real(r8), optional, pointer      :: dgnum_m(:,:,:) ! interstital aerosol dry number mode radius (m)

   ! local
   integer  :: i, k, l1, n
   integer  :: lchnk, ncol
   integer  :: list_idx, stat
   integer  :: nmodes
   integer  :: nspec

   real(r8), pointer :: dgncur_a(:,:) ! (pcols,pver)

   real(r8), pointer :: mode_num(:,:) ! mode number mixing ratio
   real(r8), pointer :: specmmr(:,:)  ! specie mmr
   real(r8)          :: specdens      ! specie density

   real(r8) :: dryvol_a(pcols,pver)   ! interstital aerosol dry volume (cm^3/mol_air)

   real(r8) :: dgnum, dgnumhi, dgnumlo
   real(r8) :: dgnyy, dgnxx           ! dgnumlo/hi of current mode
   real(r8) :: drv_a                  ! dry volume (cm3/mol_air)
   real(r8) :: dumfac, dummwdens      ! work variables
   real(r8) :: num_a0                 ! initial number (#/mol_air)
   real(r8) :: num_a                  ! final number (#/mol_air)
   real(r8) :: voltonumbhi, voltonumblo
   real(r8) :: v2nyy, v2nxx           ! voltonumblo/hi of current mode
   real(r8) :: sigmag, alnsg
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   list_idx = 0  ! climate list by default
   if (present(list_idx_in)) list_idx = list_idx_in

   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx /= 0) then
      if (.not. present(dgnum_m)) then
         call endrun('modal_aero_calcsize_diag called for'// &
                     'diagnostic list but dgnum_m pointer not present')
      end if
      allocate(dgnum_m(pcols,pver,nmodes), stat=stat)
      if (stat > 0) then
         call endrun('modal_aero_calcsize_diag: allocation FAILURE: dgnum_m')
      end if
   end if

   do n = 1, nmodes

      if (list_idx == 0) then
         call pbuf_get_field(pbuf, dgnum_idx, dgncur_a, start=(/1,1,n/), kount=(/pcols,pver,1/))
      else
         dgncur_a => dgnum_m(:,:,n)
      end if

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, n, dgnum=dgnum, dgnumhi=dgnumhi, dgnumlo=dgnumlo, &
                                   sigmag=sigmag)

      ! get mode number mixing ratio
      call rad_cnst_get_mode_num(list_idx, n, 'a', state, pbuf, mode_num)

      dgncur_a(:,:) = dgnum
      dryvol_a(:,:) = 0.0_r8

      ! compute dry volume mixrats =
      !      sum_over_components{ component_mass mixrat / density }
      call rad_cnst_get_info(list_idx, n, nspec=nspec)
      do l1 = 1, nspec

         call rad_cnst_get_aer_mmr(list_idx, n, l1, 'a', state, pbuf, specmmr)
         call rad_cnst_get_aer_props(list_idx, n, l1, density_aer=specdens)

         ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
         dummwdens = 1.0_r8 / specdens

         do k=top_lev,pver
            do i=1,ncol
               dryvol_a(i,k) = dryvol_a(i,k)    &
                  + max(0.0_r8, specmmr(i,k))*dummwdens
            end do
         end do
      end do

      alnsg  = log( sigmag )
      dumfac = exp(4.5_r8*alnsg**2)*pi/6.0_r8
      voltonumblo = 1._r8 / ( (pi/6._r8)*(dgnumlo**3)*exp(4.5_r8*alnsg**2) )
      voltonumbhi = 1._r8 / ( (pi/6._r8)*(dgnumhi**3)*exp(4.5_r8*alnsg**2) )
      v2nxx = voltonumbhi
      v2nyy = voltonumblo
      dgnxx = dgnumhi
      dgnyy = dgnumlo

      do k = top_lev, pver
         do i = 1, ncol

            drv_a = dryvol_a(i,k)
            num_a0 = mode_num(i,k)
            num_a = max( 0.0_r8, num_a0 )

            if (drv_a > 0.0_r8) then
               if (num_a <= drv_a*v2nxx) then
                  dgncur_a(i,k) = dgnxx
               else if (num_a >= drv_a*v2nyy) then
                  dgncur_a(i,k) = dgnyy
               else
                  dgncur_a(i,k) = (drv_a/(dumfac*num_a))**third
               end if
            end if

         end do
      end do

   end do ! nmodes

end subroutine modal_aero_calcsize_diag

!----------------------------------------------------------------------

end module modal_aero_calcsize
