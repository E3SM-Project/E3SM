module modal_aero_calcsize

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8, cs => shr_kind_cs
use spmd_utils,       only: masterproc
use physconst,        only: pi, gravit

use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state, physics_ptend
use physics_buffer,   only: physics_buffer_desc, pbuf_get_field

use phys_control,     only: phys_getopts
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num, N_DIAG, rad_cnst_get_mode_idx

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use shr_log_mod ,     only: errMsg => shr_log_errMsg
use cam_history,      only: addfld, horiz_only, add_default, fieldname_len, outfld
use constituents,     only: pcnst, cnst_name, cnst_get_ind

use ref_pres,         only: top_lev => clim_modal_aero_top_lev

#ifdef MODAL_AERO

use modal_aero_data, only: ntot_amode, numptr_amode, modename_amode
use modal_aero_data,  only: numptrcw_amode, mprognum_amode, qqcw_get_field, &
     modeptr_accum, modeptr_aitken, ntot_aspectype, cnst_name_cw

#endif
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!BSINGH: THINGS TO DO (TODO):
!1. Find out if we can remove the "if (masterproc)" from the subroutines called during time stepping
!2. cloud borne and interstitial species processes can be combined and refactored into one subroutine.(aitken_accum_exchange)
!3. init and register routine should be re-worked/streamlined to remove redundant logics
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------

implicit none
private
save

public :: modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
public :: modal_aero_calcsize_reg
!Mimic enumerators for aerosol types
integer, parameter:: inter_aero   = 1 !interstitial aerosols
integer, parameter:: cld_brn_aero = 2 !cloud borne species

integer :: dgnum_idx = -1 !pbuf id for dgnum

integer, parameter :: maxpair_csizxf = N_DIAG
#ifdef MODAL_AERO
integer, parameter :: maxspec_csizxf = ntot_aspectype
#else
! TODO: this is a kludge.  This value should probably be assigned
! elsewhere for the non-modal case.  S.M. Burrows.
integer, parameter :: maxspec_csizxf = 8
#endif

real(r8), parameter :: third   = 1.0_r8/3.0_r8
real(r8), parameter :: r8_huge = huge(1.0_r8)

!---------------------------------------------------------------------------------
!Note: For diagnostic calls, users can ask for any diagnostics like rad_diag_1
!      and rad_diag_3 (e.g. skipping rad_diag_2). Therefore the arrays should
!      have a length of N_DIAG (unless we define another array which maps info
!      such as to include info about the missing diagnostics)
!---------------------------------------------------------------------------------
integer, parameter :: npair_csizxf = N_DIAG           !total number of possible diagnostic calls

logical :: do_adjust_allowed                          !flag to turn on/off  aerosol size adjustment process
!flag to turn on/off  aerosol aitken<->accumulation transfer process
logical :: do_aitacc_transfer_allowed(0:npair_csizxf) ! This is an array as it can be different for each radiatio diagnostic call

!--------------------------------------------------------------------------------
!Naming convention:
!"*_csizxf": variables for transfering ("xf") species from one mode to another
!"*frm*"   : "From" mode, the mode from where the species are coming from
!*to*      : "To" mode, the mode which is receiving the species
!
! All ararys will start with a zero index as 0 index is reserved for prognostic
! radiation calls
!--------------------------------------------------------------------------------

!"modefrm_csizxf" stores mode number "from" which species will be moved
!to a mode stored in "modetoo_csizxf".
![E.g. if modefrm_csizxf(3)=2 and modetoo_csizxf(3)=1, for rad_diag_3 (notice that
!these arrays are indexed "3" as rad_diag is rad_diag_3), species will be moved
!from 2nd mode to the 1st mode.]
integer :: modefrm_csizxf(0:maxpair_csizxf)
integer :: modetoo_csizxf(0:maxpair_csizxf)

!Total number of species in the "from" mode.
![ E.g. for rad_diag_3, if nspecfrm_csizxf(3)=7,
!it means 7 species exist in modefrm_csizxf(3).
!In this case, if modefrm_csizxf(3)=2,
!it means 2nd mode has 7 species]
integer :: nspecfrm_csizxf(0:maxpair_csizxf)

!"lspec*" are the indices in constituent(cnst_name) and state Q array
! (state%q) arrays.
!These arrays hold species index for "from" and "to" for each pair [to-from pair
!(from-to pair is same as to-from) for transfering species between modes]

![E.g. for rad_diag_3:
!"lspecfrm*_csizxf" will hold index of a specie of the mode modefrm_csizxf(3)
!"lspecto*_csizxf" will hold index of a specie of the mode modetoo_csizxf(3)
!"*a_*" and "*c_* refers to interstitial and cloud borne species
integer :: lspecfrmc_csizxf(maxspec_csizxf, 0:maxpair_csizxf)
integer :: lspecfrma_csizxf(maxspec_csizxf, 0:maxpair_csizxf)
integer :: lspectooc_csizxf(maxspec_csizxf, 0:maxpair_csizxf)
integer :: lspectooa_csizxf(maxspec_csizxf, 0:maxpair_csizxf)

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_calcsize_reg()
  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

  integer :: nmodes

  !find out number of modes (nmodes) in radiation climate list
  call rad_cnst_get_info(0, nmodes=nmodes)

  !register dgnum field
  call pbuf_add_field('DGNUM', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnum_idx)
end subroutine modal_aero_calcsize_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_init( pbuf2d, species_class)
   use time_manager,   only: is_first_step
   use physics_buffer, only: pbuf_set_field

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   !    set do_adjust_allowed and do_aitacc_transfer_allowed flags
   !    create history fields for column tendencies associated with
   !       modal_aero_calcsize
   !
   ! Author: R. Easter (refactored by Balwinder Singh)
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer, intent(in) :: species_class(:)

   ! local
   integer  :: iq, ilist, iacc, iait, imode, icnst
   integer  :: nspec_ait, nspec_acc, imode_ait, imode_acc, ispec_acc, ispec_ait
   integer  :: aer_type, ind_ait, ind_acc, icnt
   integer  :: lsfrm, lstoo
   integer  :: nacc, nait
   logical  :: history_aerosol
   logical  :: history_verbose
   logical  :: accum_exists, aitken_exists

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
   character(len=32)              :: spec_name_ait, spec_name_acc
   character(len=32)              :: spec_name_ait_cw, spec_name_acc_cw
   character(len=32)              :: spec_cnst_name_acc, spec_cnst_name_ait
   character(len=32)              :: num_name_acc, num_name_ait
   character(len=32)              :: num_name_acc_cw, num_name_ait_cw
   character(len=32)              :: mass_name_ait_cw, mass_name_acc_cw
   !-----------------------------------------------------------------------

   call phys_getopts(history_aerosol_out=history_aerosol, &
                     history_verbose_out=history_verbose)

   ! init entities required for both prescribed and prognostic modes

   if (is_first_step()) then ! bsingh- Do we need this conditional in a init routine?
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnum_idx, 0.0_r8)
   endif

   !initialize
   !(0 index 0f modefrm_csizxf and  modetoo_csizxf is reserved for the prognostic call)
   do ilist = 0, npair_csizxf
      modefrm_csizxf(ilist) = -1
      modetoo_csizxf(ilist) = -1
   enddo

#ifndef MODAL_AERO
   do_adjust_allowed          = .false.
   do_aitacc_transfer_allowed = .false.
#else
   !  do_adjust_allowed allows aerosol size  adjustment to be turned on/off
   do_adjust_allowed = .true.

   !------------------------------------------------------------------------------------------
   !  "do_aitacc_transfer_allowed" allows aitken <--> accum mode transfer to be turned on/off
   !  NOTE: it can only be true when aitken & accum modes are both present
   !      and have prognosed number and diagnosed surface/sigmag
   !------------------------------------------------------------------------------------------

   !find out mode index for accum and aitken modes in the prognostic radiation list (rad_climate)
   !These are used to get strings stored for these modes in modename_amode array
   nait = modeptr_aitken !mode number of accumulation mode
   nacc = modeptr_accum  !mode number of aitken mode

   !Go through the radiation list and decide on do_aitacc_transfer_allowed value(true/false)
   do ilist = 0, npair_csizxf

      !find out accumulation and aitken modes in the radiation list
      iacc = rad_cnst_get_mode_idx(ilist, modename_amode(nacc))
      iait = rad_cnst_get_mode_idx(ilist, modename_amode(nait))

      !find out if aitken or accumulation modes exist in the radiation list
      !(a positive value means that the mode exists)
      accum_exists  = ( iacc > 0)
      aitken_exists = ( iait > 0)

      do_aitacc_transfer_allowed(ilist)=.false.
      if(accum_exists .and. aitken_exists .and. iacc .ne. iait) then
         do_aitacc_transfer_allowed(ilist)=.true.
         modefrm_csizxf(ilist) = iait
         modetoo_csizxf(ilist) = nacc
      endif
   enddo

   !if num mixing ratio is set *not* to be prognosed, turn aitken<->accumm transfer off for all radiation lists
   if (mprognum_amode(nait) <= 0) do_aitacc_transfer_allowed(:) = .false.
   if (mprognum_amode(nacc) <= 0) do_aitacc_transfer_allowed(:) = .false.

   !-------------------------------------------------------------------------------
   !Find mapping between the corresponding species of aitken and accumulation nodes
   !
   ! For example, soa may exist in accumulation and aitken modes. Here were are
   ! finding mapping between their indices so that if we need to move soa from
   ! accumulation to aitken, we can do that using this mapping
   !-------------------------------------------------------------------------------


   ! Only find the mapping if atleast one of the do_aitacc_transfer_allowed is true
   if (any(do_aitacc_transfer_allowed(:))) then

      !Go through the radiation list and find "aitken<-->accumulation" mapping for each list
      do ilist = 0, npair_csizxf
         if(do_aitacc_transfer_allowed(ilist)) then
            icnt = 0
            imode_ait = modefrm_csizxf(ilist) !aitken  mode of this list
            imode_acc = modetoo_csizxf(ilist) !accumulation mode of this list

            !--------------------------------------------------------------------------------------
            !find aerosol *number* indices mapping between aitken and accumulation modes in this list
            !--------------------------------------------------------------------------------------

            !For aitken mode
            call rad_cnst_get_info(ilist, imode_ait, num_name=num_name_ait)
            call cnst_get_ind(num_name_ait, ind_ait)

            !For accumulation mode
            call rad_cnst_get_info(ilist, imode_acc, num_name=num_name_acc)
            call cnst_get_ind(num_name_acc, ind_acc)

            icnt = icnt + 1
            lspecfrma_csizxf(icnt,ilist) = ind_ait
            lspectooa_csizxf(icnt,ilist) = ind_acc

            !index for cloudborne species (qqcw_get_field array) is same as interstitial species
            lspecfrmc_csizxf(icnt,ilist) = ind_ait
            lspectooc_csizxf(icnt,ilist) = ind_acc

            !--------------------------------------------------------------------------------------
            !find aerosol *mass* indices mapping between aitken and accumulation modes in this list
            !--------------------------------------------------------------------------------------
            !find number of species in the aitken mode of this ilist
            call rad_cnst_get_info(ilist, imode_ait, nspec = nspec_ait) !output:nspec_ait
            !find number of species in the accumulation mode of this list
            call rad_cnst_get_info(ilist, imode_acc, nspec = nspec_acc) !output:nspec_acc

            !loop through speices of aitken mode and find corresponding species in the accumulation mode
            do ispec_ait = 1, nspec_ait
               !find specie name
               call rad_cnst_get_info(ilist, imode_ait, ispec_ait,spec_name=spec_name_ait) !output:spec_name_ait

               !extract constituent name from the specie name (e.g for soa_a1, constituent name would be "soa")
               spec_cnst_name_ait  = extract_cnst_name(spec_name_ait)

               !Now find this specie index in the "to" mode

               !loop through speices of the "to" mode
               do ispec_acc = 1, nspec_acc
                  !find species name
                  call rad_cnst_get_info(ilist, imode_acc, ispec_acc,spec_name=spec_name_acc) !output:spec_name_acc

                  !extract constituent name from the specie name (e.g for soa_a1, constituent name would be "soa")
                  spec_cnst_name_acc  = extract_cnst_name(spec_name_acc)

                  !find if specie in acc mode is same as ait or not
                  if(trim(spec_cnst_name_acc) == trim(spec_cnst_name_ait)) then
                     !if there is a match, find indices of species in cnst array
                     call cnst_get_ind(spec_name_ait, ind_ait)
                     call cnst_get_ind(spec_name_acc, ind_acc)

                     icnt = icnt + 1
                     lspecfrma_csizxf(icnt,ilist) = ind_ait
                     lspectooa_csizxf(icnt,ilist) = ind_acc

                     !index for cloudborne species (qqcw_get_field array) is same as interstitial species
                     lspecfrmc_csizxf(icnt,ilist) = ind_ait
                     lspectooc_csizxf(icnt,ilist) = ind_acc
                  endif
               enddo !ispec_acc
            enddo!ispec_ait
            nspecfrm_csizxf(ilist) = icnt !count total number of matching species in each list
         endif!do_aitacc_transfer_allowed
      enddo!ilist
   endif !any(do_aitacc_transfer_allowed)

   !--------------------------------------------------------------------------------
   ! Define history fields for number-adjust source-sink for all modes
   ! NOTE: This is only done for prognostic radiation list (ilist = 0)
   !
   ! BSINGH: Did not modify the following code, it should be refactored
   ! by removing the do-loops
   !--------------------------------------------------------------------------------
9310  format( / 'subr. modal_aero_calcsize_init' / 'do_adjust_allowed, do_aitacc_transfer_allowed = ', 2l10 )
9320  format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330  format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340  format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350  format( 5x, 'no corresponding activated species' )


!  define history fields for number-adjust source-sink for all modes
do_adjust_if_block2: &
      if ( do_adjust_allowed ) then

   do imode = 1, ntot_amode
      if (mprognum_amode(imode) <= 0) cycle

      do aer_type = 1, 2
         if (aer_type == inter_aero) then
            tmpnamea = cnst_name(numptr_amode(imode))
         else
            tmpnamea = cnst_name_cw(numptrcw_amode(imode))
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


!  define history fields for aitken-accum transfer [Only for list_idx = 0 (radiation climate list)]
ilist = 0
do_aitacc_transfer_if_block2: &
      if ( do_aitacc_transfer_allowed(ilist) ) then

! check that calcsize transfer ilist=0 ("radiation_climate" list) is aitken-->accum
      if ((modefrm_csizxf(ilist) .ne. nait) .or.   &
          (modetoo_csizxf(ilist) .ne. nacc)) then
         write( iulog, '(//2a//)' ) '*** modal_aero_calcaersize_init error -- modefrm/too_csizxf(1) are wrong'
         call endrun( 'modal_aero_calcaersize_init error' )
      end if

      do iq = 1, nspecfrm_csizxf(ilist)

! aer_type=1 does interstitial ("_a"); aer_type=2 does activated ("_c");
         do aer_type = 1, 2

! the lspecfrma_csizxf (and lspecfrmc_csizxf) are aitken species
! the lspectooa_csizxf (and lspectooc_csizxf) are accum  species
            if (aer_type .eq. inter_aero) then
               lsfrm = lspecfrma_csizxf(iq,ilist)
               lstoo = lspectooa_csizxf(iq,ilist)
            else
               lsfrm = lspecfrmc_csizxf(iq,ilist)
               lstoo = lspectooc_csizxf(iq,ilist)
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
         do icnst = 1, pcnst
            write(iulog,'(2i4,2x,a)') icnst, species_class(icnst), cnst_name(icnst)
         end do
      end if
   if ( masterproc ) write(iulog,'(a)') 'modal_aero_calcsize_init ALL DONE'

#endif

end subroutine modal_aero_calcsize_init


function extract_cnst_name(spec_name) result(spec_cnst_name)

  !------------------------------------------------------------------------
  !extract constituent name from the specie name (e.g for soa_a1,
  !constituent name would be "soa")
  !------------------------------------------------------------------------
  character(len=*), intent(in) :: spec_name

  !function output
  character(len=cs)  :: spec_cnst_name

  !local variables
  integer            :: ind
  character(len=200) :: err_msg

  !ASSUMPTION: species name will have an underscore("_")

  ind = index(trim(adjustl(spec_name)),'_')
  if(ind >= 1 ) then ! if "_" exists
     spec_cnst_name = spec_name(1:ind-1)
  else
     err_msg = "Species name should have an underscore in it, "// &
          " species name is:"//trim(spec_name)//" "//errmsg(__FILE__,__LINE__)
     call endrun(err_msg)
  endif

end function extract_cnst_name

!===============================================================================

subroutine modal_aero_calcsize_sub(state, deltat, pbuf, ptend, do_adjust_in, &
   do_aitacc_transfer_in, list_idx_in, update_mmr_in, dgnumdry_m, caller)

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
   ! 09/2020: Refacotred by Balwinder Singh
   !
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), target, intent(in)    :: state       ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer
   real(r8),                    intent(in)    :: deltat      ! model time-step size (s)
   type(physics_ptend), target, optional, intent(inout) :: ptend       ! indivdual parameterization tendencies

   logical,  optional, intent(in) :: do_adjust_in
   logical,  optional, intent(in) :: do_aitacc_transfer_in
   logical,  optional, intent(in) :: update_mmr_in
   integer,  optional, intent(in) :: list_idx_in       ! diagnostic list index
   real(r8), optional, intent(inout), allocatable, target :: dgnumdry_m(:,:,:) ! interstital aerosol dry number mode radius (m)

   !This subroutine is called from various places in the code
   !"caller" optional variable can hold the name of the subroutine
   !which called this subroutine (for debugging only)
   character(len=*), optional     :: caller

#ifdef MODAL_AERO

   ! local

   logical :: do_adjust
   logical :: do_aitacc_transfer
   logical :: update_mmr

   integer  :: lchnk                ! chunk identifier
   integer  :: ncol                 ! number of columns
   integer  :: list_idx_local       !list idx local to this subroutine
   integer  :: ilist

   real(r8), parameter :: close_to_one = 1.0_r8 + 1.0e-15_r8
   real(r8), parameter :: seconds_in_a_day = 86400.0_r8

   real(r8), pointer :: state_q(:,:,:), dqdt(:,:,:)
   real(r8), pointer :: pdel(:,:)   ! pressure thickness of levels

   real(r8), pointer :: dgncur_a(:,:,:)

   real(r8), target :: fake_dqdt(1,1,1)

   integer  :: nspec, imode, klev, icol, idx, iq
   integer  :: nmodes, num_idx
   integer  :: num_mode_idx, num_cldbrn_mode_idx, lsfrm, lstoo
   integer  :: stat

   logical, pointer  :: dotend(:)
   logical, target   :: fake_dotend(1)
   logical  :: dotendqqcw(pcnst)

   character(len=fieldname_len)   :: tmpnamea, tmpnameb, name
   character(len=fieldname_len+3) :: fieldname
   character(len = 1000) :: err_msg
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
   !    3=aitken-->accum transfer
   !    4=accum -->aitken transfer
   ! 4th index --
   !    1="a" species (interstitial); 2="c" species(cloud borne)
   !-----------------------------------------------------------------------
   integer, parameter :: nsrflx = 4   ! last dimension of qsrflx
   real(r8) :: qsrflx(pcols,pcnst,nsrflx,2)

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

   !Set list_idx_local and other sanity checks
   list_idx_local  = 0
   if(present(list_idx_in)) then
      list_idx_local  = list_idx_in
      if (.not. present(dgnumdry_m)) &
           call endrun('list_idx_in is present but dgnumdry_m is missing'//errmsg(__FILE__,__LINE__))

      if(.not. allocated(dgnumdry_m)) &
           call endrun('list_idx_in is present but dgnumdry_m is not allocated'//errmsg(__FILE__,__LINE__))
      dgncur_a => dgnumdry_m(:,:,:)

      !update_mmr should be true ONLY for list_idx = 0
      if(list_idx_local > 0 .and. update_mmr) then
         call endrun('update_mmr should be FLASE for list_idx>0'//errmsg(__FILE__,__LINE__))
      endif
   else
      call pbuf_get_field(pbuf, dgnum_idx, dgncur_a)
   endif


   !For adjusting aerosol sizes
   do_adjust = do_adjust_allowed
   if (present(do_adjust_in)) then
      if(do_adjust_in .and. .not.do_adjust_allowed) then
         write(err_msg,*)'Cannot adjust aerosol size if size adjustment is not allowed', &
              ' do_adjust_allowed is:',do_adjust_allowed,' and do_adjust_in is:',do_adjust_in, &
              ';',errmsg(__FILE__,__LINE__)
         call endrun(trim(err_msg))
      endif
      do_adjust = do_adjust_in
   endif

   !For transerfering aerosols between modes (Aitkem<-->Accumulation) based on their new size
   do_aitacc_transfer = do_aitacc_transfer_allowed(list_idx_local)
   if (present(do_aitacc_transfer_in)) then
      if(do_aitacc_transfer_in .and. .not.do_aitacc_transfer_allowed(list_idx_local)) then
         write(err_msg,*)'Cannot transfer species between modes if transfer is not allowed for radiation', &
              ' list:', list_idx_local,' do_aitacc_transfer_allowed is:',do_aitacc_transfer_allowed(list_idx_local), &
              ' and do_aitacc_transfer_in is:',do_aitacc_transfer_in,';', errmsg(__FILE__,__LINE__)
         call endrun(trim(err_msg))
      endif
      do_aitacc_transfer = do_aitacc_transfer_in
   endif

   !Set tendency variable based on the presence of ptend
   if(update_mmr) then
      if(.not. present(ptend)) then
         write(err_msg,*)'ptend should be present if update_mmr is true (update_mmr ', &
              ' is true by default); '//errmsg(__FILE__,__LINE__)
         call endrun(trim(err_msg))
      endif
      dotend => ptend%lq
      dqdt   => ptend%q
      dotendqqcw(:)   = .false.
      dqqcwdt(:,:,:)  = 0.0_r8
      qsrflx(:,:,:,:) = 0.0_r8
   else
      dotend => fake_dotend
      dqdt   => fake_dqdt
      dqqcwdt(:,:,:)  = r8_huge
      qsrflx(:,:,:,:) = r8_huge
   endif

   !if(present(caller) .and. masterproc) write(iulog,*)'modal_aero_calcsize_sub has been called by ', trim(caller)

   pdel     => state%pdel !Only required if update_mmr = .true.
   state_q  => state%q    !BSINGH - it is okay to use state_q for num mmr but not for specie mmr (as diagnostic call may miss some species)

   !----------------------------------------------------------------------------
   ! tadj = adjustment time scale for number, surface when they are prognosed
   !----------------------------------------------------------------------------
   tadj    = max( seconds_in_a_day, deltat )
   tadjinv = 1.0_r8/(tadj*close_to_one)
   fracadj = max( 0.0_r8, min( 1.0_r8, deltat*tadjinv ) )

   !grid parameters
   ncol  = state%ncol  !# of columns
   lchnk = state%lchnk !chunk #

   !inverse of time step
   deltatinv = 1.0_r8/(deltat*close_to_one)

   call rad_cnst_get_info(list_idx_local, nmodes=nmodes)

   !Now compute dry diameter for both interstitial and cloud borne aerosols
   do imode = 1, nmodes

      !----------------------------------------------------------------------
      !Initialize all parameters to the default values for the mode
      !----------------------------------------------------------------------
      !interstitial
      call set_initial_sz_and_volumes(list_idx_local, top_lev, ncol, imode, & !input
           dgncur_a, v2ncur_a, dryvol_a)                                      !output

      !cloud borne
      call set_initial_sz_and_volumes(list_idx_local, top_lev, ncol, imode, & !input
           dgncur_c, v2ncur_c, dryvol_c)                                      !output

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
      call compute_dry_volume(top_lev, ncol, imode, nspec, state, pbuf, dryvol_a, dryvol_c, list_idx_local)


      ! do size adjustment based on computed dry diameter values and update the diameters
      call size_adjustment(list_idx_local, top_lev, ncol, lchnk, imode, dryvol_a, state_q, & !input
           dryvol_c, pdel, do_adjust, update_mmr, do_aitacc_transfer, deltatinv, fracadj, pbuf,  & !input
           dgncur_a, dgncur_c, v2ncur_a, v2ncur_c,                                               & !output
           drv_a_accsv, drv_c_accsv, drv_a_aitsv, drv_c_aitsv, drv_a_sv, drv_c_sv,               & !output
           num_a_accsv, num_c_accsv, num_a_aitsv, num_c_aitsv, num_a_sv, num_c_sv,               & !output
           dotend, dotendqqcw, dqdt, dqqcwdt, qsrflx)                                              !output

   end do  ! do imode = 1, ntot_amode


   !------------------------------------------------------------------------------
   ! when the aitken mode mean size is too big, the largest
   !    aitken particles are transferred into the accum mode
   !    to reduce the aitken mode mean size
   ! when the accum mode mean size is too small, the smallest
   !    accum particles are transferred into the aitken mode
   !    to increase the accum mode mean size
   !------------------------------------------------------------------------------

   if ( do_aitacc_transfer ) then
      call aitken_accum_exchange(ncol, lchnk, list_idx_local, update_mmr, tadjinv, &
           deltat, pdel, state_q, state, &
           pbuf, &
           drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv, &
           drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,  &
           dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
           dqdt, dqqcwdt, qsrflx)
   end if

   lsfrm = -huge(lsfrm) !initialize


   !----------------------------------------------------------------------
   ! apply tendencies to cloud-borne species MRs
   ! Only if "update_mmr" is TRUE and only for list_idx_local=0
   !----------------------------------------------------------------------
   if(update_mmr .and. list_idx_local == 0) then

      !update cld brn aerosols
      call update_cld_brn_mmr(list_idx_local, top_lev, ncol, lchnk, pcnst, deltat, pbuf, dotendqqcw, dqqcwdt)

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
      if ( .not. do_aitacc_transfer) return ! return if transfer of species not required
      ilist = 0 ! write output only for list_idx_local = 0
      do iq = 1, nspecfrm_csizxf(ilist)

         ! aer_type=1 does interstitial ("_a"); aer_type=2 does activated ("_c");

         lsfrm = lspecfrma_csizxf(iq,ilist)
         lstoo = lspectooa_csizxf(iq,ilist)

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

         lsfrm = lspecfrmc_csizxf(iq,ilist)
         lstoo = lspectooc_csizxf(iq,ilist)
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
#ifdef MODAL_AERO
subroutine set_initial_sz_and_volumes(list_idx, top_lev, ncol, imode, & !input
     dgncur, v2ncur, dryvol )                                                 !output

  !-----------------------------------------------------------------------------
  !Purpose: Set initial defaults for the dry diameter, volume to num
  ! and dry volume
  !
  !Called by: modal_aero_calcsize_sub
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------
  implicit none

  !inputs
  integer, intent(in) :: list_idx, top_lev !for model level loop
  integer, intent(in) :: ncol          !# of columns
  integer, intent(in) :: imode         !mode index

  !outputs
  real(r8), intent(out) :: dgncur(:,:,:) !diameter
  real(r8), intent(out) :: v2ncur(:,:,:) !volume to number
  real(r8), intent(out) :: dryvol(:,:)   !dry volume

  !local variables
  integer  :: icol, klev
  real(r8) :: dgnum, sigmag, voltonumb

  call rad_cnst_get_mode_props(list_idx, imode, dgnum=dgnum, sigmag=sigmag)
  voltonumb   = 1._r8 / ( (pi/6._r8)*(dgnum**3.0_r8)*exp(4.5_r8*log(sigmag)**2.0_r8) )

  do klev = top_lev, pver
     do icol = 1, ncol
        dgncur(icol,klev,imode) = dgnum     !diameter
        v2ncur(icol,klev,imode) = voltonumb !volume to number
        dryvol(icol,klev)       = 0.0_r8    !initialize dry vol
     end do
  end do

  return

end subroutine set_initial_sz_and_volumes

!---------------------------------------------------------------------------------------------

subroutine compute_dry_volume(top_lev, ncol, imode, nspec, state, pbuf, &
     dryvol_a, dryvol_c, list_idx_in)

  !-----------------------------------------------------------------------------
  !Purpose: Compute initial dry volume based on mmr and specie density
  ! volume = mmr/density
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : rad_cnst_get_aer_mmr, rad_cnst_get_aer_props
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: top_lev !for model level loop
  integer,  intent(in) :: ncol          !# of columns
  integer,  intent(in) :: imode         !mode index
  integer,  intent(in) :: nspec         !# of species in mode "imode"
  integer,  intent(in) :: list_idx_in
  type(physics_state), target, intent(in)    :: state       ! Physics state variables
  type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer

  !in-outs
  real(r8), intent(inout) :: dryvol_a(:,:)                    ! interstital aerosol dry
  real(r8), intent(inout) :: dryvol_c(:,:)                    ! interstital aerosol dry

  !local vars
  integer  :: ispec, icol, klev, lchnk, idx_cw
  real(r8) :: specdens  !specie density
  real(r8) :: dummwdens !density inverse
  real(r8), pointer :: specmmr(:,:)  !specie mmr (interstitial)
  real(r8), pointer :: specmmr_cld(:,:)  !specie mmr (cloud borne)

  character(len=32) :: spec_name

  lchnk = state%lchnk !get chunk info for retrieving cloud borne aerosols mmr

  do ispec = 1, nspec
     ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air

     !Get mmr for the interstitial specie
     call rad_cnst_get_aer_mmr(list_idx_in, imode, ispec, 'a', state, pbuf, specmmr)
     !Get mmr for the cloud borne species
     call rad_cnst_get_aer_mmr(list_idx_in, imode, ispec, 'c', state, pbuf, specmmr_cld)

     !get density of the aerosol specie
     call rad_cnst_get_aer_props(list_idx_in, imode, ispec, density_aer=specdens)
     dummwdens = 1.0_r8 / specdens !inverse of density

     !compute dry volume as a function of space (i,k)
     do klev = top_lev, pver
        do icol = 1, ncol
           dryvol_a(icol,klev) = dryvol_a(icol,klev) + max(0.0_r8,specmmr(icol,klev))*dummwdens
           dryvol_c(icol,klev) = dryvol_c(icol,klev) + max(0.0_r8,specmmr_cld(icol,klev))*dummwdens
        end do
     end do

  end do ! nspec loop

  return
end subroutine compute_dry_volume

!---------------------------------------------------------------------------------------------

subroutine size_adjustment(list_idx, top_lev, ncol, lchnk, imode, dryvol_a, state_q, & !input
     dryvol_c, pdel, do_adjust, update_mmr, do_aitacc_transfer, deltatinv, fracadj, pbuf,  & !input
     dgncur_a, dgncur_c, v2ncur_a, v2ncur_c, &
     drv_a_accsv, drv_c_accsv, drv_a_aitsv, drv_c_aitsv, drv_a_sv, drv_c_sv, &
     num_a_accsv, num_c_accsv, num_a_aitsv, num_c_aitsv, num_a_sv, num_c_sv, &
     dotend, dotendqqcw, dqdt, dqqcwdt, qsrflx)

  !-----------------------------------------------------------------------------
  !Purpose: Do the aerosol size adjustment if needed
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : compute_dgn_vol_limits, adjust_num_sizes, update_dgn_voltonum
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: list_idx
  integer,  intent(in) :: top_lev!for model level loop
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

  !outputs
  real(r8), intent(inout) :: dgncur_a(:,:,:) ! TODO: Add comments here!
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
  integer :: num_mode_idx, num_cldbrn_mode_idx, mam_ait, mam_acc, nait, nacc
  real(r8) :: v2nmax, v2nmin                  ! voltonumblo/hi of current mode
  real(r8) :: v2nmaxrl, v2nminrl              ! relaxed voltonumblo/hi
  real(r8) :: dgnyy, dgnxx                  ! dgnumlo/hi of current mode
  real(r8) :: drv_a, drv_c                  ! dry volume (cm3/mol_air)
  real(r8) :: num_a0, num_c0                ! initial number (#/mol_air)
  real(r8) :: num_a, num_c                  ! working number (#/mol_air)
  real(r8) :: sigmag, cmn_factor, dgnumlo, dgnumhi
  real(r8), pointer :: fldcw(:,:)           !specie mmr (cloud borne)

  !find state q array number mode indices for interstitial and cloud borne aerosols
  !Both num_mode_idx and num_cldbrn_mode_idx should be exactly same and should be same
  !for both prognostic and diagnostic radiation lists
  num_mode_idx        = numptr_amode(imode)
  num_cldbrn_mode_idx = numptrcw_amode(imode)

  mam_ait = modeptr_aitken !aitken mode number in this mam package
  mam_acc = modeptr_accum  !accumulation mode number in this mam package

  !find out accumulation and aitken modes in the radiation list
  nacc = rad_cnst_get_mode_idx(list_idx, modename_amode(mam_acc))
  nait = rad_cnst_get_mode_idx(list_idx, modename_amode(mam_ait))

  !set hi/lo limits (min, max and their relaxed counterparts) for volumes and diameters
  call compute_dgn_vol_limits(list_idx, imode, nait, nacc, do_aitacc_transfer, & !input
       v2nmin, v2nmax, v2nminrl, v2nmaxrl, dgnxx, dgnyy) !output

  !set tendency logicals to true if we need to update mmr
  if (update_mmr) then
     dotend(num_mode_idx)            = .true.
     dotendqqcw(num_cldbrn_mode_idx) = .true.
  end if

  !pointer to cloud borne number mmr for mode imode
  !(it is okay to use qqcw_get_field to get number mixing ratios for diagnostic calls)
  fldcw => qqcw_get_field(pbuf,num_cldbrn_mode_idx,lchnk,.true.)

  !Compute a common factor for size computations
  call rad_cnst_get_mode_props(list_idx, imode, sigmag=sigmag, dgnumhi=dgnumhi, dgnumlo=dgnumlo)
  cmn_factor = exp(4.5_r8*log(sigmag)**2.0_r8)*pi/6.0_r8

  do  klev = top_lev, pver
     do  icol = 1, ncol

        drv_a  = dryvol_a(icol,klev)

        !NOTE: number mixing ratios are always present for diagnostic
        !calls, so it is okay to use "state_q" instead of rad_cnst calls
        num_a0 = state_q(icol,klev,num_mode_idx)
        num_a  = max( 0.0_r8, num_a0 )
        drv_c  = dryvol_c(icol,klev)
        num_c0 = fldcw(icol,klev)
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
           call adjust_num_sizes(icol, klev, update_mmr, num_mode_idx, num_cldbrn_mode_idx, &                   !input
                   drv_a, num_a0, drv_c, num_c0, deltatinv, v2nmin, v2nminrl, v2nmax, v2nmaxrl, fracadj, & !input
                   num_a, num_c, dqdt, dqqcwdt)                                                        !output

        endif !do_adjust

        !
        ! now compute current dgn and v2n
        !
        call update_dgn_voltonum(icol, klev, imode, update_mmr, num_mode_idx, inter_aero, gravit, cmn_factor, drv_a, num_a, &
             v2nmin, v2nmax, dgnxx, dgnyy, pdel, &
             dgncur_a, v2ncur_a, qsrflx, dqdt)

        !for cloud borne aerosols
        call update_dgn_voltonum(icol, klev, imode, update_mmr, num_cldbrn_mode_idx, cld_brn_aero, gravit, cmn_factor, drv_c, num_c, &
             v2nmin, v2nmax, dgnumhi, dgnumlo, pdel, &
             dgncur_c, v2ncur_c, qsrflx, dqqcwdt)

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

subroutine compute_dgn_vol_limits(list_idx, imode, nait, nacc, do_aitacc_transfer, & !input
           v2nmin, v2nmax, v2nminrl, v2nmaxrl, dgnxx, dgnyy) !output

  !-----------------------------------------------------------------------------
  !Purpose: Compute hi and lo limits for diameter and volume based on the
  ! relaxation factor
  !
  !Called by: size_adjustment
  !Calls    : None
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

implicit none

!inputs
integer,  intent(in) :: list_idx, imode
integer,  intent(in) :: nait, nacc
logical,  intent(in) :: do_aitacc_transfer

!outputs
real(r8), intent(out) :: v2nmin
real(r8), intent(out) :: v2nmax
real(r8), intent(out) :: v2nminrl
real(r8), intent(out) :: v2nmaxrl
real(r8), intent(out) :: dgnxx
real(r8), intent(out) :: dgnyy

!local
real(r8), parameter :: relax_factor = 27.0_r8 !relax_factor=3**3=27,
                                              !i.e. dgnumlo_relaxed = dgnumlo/3 and dgnumhi_relaxed = dgnumhi*3
real(r8), parameter :: szadj_block_fac = 1.0e6_r8
real(r8) :: dgnumlo, dgnumhi, sigmag, cmn_factor


call rad_cnst_get_mode_props(list_idx, imode, dgnumlo=dgnumlo, dgnumhi=dgnumhi, sigmag=sigmag)
cmn_factor = exp(4.5_r8*log(sigmag)**2.0_r8)*pi/6.0_r8
!v2nmin = voltonumbhi is proportional to dgnumhi**(-3),
!        and produces the minimum allowed number for a given volume
v2nmin   = 1._r8 / ( (pi/6._r8)*(dgnumhi**3.0_r8)*exp(4.5_r8*log(sigmag)**2.0_r8) )

!v2nmax = voltonumblo is proportional to dgnumlo**(-3),
!        and produces the maximum allowed number for a given volume
v2nmax   = 1._r8 / ( (pi/6._r8)*(dgnumlo**3.0_r8)*exp(4.5_r8*log(sigmag)**2.0_r8) )

!v2nminrl and v2nmaxrl are their "relaxed" equivalents.
v2nminrl = v2nmin/relax_factor
v2nmaxrl = v2nmax*relax_factor
dgnxx   = dgnumhi
dgnyy   = dgnumlo


!if do_aitacc_transfer is turned on, we will do the ait<->acc tranfer separately in
!aitken_accum_exchange subroutine, so we are turning the size adjustment for these
!two modes here.
if ( do_aitacc_transfer ) then
   !for n=nait, divide v2nmin by 1.0e6 to effectively turn off the
   !         adjustment when number is too small (size is too big)
   if (imode == nait) v2nmin = v2nmin/szadj_block_fac

   !for n=nacc, multiply v2nmax by 1.0e6 to effectively turn off the
   !         adjustment when number is too big (size is too small)
   if (imode == nacc) v2nmax = v2nmax*szadj_block_fac

   !Also change the v2nmaxrl/v2nminrl so that
   !the interstitial<-->activated adjustment is turned off
   v2nminrl = v2nmin/relax_factor
   v2nmaxrl = v2nmax*relax_factor
end if

return
end subroutine compute_dgn_vol_limits

!---------------------------------------------------------------------------------------------

subroutine adjust_num_sizes(icol, klev, update_mmr, num_mode_idx, num_cldbrn_mode_idx, &
     drv_a, num_a0, drv_c, num_c0, deltatinv, v2nmin, v2nminrl, v2nmax, v2nmaxrl, fracadj, &
     num_a, num_c, dqdt, dqqcwdt)

  !-----------------------------------------------------------------------------
  !Purpose: Adjust num sizes if needed
  !
  !Called by: size_adjustment
  !Calls    : None
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: icol, klev
  logical,  intent(in) :: update_mmr
  real(r8), intent(in) :: drv_a, num_a0, drv_c, num_c0
  real(r8), intent(in) :: deltatinv, v2nmin, v2nminrl, v2nmax, v2nmaxrl, fracadj
  integer,  intent(in) :: num_mode_idx, num_cldbrn_mode_idx

  !outputs
  real(r8), intent(inout) :: num_a, num_c
  real(r8), intent(inout) :: dqdt(:,:,:), dqqcwdt(:,:,:)

  !local
  real(r8) :: num_a1, num_c1, numbnd
  real(r8) :: num_a2, num_c2, delnum_a2, delnum_c2
  real(r8) :: delnum_a3, delnum_c3
  real(r8) :: num_t2, drv_t, delnum_t3


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
     numbnd = max( drv_a*v2nmin, min( drv_a*v2nmax, num_a1 ) )
     num_a  = num_a1 + (numbnd - num_a1)*fracadj
     if(update_mmr) then
        dqdt(icol,klev,num_mode_idx)           = (num_a - num_a0)*deltatinv
        dqqcwdt(icol,klev,num_cldbrn_mode_idx) = -num_c0*deltatinv
     endif
  else if (drv_a <= 0.0_r8) then
     ! interstitial volume is zero, treat similar to above
     num_a = 0.0_r8
     num_c1 = num_c
     numbnd = max( drv_c*v2nmin, min( drv_c*v2nmax, num_c1 ) )
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
     numbnd = max( drv_a*v2nminrl, min( drv_a*v2nmaxrl, num_a1 ) )
     delnum_a2 = (numbnd - num_a1)*fracadj
     num_a2 = num_a1 + delnum_a2
     numbnd = max( drv_c*v2nminrl, min( drv_c*v2nmaxrl, num_c1 ) )
     delnum_c2 = (numbnd - num_c1)*fracadj
     num_c2 = num_c1 + delnum_c2
     if ((delnum_a2 == 0.0_r8) .and. (delnum_c2 /= 0.0_r8)) then
        num_a2 = max( drv_a*v2nminrl, min( drv_a*v2nmaxrl,   &
             num_a1-delnum_c2 ) )
     else if ((delnum_a2 /= 0.0_r8) .and. (delnum_c2 == 0.0_r8)) then
        num_c2 = max( drv_c*v2nminrl, min( drv_c*v2nmaxrl,   &
             num_c1-delnum_a2 ) )
     end if
     ! step3:  num_a,c2 --> num_a,c3 applies stricter bounds to the
     !    combined/total number
     drv_t = drv_a + drv_c
     num_t2 = num_a2 + num_c2
     delnum_a3 = 0.0_r8
     delnum_c3 = 0.0_r8
     if (num_t2 < drv_t*v2nmin) then
        delnum_t3 = (drv_t*v2nmin - num_t2)*fracadj
        ! if you are here then (num_a2 < drv_a*v2nmin) and/or
        !                      (num_c2 < drv_c*v2nmin) must be true
        if ((num_a2 < drv_a*v2nmin) .and. (num_c2 < drv_c*v2nmin)) then
           delnum_a3 = delnum_t3*(num_a2/num_t2)
           delnum_c3 = delnum_t3*(num_c2/num_t2)
        else if (num_c2 < drv_c*v2nmin) then
           delnum_c3 = delnum_t3
        else if (num_a2 < drv_a*v2nmin) then
           delnum_a3 = delnum_t3
        end if
     else if (num_t2 > drv_t*v2nmax) then
        delnum_t3 = (drv_t*v2nmax - num_t2)*fracadj
        ! if you are here then (num_a2 > drv_a*v2nmax) and/or
        !                      (num_c2 > drv_c*v2nmax) must be true
        if ((num_a2 > drv_a*v2nmax) .and. (num_c2 > drv_c*v2nmax)) then
           delnum_a3 = delnum_t3*(num_a2/num_t2)
           delnum_c3 = delnum_t3*(num_c2/num_t2)
        else if (num_c2 > drv_c*v2nmax) then
           delnum_c3 = delnum_t3
        else if (num_a2 > drv_a*v2nmax) then
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

subroutine update_dgn_voltonum(icol, klev, imode, update_mmr, num_idx, aer_type, gravit, cmn_factor, drv, num, &
     v2nmin, v2nmax, dgnxx, dgnyy, pdel, &
     dgncur, v2ncur, qsrflx, dqdt)

  !-----------------------------------------------------------------------------
  !Purpose: updates diameter and volume to num based on limits
  !
  !Called by: size_adjustment
  !Calls    : None
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: icol, klev, imode, num_idx, aer_type
  logical,  intent(in) :: update_mmr
  real(r8), intent(in) :: gravit, cmn_factor
  real(r8), intent(in) :: drv, num
  real(r8), intent(in) :: v2nmin, v2nmax, dgnxx, dgnyy
  real(r8), intent(in) :: pdel(:,:)
  real(r8), intent(in) :: dqdt(:,:,:)
  !output
  real(r8), intent(inout) :: dgncur(:,:,:), v2ncur(:,:,:)
  real(r8), intent(inout) :: qsrflx(:,:,:,:)

  !local
  real(r8) :: pdel_fac

  if (drv > 0.0_r8) then
     if (num <= drv*v2nmin) then
        dgncur(icol,klev,imode) = dgnxx
        v2ncur(icol,klev,imode) = v2nmin
     else if (num >= drv*v2nmax) then
        dgncur(icol,klev,imode) = dgnyy
        v2ncur(icol,klev,imode) = v2nmax
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

subroutine  aitken_accum_exchange(ncol, lchnk, list_idx, update_mmr, tadjinv, &
     deltat, pdel, state_q, state, &
     pbuf, &
     drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv,     &
     drv_a_accsv,num_a_accsv, drv_c_accsv, num_c_accsv,      &
     dgncur_a, v2ncur_a, dgncur_c, v2ncur_c, dotend, dotendqqcw, &
     dqdt, dqqcwdt, qsrflx)

  !-----------------------------------------------------------------------------
  !Purpose: Exchange aerosols between aitken and accumulation modes based on new
  ! sizes
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : endrun
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

  implicit none

  !inputs
  integer,  intent(in) :: ncol, lchnk, list_idx
  logical,  intent(in) :: update_mmr
  real(r8), intent(in) :: tadjinv
  real(r8), intent(in) :: deltat
  real(r8), intent(in) :: pdel(:,:)
  real(r8), intent(in) :: state_q(:,:,:) !state_q should only be used for list_idx==0
  real(r8), intent(in) :: drv_a_aitsv(:,:), num_a_aitsv(:,:)
  real(r8), intent(in) :: drv_a_accsv(:,:), num_a_accsv(:,:)
  real(r8), intent(in) :: drv_c_aitsv(:,:), num_c_aitsv(:,:)
  real(r8), intent(in) :: drv_c_accsv(:,:), num_c_accsv(:,:)
  type(physics_state), intent(in)    :: state       ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer

  !outputs
  real(r8), intent(inout) :: dgncur_a(:,:,:) !TODO: Add comments here!
  real(r8), intent(inout) :: dgncur_c(:,:,:)
  real(r8), intent(inout) :: v2ncur_a(:,:,:)
  real(r8), intent(inout) :: v2ncur_c(:,:,:)
  logical,  intent(inout) :: dotend(:), dotendqqcw(:)
  real(r8), intent(inout) :: dqdt(:,:,:), dqqcwdt(:,:,:), qsrflx(:,:,:,:)

  !local
  integer  :: imode, jmode, aer_type, jsrflx, icol, klev, iq
  integer  :: lsfrm, lstoo, ispec_acc, idx
  integer  :: ixfer_ait2acc, ixfer_acc2ait, mam_ait, mam_acc
  integer  :: iait, iacc, idx_cw, nspec_acc
  integer, save  :: idiagaa = 1
  real(r8) :: num_a, drv_a, num_c, drv_c
  real(r8) :: num_a_acc, num_c_acc
  real(r8) :: drv_a_acc, drv_c_acc
  real(r8) :: v2n_geomean, pdel_fac, num_t, drv_t
  real(r8) :: drv_a_noxf, drv_c_noxf, drv_t_noxf, dummwdens
  real(r8) :: num_t0, num_t_noxf
  real(r8) :: duma, cmn_factor, dgnum, sigmag, voltonumb_ait, voltonumb_acc
  real(r8) :: xfercoef
  real(r8) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
  real(r8) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
  real(r8) :: xfertend, xfertend_num(2,2)
  logical  :: noxf_acc2ait(ntot_aspectype), accum_exists, aitken_exists

  character(len=32)  :: spec_name
  character(len=800) :: err_msg

  if (npair_csizxf .le. 0)call endrun('npair_csizxf <= 0'//errmsg(__FILE__,__LINE__))

  ! check that calcsize transfer is aitken-->accum
  mam_ait = modeptr_aitken !aitken mode number in this mam package
  mam_acc = modeptr_accum  !accumulation mode number in this mam package

  !find out accumulation and aitken modes in the radiation list
  iacc = rad_cnst_get_mode_idx(list_idx, modename_amode(mam_acc))
  iait = rad_cnst_get_mode_idx(list_idx, modename_amode(mam_ait))

  !find out if aitken or accumulation modes exist in the radiation list
  !(a positive value means that the mode exists)
  accum_exists  = ( iacc > 0)
  aitken_exists = ( iait > 0)

  if(.not.(accum_exists .and. aitken_exists)) then
     write(err_msg,*)'Accumulation or the Aitken mode do not exist in list:', &
          list_idx,', Accu mode:',iacc,', Aitken mode:',iait,', a negative mode is', &
          ' the non-existent mode ',errmsg(__FILE__,__LINE__)
     call endrun(trim(err_msg))
  endif

  if (modefrm_csizxf(list_idx) .ne. iait .or.modetoo_csizxf(list_idx) .ne. iacc) then
     write(err_msg,*)'modefrm/too_csizxf are wrong for radiation list:',list_idx,' ',errmsg(__FILE__,__LINE__)
     call endrun(trim(err_msg))
  endif

  ! set dotend() for species that will be transferred
  ! for both mass and number
  if(update_mmr) then
     do iq = 1, nspecfrm_csizxf(list_idx)
        lsfrm = lspecfrma_csizxf(iq,list_idx)
        lstoo = lspectooa_csizxf(iq,list_idx)
        if ((lsfrm > 0) .and. (lstoo > 0)) then
           dotend(lsfrm) = .true.
           dotend(lstoo) = .true.
        end if
        !for cloud borne aerosols
        lsfrm = lspecfrmc_csizxf(iq,list_idx)
        lstoo = lspectooc_csizxf(iq,list_idx)
        if ((lsfrm > 0) .and. (lstoo > 0)) then
           dotendqqcw(lsfrm) = .true.
           dotendqqcw(lstoo) = .true.
        end if
     end do
  endif

  !------------------------------------------------------------------------
  ! Identify accum species cannot be transferred to aitken mode
  !
  ! Accumulation mode have more species than Aitken mode. Therefore, there
  ! will be some species which cannot be transferred from accumulation to
  ! Aitken mode as they don't exist in the Aitken mode
  !------------------------------------------------------------------------

  noxf_acc2ait(:) = .true. ! let us assume we are not going to move any species

  !Now compute species which can be transfered
  !Get # of species in accumulation mode
  call rad_cnst_get_info(list_idx, iacc, nspec = nspec_acc) !output:nspec_acc

  do ispec_acc = 1, nspec_acc     !Go through all species within accumulation mode (only species, not number (e.g. num_a1))
     call rad_cnst_get_info(list_idx, iacc, ispec_acc,spec_name=spec_name) !output:spec_name
     call cnst_get_ind(spec_name, idx)
     do iq = 1, nspecfrm_csizxf(list_idx) !Go through all mapped species (and number) for the accumulation mode
        if (lspectooa_csizxf(iq,list_idx) == idx) then !compare idx with mapped species in the accumulation mode
           noxf_acc2ait(ispec_acc) = .false. ! species which can be tranferred
        end if
     end do
  end do


  ! v2n_geomean is voltonumb at the "geometrically-defined" mid-point
  ! between the aitken and accum modes

  call rad_cnst_get_mode_props(list_idx, iait, dgnum=dgnum, sigmag=sigmag)
  voltonumb_ait   = 1._r8 / ( (pi/6._r8)*(dgnum**3.0_r8)*exp(4.5_r8*log(sigmag)**2.0_r8) )

  call rad_cnst_get_mode_props(list_idx, iacc, dgnum=dgnum, sigmag=sigmag)
  voltonumb_acc   = 1._r8 / ( (pi/6._r8)*(dgnum**3._r8)*exp(4.5_r8*log(sigmag)**2._r8) )

  v2n_geomean = sqrt(voltonumb_ait*voltonumb_acc)

  ! loop over columns and levels
  do  klev = top_lev, pver
     do  icol = 1, ncol

        pdel_fac = pdel(icol,klev)/gravit   ! = rho*dz

        !Compute aitken->accumulation transfer
        call compute_coef_ait_acc_transfer(iacc, v2n_geomean, tadjinv, drv_a_aitsv(icol,klev),       & !input
             drv_c_aitsv (icol,klev), num_a_aitsv(icol,klev), num_c_aitsv(icol,klev), voltonumb_acc, & !input
             ixfer_ait2acc, xfercoef_num_ait2acc, xfercoef_vol_ait2acc, xfertend_num)                  !output


        !----------------------------------------------------------------------------------------
        ! compute accum --> aitken transfer rates
        !
        ! accum may have some species (seasalt, dust, poa etc.) that are
        !    not in aitken mode
        ! so first divide the accum drv & num into not-transferred (noxf) species
        !    and transferred species, and use the transferred-species
        !    portion in what follows
        !----------------------------------------------------------------------------------------
        call compute_coef_acc_ait_transfer(iait, iacc, icol, klev, list_idx, lchnk, v2n_geomean, & !input
             tadjinv, state, pbuf, drv_a_accsv(icol,klev), drv_c_accsv (icol, klev),       & !input
             num_a_accsv(icol,klev), num_c_accsv(icol,klev), noxf_acc2ait, voltonumb_ait,  & !input
             drv_a_noxf, drv_c_noxf, ixfer_acc2ait, xfercoef_num_acc2ait,                  & !output
             xfercoef_vol_acc2ait, xfertend_num)                                             !output


        ! jump to end-of-loop if no transfer is needed at current icol,klev
        if (ixfer_ait2acc+ixfer_acc2ait > 0) then
           !
           ! compute new dgncur & v2ncur for aitken & accum modes
           !
           ! currently inactive (??? BSINGH: Not sure what this comment refers to...)

           !interstitial species
           duma = (xfertend_num(1,1) - xfertend_num(2,1))*deltat    !diff in num from  ait->accum and accum->ait transfer
           num_a     = max( 0.0_r8, num_a_aitsv(icol,klev) - duma ) !num removed/added from aitken mode
           num_a_acc = max( 0.0_r8, num_a_accsv(icol,klev) + duma ) !num added/removed to accumulation mode

           duma = (drv_a_aitsv(icol,klev)*xfercoef_vol_ait2acc -   &
                (drv_a_accsv(icol,klev)-drv_a_noxf)*xfercoef_vol_acc2ait)*deltat ! diff in volume transfer fomr ait->accum and accum->ait transfer
           drv_a     = max( 0.0_r8, drv_a_aitsv(icol,klev) - duma ) !drv removed/added from aitken mode
           drv_a_acc = max( 0.0_r8, drv_a_accsv(icol,klev) + duma ) !drv added/removed to accumulation mode

           !cloud borne species
           duma = (xfertend_num(1,2) - xfertend_num(2,2))*deltat    !same as above for cloud borne aerosols

           num_c     = max( 0.0_r8, num_c_aitsv(icol,klev) - duma )
           num_c_acc = max( 0.0_r8, num_c_accsv(icol,klev) + duma )
           duma = (drv_c_aitsv(icol,klev)*xfercoef_vol_ait2acc -   &
                (drv_c_accsv(icol,klev)-drv_c_noxf)*xfercoef_vol_acc2ait)*deltat
           drv_c     = max( 0.0_r8, drv_c_aitsv(icol,klev) - duma )
           drv_c_acc = max( 0.0_r8, drv_c_accsv(icol,klev) + duma )

           !interstitial species (aitken mode)
           call compute_new_sz_after_transfer(list_idx, iait, drv_a, num_a, &
                dgncur_a(icol,klev,iait), v2ncur_a(icol,klev,iait))

           !cloud borne species (aitken mode)
           call compute_new_sz_after_transfer(list_idx, iait, drv_c, num_c, &
                dgncur_c(icol,klev,iait), v2ncur_c(icol,klev,iait))

           num_a = num_a_acc
           drv_a = drv_a_acc
           num_c = num_c_acc
           drv_c = drv_c_acc

           !interstitial species (accumulation mode)
           call compute_new_sz_after_transfer(list_idx, iacc, drv_a, num_a, &
                dgncur_a(icol,klev,iacc), v2ncur_a(icol,klev,iacc))

           !cloud borne species (accumulation mode)
           call compute_new_sz_after_transfer(list_idx, iacc, drv_c, num_c, &
                dgncur_c(icol,klev,iacc), v2ncur_c(icol,klev,iacc))

           !Printout just once in hte log file (idiagaa is set to negaitive after this if condition)
           if ( masterproc ) then
              if (idiagaa > 0 .and. list_idx == 0) then
                 do jmode = 1, 2
                    do iq = 1, nspecfrm_csizxf(list_idx)
                       do aer_type = 1, 2
                          if (jmode .eq. 1) then
                             if (aer_type .eq. inter_aero) then
                                lsfrm = lspecfrma_csizxf(iq,list_idx)
                                lstoo = lspectooa_csizxf(iq,list_idx)
                             else
                                lsfrm = lspecfrmc_csizxf(iq,list_idx)
                                lstoo = lspectooc_csizxf(iq,list_idx)
                             end if
                          else
                             if (aer_type .eq. inter_aero) then
                                lsfrm = lspectooa_csizxf(iq,list_idx)
                                lstoo = lspecfrma_csizxf(iq,list_idx)
                             else
                                lsfrm = lspectooc_csizxf(iq,list_idx)
                                lstoo = lspecfrmc_csizxf(iq,list_idx)
                             end if
                          end if
                          write( iulog, '(a,3i3,2i4)' ) 'calcsize jmode,iq,aer_type, lsfrm,lstoo',   &
                               jmode,iq,aer_type, lsfrm,lstoo
                       end do
                    end do
                 end do
              end if
           end if
           if(list_idx == 0)idiagaa = -1

           !------------------------------------------------------------------
           ! compute tendency amounts for aitken <--> accum transfer
           !------------------------------------------------------------------

           !ASSUMPTION: "update_mmr" will only be true for the prognostic radiation list(i.e. list_idx=0, "radiation_climate")
           !If list_idx=0, it is okay to get specie mmr from state_q array. Therefore, state_q is used in update_tends_flx calls
           !below
           if(update_mmr) then
              ! jmode=1 does aitken-->accum
              if(ixfer_ait2acc > 0) then
                 jmode = 1
                 call update_tends_flx(icol, klev, jmode, list_idx, lchnk , lspecfrma_csizxf, lspectooa_csizxf, &
                      lspecfrmc_csizxf, lspectooc_csizxf, xfertend_num, xfercoef_vol_ait2acc, state_q, pbuf, &
                      pdel_fac, dqdt, dqqcwdt, qsrflx)
              endif

              !jmode=2 does accum-->aitken
              if(ixfer_acc2ait > 0) then
                 jmode = 2
                 !suboutine (update_tends_flx) is called but lspectooa and lspecfrma are switched and
                 !xfercoef_vol_acc2ait is also used instead xfercoef_vol_ait2acc for accum->aitken transfer
                 call update_tends_flx(icol, klev, jmode, list_idx, lchnk , lspectooa_csizxf, lspecfrma_csizxf, &
                      lspectooc_csizxf, lspecfrmc_csizxf, xfertend_num, xfercoef_vol_acc2ait, state_q, pbuf, &
                      pdel_fac, dqdt, dqqcwdt, qsrflx)
              endif
           endif !update_mmr

        end if !ixfer_ait2acc+ixfer_acc2ait > 0
     end do !ncol
  end do !pver

  return
end subroutine aitken_accum_exchange

!---------------------------------------------------------------------------------------------

subroutine compute_coef_ait_acc_transfer(iacc, v2n_geomean, tadjinv, drv_a_aitsv, &
     drv_c_aitsv, num_a_aitsv, num_c_aitsv,  voltonumb_acc, &
     ixfer_ait2acc, xfercoef_num_ait2acc, xfercoef_vol_ait2acc, xfertend_num)

  !------------------------------------------------------------
  ! Purpose: Computes coefficients for transfer from aitken to accumulation mode
  !
  ! Author: Richard Easter (Refactored by Balwinder Singh)
  !------------------------------------------------------------

  !intent ins
  integer,  intent(in) :: iacc
  real(r8), intent(in) :: v2n_geomean
  real(r8), intent(in) :: tadjinv
  real(r8), intent(in) :: drv_a_aitsv, drv_c_aitsv
  real(r8), intent(in) :: num_a_aitsv, num_c_aitsv, voltonumb_acc

  !intent outs
  integer,  intent(inout) :: ixfer_ait2acc
  real(r8), intent(inout) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
  real(r8), intent(inout) :: xfertend_num(2,2)

  !local
  real(r8) :: drv_t, num_t
  real(r8) :: xferfrac_num_ait2acc, xferfrac_vol_ait2acc

  !initialize
  ixfer_ait2acc        = 0
  xfercoef_num_ait2acc = 0.0_r8
  xfercoef_vol_ait2acc = 0.0_r8
  xfertend_num(:,:)    = 0.0_r8

  ! compute aitken --> accum transfer rates

  drv_t = drv_a_aitsv + drv_c_aitsv
  num_t = num_a_aitsv + num_c_aitsv
  if (drv_t > 0.0_r8) then
     !if num is less than the mean value, we have large particles (keeping volume constant drv_t)
     !which needs to be moved to accumulation mode
     if (num_t < drv_t*v2n_geomean) then
        ixfer_ait2acc = 1
        if (num_t < drv_t*voltonumb_acc) then ! move all particles if number is smaller than the acc mean
           xferfrac_num_ait2acc = 1.0_r8
           xferfrac_vol_ait2acc = 1.0_r8
        else !otherwise scale the transfer
           xferfrac_vol_ait2acc = ((num_t/drv_t) - v2n_geomean)/   &
                (voltonumb_acc - v2n_geomean)
           xferfrac_num_ait2acc = xferfrac_vol_ait2acc*   &
                (drv_t*voltonumb_acc/num_t)
           !bound the transfer coefficients between 0 and 1
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
        xfertend_num(1,1) = num_a_aitsv*xfercoef_num_ait2acc
        xfertend_num(1,2) = num_c_aitsv*xfercoef_num_ait2acc
     end if
  end if

end subroutine compute_coef_ait_acc_transfer

!---------------------------------------------------------------------------------------------

subroutine compute_coef_acc_ait_transfer( iait, iacc, icol, klev, list_idx, lchnk,     &
     v2n_geomean, tadjinv, state, pbuf, drv_a_accsv, drv_c_accsv, num_a_accsv,      &
     num_c_accsv, noxf_acc2ait, voltonumb_ait,                                &
     drv_a_noxf, drv_c_noxf, ixfer_acc2ait, xfercoef_num_acc2ait, &
     xfercoef_vol_acc2ait, xfertend_num)

  !intent -ins
  integer,  intent(in) :: iait, iacc, icol, klev, list_idx, lchnk
  real(r8), intent(in) :: v2n_geomean
  real(r8), intent(in) :: tadjinv
  real(r8), intent(in) :: drv_a_accsv, drv_c_accsv
  real(r8), intent(in) :: num_a_accsv, num_c_accsv, voltonumb_ait
  logical,  intent(in) :: noxf_acc2ait(:)

  type(physics_state), intent(in)    :: state       ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer

  !intent - outs
  integer,  intent(inout) :: ixfer_acc2ait
  real(r8), intent(inout) :: drv_a_noxf, drv_c_noxf
  real(r8), intent(inout) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
  real(r8), intent(inout) :: xfertend_num(2,2)

  !local
  integer  :: ispec_acc, idx, nspec_acc
  real(r8) :: drv_t, num_t, drv_t_noxf, num_t0
  real(r8) :: num_t_noxf
  real(r8) :: dummwdens, specdens, dgnumlo, sigmag, voltonumblo
  real(r8) :: xferfrac_num_acc2ait, xferfrac_vol_acc2ait
  real(r8), pointer :: specmmr(:,:)  !specie mmr (interstitial)
  real(r8), parameter :: zero_div_fac = 1.0e-37_r8

  ixfer_acc2ait = 0
  xfercoef_num_acc2ait = 0.0_r8
  xfercoef_vol_acc2ait = 0.0_r8

  drv_t = drv_a_accsv + drv_c_accsv
  num_t = num_a_accsv + num_c_accsv
  drv_a_noxf = 0.0_r8
  drv_c_noxf = 0.0_r8
  if (drv_t > 0.0_r8) then
     !if number is larger than the mean, it means we have small particles (keeping volume constant drv_t),
     !we need to move particles to aitken mode
     if (num_t > drv_t*v2n_geomean) then
        !As there may be more species in the accumulation mode which are not present in the aitken mode,
        !we need to compute the num and volume only for the species which can be transferred
        call rad_cnst_get_info(list_idx, iacc, nspec = nspec_acc)
        do ispec_acc = 1, nspec_acc

           if ( noxf_acc2ait(ispec_acc) ) then !species which can't be transferred
              call rad_cnst_get_aer_props(list_idx, iacc, ispec_acc, density_aer=specdens)    !get density
              ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
              dummwdens = 1.0_r8 / specdens
              call rad_cnst_get_aer_mmr(list_idx, iacc, ispec_acc, 'a', state, pbuf, specmmr) !get mmr
              drv_a_noxf = drv_a_noxf + max(0.0_r8,specmmr(icol,klev))*dummwdens

              call rad_cnst_get_aer_mmr(list_idx, iacc, ispec_acc, 'c', state, pbuf, specmmr) !get mmr
              drv_c_noxf = drv_c_noxf + max(0.0_r8,specmmr(icol,klev))*dummwdens
           end if
        end do
        drv_t_noxf = drv_a_noxf + drv_c_noxf !total volume which can't be moved to the aitken mode

        !Compute voltonumlo
        call rad_cnst_get_mode_props(list_idx, iacc, dgnumlo=dgnumlo, sigmag=sigmag)
        voltonumblo = 1._r8 / ( (pi/6._r8)*(dgnumlo**3.0_r8)*exp(4.5_r8*log(sigmag)**2.0_r8) )

        num_t_noxf = drv_t_noxf*voltonumblo !total number which can't be moved to the aitken mode
        num_t0 = num_t
        num_t = max( 0.0_r8, num_t - num_t_noxf )
        drv_t = max( 0.0_r8, drv_t - drv_t_noxf )
     end if
  end if

  if (drv_t > 0.0_r8) then
     !Find out if we need to transfer based on the new num_t
     if (num_t > drv_t*v2n_geomean) then
        ixfer_acc2ait = 1
        if (num_t > drv_t*voltonumb_ait) then! if number of larger than the aitken mean, move all particles
           xferfrac_num_acc2ait = 1.0_r8
           xferfrac_vol_acc2ait = 1.0_r8
        else ! scale the transfer
           xferfrac_vol_acc2ait = ((num_t/drv_t) - v2n_geomean)/   &
                (voltonumb_ait - v2n_geomean)
           xferfrac_num_acc2ait = xferfrac_vol_acc2ait*   &
                (drv_t*voltonumb_ait/num_t)
           !bound the transfer coefficients between 0 and 1
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
        xferfrac_num_acc2ait = xferfrac_num_acc2ait*   &
             num_t/max( zero_div_fac, num_t0 )
        xfercoef_num_acc2ait = xferfrac_num_acc2ait*tadjinv
        xfercoef_vol_acc2ait = xferfrac_vol_acc2ait*tadjinv
        xfertend_num(2,1) = num_a_accsv*xfercoef_num_acc2ait
        xfertend_num(2,2) = num_c_accsv*xfercoef_num_acc2ait
     end if
  end if

end subroutine compute_coef_acc_ait_transfer


subroutine compute_new_sz_after_transfer(list_idx, imode, drv, num, &
              dgncur, v2ncur)

implicit none

!intent-ins
integer,  intent(in) :: imode, list_idx
real(r8), intent(in) :: drv, num

!intent-outs
real(r8), intent(inout) :: dgncur, v2ncur

!local
real(r8) :: cmn_factor, sigmag, dgnum, dgnumlo, dgnumhi, voltonumb, voltonumbhi, voltonumblo

!Compute a common factor for size computations
call rad_cnst_get_mode_props(list_idx, imode, dgnum=dgnum, dgnumlo=dgnumlo, dgnumhi=dgnumhi, sigmag=sigmag)
cmn_factor = exp(4.5_r8*log(sigmag)**2.0_r8)*pi/6.0_r8

voltonumbhi = 1._r8 / (cmn_factor*dgnumhi**3.0_r8)
voltonumblo = 1._r8 / (cmn_factor*dgnumlo**3.0_r8)
voltonumb   = 1._r8 / (cmn_factor*dgnum**3.0_r8)

if (drv > 0.0_r8) then
   if (num <= drv*voltonumbhi) then
      dgncur = dgnumhi
      v2ncur = voltonumbhi
   else if (num >= drv*voltonumblo) then
      dgncur = dgnumlo
      v2ncur = voltonumblo
   else
      dgncur = (drv/(cmn_factor*num))**third
      v2ncur = num/drv
   end if
else
   dgncur = dgnum
   v2ncur = voltonumb
end if

end subroutine compute_new_sz_after_transfer


subroutine update_tends_flx(icol, klev, jmode, list_idx, lchnk , frm_spec_a, to_spec_a, &
     frm_spec_c, to_spec_c, xfertend_num, xfercoef, state_q, pbuf, &
     pdel_fac, dqdt, dqqcwdt, qsrflx)

  implicit none

  !intent - ins
  integer,  intent(in) :: icol, klev, jmode, list_idx, lchnk
  integer,  intent(in) :: frm_spec_a(maxspec_csizxf, 0:maxpair_csizxf)
  integer,  intent(in) :: to_spec_a(maxspec_csizxf, 0:maxpair_csizxf)
  integer,  intent(in) :: frm_spec_c(maxspec_csizxf, 0:maxpair_csizxf)
  integer,  intent(in) :: to_spec_c(maxspec_csizxf, 0:maxpair_csizxf)

  real(r8), intent(in) :: xfertend_num(2,2)
  real(r8), intent(in) :: xfercoef
  real(r8), intent(in) :: state_q(:,:,:)
  real(r8), intent(in) :: pdel_fac

  type(physics_buffer_desc), pointer :: pbuf(:)     ! physics buffer

  !intent -inout
  real(r8), intent(inout) :: dqdt(:,:,:), dqqcwdt(:,:,:), qsrflx(:,:,:,:)

  !local
  integer  :: iq, lsfrm, lstoo, jsrflx
  real(r8) :: xfertend
  real(r8), pointer :: fldcw(:,:)    !specie mmr (cloud borne)

  character(len=800) :: err_msg

  !check if list_idx > 0
  !This subroutine should ONLY be called for list_idx=0 as we are using state_q for specie mmr
  if(list_idx > 0) then
     write(err_msg,*)'updating mmrs is not supported for list_idx>0, current list_idx is:',list_idx,' ',errmsg(__FILE__,__LINE__)
     call endrun(trim(err_msg))
  endif

  jsrflx = jmode + 2

  !interstiatial species
  iq = 1 !iq = 1 is for num_* species
  lsfrm = frm_spec_a(iq,list_idx)
  lstoo = to_spec_a(iq,list_idx)
  if((lsfrm > 0) .and. (lstoo > 0)) then
     call update_num_tends(icol, klev, jmode, jsrflx, lsfrm, lstoo, inter_aero, pdel_fac, xfertend_num, dqdt, qsrflx)
  endif


  do iq = 2, nspecfrm_csizxf(list_idx)
     lsfrm = frm_spec_a(iq,list_idx)
     lstoo = to_spec_a(iq,list_idx)
     if((lsfrm > 0) .and. (lstoo > 0)) then
        xfertend = max(0.0_r8,state_q(icol,klev,lsfrm))*xfercoef
        dqdt(icol,klev,lsfrm) = dqdt(icol,klev,lsfrm) - xfertend
        dqdt(icol,klev,lstoo) = dqdt(icol,klev,lstoo) + xfertend
        qsrflx(icol,lsfrm,jsrflx,inter_aero) = qsrflx(icol,lsfrm,jsrflx,inter_aero) - xfertend*pdel_fac
        qsrflx(icol,lstoo,jsrflx,inter_aero) = qsrflx(icol,lstoo,jsrflx,inter_aero) + xfertend*pdel_fac
     endif
  enddo

  !cloud borne apecies
  iq = 1 !number species
  lsfrm = frm_spec_c(iq,list_idx)
  lstoo = to_spec_c(iq,list_idx)

  if((lsfrm > 0) .and. (lstoo > 0)) then
     call update_num_tends(icol, klev, jmode, jsrflx, lsfrm, lstoo, cld_brn_aero, pdel_fac, xfertend_num, dqqcwdt, qsrflx)
  endif

  !mass species
  do iq = 2, nspecfrm_csizxf(list_idx)
     lsfrm = frm_spec_c(iq,list_idx)
     lstoo = to_spec_c(iq,list_idx)
     if((lsfrm > 0) .and. (lstoo > 0)) then
        fldcw => qqcw_get_field(pbuf,lsfrm,lchnk)
        xfertend = max(0.0_r8,fldcw(icol,klev))*xfercoef
        dqqcwdt(icol,klev,lsfrm) = dqqcwdt(icol,klev,lsfrm) - xfertend
        dqqcwdt(icol,klev,lstoo) = dqqcwdt(icol,klev,lstoo) + xfertend
        qsrflx(icol,lsfrm,jsrflx,cld_brn_aero) = qsrflx(icol,lsfrm,jsrflx,cld_brn_aero) - xfertend*pdel_fac
        qsrflx(icol,lstoo,jsrflx,cld_brn_aero) = qsrflx(icol,lstoo,jsrflx,cld_brn_aero) + xfertend*pdel_fac
     end if
  enddo

end subroutine update_tends_flx

subroutine update_num_tends(icol, klev, jmode, jsrflx, lsfrm, lstoo, aer_type, pdel_fac, xfertend_num, dqdt, qsrflx)
  !intent ins
  integer,  intent(in) :: icol, klev, jmode, jsrflx, lsfrm, lstoo, aer_type
  real(r8), intent(in) :: pdel_fac
  real(r8), intent(in) :: xfertend_num(:,:)

  !intent inouts
  real(r8), intent(inout) :: dqdt(:,:,:)
  real(r8), intent(inout) :: qsrflx(:,:,:,:)

  !local
  real(r8) :: xfertend

  xfertend = xfertend_num(jmode,aer_type)
  dqdt(icol,klev,lsfrm) = dqdt(icol,klev,lsfrm) - xfertend
  dqdt(icol,klev,lstoo) = dqdt(icol,klev,lstoo) + xfertend
  qsrflx(icol,lsfrm,jsrflx,aer_type) = qsrflx(icol,lsfrm,jsrflx,aer_type) - xfertend*pdel_fac
  qsrflx(icol,lstoo,jsrflx,aer_type) = qsrflx(icol,lstoo,jsrflx,aer_type) + xfertend*pdel_fac

end subroutine update_num_tends


!---------------------------------------------------------------------------------------------

subroutine update_cld_brn_mmr(list_idx, top_lev, ncol, lchnk, pcnst, deltat, pbuf, dotendqqcw, dqqcwdt)

  !-----------------------------------------------------------------------------
  !Purpose: updates mmr of cloud borne aerosols
  !
  !Called by: modal_aero_calcsize_sub
  !Calls    : None
  !
  !Author: Richard Easter (Refactored by Balwinder Singh)
  !-----------------------------------------------------------------------------

  implicit none

  !input
  integer,  intent(in) :: list_idx, top_lev !for model level loop
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
  character(len = 1000) :: err_msg

  !check if list_idx > 0
  !This subroutine should ONLY be called for list_idx=0 as we are using qqcw_get_field
  if(list_idx > 0) then
     write(err_msg,*)'updating mmrs is not supported for list_idx>0, current list_idx is:',list_idx,' ',errmsg(__FILE__,__LINE__)
     call endrun(trim(err_msg))
  endif

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
  !Author: Richard Easter (Refactored by Balwinder Singh)
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

#endif

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
   real(r8), optional, target, allocatable, intent(inout)  :: dgnum_m(:,:,:) ! interstital aerosol dry number mode radius (m)

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
   real(r8) :: v2nmax, v2nmin           ! voltonumblo/hi of current mode
   real(r8) :: sigmag, alnsg
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   list_idx = 0  ! climate list by default
   if (present(list_idx_in)) list_idx = list_idx_in

   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (present(list_idx_in)) then
      if (.not. present(dgnum_m)) then
         call endrun('modal_aero_calcsize_diag called with'// &
              'list_idx_in but dgnum_m is not present '//errmsg(__FILE__,__LINE__))
      end if
      if (.not. allocated(dgnum_m)) then
         call endrun('modal_aero_calcsize_diag called with'// &
              'list_idx_in but dgnum_m is not allocated '//errmsg(__FILE__,__LINE__))
      endif
   end if

   do n = 1, nmodes

      if (.not.present(dgnum_m)) then
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
      v2nmin = voltonumbhi
      v2nmax = voltonumblo
      dgnxx = dgnumhi
      dgnyy = dgnumlo

      do k = top_lev, pver
         do i = 1, ncol

            drv_a = dryvol_a(i,k)
            num_a0 = mode_num(i,k)
            num_a = max( 0.0_r8, num_a0 )

            if (drv_a > 0.0_r8) then
               if (num_a <= drv_a*v2nmin) then
                  dgncur_a(i,k) = dgnxx
               else if (num_a >= drv_a*v2nmax) then
                  dgncur_a(i,k) = dgnyy
               else
                  dgncur_a(i,k) = (drv_a/(dumfac*num_a))**third
               end if
            end if

         end do
      end do

   end do ! nmodes

end subroutine modal_aero_calcsize_diag


end module modal_aero_calcsize
