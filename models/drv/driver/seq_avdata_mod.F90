!===============================================================================
! SVN $Id: seq_avdata_mod.F90 45286 2013-03-26 18:17:04Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq4_2_33/driver/seq_avdata_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_avdata_mod -- provides use access to public cpl7 aVect, domain,
!     and fraction data.
!
! !DESCRIPTION:
!
!    provides use access to public cpl7 aVect, domain, and fraction data.
!
! !REMARKS:
!
!    use access to public cpl7 aVect, domain, and fraction info is to avoid 
!    excessively long routine arg lists, eg. for history & restart modules.
!    Note: while cpl7's non-main program ("driver") routines CAN access this
!    data by use'ing this module, they SHOULD access it via agrument lists 
!    if it is reasonable to do so.  Do the right thing.
!
! !REVISION HISTORY:
!     2009-Sep-25 - B. Kauffman - initial version
!
! !INTERFACE: ------------------------------------------------------------------

module seq_avdata_mod

! !USES:

   use shr_kind_mod  ,only: IN => SHR_KIND_IN
   use mct_mod           ! mct_ wrappers for mct lib

   use seq_cdata_mod     ! "cdata" type & methods (domain + decomp + infodata in one datatype)
   use seq_infodata_mod  ! "infodata" gathers various control flags into one datatype
   use seq_comm_mct, only : num_inst_atm, num_inst_lnd, num_inst_ocn, &
                            num_inst_ice, num_inst_glc, num_inst_rof, &
                            num_inst_wav

   implicit none

   public  ! default is public
   save

! !PUBLIC DATA MEMBERS:

   !----------------------------------------------------------------------------
   ! Infodata: inter-model control flags, domain info
   !----------------------------------------------------------------------------

   type (seq_infodata_type), target :: infodata ! single instance for cpl and all comps

   !----------------------------------------------------------------------------
   ! cdata types: contains pointers to domain info + component ID + infobuffer
   !----------------------------------------------------------------------------

   type (seq_cdata) :: cdata_aa(num_inst_atm)    ! on component pes
   type (seq_cdata) :: cdata_ll(num_inst_lnd)
   type (seq_cdata) :: cdata_oo(num_inst_ocn)
   type (seq_cdata) :: cdata_ii(num_inst_ice)
   type (seq_cdata) :: cdata_rr(num_inst_rof)
   type (seq_cdata) :: cdata_gg(num_inst_glc)
   type (seq_cdata) :: cdata_ss(num_inst_lnd)
   type (seq_cdata) :: cdata_ww(num_inst_wav)

   type (seq_cdata) :: cdata_ax    ! on cpl pes
   type (seq_cdata) :: cdata_lx
   type (seq_cdata) :: cdata_ox
   type (seq_cdata) :: cdata_ix
   type (seq_cdata) :: cdata_rx
   type (seq_cdata) :: cdata_gx
   type (seq_cdata) :: cdata_sx
   type (seq_cdata) :: cdata_wx

   !----------------------------------------------------------------------------
   ! domain info: coords, fractions, decomps, area correction factors
   !----------------------------------------------------------------------------

   !--- domain coords, area, mask  (MCT General Grids) --

   type(mct_gGrid), target :: dom_aa(num_inst_atm)      ! atm domain on atm pes
   type(mct_gGrid), target :: dom_ll(num_inst_lnd)      ! lnd domain
   type(mct_gGrid), target :: dom_ii(num_inst_ice)      ! ice domain
   type(mct_gGrid), target :: dom_oo(num_inst_ocn)      ! ocn domain
   type(mct_gGrid), target :: dom_rr(num_inst_rof)      ! runoff domain
   type(mct_gGrid), target :: dom_gg(num_inst_glc)      ! glc domain
   type(mct_gGrid), target :: dom_ss(num_inst_lnd)      ! sno domain
   type(mct_gGrid), target :: dom_ww(num_inst_wav)      ! wav domain

   type(mct_gGrid), target :: dom_ax      ! atm domain on cpl pes
   type(mct_gGrid), target :: dom_lx      ! lnd domain
   type(mct_gGrid), target :: dom_ix      ! ice domain
   type(mct_gGrid), target :: dom_ox      ! ocn domain
   type(mct_gGrid), target :: dom_rx      ! runoff domain
   type(mct_gGrid), target :: dom_gx      ! glc domain
   type(mct_gGrid), target :: dom_sx      ! sno domain
   type(mct_gGrid), target :: dom_wx      ! wav domain

   !--- domain fractions (only defined on cpl pes) ---

   type(mct_aVect),pointer :: fractions_ax(:)   ! Fractions on atm grid
   type(mct_aVect),pointer :: fractions_lx(:)   ! Fractions on lnd grid 
   type(mct_aVect),pointer :: fractions_ix(:)   ! Fractions on ice grid
   type(mct_aVect),pointer :: fractions_ox(:)   ! Fractions on ocn grid
   type(mct_aVect),pointer :: fractions_gx(:)   ! Fractions on glc grid
   type(mct_aVect),pointer :: fractions_rx(:)   ! Fractions on rof grid
   type(mct_aVect),pointer :: fractions_wx(:)   ! Fractions on wav grid

   !----------------------------------------------------------------------------
   ! State/flux field bundles (MCT attribute vectors)
   !----------------------------------------------------------------------------

   type(mct_aVect) :: x2a_aa(num_inst_atm)    ! Atm import, atm grid, atm pes - allocated in atm gc
   type(mct_aVect) :: a2x_aa(num_inst_atm)    ! Atm export, atm grid, atm pes - allocated in atm gc

   type(mct_aVect) :: x2a_ax(num_inst_atm)    ! Atm import, atm grid, cpl pes - allocated in driver
   type(mct_aVect) :: a2x_ax(num_inst_atm)    ! Atm export, atm grid, cpl pes - allocated in driver
   type(mct_aVect) :: a2x_lx(num_inst_atm)    ! Atm export, lnd grid, cpl pes - allocated in driver
   type(mct_aVect) :: a2x_ix(num_inst_atm)    ! Atm export, ice grid, cpl pes - allocated in driver
   type(mct_aVect) :: a2x_ox(num_inst_atm)    ! Atm export, ocn grid, cpl pes - allocated in driver
   type(mct_aVect) :: a2x_wx(num_inst_atm)    ! Atm export, wav grid, cpl pes - allocated in driver

   type(mct_aVect) :: x2l_ll(num_inst_lnd)    ! Lnd import, lnd grid, lnd pes - allocated in lnd gc
   type(mct_aVect) :: l2x_ll(num_inst_lnd)    ! Lnd export, lnd grid, lnd pes - allocated in lnd gc

   type(mct_aVect) :: x2l_lx(num_inst_lnd)    ! Lnd import, lnd grid, cpl pes - allocated in driver
   type(mct_aVect) :: l2x_lx(num_inst_lnd)    ! Lnd export, lnd grid, cpl pes - allocated in driver
   type(mct_aVect) :: l2x_ax(num_inst_lnd)    ! Lnd export, atm grid, cpl pes - allocated in driver
   type(mct_aVect) :: x2racc_lx(num_inst_lnd) ! Lnd export, lnd grid, cpl pes - allocated in driver

   type(mct_aVect) :: x2r_rr(num_inst_rof)    ! Rof import, rof grid, lnd pes - allocated in lnd gc
   type(mct_aVect) :: r2x_rr(num_inst_rof)    ! Rof export, rof grid, lnd pes - allocated in lnd gc

   type(mct_aVect) :: x2r_rx(num_inst_rof)    ! Rof import, rof grid, lnd pes - allocated in lnd gc
   type(mct_aVect) :: r2x_rx(num_inst_rof)    ! Rof export, rof grid, cpl pes - allocated in driver
   type(mct_aVect) :: r2xacc_rx(num_inst_rof) ! Rof export, rof grid, cpl pes - allocated in driver
   type(mct_aVect) :: r2x_lx(num_inst_rof)    ! Rof export, lnd grid, lnd pes - allocated in lnd gc
   type(mct_aVect) :: r2x_ox(num_inst_rof)    ! Rof export, ocn grid, cpl pes - allocated in driver

   type(mct_aVect) :: x2s_ss(num_inst_lnd)    ! Sno import, sno grid, sno pes - allocated in lnd gc
   type(mct_aVect) :: s2x_ss(num_inst_lnd)    ! Sno export, sno grid, sno pes - allocated in lnd gc

   type(mct_aVect) :: x2s_sx(num_inst_lnd)    ! Sno import, sno grid, cpl pes - allocated in driver
   type(mct_aVect) :: s2x_sx(num_inst_lnd)    ! Sno export, sno grid, cpl pes - allocated in driver
   type(mct_aVect) :: s2x_gx(num_inst_lnd)    ! Sno export, glc grid, cpl pes - allocated in driver

   type(mct_aVect) :: x2i_ii(num_inst_ice)    ! Ice import, ice grid, ice pes - allocated in ice gc
   type(mct_aVect) :: i2x_ii(num_inst_ice)    ! Ice export, ice grid, ice pes - allocated in ice gc

   type(mct_aVect) :: x2i_ix(num_inst_ice)    ! Ice import, ice grid, cpl pes - allocated in driver
   type(mct_aVect) :: i2x_ix(num_inst_ice)    ! Ice export, ice grid, cpl pes - allocated in driver
   type(mct_aVect) :: i2x_ax(num_inst_ice)    ! Ice export, atm grid, cpl pes - allocated in driver
   type(mct_aVect) :: i2x_ox(num_inst_ice)    ! Ice export, ocn grid, cpl pes - allocated in driver
   type(mct_aVect) :: i2x_wx(num_inst_ice)    ! Ice export, wav grid, cpl pes - allocated in driver

   type(mct_aVect) :: x2o_oo(num_inst_ocn)    ! Ocn import, ocn grid, ocn pes - allocated in ocn gc
   type(mct_aVect) :: o2x_oo(num_inst_ocn)    ! Ocn export, ocn grid, ocn pes - allocated in ocn gc

   type(mct_aVect) :: x2o_ox(num_inst_ocn)    ! Ocn import, ocn grid, cpl pes - allocated in driver
   type(mct_aVect) :: x2oacc_ox(num_inst_ocn) ! Ocn import, ocn grid, cpl pes - allocated in driver
   type(mct_aVect) :: o2x_ox(num_inst_ocn)    ! Ocn export, ocn grid, cpl pes - allocated in driver
   type(mct_aVect) :: o2x_ax(num_inst_ocn)    ! Ocn export, atm grid, cpl pes - allocated in driver
   type(mct_aVect) :: o2x_ix(num_inst_ocn)    ! Ocn export, ice grid, cpl pes - allocated in driver
   type(mct_aVect) :: o2x_wx(num_inst_ocn)    ! Ocn export, wav grid, cpl pes - allocated in driver

   type(mct_aVect),pointer :: xao_ox(:)       ! Atm-ocn fluxes, ocn grid, cpl pes - allocated in flux_ao gc 
   type(mct_aVect),pointer :: xao_ax(:)       ! Atm-ocn fluxes, atm grid, cpl pes - allocated in flux_ao gc 

   type(mct_aVect) :: x2g_gg(num_inst_glc)    ! Glc import, glc grid, ice pes - allocated in glc gc
   type(mct_aVect) :: g2x_gg(num_inst_glc)    ! Glc export, glc grid, ice pes - allocated in glc gc

   type(mct_aVect) :: x2g_gx(num_inst_glc)    ! Glc import, glc grid, cpl pes - allocated in driver
   type(mct_aVect) :: g2x_gx(num_inst_glc)    ! Glc export, glc grid, cpl pes - allocated in driver
   type(mct_aVect) :: g2x_sx(num_inst_glc)    ! Glc export, sno grid, cpl pes - allocated in driver

   type(mct_aVect) :: x2w_ww(num_inst_wav)    ! Wav import, wav grid, ice pes - allocated in wav gc
   type(mct_aVect) :: w2x_ww(num_inst_wav)    ! Wav export, wav grid, ice pes - allocated in wav gc

   type(mct_aVect) :: x2w_wx(num_inst_wav)    ! Wav import, wav grid, cpl pes - allocated in driver
   type(mct_aVect) :: w2x_wx(num_inst_wav)    ! Wav export, wav grid, cpl pes - allocated in driver
   type(mct_aVect) :: w2x_ox(num_inst_wav)    ! Wav export, ocn grid, cpl pes - allocated in driver

   integer(IN)     :: r2xacc_rx_cnt ! r2xacc_rx: number of time samples accumulated
   integer(IN)     :: x2oacc_ox_cnt ! x2oacc_ox: number of time samples accumulated
   integer(IN)     :: x2racc_lx_cnt ! x2racc_lx: number of time samples accumulated

! !PUBLIC MEMBER FUNCTIONS

   ! no public routines

! !PUBLIC TYPES:

   ! no public types

end module seq_avdata_mod

!===============================================================================
