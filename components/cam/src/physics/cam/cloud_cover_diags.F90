
!===============================================================================
! cloud cover output
!===============================================================================
module cloud_cover_diags

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver,pverp
  use cam_history,   only: addfld, horiz_only, add_default, outfld
  use phys_control,  only: phys_getopts

  implicit none

  private

  public :: cloud_cover_diags_init
  public :: cloud_cover_diags_out

contains

!===============================================================================
!===============================================================================
subroutine cloud_cover_diags_init(sampling_seq)

  character(len=*), intent(in) :: sampling_seq
  logical :: history_amwg         ! output the variables used by the AMWG diag package

  call addfld ('CLOUD',(/ 'lev' /), 'A','fraction','Cloud fraction'                        , sampling_seq=sampling_seq)
  call addfld ('CLDTOT',horiz_only,    'A','fraction','Vertically-integrated total cloud'     , sampling_seq=sampling_seq)
  call addfld ('CLDLOW',horiz_only,    'A','fraction','Vertically-integrated low cloud'       , sampling_seq=sampling_seq)
  call addfld ('CLDMED',horiz_only,    'A','fraction','Vertically-integrated mid-level cloud' , sampling_seq=sampling_seq)
  call addfld ('CLDHGH',horiz_only,    'A','fraction','Vertically-integrated high cloud'      , sampling_seq=sampling_seq)

  ! determine the add_default fields
  call phys_getopts(history_amwg_out           = history_amwg  )
 
  if (history_amwg) then
      call add_default ('CLOUD   ', 1, ' ')
      call add_default ('CLDTOT  ', 1, ' ')
      call add_default ('CLDLOW  ', 1, ' ')
      call add_default ('CLDMED  ', 1, ' ')
      call add_default ('CLDHGH  ', 1, ' ')
  endif
    

end subroutine cloud_cover_diags_init

!===============================================================================
!===============================================================================
subroutine cloud_cover_diags_out(lchnk, ncol, cld, pmid, nmxrgn, pmxrgn )

  integer,  intent(in) :: lchnk, ncol
  real(r8), intent(in) :: cld(pcols,pver)
  real(r8), intent(in) :: pmid(pcols,pver)
  integer,  intent(in) :: nmxrgn(pcols)
  real(r8), intent(in) :: pmxrgn(pcols,pverp)

  real(r8) :: cltot(pcols)            ! Diagnostic total cloud cover
  real(r8) :: cllow(pcols)            !       "     low  cloud cover
  real(r8) :: clmed(pcols)            !       "     mid  cloud cover
  real(r8) :: clhgh(pcols)            !       "     hgh  cloud cover

  call cldsav (lchnk, ncol, cld, pmid, cltot, cllow, clmed, clhgh, nmxrgn, pmxrgn)

  !
  ! Dump cloud field information to history tape buffer (diagnostics)
  !
  call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
  call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
  call outfld('CLDMED  ',clmed  ,pcols,lchnk)
  call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)

  call outfld('CLOUD   ',cld    ,pcols,lchnk) 

end subroutine cloud_cover_diags_out

!===============================================================================
!===============================================================================
subroutine cldsav(lchnk   ,ncol    , &
                  cld     ,pmid    ,cldtot  ,cldlow  ,cldmed  , &
                  cldhgh  ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
! Pressure ranges for the 3 cloud levels are specified.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

   implicit none
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: cld(pcols,pver)     ! Cloud fraction
   real(r8), intent(in) :: pmid(pcols,pver)    ! Level pressures
   real(r8), intent(in) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc

   integer, intent(in) :: nmxrgn(pcols)        ! Number of maximally overlapped regions
!
! Output arguments
!
   real(r8), intent(out) :: cldtot(pcols)       ! Total random overlap cloud cover
   real(r8), intent(out) :: cldlow(pcols)       ! Low random overlap cloud cover
   real(r8), intent(out) :: cldmed(pcols)       ! Middle random overlap cloud cover
   real(r8), intent(out) :: cldhgh(pcols)       ! High random overlap cloud cover

!
!---------------------------Local workspace-----------------------------
!
   integer i,k                  ! Longitude,level indices
   integer irgn(pcols)          ! Max-overlap region index
   integer max_nmxrgn           ! maximum value of nmxrgn over columns
   integer ityp                 ! Type counter
   real(r8) clrsky(pcols)       ! Max-random clear sky fraction
   real(r8) clrskymax(pcols)    ! Maximum overlap clear sky fraction
!------------------------------Parameters-------------------------------
   real(r8) plowmax             ! Max prs for low cloud cover range
   real(r8) plowmin             ! Min prs for low cloud cover range
   real(r8) pmedmax             ! Max prs for mid cloud cover range
   real(r8) pmedmin             ! Min prs for mid cloud cover range
   real(r8) phghmax             ! Max prs for hgh cloud cover range
   real(r8) phghmin             ! Min prs for hgh cloud cover range
!
   parameter (plowmax = 120000._r8,plowmin = 70000._r8, &
              pmedmax =  70000._r8,pmedmin = 40000._r8, &
              phghmax =  40000._r8,phghmin =  5000._r8)

   real(r8) ptypmin(4)
   real(r8) ptypmax(4)

   data ptypmin /phghmin, plowmin, pmedmin, phghmin/
   data ptypmax /plowmax, plowmax, pmedmax, phghmax/
!
!-----------------------------------------------------------------------
!
! Initialize region number
!
   max_nmxrgn = -1
   do i=1,ncol
      max_nmxrgn = max(max_nmxrgn,nmxrgn(i))
   end do

   do ityp = 1, 4
      irgn(1:ncol) = 1
      do k =1,max_nmxrgn-1
         do i=1,ncol
            if (pmxrgn(i,irgn(i)) < ptypmin(ityp) .and. irgn(i) < nmxrgn(i)) then
               irgn(i) = irgn(i) + 1
            end if
         end do
      end do
!
! Compute cloud amount by estimating clear-sky amounts
!
      clrsky(1:ncol)    = 1.0_r8
      clrskymax(1:ncol) = 1.0_r8
      do k = 1, pver
         do i=1,ncol
            if (pmid(i,k) >= ptypmin(ityp) .and. pmid(i,k) <= ptypmax(ityp)) then
               if (pmxrgn(i,irgn(i)) < pmid(i,k) .and. irgn(i) < nmxrgn(i)) then
                  irgn(i) = irgn(i) + 1
                  clrsky(i) = clrsky(i) * clrskymax(i)
                  clrskymax(i) = 1.0_r8
               endif
               clrskymax(i) = min(clrskymax(i),1.0_r8-cld(i,k))
            endif
         end do
      end do
      if (ityp == 1) cldtot(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 2) cldlow(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 3) cldmed(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 4) cldhgh(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
   end do

   return
end subroutine cldsav

end module cloud_cover_diags
