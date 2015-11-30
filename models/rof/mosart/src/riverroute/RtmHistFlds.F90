module RtmHistFlds

!-----------------------------------------------------------------------
! !DESCRIPTION:
! Module containing initialization of RTM history fields and files
! This is the module that the user must modify in order to add new
! history fields or modify defaults associated with existing history
! fields.
!
! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use RunoffMod      , only : rtmCTL
  use RtmHistFile    , only : RtmHistAddfld, RtmHistPrintflds
  use rof_cpl_indices, only : nt_rtm, rtm_tracers  

  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: RtmHistFldsInit 
  public :: RtmHistFldsSet
!
!------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

  subroutine RtmHistFldsInit()

    !-------------------------------------------------------
    ! DESCRIPTION:
    ! Build master field list of all possible fields in a history file.
    ! Each field has associated with it a ``long\_name'' netcdf attribute that
    ! describes what the field is, and a ``units'' attribute. A subroutine is
    ! called to add each field to the masterlist.
    !
    ! ARGUMENTS:
    implicit none
    !-------------------------------------------------------

    call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART river flow: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%runofflnd_nt1, default='active')

    call RtmHistAddfld (fname='QCHANR'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART river flow: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%runofflnd_nt2, default='active')

    call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(1)), units='m3/s', &
         avgflag='A', long_name='MOSART river discharge into ocean: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%runoffocn_nt1, default='active')

    call RtmHistAddfld (fname='QCHOCNR'//'_'//trim(rtm_tracers(2)), units='m3/s', &
         avgflag='A', long_name='MOSART river discharge into ocean: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%runoffocn_nt2, default='active')

    call RtmHistAddfld (fname='VOLR'//'_'//trim(rtm_tracers(1)), units='m3',  &
         avgflag='A', long_name='MOSART storage: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%volr_nt1, default='active')

    call RtmHistAddfld (fname='VOLR'//'_'//trim(rtm_tracers(2)), units='m3',  &
         avgflag='A', long_name='MOSART storage: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%volr_nt2, default='active')

    call RtmHistAddfld (fname='DVOLRDT_LND'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART land change in storage: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%dvolrdtlnd_nt1, default='inactive')

    call RtmHistAddfld (fname='DVOLRDT_LND'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART land change in storage: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%dvolrdtlnd_nt2, default='inactive')

    call RtmHistAddfld (fname='DVOLRDT_OCN'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART ocean change of storage: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%dvolrdtocn_nt1, default='inactive')

    call RtmHistAddfld (fname='DVOLRDT_OCN'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART ocean change of storage: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%dvolrdtocn_nt2, default='inactive')

    call RtmHistAddfld (fname='QSUR'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART input surface runoff: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%qsur_nt1, default='active')

    call RtmHistAddfld (fname='QSUR'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART input surface runoff: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%qsur_nt2, default='active')

    call RtmHistAddfld (fname='QSUB'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART input subsurface runoff: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%qsub_nt1, default='active')

    call RtmHistAddfld (fname='QSUB'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART input subsurface runoff: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%qsub_nt2, default='active')

    call RtmHistAddfld (fname='QGWL'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART input GWL runoff: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%qgwl_nt1, default='active')

    call RtmHistAddfld (fname='QGWL'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART input GWL runoff: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%qgwl_nt2, default='active')

    call RtmHistAddfld (fname='QDTO'//'_'//trim(rtm_tracers(1)), units='m3/s',  &
         avgflag='A', long_name='MOSART input direct to ocean runoff: '//trim(rtm_tracers(1)), &
         ptr_rof=rtmCTL%qdto_nt1, default='active')

    call RtmHistAddfld (fname='QDTO'//'_'//trim(rtm_tracers(2)), units='m3/s',  &
         avgflag='A', long_name='MOSART input direct to ocean runoff: '//trim(rtm_tracers(2)), &
         ptr_rof=rtmCTL%qdto_nt2, default='active')

    ! Print masterlist of history fields

    call RtmHistPrintflds()

  end subroutine RtmHistFldsInit

!-----------------------------------------------------------------------

  subroutine RtmHistFldsSet()

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set mosart history fields as 1d poitner arrays
    !
    implicit none
    !-----------------------------------------------------------------------

    ! Currently only have two tracers

    rtmCTL%runofflnd_nt1(:)  = rtmCTL%runofflnd(:,1)
    rtmCTL%runofflnd_nt2(:)  = rtmCTL%runofflnd(:,2)

    rtmCTL%runoffocn_nt1(:)  = rtmCTL%runoffocn(:,1)
    rtmCTL%runoffocn_nt2(:)  = rtmCTL%runoffocn(:,2)

    rtmCTL%dvolrdtlnd_nt1(:) = rtmCTL%dvolrdtlnd(:,1)
    rtmCTL%dvolrdtlnd_nt2(:) = rtmCTL%dvolrdtlnd(:,2)

    rtmCTL%dvolrdtocn_nt1(:) = rtmCTL%dvolrdtocn(:,1)
    rtmCTL%dvolrdtocn_nt2(:) = rtmCTL%dvolrdtocn(:,2)

    rtmCTL%volr_nt1(:)       = rtmCTL%volr(:,1)
    rtmCTL%volr_nt2(:)       = rtmCTL%volr(:,2)

    rtmCTL%qsub_nt1(:)       = rtmCTL%qsub(:,1)
    rtmCTL%qsub_nt2(:)       = rtmCTL%qsub(:,2)

    rtmCTL%qsur_nt1(:)       = rtmCTL%qsur(:,1)
    rtmCTL%qsur_nt2(:)       = rtmCTL%qsur(:,2)

    rtmCTL%qgwl_nt1(:)       = rtmCTL%qgwl(:,1)
    rtmCTL%qgwl_nt2(:)       = rtmCTL%qgwl(:,2)

    rtmCTL%qdto_nt1(:)       = rtmCTL%qdto(:,1)
    rtmCTL%qdto_nt2(:)       = rtmCTL%qdto(:,2)

  end subroutine RtmHistFldsSet


end module RtmHistFlds
