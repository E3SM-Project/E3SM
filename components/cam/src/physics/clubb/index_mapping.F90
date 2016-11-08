!---------------------------------------------------------------------------
! $Id: index_mapping.F90 7118 2014-07-25 00:12:15Z raut@uwm.edu $
!===============================================================================
module index_mapping

  ! Description:
  ! Functions to map back and forth between the PDF arrays and the hydrometeor
  ! arrays.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  ! Hydrometeor array indices
  use array_index, only: &
      iirrm, & ! Hydrometeor array index for rain water mixing ratio, rr
      iirsm, & ! Hydrometeor array index for snow mixing ratio, rs
      iirim, & ! Hydrometeor array index for ice mixing ratio, ri
      iirgm, & ! Hydrometeor array index for graupel mixing ratio, rg
      iiNrm, & ! Hydrometeor array index for rain drop concentration, Nr
      iiNsm, & ! Hydrometeor array index for snow concentration, Ns
      iiNim, & ! Hydrometeor array index for ice concentration, Ni
      iiNgm    ! Hydrometeor array index for graupel concentration, Ng

  ! PDF array indices
  use corr_varnce_module, only: &
      iiPDF_rr, & ! PDF array index for rain water mixing ratio, rr
      iiPDF_rs, & ! PDF array index for snow mixing ratio, rs
      iiPDF_ri, & ! PDF array index for ice mixing ratio, ri
      iiPDF_rg, & ! PDF array index for graupel mixing ratio, rg
      iiPDF_Nr, & ! PDF array index for rain drop concentration, Nr
      iiPDF_Ns, & ! PDF array index for snow concentration, Ns
      iiPDF_Ni, & ! PDF array index for ice concentration, Ni
      iiPDF_Ng    ! PDF array index for graupel concentration, Ng

  implicit none

  private ! Default Scope

  public :: pdf2hydromet_idx, &
            hydromet2pdf_idx, &
            rx2Nx_hm_idx,     &
            Nx2rx_hm_idx,     &
            mvr_hm_max

contains

  !=============================================================================
  function pdf2hydromet_idx( pdf_idx ) result( hydromet_idx )

    ! Description:
    ! Returns the position of a specific precipitating hydrometeor corresponding
    ! to the PDF index (pdf_idx) in the precipitating hydrometeor array
    ! (hydromet_idx).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      pdf_idx    ! Index of a hydrometeor in the PDF array.

    ! Return Variable
    integer :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.


    ! Initialize hydromet_idx
    hydromet_idx = 0

    if ( pdf_idx == iiPDF_rr ) then

       ! Index for rain water mixing ratio, rr.
       hydromet_idx = iirrm

    elseif ( pdf_idx == iiPDF_Nr ) then

       ! Index for rain drop concentration, Nr.
       hydromet_idx = iiNrm

    elseif ( pdf_idx == iiPDF_rs ) then

       ! Index for snow mixing ratio, rs.
       hydromet_idx = iirsm

    elseif ( pdf_idx == iiPDF_Ns ) then

       ! Index for snow flake concentration, Ns.
       hydromet_idx = iiNsm

    elseif ( pdf_idx == iiPDF_rg ) then

       ! Index for graupel mixing ratio, rg.
       hydromet_idx = iirgm

    elseif ( pdf_idx == iiPDF_Ng ) then

       ! Index for graupel concentration, Ng.
       hydromet_idx = iiNgm

    elseif ( pdf_idx == iiPDF_ri ) then

       ! Index for ice mixing ratio, ri.
       hydromet_idx = iirim

    elseif ( pdf_idx == iiPDF_Ni ) then

       ! Index for ice concentration, Ni.
       hydromet_idx = iiNim

    endif


    return

  end function pdf2hydromet_idx

  !=============================================================================
  function hydromet2pdf_idx( hydromet_idx ) result( pdf_idx )

    ! Description:
    ! Returns the position of a specific precipitating hydrometeor corresponding
    ! to the precipitating hydrometeor index (hydromet_idx) in the PDF array
    ! (pdf_idx).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.

    ! Return Variable
    integer :: &
      pdf_idx    ! Index of a hydrometeor in the PDF array.


    ! Initialize pdf_idx.
    pdf_idx = 0

    if ( hydromet_idx == iirrm ) then

       ! Index for rain water mixing ratio, rr.
       pdf_idx = iiPDF_rr

    elseif ( hydromet_idx == iiNrm ) then

       ! Index for rain drop concentration, Nr.
       pdf_idx = iiPDF_Nr

    elseif ( hydromet_idx == iirim ) then

       ! Index for ice mixing ratio, ri.
       pdf_idx = iiPDF_ri

    elseif ( hydromet_idx == iiNim ) then

       ! Index for ice concentration, Ni.
       pdf_idx = iiPDF_Ni

    elseif ( hydromet_idx == iirsm ) then

       ! Index for snow mixing ratio, rs.
       pdf_idx = iiPDF_rs

    elseif ( hydromet_idx == iiNsm ) then

       ! Index for snow flake concentration, Ns.
       pdf_idx = iiPDF_Ns

    elseif ( hydromet_idx == iirgm ) then

       ! Index for graupel mixing ratio, rg.
       pdf_idx = iiPDF_rg

    elseif ( hydromet_idx == iiNgm ) then

       ! Index for graupel concentration, Ng.
       pdf_idx = iiPDF_Ng

    endif


    return

  end function hydromet2pdf_idx

  !=============================================================================
  function rx2Nx_hm_idx( rx_idx ) result( Nx_idx )

    ! Description:
    ! Returns the position in the hydrometeor array of the specific
    ! precipitating hydrometeor concentration (Nx_idx) corresponding to the
    ! precipitating hydrometeor mixing ratio (rx_idx) of the same species of
    ! precipitating hydrometeor (rain, ice, snow, or graupel).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      rx_idx    ! Index of the mixing ratio in the hydrometeor array.

    ! Return Variable
    integer :: &
      Nx_idx    ! Index of the concentration in the hydrometeor array.


    ! Initialize Nx_idx.
    Nx_idx = 0

    if ( rx_idx == iirrm ) then

       ! Index for rain drop concentration, Nr.
       Nx_idx = iiNrm

    elseif ( rx_idx == iirim ) then

       ! Index for ice crystal concentration, Ni.
       Nx_idx = iiNim

    elseif ( rx_idx == iirsm ) then

       ! Index for snow flake concentration, Ns.
       Nx_idx = iiNsm

    elseif ( rx_idx == iirgm ) then

       ! Index for graupel concentration, Ng.
       Nx_idx = iiNgm

    endif


    return

  end function rx2Nx_hm_idx

  !=============================================================================
  function Nx2rx_hm_idx( Nx_idx ) result( rx_idx )

    ! Description:
    ! Returns the position in the hydrometeor array of the specific
    ! precipitating hydrometeor mixing ratio (rx_idx) corresponding to the
    ! precipitating hydrometeor concentration (Nx_idx) of the same species of
    ! precipitating hydrometeor (rain, ice, snow, or graupel).

    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      Nx_idx    ! Index of the concentration in the hydrometeor array.

    ! Return Variable
    integer :: &
      rx_idx    ! Index of the mixing ratio in the hydrometeor array.


    ! Initialize rx_idx.
    rx_idx = 0

    if ( Nx_idx == iiNrm ) then

       ! Index for rain water mixing ratio, rr.
       rx_idx = iirrm

    elseif ( Nx_idx == iiNim ) then

       ! Index for ice mixing ratio, ri.
       rx_idx = iirim

    elseif ( Nx_idx == iiNsm ) then

       ! Index for snow mixing ratio, rs.
       rx_idx = iirsm

    elseif ( Nx_idx == iiNgm ) then

       ! Index for graupel mixing ratio, rg.
       rx_idx = iirgm

    endif


    return

  end function Nx2rx_hm_idx

  !=============================================================================
  function mvr_hm_max( hydromet_idx ) result( mvr_hydromet_max )

    ! Description:
    ! Returns the maximum allowable mean volume radius of a specific
    ! precipitating hydrometeor type (rain, ice, snow, or graupel) corresponding
    ! to the precipitating hydrometeor index, whether that index is for the
    ! mixing ratio or concentration associated with that hydrometeor type.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        mvr_rain_max,    & ! Constant(s)
        mvr_ice_max,     &
        mvr_snow_max,    &
        mvr_graupel_max, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      hydromet_idx    ! Index of a hydrometeor in the hydromet array.

    ! Return Variable
    real( kind = core_rknd ) :: &
      mvr_hydromet_max    ! Maximum allowable mean volume radius    [m]


    ! Initialize mvr_hydromet_max.
    mvr_hydromet_max = zero

    if ( hydromet_idx == iirrm .or. hydromet_idx == iiNrm ) then

       ! Maximum allowable mean volume radius for rain drops.
       mvr_hydromet_max = mvr_rain_max

    elseif ( hydromet_idx == iirim .or. hydromet_idx == iiNim ) then

       ! Maximum allowable mean volume radius for ice crystals.
       mvr_hydromet_max = mvr_ice_max

    elseif ( hydromet_idx == iirsm .or. hydromet_idx == iiNsm ) then

       ! Maximum allowable mean volume radius for snow flakes.
       mvr_hydromet_max = mvr_snow_max

    elseif ( hydromet_idx == iirgm .or. hydromet_idx == iiNgm ) then

       ! Maximum allowable mean volume radius for graupel.
       mvr_hydromet_max = mvr_graupel_max

    endif


    return

  end function mvr_hm_max

!===============================================================================

end module index_mapping
