!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module index_mapping

  ! Description:
  ! Functions to map back and forth between the PDF arrays and the hydrometeor
  ! arrays.

  ! The “iiPDF” indices are used to index all PDF variates, including all
  ! hydrometeor variates.  
  ! The “ii” indices are used to index hydrometeor arrays.  
  ! The “ii” variates are a subset of the “iiPDF” variates.  
  ! Conversions between the two sets of indices are done by the 
  ! functions pdf2hydromet_idx and hydromet2pdf_idx below.
  ! 
  ! ------------------------------------------------------------------------
  ! 
  ! iiPDF indices:
  ! 
  ! Included indices:  
  ! iiPDF_chi, iiPDF_eta, iiPDF_w, iiPDF_Ncn, iiPDF_rr, & all other hydrometeors
  ! 
  ! Number of indices:  pdf_dim
  ! 
  ! Examples of arrays dimensioned by pdf_dim:
  ! mu_x_1_n, corr_array_n_cloud, . . .
  ! 
  ! Declared as module variables in module array_index
  ! 
  ! Initialized in subroutine setup_pdf_indices
  ! 
  ! ----------------------------------------------------------------------
  ! 
  ! ii indices:
  ! 
  ! Included indices:  
  ! iirr, iiNr, iiri, iiNi, iirs, iiNs, iirg, iiNg
  ! 
  ! Number of indices:  hydromet_dim
  ! 
  ! Examples of arrays dimensioned by hydromet_dim: 
  ! hydromet, wphydrometp, . . .
  ! 
  ! Declared as module variables in module array_index.
  ! 
  ! Initialized in subroutine init_microphys
  ! 
  ! -----------------------------------------------------------------------
  !
  ! References:
  !   None
  !-------------------------------------------------------------------------

  ! Hydrometeor array indices
  use array_index, only: &
      iirr, & ! Hydrometeor array index for rain water mixing ratio, rr
      iirs, & ! Hydrometeor array index for snow mixing ratio, rs
      iiri, & ! Hydrometeor array index for ice mixing ratio, ri
      iirg, & ! Hydrometeor array index for graupel mixing ratio, rg
      iiNr, & ! Hydrometeor array index for rain drop concentration, Nr
      iiNs, & ! Hydrometeor array index for snow concentration, Ns
      iiNi, & ! Hydrometeor array index for ice concentration, Ni
      iiNg, &    ! Hydrometeor array index for graupel concentration, Ng
  ! PDF array indices
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
       hydromet_idx = iirr

    elseif ( pdf_idx == iiPDF_Nr ) then

       ! Index for rain drop concentration, Nr.
       hydromet_idx = iiNr

    elseif ( pdf_idx == iiPDF_rs ) then

       ! Index for snow mixing ratio, rs.
       hydromet_idx = iirs

    elseif ( pdf_idx == iiPDF_Ns ) then

       ! Index for snow flake concentration, Ns.
       hydromet_idx = iiNs

    elseif ( pdf_idx == iiPDF_rg ) then

       ! Index for graupel mixing ratio, rg.
       hydromet_idx = iirg

    elseif ( pdf_idx == iiPDF_Ng ) then

       ! Index for graupel concentration, Ng.
       hydromet_idx = iiNg

    elseif ( pdf_idx == iiPDF_ri ) then

       ! Index for ice mixing ratio, ri.
       hydromet_idx = iiri

    elseif ( pdf_idx == iiPDF_Ni ) then

       ! Index for ice concentration, Ni.
       hydromet_idx = iiNi

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

    if ( hydromet_idx == iirr ) then

       ! Index for rain water mixing ratio, rr.
       pdf_idx = iiPDF_rr

    elseif ( hydromet_idx == iiNr ) then

       ! Index for rain drop concentration, Nr.
       pdf_idx = iiPDF_Nr

    elseif ( hydromet_idx == iiri ) then

       ! Index for ice mixing ratio, ri.
       pdf_idx = iiPDF_ri

    elseif ( hydromet_idx == iiNi ) then

       ! Index for ice concentration, Ni.
       pdf_idx = iiPDF_Ni

    elseif ( hydromet_idx == iirs ) then

       ! Index for snow mixing ratio, rs.
       pdf_idx = iiPDF_rs

    elseif ( hydromet_idx == iiNs ) then

       ! Index for snow flake concentration, Ns.
       pdf_idx = iiPDF_Ns

    elseif ( hydromet_idx == iirg ) then

       ! Index for graupel mixing ratio, rg.
       pdf_idx = iiPDF_rg

    elseif ( hydromet_idx == iiNg ) then

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

    if ( rx_idx == iirr ) then

       ! Index for rain drop concentration, Nr.
       Nx_idx = iiNr

    elseif ( rx_idx == iiri ) then

       ! Index for ice crystal concentration, Ni.
       Nx_idx = iiNi

    elseif ( rx_idx == iirs ) then

       ! Index for snow flake concentration, Ns.
       Nx_idx = iiNs

    elseif ( rx_idx == iirg ) then

       ! Index for graupel concentration, Ng.
       Nx_idx = iiNg

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

    if ( Nx_idx == iiNr ) then

       ! Index for rain water mixing ratio, rr.
       rx_idx = iirr

    elseif ( Nx_idx == iiNi ) then

       ! Index for ice mixing ratio, ri.
       rx_idx = iiri

    elseif ( Nx_idx == iiNs ) then

       ! Index for snow mixing ratio, rs.
       rx_idx = iirs

    elseif ( Nx_idx == iiNg ) then

       ! Index for graupel mixing ratio, rg.
       rx_idx = iirg

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

    if ( hydromet_idx == iirr .or. hydromet_idx == iiNr ) then

       ! Maximum allowable mean volume radius for rain drops.
       mvr_hydromet_max = mvr_rain_max

    elseif ( hydromet_idx == iiri .or. hydromet_idx == iiNi ) then

       ! Maximum allowable mean volume radius for ice crystals.
       mvr_hydromet_max = mvr_ice_max

    elseif ( hydromet_idx == iirs .or. hydromet_idx == iiNs ) then

       ! Maximum allowable mean volume radius for snow flakes.
       mvr_hydromet_max = mvr_snow_max

    elseif ( hydromet_idx == iirg .or. hydromet_idx == iiNg ) then

       ! Maximum allowable mean volume radius for graupel.
       mvr_hydromet_max = mvr_graupel_max

    endif


    return

  end function mvr_hm_max

!===============================================================================

end module index_mapping
