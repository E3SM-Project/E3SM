
MODULE MOSARTinund_PreProcs_MOD

!-----------------------------------------------------------------------------
! DESCRIPTION: Pre-processing for MOSART-Inundation.
! 
! HISTORY:
! 2014-2016: Created and improved in offline MOSART-Inundation.
! 2017: Integrated with ACME.
! ... ... ...
!
!-----------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_sys_mod, only: shr_sys_abort
  use RunoffMod, only: Tctl, TUnit, rtmCTL
  use RtmVar, only: iulog
  
  implicit none
  
  real( r8 ), public, parameter :: con1Em3 = 1.0e-3_r8
  public calc_chnlMannCoe, preprocess_elevProf, interpolation_linear
  
!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine calc_chnlMannCoe ( )
    ! DESCRIPTION: Calculate channel Manning roughness coefficients.
    
    ! HISTORY:
    ! 2014-2016: Created and improved in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.   
    ! ... ... ...
    
    implicit none
    integer :: iu
    
    ! 1 -- use channel depth (Luo et al. 2017 GMD) :
    if ( Tctl%OPT_calcNr .eq. 1 ) then    
      
      do iu = rtmCTL%begr, rtmCTL%endr
        !if ( TUnit%mask( iu ) .gt. 0 ) then
        if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
        
          if ( Tctl%rdepth_min .lt. TUnit%rdepth(iu) .and. TUnit%rdepth(iu) .lt. Tctl%rdepth_max ) then
            TUnit%nr(iu) = Tctl%nr_min + (Tctl%nr_max - Tctl%nr_min) * ( ( Tctl%rdepth_max - TUnit%rdepth(iu) ) / (Tctl%rdepth_max - Tctl%rdepth_min) )
          elseif ( TUnit%rdepth(iu) .le. Tctl%rdepth_min ) then
            TUnit%nr(iu) = Tctl%nr_max
          elseif ( Tctl%rdepth_max .le. TUnit%rdepth(iu) ) then
            TUnit%nr(iu) = Tctl%nr_min
          end if
          
        end if
      enddo

    ! 2 -- use channel depth and exponent of 1/3 (Getirana et al. 2012 JHM) :
    elseif ( Tctl%OPT_calcNr .eq. 2 ) then
      
      do iu = rtmCTL%begr, rtmCTL%endr
        !if ( TUnit%mask( iu ) .gt. 0 ) then
        if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
        
          if ( Tctl%rdepth_min .lt. TUnit%rdepth(iu) .and. TUnit%rdepth(iu) .lt. Tctl%rdepth_max ) then
            TUnit%nr(iu) = Tctl%nr_min + (Tctl%nr_max - Tctl%nr_min) * ( ( Tctl%rdepth_max - TUnit%rdepth(iu) ) / (Tctl%rdepth_max - Tctl%rdepth_min) ) ** (1.0_r8/3.0_r8)
          elseif ( TUnit%rdepth(iu) .le. Tctl%rdepth_min ) then 
            TUnit%nr(iu) = Tctl%nr_max
          elseif ( Tctl%rdepth_max .le. TUnit%rdepth(iu) ) then
            TUnit%nr(iu) = Tctl%nr_min
          end if  
        
        end if
      enddo       
      
    ! 3 -- use channel width (Decharme et al. 2010 JHM) :
    elseif ( Tctl%OPT_calcNr .eq. 3 ) then
      
      do iu = rtmCTL%begr, rtmCTL%endr
        !if ( TUnit%mask( iu ) .gt. 0 ) then
        if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean). 
      
          if ( Tctl%rwidth_min .lt. TUnit%rwidth(iu) .and. TUnit%rwidth(iu) .lt. Tctl%rwidth_max ) then
            TUnit%nr(iu) = Tctl%nr_min + (Tctl%nr_max - Tctl%nr_min) * ( ( Tctl%rwidth_max - TUnit%rwidth(iu) ) / (Tctl%rwidth_max - Tctl%rwidth_min) )
          elseif ( TUnit%rwidth(iu) .le. Tctl%rwidth_min ) then 
            TUnit%nr(iu) = Tctl%nr_max
          elseif ( Tctl%rwidth_max .le. TUnit%rwidth(iu) ) then
            TUnit%nr(iu) = Tctl%nr_min
          end if
      
        end if
      enddo       
      
    ! 4 -- use one uniform value :
    elseif ( Tctl%OPT_calcNr .eq. 4 ) then
      
      do iu = rtmCTL%begr, rtmCTL%endr
        !if ( TUnit%mask( iu ) .gt. 0 ) then
        if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
          TUnit%nr(iu) = Tctl%nr_uniform
        end if  
      enddo 
      
    endif
    
    ! ---------------------------------  
    ! References: 
    ! Decharme, B., Alkama, R., Douville, H., Becker, M. and Cazenave, A.: 
    !     Global Evaluation of the ISBA-TRIP 15 Continental Hydrological System. Part II: Uncertainties in River Routing Simulation 
    !     Related to Flow Velocity and Groundwater Storage, J. Hydrometeorol., 11(3), 601–617, doi:10.1175/2010JHM1212.1, 2010.
    ! Luo, X., Li, H.-Y., Leung, L. R., Tesfa, T. K., Getirana, A., Papa, F., and Hess, L. L.: 
    !     Modeling surface water dynamics in the Amazon Basin using MOSART-Inundation v1.0: impacts of geomorphological parameters 
    !     and river flow representation, Geosci. Model Dev., 10(3), 1233-1259, doi:10.5194/gmd-10-1233-2017, 2017.
    ! Getirana, A. C. V., Boone, A., Yamazaki, D., Decharme, B., Papa, F. and Mognard, N.: 
    !     The Hydrological Modeling and Analysis Platform (HyMAP): Evaluation in the Amazon Basin, J. Hydrometeorol., 13(6), 1641–1665,
    !     doi:10.1175/JHM-D-12-021.1, 2012.
    ! ---------------------------------  
  
  end subroutine calc_chnlMannCoe

!-----------------------------------------------------------------------------

  subroutine preprocess_elevProf ( )
    ! DESCRIPTION: Pre-process elevation profile parameters.
    
    ! HISTORY:
    ! 2014-2016: Created and improved in offline MOSART-Inundation ( X.Luo ).
    ! 2017: Integrated with ACME.
    ! ... ... ...
    
    implicit none
    integer :: j, k, iu     ! Indexes.
    real( r8 ) :: ds        ! The volume between two elevations of an adjusted elevation profile (m^3).
    character( len = * ), parameter :: subname = '(preprocess_elevProf)'
    
    
    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
	    ! Area fraction of a computation unit (grid cell or subbasin) (dimensionless) :
        TUnit%a_eprof(iu,1) = 0.0_r8
        TUnit%a_eprof(iu,2) = 0.1_r8
        TUnit%a_eprof(iu,3) = 0.2_r8
        TUnit%a_eprof(iu,4) = 0.3_r8
        TUnit%a_eprof(iu,5) = 0.4_r8
        TUnit%a_eprof(iu,6) = 0.5_r8
        TUnit%a_eprof(iu,7) = 0.6_r8
        TUnit%a_eprof(iu,8) = 0.7_r8
        TUnit%a_eprof(iu,9) = 0.8_r8
        TUnit%a_eprof(iu,10) = 0.9_r8
        TUnit%a_eprof(iu,11) = 1.0_r8
        TUnit%a_eprof(iu,12) = 1.0_r8
      end if
    end do

    ! 1 -- use real elevation profile data :
    if ( Tctl%OPT_elevProf .eq. 1 ) then      
      
      do iu = rtmCTL%begr, rtmCTL%endr
        !if ( TUnit%mask( iu ) .gt. 0 ) then
        if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
          do k = 1, 11
            TUnit%e_eprof( iu, k ) = TUnit%e_eprof_in2( k, iu )
          end do
          !TUnit%e_eprof(iu,12) = TUnit%e_eprof_std(12)      ! The last point is hypothetical.
          TUnit%e_eprof(iu,12) = Tctl%e_eprof_std(12)        ! The last point is hypothetical.
        end if
      enddo
      
    ! 2 -- use the hypothetical elevation profile :
    elseif ( Tctl%OPT_elevProf .eq. 2 ) then    
      
      do iu = rtmCTL%begr, rtmCTL%endr
        !if ( TUnit%mask( iu ) .gt. 0 ) then
        if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
          do k = 1, 12
            TUnit%e_eprof( iu, k ) = Tctl%e_eprof_std(k)
          end do        
        end if  
      enddo
      
    end if
      
    ! ---------------------------------  
    ! Alleviate the effect of DEM pits on elevation profiles: DEM pits could make the lowest section of a elevation profile much steeper than the "true" section.
    ! Therefore, when the lowest section is much steeper than the adjacent section, it is deemed that there exist DEM pits. For this situation, the lowest section 
    ! is forced to be milder.
    ! ---------------------------------  

    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
        ! It is deemed that there exist DEM pits :
        if ( (TUnit%e_eprof(iu, 2) - TUnit%e_eprof(iu, 1) ) .gt. Tctl%threshold_slpRatio * (TUnit%e_eprof(iu, 3) - TUnit%e_eprof(iu, 2) ) ) then    
          TUnit%e_eprof(iu, 1) = TUnit%e_eprof(iu, 2) - Tctl%threshold_slpRatio * (TUnit%e_eprof(iu, 3) - TUnit%e_eprof(iu, 2) )
        end if
      end if
    enddo

    ! ---------------------------------  
    ! Slightly modify a elevation profile if it does not rise monotonously :
    ! ---------------------------------  

    do iu = rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
        do k = 2, 11
          if (TUnit%e_eprof(iu, k-1) .ge. TUnit%e_eprof(iu, k)) then
            TUnit%e_eprof(iu, k) = TUnit%e_eprof(iu, k-1) + 0.01_r8    ! (Unit: m)
          end if
        enddo
      end if
    enddo
      
    ! --------------------------------- 
    ! Calculate the elevation of channel banktop :
    ! --------------------------------- 

    do iu=rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
    
        ! Fraction of channel area ( = channel width * channel length / unit area ) :

        if ( TUnit%area(iu) .gt. 0._r8 .and. TUnit%frac(iu) .gt. 0._r8 ) then
          !TUnit%a_chnl(iu) = TUnit%rwidth(iu) * TUnit%rlen(iu) / TUnit%area(iu)   
          TUnit%a_chnl(iu) = TUnit%rwidth(iu) * TUnit%rlen(iu) / ( TUnit%area(iu) * TUnit%frac(iu) )
        else
          TUnit%a_chnl(iu) = 0._r8
          write( iulog, * ) trim( subname ) // 'TUnit%area(iu) and TUnit%frac(iu) are:', TUnit%area(iu), TUnit%frac(iu)
        end if

        do j=1, Tctl%npt_elevProf -1
          !if (TUnit%a_eprof(iu,j) < TUnit%a_chnl(iu) .and. TUnit%a_chnl(iu) <= TUnit%a_eprof(iu,j+1)) then
          if (TUnit%a_eprof(iu,j) <= TUnit%a_chnl(iu) .and. TUnit%a_chnl(iu) < TUnit%a_eprof(iu,j+1)) then          
            ! Channel banktop elevation :
            TUnit%e_chnl(iu) = interpolation_linear (TUnit%a_chnl(iu), TUnit%a_eprof(iu,j),  TUnit%a_eprof(iu,j+1), TUnit%e_eprof(iu,j), TUnit%e_eprof(iu,j+1))
          
            ! The index of the point right below the banktop in the elevation profile :
            TUnit%ipt_bl_bktp(iu) = j   
            exit    ! Jump out of the loop.
          endif
        enddo         
    
      end if
    enddo
      
    ! --------------------------------- 
    ! Create the adjusted elevation profile (i.e., replace the section below banktop with a level line), and the related coefficients :
    ! --------------------------------- 

    do iu=rtmCTL%begr, rtmCTL%endr
      !if ( TUnit%mask( iu ) .gt. 0 ) then
      if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
    
        ! ---------------------------------  
        ! Create the adjusted elevation profile :
        ! ---------------------------------  

        ! The 1st point :
        TUnit%a_eprof3(iu, 1) = 0.0_r8
        TUnit%e_eprof3(iu, 1) = 0.0_r8
      
        ! The 2nd point :
        TUnit%a_eprof3(iu, 2) = TUnit%a_chnl(iu) 
        TUnit%e_eprof3(iu, 2) = 0.0_r8
            
        ! Other upper points :
        do j=TUnit%ipt_bl_bktp(iu)+1, Tctl%npt_elevProf       
          ! Index of a point in the adjusted elevation profile :
          k = j - TUnit%ipt_bl_bktp(iu) + 2                             

          if (k > size(Tunit%a_eprof3,dim=2)) then
            write( iulog, * ) trim( subname ) // ' ERROR: k > size(Tunit%a_eprof3,dim=2) !'
            call shr_sys_abort( trim( subname ) // ' ERROR: k > size(Tunit%a_eprof3,dim=2) !' )
          end if

          ! Area fraction does not change :
          TUnit%a_eprof3(iu, k) = TUnit%a_eprof(iu, j)                    

          ! Elevation is lowered by the channel banktop elevation :
          TUnit%e_eprof3(iu, k) = TUnit%e_eprof(iu, j) - TUnit%e_chnl(iu)   
        enddo

        if (k+1 > size(Tunit%a_eprof3,dim=2)) then
          write( iulog, * ) trim( subname ) // ' ERROR: k+1 > size(Tunit%a_eprof3,dim=2) !'
          call shr_sys_abort( trim( subname ) // ' ERROR: k+1 > size(Tunit%a_eprof3,dim=2) !' )
        end if

        ! The upmost point (which is a virtual point) :
        TUnit%a_eprof3(iu, k+1) = TUnit%a_eprof3(iu, k)   ! Actually this area fraction equals 1.0 .      
        TUnit%e_eprof3(iu, k+1) = 10000.0_r8              ! Arbitrary high elevation (assume that water cannot overflow the watershed) (m).

        ! Number of points in the adjusted elevation profile :
        TUnit%npt_eprof3(iu) = k + 1
      
        ! --------------------------------- 
        ! In the inundation calculation, a quadratic equation is solved to derive the water depth based on water volume. The following coefficients are for the quadratic equation :
        ! (Please find more information in "MOSARTinund_Core_MOD.F90".)
        ! ---------------------------------  

        do j = 2, TUnit%npt_eprof3(iu) - 1    ! NOTE: j start from 2 .
          if ( abs( TUnit%e_eprof3(iu, j+1) - TUnit%e_eprof3(iu, j) ) > 1.0e-10_r8 ) then      ! Precision is 15-16 decimal digits.
            TUnit%alfa3(iu, j) = (TUnit%a_eprof3(iu, j+1) - TUnit%a_eprof3(iu, j)) / (TUnit%e_eprof3(iu, j+1) - TUnit%e_eprof3(iu, j))    ! (1/m)
          else
            write( iulog, * ) trim( subname ) // ' ERROR: Divided by zero !'
            call shr_sys_abort( trim( subname ) // ' ERROR: Divided by zero !' )
          end if

          !TUnit%p3(iu, j) = TUnit%alfa3(iu, j) * TUnit%area(iu) / 2.0_r8                 ! (m)
          TUnit%p3(iu, j) = TUnit%alfa3(iu, j) * TUnit%area(iu) * TUnit%frac(iu) / 2.0_r8 ! (m)

          !TUnit%q3(iu, j) = TUnit%a_eprof3(iu, j) * TUnit%area(iu)                       ! (m^2)
          TUnit%q3(iu, j) = TUnit%a_eprof3(iu, j) * TUnit%area(iu) * TUnit%frac(iu)       ! (m^2)
        enddo
            
        ! ---------------------------------  
        ! Calculate total volume below the level through a point of the adjusted elevation profile (i.e., the volume between channel banktop and an elevation of the adjusted elevation profile) :
        ! ---------------------------------  

        ! The volume below the 1st elevation :
        TUnit%s_eprof3(iu, 1) = 0.0_r8

        ! The volume below the 2nd elevation :
        TUnit%s_eprof3(iu, 2) = 0.0_r8

        ! The volumes below other upper elevations :
        do j = 3, TUnit%npt_eprof3(iu)
          ! The volume between the two levels through two points in the adjusted elevation profile (i.e., the volume between two elevations of the adjusted elevation profile) (m^3) :
          !ds = ( TUnit%a_eprof3(iu, j-1) + TUnit%a_eprof3(iu, j) ) * ( TUnit%e_eprof3(iu, j) - TUnit%e_eprof3(iu, j-1) ) * TUnit%area(iu) / 2.0_r8 
          ds = ( TUnit%a_eprof3(iu, j-1) + TUnit%a_eprof3(iu, j) ) * ( TUnit%e_eprof3(iu, j) - TUnit%e_eprof3(iu, j-1) ) * TUnit%area(iu) * TUnit%frac(iu) / 2.0_r8
        
          ! Total volume below the level through the current point in the adjusted elevation profile (m^3) :
          TUnit%s_eprof3(iu, j) = TUnit%s_eprof3(iu, j-1) + ds
        enddo
    
      end if    ! end of if ( rtmCTL%mask(iu) .eq. 1 .or. rtmCTL%mask(iu) .eq. 3 )
    enddo   ! do iu=rtmCTL%begr, rtmCTL%endr

  end subroutine preprocess_elevProf
  
!-----------------------------------------------------------------------------

  !function interpolation_linear (x,x1,x2,y1,y2) result (y)
  real( r8 ) function interpolation_linear (x, x1, x2, y1, y2)
    ! DESCRIPTION: For linear interpolation.

    implicit none
    real( r8 ), intent(in) :: x, x1, x2, y1, y2
    !real( r8 ), intent(out) :: y
    character( len = * ), parameter :: subname = '(interpolation_linear)'

    if (abs(x1-x2) > 1.0e-10_r8) then              ! Precision is 15-16 decimal digits.

      !y = (y2-y1)*(x-x1)/(x2-x1) + y1
      interpolation_linear = (y2 - y1) * (x - x1) / (x2 - x1) + y1

    else

      write( iulog, * ) trim( subname ) // ' ERROR: Divided by zero !'
      call shr_sys_abort( trim( subname ) // ' ERROR: Divided by zero !' )

    endif

    return
  end function interpolation_linear

!-----------------------------------------------------------------------------
  
end MODULE MOSARTinund_PreProcs_MOD
