
  subroutine diag_dynvar_ic(phis, ps, t3, u3, v3, q3)
!
!----------------------------------------------------------------------- 
! 
! Purpose: record state variables to IC file
!
!-----------------------------------------------------------------------
!
    use shr_kind_mod, only: r8 => shr_kind_r8
    use pmgrid
    use cam_history , only: outfld, write_inithist
    use constituents, only: pcnst, cnst_name

    implicit none
!
!-----------------------------------------------------------------------
!
! Arguments
!
    real(r8), intent(in) :: phis(plon,               beglat:endlat) ! Surface geopotential
    real(r8), intent(in) :: ps  (plon,               beglat:endlat) ! surface pressure
    real(r8), intent(in) :: t3  (plon, plev,         beglat:endlat) ! temperature
    real(r8), intent(in) :: u3  (plon, plev,         beglat:endlat) ! u-wind component
    real(r8), intent(in) :: v3  (plon, plev,         beglat:endlat) ! v-wind component
    real(r8), intent(in) :: q3  (plon, plev, pcnst, beglat:endlat)  ! constituents
!
!---------------------------Local workspace-----------------------------
!
    integer lat, m   ! indices
!
!-----------------------------------------------------------------------
!
    if( write_inithist() ) then

!$OMP PARALLEL DO PRIVATE (LAT, M)
       do lat=beglat,endlat

          call outfld('PS&IC      ' , ps  (1  ,lat), plon, lat)
          call outfld('T&IC       ' , t3  (1,1,lat), plon, lat)
          call outfld('U&IC       ' , u3  (1,1,lat), plon, lat)
          call outfld('V&IC       ' , v3  (1,1,lat), plon, lat)

          do m=1,pcnst
             call outfld(trim(cnst_name(m))//'&IC', q3(1,1,m,lat), plon, lat)
          end do

       end do

    end if

    return
  end subroutine diag_dynvar_ic
