
module mo_rtlink

  use shr_kind_mod, only : r8 => shr_kind_r8

  private
  public :: rtlink

contains

  subroutine rtlink( z, nw, albedo, zen, dsdh, &
                     nid, &
                     dtrl, &
                     dto3, & 
                     dto2, & 
                     dtcld, omcld, gcld, &
                     dtcbs1, omcbs1, gcbs1, &
                     dtcbs2, omcbs2, gcbs2, &
                     dtocs1, omocs1, gocs1, &
                     dtocs2, omocs2, gocs2, &
                     dtant, omant, gant, &
                     dtso4, omso4, gso4, &
                     dtsal, omsal, gsal, &
                     dtds1, omds1, gds1, &
                     dtds2, omds2, gds2, &
                     dtds3, omds3, gds3, &
                     dtds4, omds4, gds4, &
                     radfld )

    !-----------------------------------------------------------------------
    !
    ! Rewritten by P. Hess, April 2005 to account for new mie lookup table
    !
    !
    ! prefix dt = optical depth
    ! prefix om = single scattering albedo
    ! prefix g  = asymmetery parameter
    ! prefix ds = optical depth x single scattering albedo
    ! prefix da = optical depth x (1-single scattering albedo)
    ! suffix 1 = dry
    ! suffix 2 = wet
    ! cgs = soot, 
    ! ocs = organic carbon + soa (all soluble)
    ! ant = ammonia nitrate
    ! sal = sea-salt (4 bins) 
    ! ds1 - ds4 = dust
    !-----------------------------------------------------------------------

    use mo_params,  only : smallest
    use mo_ps2str,  only : ps2str
    use ppgrid,     only : pver, pverp

    implicit none

    !-----------------------------------------------------------------------
    !	... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) :: nw
    integer, intent(in) :: nid(0:pver)
    real(r8), intent(in)  :: z(pverp)
    real(r8), intent(in)  :: albedo(nw)
    real(r8), intent(in)  :: zen
    real(r8), intent(in)  :: dtrl(pver,nw)
    real(r8), intent(in)  :: dto3(pver,nw)
    real(r8), intent(in)  :: dto2(pver,nw)
    real(r8), intent(in)  :: dtcld(pver,nw)
    real(r8), intent(in)  :: omcld(pver,nw)
    real(r8), intent(in)  :: gcld(pver,nw)

    real(r8), intent(in)  :: dtcbs1(pver,nw)
    real(r8), intent(in)  :: omcbs1(pver,nw)
    real(r8), intent(in)  :: gcbs1(pver,nw)

    real(r8), intent(in)  :: dtcbs2(pver,nw)
    real(r8), intent(in)  :: omcbs2(pver,nw)
    real(r8), intent(in)  :: gcbs2(pver,nw)

    real(r8), intent(in)  :: dtocs1(pver,nw)
    real(r8), intent(in)  :: omocs1(pver,nw)
    real(r8), intent(in)  :: gocs1(pver,nw)

    real(r8), intent(in)  :: dtocs2(pver,nw)
    real(r8), intent(in)  :: omocs2(pver,nw)
    real(r8), intent(in)  :: gocs2(pver,nw)

    real(r8), intent(in)  :: dtant(pver,nw)
    real(r8), intent(in)  :: omant(pver,nw)
    real(r8), intent(in)  :: gant(pver,nw)
    real(r8), intent(in)  :: dtso4(pver,nw)
    real(r8), intent(in)  :: omso4(pver,nw)
    real(r8), intent(in)  :: gso4(pver,nw)
    real(r8), intent(in)  :: dtsal(pver,nw,4)
    real(r8), intent(in)  :: omsal(pver,nw,4)
    real(r8), intent(in)  :: gsal(pver,nw,4)
    real(r8), intent(in)  :: dtds1(pver,nw)
    real(r8), intent(in)  :: omds1(pver,nw)
    real(r8), intent(in)  :: gds1(pver,nw)
    real(r8), intent(in)  :: dtds2(pver,nw)
    real(r8), intent(in)  :: omds2(pver,nw)
    real(r8), intent(in)  :: gds2(pver,nw)
    real(r8), intent(in)  :: dtds3(pver,nw)
    real(r8), intent(in)  :: omds3(pver,nw)
    real(r8), intent(in)  :: gds3(pver,nw)
    real(r8), intent(in)  :: dtds4(pver,nw)
    real(r8), intent(in)  :: omds4(pver,nw)
    real(r8), intent(in)  :: gds4(pver,nw)

    real(r8), intent(in)  :: dsdh(0:pver,pver)
    real(r8), intent(out) :: radfld(pverp,nw)

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer     :: k, kk, wn
    real(r8)    :: daaer, dtsct, dtabs, dsaer, dscld, dacld
    real(r8)    :: dscbs1, dacbs1                            
    real(r8)    :: dscbs2, dacbs2                            
    real(r8)    :: dsocs1, daocs1  
    real(r8)    :: dsocs2, daocs2  
    real(r8)    :: dsant, daant                            
    real(r8)    :: dsso4, daso4  
    real(r8)    :: dssal1, dssal2, dssal3, dssal4
    real(r8)    :: dasal1, dasal2, dasal3, dasal4
    real(r8)    :: dsds1, dads1
    real(r8)    :: dsds2, dads2
    real(r8)    :: dsds3, dads3
    real(r8)    :: dsds4, dads4
    real(r8)    :: wrk
    real(r8)    :: dt(pver,nw)
    real(r8)    :: om(pver,nw)
    real(r8)    :: g(pver,nw)

    !-----------------------------------------------------------------------
    !  	... set any coefficients specific to rt scheme
    !-----------------------------------------------------------------------
    wave_loop : do wn = 1,nw
       level_loop : do k = 1,pver
          kk = pverp - k
          !-----------------------------------------------------------------------
          ! scattering and absorbing optical depths
          !-----------------------------------------------------------------------
          dscld = dtcld(k,wn)*omcld(k,wn)
          dacld = dtcld(k,wn)*abs( 1._r8 - omcld(k,wn) )
          !-----------------------------------------------------------------------
          ! black carbon
          !-----------------------------------------------------------------------
          wrk     = max( min( omcbs1(k,wn),1._r8 ),smallest )
          dscbs1 = dtcbs1(k,wn)*wrk
          dacbs1 = dtcbs1(k,wn)*(1._r8 - wrk)

          wrk     = max( min( omcbs2(k,wn),1._r8 ),smallest )
          dscbs2 = dtcbs2(k,wn)*wrk
          dacbs2 = dtcbs2(k,wn)*(1._r8 - wrk)
          !-----------------------------------------------------------------------
          ! organic carbon and soa
          !-----------------------------------------------------------------------
          wrk     = max( min( omocs1(k,wn),1._r8 ),smallest )
          dsocs1 = dtocs1(k,wn)*wrk
          daocs1 = dtocs1(k,wn)*(1._r8 - wrk)

          wrk     = max( min( omocs2(k,wn),1._r8 ),smallest )
          dsocs2 = dtocs2(k,wn)*wrk
          daocs2 = dtocs2(k,wn)*(1._r8 - wrk)
          !-----------------------------------------------------------------------
          ! ammonia sulfate
          !-----------------------------------------------------------------------
          wrk     = max( min( omant(k,wn),1._r8 ),smallest )
          dsant   = dtant(k,wn)*wrk
          daant =   dtant(k,wn)*(1._r8 - wrk)
          !-----------------------------------------------------------------------
          ! ammonia sulfate
          !-----------------------------------------------------------------------
          wrk     = max( min( omso4(k,wn),1._r8 ),smallest )
          dsso4   = dtso4(k,wn)*wrk
          daso4 =   dtso4(k,wn)*(1._r8 - wrk)
          !-----------------------------------------------------------------------
          ! summation to this point
          !-----------------------------------------------------------------------
          dtsct = dtrl(k,wn)  + dscld  &
               + dscbs2 + dscbs1 + dsocs2 + dsocs1 + dsant + dsso4 
          dtabs = dto3(k,wn) + dto2(k,wn) + dacld &
               + dacbs2 + dacbs1 + daocs2 + daocs1 + daant + daso4 
          !-----------------------------------------------------------------------
          ! sea salt
          !-----------------------------------------------------------------------
          wrk     = max( min( omsal(k,wn,1),1._r8 ),smallest )
          dssal1 = dtsal(k,wn,1)*wrk
          dasal1 = dtsal(k,wn,1)*(1._r8 - wrk)

          wrk     = max( min( omsal(k,wn,2),1._r8 ),smallest )
          dssal2 = dtsal(k,wn,2)*wrk
          dasal2 = dtsal(k,wn,2)*(1._r8 - wrk)

          wrk     = max( min( omsal(k,wn,3),1._r8 ),smallest )
          dssal3 = dtsal(k,wn,3)*wrk
          dasal3 = dtsal(k,wn,3)*(1._r8 - wrk)

          wrk     = max( min( omsal(k,wn,4),1._r8 ),smallest )
          dssal4 = dtsal(k,wn,4)*wrk
          dasal4 = dtsal(k,wn,4)*(1._r8 - wrk)
          !-----------------------------------------------------------------------
          ! summation
          !-----------------------------------------------------------------------
          dtsct = dtsct + dssal1 + dssal2 + dssal3 + dssal4
          dtabs = dtabs + dasal1 + dasal2 + dasal3 + dasal4

          wrk     = max( min( omds1(k,wn),1._r8 ),smallest )
          dsds1 = dtds1(k,wn)*wrk
          dads1 = dsds1*(1._r8 - wrk)

          wrk     = max( min( omds2(k,wn),1._r8 ),smallest )
          dsds2 = dtds2(k,wn)*wrk
          dads2 = dsds2*(1._r8 - wrk)

          wrk     = max( min( omds3(k,wn),1._r8 ),smallest )
          dsds3 = dtds3(k,wn)*wrk
          dads3 = dsds3*(1._r8 - wrk)

          wrk     = max( min( omds4(k,wn),1._r8 ),smallest )
          dsds4 = dtds4(k,wn)*wrk
          dads4 = dsds4*(1._r8 - wrk)

          dtsct = dtsct + dsds1 + dsds2 + dsds3 + dsds4
          dtabs = dtabs + dads1 + dads2 + dads3 + dads4

          dtabs = max( dtabs,smallest )
          dtsct = max( dtsct,smallest )
          !-----------------------------------------------------------------------
          ! 	... invert z-coordinate
          !-----------------------------------------------------------------------
          dt(kk,wn) = dtsct + dtabs
          if( dtsct /= smallest ) then
             om(kk,wn) = dtsct/(dtsct + dtabs)
             g(kk,wn) = gcld(k,wn)*dscld     &
                      + gcbs1(k,wn)*dscbs1 &
                      + gcbs2(k,wn)*dscbs2 &
                      + gocs2(k,wn)*dsocs2 &
                      + gocs1(k,wn)*dsocs1 &
                      + gant(k,wn)*dsant     &
                      + gso4(k,wn)*dsso4     &            
                      + gsal(k,wn,1)*dssal1  &
                      + gsal(k,wn,2)*dssal2  &
                      + gsal(k,wn,3)*dssal3  &
                      + gsal(k,wn,4)*dssal4

             g(kk,wn) = g(kk,wn) &
                  + gds1(k,wn)*dsds1 + gds2(k,wn)*dsds2 &
                  + gds3(k,wn)*dsds3 + gds4(k,wn)*dsds4
             g(kk,wn) = min( max( g(kk,wn)/dtsct,smallest ),1._r8 )
          else
             om(kk,wn) = smallest
             g(kk,wn)  = smallest
          end if
       end do level_loop
    end do wave_loop

    !-----------------------------------------------------------------------
    !  	... call rt routine
    !-----------------------------------------------------------------------
    call ps2str( nw, zen, albedo, dt, om, &
                 g, dsdh, nid, radfld )

  end subroutine rtlink

end module mo_rtlink
