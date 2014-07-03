
      module mo_setz

      use shr_kind_mod, only : r8 => shr_kind_r8
      use abortutils,   only : endrun
      use cam_logfile,  only : iulog

      private
      public :: setz

      contains

      subroutine setz( cz, tlev, c, zen, adjcoe, pht_tag )
!-----------------------------------------------------------------------------
!   adjcoe - adjust cross section coefficients                        
!-----------------------------------------------------------------------------

      use mo_params,   only : kj
      use mo_calcoe,   only : calcoe
      use ppgrid,      only : pverp
      use chem_mods,   only : phtcnt
      use mo_tuv_inti, only : nlng, ncof

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)    :: zen                      ! zenith angle (degrees)
      real(r8), intent(in)    :: cz(pverp)
      real(r8), intent(in)    :: tlev(pverp) 
      real(r8), intent(in)    :: c(:,:,:)
      real(r8), intent(out) :: adjcoe(:,:)
      character(len=32), intent(in) :: pht_tag(:)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer, parameter  :: nzen = 4
      real(r8), parameter :: zen_angles(nzen) = (/ 20.5_r8, 40.5_r8, 60.5_r8, 80._r8 /)

      integer :: astat
      integer :: ndx
      integer :: m
      real(r8)    :: tt
      real(r8)    :: adj_fac
      real(r8)    :: interp_factor
      real(r8)    :: c0, c1, c2
      real(r8)    :: xz(pverp)
      real(r8), allocatable :: wrk(:,:)
      character(len=32) :: jname

!-----------------------------------------------------------------------------
! 1 o2 + hv -> o + o                        
! 2 o3 -> o2 + o(1d)                        
! 3 o3 -> o2 + o(3p)                        
! 4 no2 -> no + o(3p)                       
! 5 no3 -> no + o2                          
! 6 no3 -> no2 + o(3p)                      
! 7 n2o5 -> no3 + no + o(3p)                
! 8 n2o5 -> no3 + no2                       
! 9 n2o + hv -> n2 + o(1d)                  
! 10 ho2 + hv -> oh + o                      
! 11 h2o2 -> 2 oh                            
! 12 hno2 -> oh + no                         
! 13 hno3 -> oh + no2                        
! 14 hno4 -> ho2 + no2                       
! 15 ch2o -> h + hco                         
! 16 ch2o -> h2 + co                         
! 17 ch3cho -> ch3 + hco                     
! 18 ch3cho -> ch4 + co                      
! 19 ch3cho -> ch3co + h                     
! 20 c2h5cho -> c2h5 + hco                   
! 21 chocho -> products                      
! 22 ch3cocho -> products                    
! 23 ch3coch3                                
! 24 ch3ooh -> ch3o + oh                     
! 25 ch3ono2 -> ch3o+no2                     
! 26 pan + hv -> products                    
!-----------------------------------------------------------------------------

      xz(1:pverp) = cz(1:pverp)*1.e-18_r8
      do m = 1,nlng
         adjcoe(1:pverp,m) = 1._r8
      end do
      if( zen < zen_angles(1) ) then
         ndx = 1
         interp_factor = 0._r8
      else if( zen >= zen_angles(nzen) ) then
         ndx = nzen
         interp_factor = 0._r8
      else
         do ndx = 1,nzen-1
            if( zen >= zen_angles(ndx) .and. zen < zen_angles(ndx+1) ) then
!!$               interp_factor = (zen - zen_angles(ndx-1))/(zen_angles(ndx) - zen_angles(ndx-1))
               interp_factor = (zen - zen_angles(ndx))/(zen_angles(ndx+1) - zen_angles(ndx))
               exit
            end if
         end do
      end if

      allocate( wrk(pverp,2), stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'setz: failed to all wrk; error = ',astat
         call endrun
      end if

      tt = tlev(1)/281._r8
rate_loop : &
      do m = 1,nlng
         jname = trim(pht_tag(m))
	 if( jname /= 'jo1d' .and. jname /= 'j2oh' .and. jname /= 'jh2o2' ) then
	    adj_fac = 1._r8
	 else if( (jname == 'jo1d') .or. (jname == 'j2oh') ) then
!----------------------------------------------------------------------
!    	... temperature modification
!           t0.9 (1.05) t0.95(1.025)  t1.0(1.0)  t1.15(1.02)  t1.1(1.04)  
!----------------------------------------------------------------------
            select case( ndx )
	       case( 1 )
                  c0 = 4.52372_r8 ; c1 = -5.94317_r8 ; c2 = 2.63156_r8
	       case( 2 )
                  c0 = 4.99378_r8 ; c1 = -7.92752_r8 ; c2 = 3.94715_r8
	       case( 3 )
                  c0 = .969867_r8 ; c1 = -.841035_r8 ; c2 = .878835_r8
	       case( 4 )
                  c0 = 1.07801_r8 ; c1 = -2.39580_r8 ; c2 = 2.32632_r8
	    end select
            adj_fac = c0 + tt*(c1 + c2*tt)
	 else if( jname == 'jh2o2' ) then
!----------------------------------------------------------------------
!      	... temperature modification
!           t0.9 (1.05) t0.95(1.025)  t1.0(1.0)  t1.15(1.02)  t1.1(1.04)  
!----------------------------------------------------------------------
            select case( ndx )
	       case( 1 )
                  c0 = 2.43360_r8 ; c1 = -3.61363_r8 ; c2 = 2.19018_r8
	       case( 2 )
                  c0 = 3.98265_r8 ; c1 = -6.90516_r8 ; c2 = 3.93602_r8
	       case( 3 )
                  c0 = 3.49843_r8 ; c1 = -5.98839_r8 ; c2 = 3.50262_r8
	       case( 4 )
                  c0 = 3.06312_r8 ; c1 = -5.26281_r8 ; c2 = 3.20980_r8
	    end select
            adj_fac = c0 + tt*(c1 + c2*tt)
	 end if
         call calcoe( c(:,m,ndx), xz, tt, adj_fac, wrk(:,1) )
         if( interp_factor /= 0._r8 ) then
            call calcoe( c(:,m,ndx+1), xz, tt, adj_fac, wrk(:,2) )
            adjcoe(:,m) = wrk(:,1) + interp_factor * (wrk(:,2) - wrk(:,1))
         else
            adjcoe(:,m) = wrk(:,1)
         end if
      end do rate_loop

      deallocate( wrk )

      end subroutine setz

      end module mo_setz
