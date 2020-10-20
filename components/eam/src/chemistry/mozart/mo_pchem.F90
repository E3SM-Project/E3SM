
      module mo_pchem

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: pchem

      contains

      subroutine pchem( nw, wl, wc, tlev, &
                        airlev, nlng, pht_tag, xs )
!-----------------------------------------------------------------------------
!   purpose:                                                                 
!   load various "weighting functions" (products of cross section and        
!   quantum yield at each altitude and for wavelength).  the altitude        
!   dependence is necessary to ensure the consideration of pressure and      
!   temperature dependence of the cross sections or quantum yields.          
!   the actual reading, evaluation and interpolation is done is separate     
!   subroutines for ease of management and manipulation.  please refer to    
!   the inline documentation of the specific subroutines for detail information.
!-----------------------------------------------------------------------------
!   parameters:                                                              
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid                                                 
!   wl     - real(r8), vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid                                         
!   nz     - integer, number of altitude levels in working altitude grid  (i)
!   tlev   - real(r8), temperature (k) at each specified altitude level       (i)
!   airlev - real(r8), air density (molec/cc) at each altitude level          (i)
!   j      - integer, counter for number of weighting functions defined  (io)
!   sq     - real(r8), cross section x quantum yield (cm^2) for each          (o)
!            photolysis reaction defined, at each defined wavelength and     
!            at each defined altitude level                                  
!   jlabel - character*40, string identifier for each photolysis reaction (o)
!            defined                                                         
!-----------------------------------------------------------------------------

      use mo_params,   only : kw
      use ppgrid,      only : pverp
      use cam_logfile, only : iulog
      use mo_xsections

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)           :: nw
      integer, intent(in)           :: nlng
      real(r8), intent(in)          :: wl(kw)
      real(r8), intent(in)          :: wc(kw)
      real(r8), intent(in)          :: tlev(pverp)
      real(r8), intent(in)          :: airlev(pverp)
      real(r8), intent(inout)       :: xs(:,:,:)
      character(len=32), intent(in) :: pht_tag(:)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer     :: m
      integer     :: astat
      character(len=32) :: jtag

rate_loop : &
      do m = 1,nlng
         jtag = trim( pht_tag(m))
         select case( jtag )
            case( 'jo3p', 'jo1d', 'j2oh' )
!-----------------------------------------------------------------------------
! 	... o3 + hv ->  (both channels) (tmp dep)
!-----------------------------------------------------------------------------
              call r01( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jn2o5' )
!-----------------------------------------------------------------------------
! 	... n2o5 + hv -> (both channels) (tmp dep)
!-----------------------------------------------------------------------------
              call r04( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jn2o' )
!-----------------------------------------------------------------------------
! 	... n2o + hv -> n2 + o(1d)         (tmp dep)
!-----------------------------------------------------------------------------
              call r44( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jh2o2' )
!-----------------------------------------------------------------------------
! 	... h2o2 + hv -> 2 oh
!-----------------------------------------------------------------------------
              call r08( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jhno3' )
!----------------------------------------------------------------------------
! 	... hno3 + hv -> oh + no2
!-----------------------------------------------------------------------------
              call r06( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jch2o_a', 'jch2o_b' )
!-----------------------------------------------------------------------------
! 	... ch2o + hv -> (both channels)
!-----------------------------------------------------------------------------
              call r10( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jch3cocho','jmgly' )
!-----------------------------------------------------------------------------
! 	... ch3cocho + hv -> products
!-----------------------------------------------------------------------------
              call r14( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jch3coch3', 'jacet' )
!-----------------------------------------------------------------------------
! 	... ch3coch3 + hv -> products
!-----------------------------------------------------------------------------
              call r15( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jch3ono2' )
!-----------------------------------------------------------------------------
! 	... ch3ono2 + hv -> ch3o + no2
!-----------------------------------------------------------------------------
              call r17( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jpan' )
!-----------------------------------------------------------------------------
! 	... pan + hv -> products
!-----------------------------------------------------------------------------
              call r18( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
            case( 'jmvk' )
!-----------------------------------------------------------------------------
! 	... mvk + hv -> products
!-----------------------------------------------------------------------------
              call xs_mvk( nw, wl, wc, tlev, airlev, xs(:,:,m) )
            case( 'jch3cho_a', 'jch3cho_b', 'jch3cho_c' )
!-----------------------------------------------------------------------------
! 	... ch3cho + hv -> products
!       (a)  ch3cho + hv -> ch3 + hco
!       (b)  ch3cho + hv -> ch4 + co
!       (c)  ch3cho + hv -> ch3co + h
!-----------------------------------------------------------------------------
              call r11( nw, wl, wc, tlev, airlev, jtag, xs(:,:,m) )
         end select
      end do rate_loop

      end subroutine pchem

      end module mo_pchem
