
module mo_photoin

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use abortutils,   only : endrun

  implicit none

  save

  public :: photoin_inti
  public :: photoin
  private

  integer               :: jo2_ndx = 0

  logical, allocatable  :: z_dep(:)
  character(len=32), allocatable :: pht_tag(:)

contains

  subroutine photoin_inti( nlng, lng_indexer )
    !-------------------------------------------------------------
    ! 	... assign use masks
    !-------------------------------------------------------------

    use mo_params,   only : largest
    use mo_setcld,   only : setcld_inti
    use chem_mods,   only : phtcnt, rxt_tag_lst

    implicit none

    !-------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------

    integer, intent(in) :: nlng
    integer, intent(in) :: lng_indexer(:)

    !-------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------
    integer :: astat
    integer :: m
    integer :: ndx
    character(len=32) :: jname

    !-------------------------------------------------------------
    ! 	... allocate module arrays
    !-------------------------------------------------------------
    has_photorates : if( nlng > 0 ) then
       allocate( z_dep(nlng), pht_tag(nlng), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'photoin_inti: failed to allocate z_dep; error = ',astat
          call endrun
       end if
       ndx = 0
       do m = 1,phtcnt
          if( lng_indexer(m) > 0 ) then
             if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
                cycle
             end if
             ndx = ndx + 1
             pht_tag(ndx) = trim( rxt_tag_lst(m))
          end if
       end do
       if( ndx /= nlng ) then
          write(iulog,*) 'photoin_inti: corrupted lng_indexer'
          call endrun
       end if
       write(iulog,*) ' '
       write(iulog,*) 'photoin_inti: lng_indexer name mapping'
       write(iulog,'(5a)') pht_tag(:)
       write(iulog,*) ' '
       !-------------------------------------------------------------
       ! 	... search for jo2
       !-------------------------------------------------------------
       do m = 1,nlng
          if( pht_tag(m) == 'jo2' ) then
             jo2_ndx = m
             exit
          end if
       end do
       write(iulog,*) ' '
       write(iulog,*) 'photoin_inti: jo2 index = ',jo2_ndx
       write(iulog,*) ' '
       !-------------------------------------------------------------
       ! 	... set altitude dependence logical array
       !-------------------------------------------------------------
       z_dep(:) = .true.
       do m = 1,nlng
          jname = pht_tag(m)
          select case( jname )
          case( 'jno2', 'jno3', 'jho2', 'jhno2', 'jho2no2' )
             z_dep(m) = .false.
          case( 'jc2h5cho', 'jchocho', 'jch3ooh' )
             z_dep(m) = .false.
          end select
       end do

       !-------------------------------------------------------------
       ! 	... intialize cloud layer module
       !-------------------------------------------------------------
       call setcld_inti

    end if has_photorates

  end subroutine photoin_inti

  subroutine photoin( idate, alat, along, &
       ut, esfact, o3top, o2top, albedo, &
       z, tlev, tlay, xlwc, xfrc, &
       airlev, aocs1, aocs2, acbs1, acbs2, &
       asoa, aant, aso4, asal, ads, o3, rh, &
       prate,  zen, nw, dt_xdiag )
    !----------------------------------------------------------
    !     	... interactive photolysis interface routine
    !----------------------------------------------------------
    
    use mo_tuv_inti, only : nlng
    use mo_params,   only : kj, kw
    use mo_wavelen,  only : deltaw, sflx, wc, wl, wu
    use mo_wavelab , only : sj
    use mo_zadj,     only : adj_coeffs
    use mo_setair,   only : setair
    use mo_setozo,   only : setozo
    use mo_pchem,    only : pchem
    use mo_sphers,   only : sphers
    use mo_airmas,   only : airmas
    use mo_setz,     only : setz
    use mo_seto2,    only : set_o2_xsect
    use mo_rtlink,   only : rtlink
    use mo_setcld,   only : setcld   !, mreg
    use mo_setaer,   only : setaer
    use ppgrid,      only : pver, pverp

    implicit none

    !----------------------------------------------------------
    !     	... dummy arguments
    !----------------------------------------------------------
    integer, intent(in)     ::  idate
    integer, intent(in)     ::  nw
    real(r8), intent(in)    ::  alat, along, o3top, o2top
    real(r8), intent(in)    ::  ut, esfact
    real(r8), intent(in)    ::  zen
    real(r8), intent(in)    ::  albedo(kw)
    real(r8), intent(in)    ::  tlay(pver)
    real(r8), intent(in)    ::  xlwc(pverp)           ! cloud water (g/m3)
    real(r8), intent(in)    ::  xfrc(pverp)           ! cloud fraction
    real(r8), intent(in)    ::  tlev(pverp)
    real(r8), intent(in)    ::  airlev(pverp)
    real(r8), intent(in)    ::  z(pverp)
    real(r8), intent(in)    ::  aocs1(pverp)
    real(r8), intent(in)    ::  aocs2(pverp)
    real(r8), intent(in)    ::  acbs1(pverp)
    real(r8), intent(in)    ::  acbs2(pverp)
    real(r8), intent(in)    ::  asoa(pverp)
    real(r8), intent(in)    ::  aant(pverp)
    real(r8), intent(in)    ::  aso4(pverp)
    real(r8), intent(in)    ::  asal(pverp,4)
    real(r8), intent(in)    ::  ads(4,pverp)
    real(r8), intent(in)    ::  rh(pverp)
    real(r8), intent(inout) ::  o3(pverp)
    real(r8), intent(out)   ::  prate(pverp,nlng)
    real(r8), intent(out)   ::  dt_xdiag(:)

    !----------------------------------------------------------
    !     	... local variables
    !----------------------------------------------------------
    integer  :: i, j, k, km, wn, n, astat
    real(r8) :: factor, delzint
    real(r8) :: wcen
    real(r8), allocatable :: xs(:,:,:)
    real(r8), allocatable :: adjcoe(:,:)    ! ftuv adjustment factor

    !----------------------------------------------------------
    ! 	... altitude grid
    !----------------------------------------------------------
    real(r8)    :: colinc(pverp)
    real(r8)    :: vcol(pverp)
    real(r8)    :: scol(pverp)
    real(r8)    :: to3(pverp)

    !----------------------------------------------------------
    ! 	... solar zenith angle
    !           slant pathlengths in spherical geometry
    !----------------------------------------------------------
    integer     :: nid(0:pver)
    real(r8)    :: dsdh(0:pver,pver)

    !----------------------------------------------------------
    ! 	... extra terrestrial solar flux and earth-sun distance ^-2
    !----------------------------------------------------------
    real(r8)    :: etf(nw)
    real(r8)    :: delw(nw)
    real(r8)    :: xsec(nw)

    !--------------------------------------------------------------
    ! 	... atmospheric optical parameters:
    !--------------------------------------------------------------
    integer, parameter :: mreg = 16   
    integer :: nreg                  ! regions at each grid
    real(r8)    :: dtrl(pver,nw)
    real(r8)    :: dto3(pver,nw)
    real(r8)    :: dto2(pver,nw)
    real(r8)    :: dtcld(pver,nw)
    real(r8)    :: omcld(pver,nw)
    real(r8)    :: gcld(pver,nw)

    real(r8)    :: dtcbs1(pver,nw)
    real(r8)    :: dtcbs2(pver,nw)
    real(r8)    :: omcbs1(pver,nw)
    real(r8)    :: omcbs2(pver,nw)
    real(r8)    :: gcbs1(pver,nw)
    real(r8)    :: gcbs2(pver,nw)

    real(r8)    :: dtocs1(pver,nw)
    real(r8)    :: dtocs2(pver,nw)
    real(r8)    :: omocs1(pver,nw)
    real(r8)    :: omocs2(pver,nw)
    real(r8)    :: gocs1(pver,nw)
    real(r8)    :: gocs2(pver,nw)

    real(r8)    :: dtant(pver,nw)
    real(r8)    :: omant(pver,nw)
    real(r8)    :: gant(pver,nw)

    real(r8)    :: dtsoa(pver,nw)
    real(r8)    :: dtso4(pver,nw)
    real(r8)    :: omso4(pver,nw)
    real(r8)    :: gso4(pver,nw)

    real(r8)    :: dtsal(pver,nw,4)
    real(r8)    :: omsal(pver,nw,4)
    real(r8)    :: gsal(pver,nw,4)

    real(r8)    :: dtds1(pver,nw)
    real(r8)    :: dtds2(pver,nw)
    real(r8)    :: dtds3(pver,nw)
    real(r8)    :: dtds4(pver,nw)
    real(r8)    :: omds1(pver,nw)
    real(r8)    :: omds2(pver,nw)
    real(r8)    :: omds3(pver,nw)
    real(r8)    :: omds4(pver,nw)
    real(r8)    :: gds1(pver,nw)
    real(r8)    :: gds2(pver,nw)
    real(r8)    :: gds3(pver,nw)
    real(r8)    :: gds4(pver,nw)

    real(r8)    :: optr(pver,mreg)         ! cld opt (z dependent) at each region
    real(r8)    :: fp(mreg)                ! probability at each region

    real(r8)    :: xso2(nw,pverp)

    !--------------------------------------------------------------
    ! 	... spectral irradiance and actinic flux (scalar irradiance):
    !--------------------------------------------------------------
    real(r8)    :: radfld(pverp,nw)
    real(r8)    :: radxx(pverp,nw)

    !-------------------------------------------------------------
    !  	... j-values:
    !-------------------------------------------------------------
    integer :: jn, m
    
    !-------------------------------------------------------------
    ! 	... location and time
    !-------------------------------------------------------------
    integer  :: iyear, imonth, iday
    real(r8) :: dtime, ut0

    !-------------------------------------------------------------
    ! 	... allocate wrking xsection array
    !-------------------------------------------------------------
    allocate( xs(nw,pverp,nlng), adjcoe(pverp,nlng), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photoin: failed to allocate xs, adjcoe; error = ',astat
       call endrun
    end if

    etf(1:nw) = sflx(1:nw) * esfact   ! earth-sun distance effect
    !-------------------------------------------------------
    !  	... air profile and rayleigh optical depths (inter-face)
    !-------------------------------------------------------
    call setair( z, nw, wc, airlev, dtrl, colinc, o2top )

    !-------------------------------------------------------------
    ! 	... ozone optical depths (must give temperature) (inter-face)
    !-------------------------------------------------------------
    call setozo( z, nw, wl, tlay, dto3, to3, o3, airlev, o3top )

    !-------------------------------------------------------------
    ! 	... cloud optical depths
    !-------------------------------------------------------------
    call setcld( z, xlwc, xfrc, nreg, fp, optr )

    !-------------------------------------------------------------
    ! 	... aerosol optical depths
    !-------------------------------------------------------------
    call setaer( z, airlev, rh, aocs1, aocs2, &
         acbs1, acbs2,&
         aant, aso4, asal, ads, asoa, &
         dtcbs1, dtcbs2, omcbs1, omcbs2, gcbs1, gcbs2, &
         dtocs1, dtocs2, omocs1, omocs2, gocs1, gocs2, &
         dtant, omant, gant, &
         dtso4, omso4, gso4, &
         dtsal, omsal, gsal, &
         dtds1, dtds2, dtds3, dtds4, &
         omds1, omds2, omds3, omds4, &
         gds1, gds2, gds3, gds4, dtsoa, nw )
    dt_xdiag(1) = sum( dtcbs1(:,16) + dtcbs2(:,16) )
    dt_xdiag(2) = sum( dtocs1(:,16) + dtocs2(:,16) )
    dt_xdiag(3) = sum( dtso4(:,16) )
    dt_xdiag(4) = sum( dtant(:,16) )
    dt_xdiag(5) = sum( dtsal(:,16,1) + dtsal(:,16,2) + dtsal(:,16,3) + dtsal(:,16,4) )
    dt_xdiag(6) = sum( dtds1(:,16) + dtds2(:,16) + dtds3(:,16) + dtds4(:,16) )
    dt_xdiag(7) = sum( dtsoa(:,16) )
    dt_xdiag(8) = sum( dt_xdiag(1:6) )

    !------------------------------------------------------------
    ! 	... photo-chemical and photo-biological weigting functions. 
    !           for pchem, need to know temperature and pressure profiles.
    !           output:
    !           from pchem:  sj(kj,kz,kw) - for each reaction
    !-------------------------------------------------------------
    xs(:,:,1:nlng) = sj(:,:,1:nlng)
    call pchem( nw, wl, wc, tlev, &
         airlev, nlng, pht_tag, xs )

    !-------------------------------------------------------------
    ! 	... slant path lengths for spherical geometry
    !-------------------------------------------------------------
    call sphers( z, zen, dsdh, nid )
    call airmas( z, zen, dsdh, nid, colinc, vcol, scol )

    !---------------------------------------------------------------
    !    	... modification of coefficent of j-vales function of to3 and zenith
    !---------------------------------------------------------------
    call setz( to3, tlev, adj_coeffs, zen, adjcoe, pht_tag )

    !------------------------------------------------------------------
    ! 	... effective o2 optical depth (sr bands, must know zenith angle!)
    !           assign o2 cross section to sj(1,*,*)
    !------------------------------------------------------------------
    call set_o2_xsect( z, nw, wl, colinc, vcol, scol, dto2, xso2 )
    if( jo2_ndx > 0 ) then
       xs(:,:,jo2_ndx) = xso2(:,:)
    end if

    delw(:nw) = deltaw(:nw) * etf(:nw)

    !---------------------------------------------------
    !  	... monochromatic radiative transfer:
    !           outputs are  fdir, fdn, fup
    !---------------------------------------------------

    ! set for cloud only

    do wn = 1,nw
       radfld(:,wn) = 0._r8 
       omcld(:,wn)  = .9999_r8
       gcld (:,wn)  = .85_r8
    end do

    Cld_reg_loop : do n = 1,nreg
       factor = fp(n)
       do wn = 1,nw
          dtcld(:,wn) = optr(:,n)
       end do

#ifdef NO_AEROSOL
       dtcbs1(:,:)  = 0._r8
       dtcbs2(:,:)  = 0._r8
       dtocs1(:,:)  = 0._r8
       dtocs2(:,:)  = 0._r8
       dtant(:,:)   = 0._r8
       dtso4(:,:)   = 0._r8
       dtsal(:,:,:) = 0._r8
       dtds1(:,:)   = 0._r8
       dtds2(:,:)   = 0._r8
       dtds3(:,:)   = 0._r8
       dtds4(:,:)   = 0._r8

       omcbs1(:,:)  = 0._r8
       omcbs2(:,:)  = 0._r8
       omocs1(:,:)  = 0._r8
       omocs2(:,:)  = 0._r8
       omant(:,:)   = 0._r8
       omso4(:,:)   = 0._r8
       omsal(:,:,:) = 0._r8
       omds1(:,:)   = 0._r8
       omds2(:,:)   = 0._r8
       omds3(:,:)   = 0._r8
       omds4(:,:)   = 0._r8

       gcbs1(:,:)  = 0._r8
       gcbs2(:,:)  = 0._r8
       gocs1(:,:)  = 0._r8
       gocs2(:,:)  = 0._r8
       gant(:,:)   = 0._r8
       gso4(:,:)   = 0._r8
       gsal(:,:,:) = 0._r8
       gds1(:,:)   = 0._r8
       gds2(:,:)   = 0._r8
       gds3(:,:)   = 0._r8
       gds4(:,:)   = 0._r8
#endif
       call rtlink( z, nw, albedo, zen, dsdh, &
            nid, dtrl, dto3, dto2, & 
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
            dtds4, omds4, gds4, radxx )
       do wn = 1,nw
          radfld(:,wn) = radfld(:,wn) + radxx(:,wn)*factor
       end do
    end do Cld_reg_loop

    !----------------------------------------------------------
    !     	... interplation at the top level
    !----------------------------------------------------------
    delzint = (z(pver-1) - z(pver-2))/(z(pver) - z(pver-1))
    do wn = 1,nw
       radfld(1,wn) = radfld(2,wn) + (radfld(2,wn) - radfld(3,wn))*delzint
       radfld(1,wn) = max( radfld(1,wn),radfld(2,wn) )
    end do

    !----------------------------------------------------------
    !   	... j-val calculation
    !           spherical irradiance (actinic flux)
    !           as a function of altitude
    !           convert to quanta s-1 nm-1 cm-2
    !           (1.e-4 * (wc*1e-9) / (hc = 6.62e-34 * 2.998e8))
    !----------------------------------------------------------
    rate_loop : do m = 1,nlng
       if( .not. z_dep(m) ) then
          xsec(:nw)     = xs(:nw,1,m) * delw(:nw)
          prate(:pverp,m) = matmul( radfld, xsec )
       else
          do k = 1,pverp
             km = pverp - k + 1
             xsec(:nw) = xs(:nw,km,m) * delw(:nw)
             prate(k,m) = dot_product( radfld(k,:nw), xsec(:nw) )
          end do
       end if
       prate(1:pverp,m) = prate(1:pverp,m) * adjcoe(pverp:1:-1,m)  
    end do rate_loop

    deallocate( xs, adjcoe )

  end subroutine photoin

end module mo_photoin
