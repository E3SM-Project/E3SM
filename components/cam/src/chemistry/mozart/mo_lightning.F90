module mo_lightning
  !----------------------------------------------------------------------
  ! ... the lightning module
  !----------------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use ppgrid,        only : begchunk, endchunk, pcols, pver
  use phys_grid,     only : ngcols_p
  use cam_abortutils,    only : endrun
  use cam_logfile,   only : iulog
  use spmd_utils,    only : masterproc, mpicom

  implicit none

  private
  public  :: lightning_inti
  public  :: lightning_no_prod
  public  :: prod_no

  save

  real(r8) :: csrf
  real(r8) :: factor = 0.1_r8              ! user-controlled scaling factor to achieve arbitrary no prod.
  real(r8) :: geo_factor                   ! grid cell area factor
  real(r8) :: vdist(16,3)                  ! vertical distribution of lightning
  real(r8), allocatable :: prod_no(:,:,:)
  real(r8), allocatable :: glob_prod_no_col(:,:)
  real(r8), allocatable :: flash_freq(:,:)
  integer :: no_ndx,xno_ndx
  logical :: has_no_lightning_prod = .false.

contains

  subroutine lightning_inti( lght_no_prd_factor )
    !----------------------------------------------------------------------
    !       ... initialize the lightning module
    !----------------------------------------------------------------------
    use mo_constants,  only : pi
    use ioFileMod,     only : getfil
    use mo_chem_utls,  only : get_spc_ndx

    use cam_history,   only : addfld, horiz_only
    use dyn_grid,      only : get_dyn_grid_parm

    implicit none

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    real(r8), intent(in) :: lght_no_prd_factor        ! lightning no production factor

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    integer  :: astat
    integer  :: ncid
    integer  :: dimid
    integer  :: vid
    integer  :: gndx
    integer  :: jl, ju
    integer  :: nlat, nlon
    integer  :: plon, plat
    real(r8), allocatable :: lats(:)
    real(r8), allocatable :: lons(:)
    real(r8), allocatable :: landmask(:,:)
    character(len=256) :: locfn

    no_ndx = get_spc_ndx('NO')
    xno_ndx = get_spc_ndx('XNO')

    has_no_lightning_prod = no_ndx>0 .or. xno_ndx>0
    if (.not.has_no_lightning_prod) return

    
    if( lght_no_prd_factor /= 1._r8 ) then
       factor = factor*lght_no_prd_factor
    end if


    if (masterproc) write(iulog,*) 'lght_inti: lightning no production scaling factor = ',factor

    !----------------------------------------------------------------------
    !       ... vdist(kk,itype) = % of lightning nox between (kk-1) and (kk)
    !           km for profile itype
    !----------------------------------------------------------------------
    vdist(:,1) = (/  3.0_r8, 3.0_r8, 3.0_r8, 3.0_r8, 3.4_r8, 3.5_r8, 3.6_r8, 4.0_r8, &       ! midlat cont
                     5.0_r8, 7.0_r8, 9.0_r8, 14.0_r8, 16.0_r8, 14.0_r8, 8.0_r8, 0.5_r8 /)
    vdist(:,2) = (/  2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 6.1_r8, &       ! trop marine
                     17.0_r8, 15.4_r8, 14.5_r8, 13.0_r8, 12.5_r8, 2.8_r8, 0.9_r8, 0.3_r8 /)
    vdist(:,3) = (/  2.0_r8, 2.0_r8, 2.0_r8, 1.5_r8, 1.5_r8, 1.5_r8, 3.0_r8, 5.8_r8, &       ! trop cont
                     7.6_r8, 9.6_r8, 11.0_r8, 14.0_r8, 14.0_r8, 14.0_r8, 8.2_r8, 2.3_r8 /)

    allocate( prod_no(pcols,pver,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'lght_inti: failed to allocate prod_no; error = ',astat
       call endrun
    end if
    allocate( flash_freq(pcols,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'lght_inti: failed to allocate flash_freq; error = ',astat
       call endrun
    end if
    allocate( glob_prod_no_col(pcols,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'lght_inti: failed to allocate glob_prod_no_col; error = ',astat
       call endrun
    end if
    prod_no(:,:,:)   = 0._r8
    flash_freq(:,:)  = 0._r8
    geo_factor = ngcols_p/(4._r8*pi)


    call addfld( 'LNO_COL_PROD', horiz_only,    'I','TG N/YR', 'lighting column NO source' )
    call addfld( 'LNO_PROD',  (/ 'lev' /), 'I',    '/cm3/s', 'lighting insitu NO source' )
    call addfld( 'FLASHFRQ',   horiz_only,    'I',    '1/MIN', 'lighting flash rate' )        ! flash frequency in grid box per minute (PPP)
    call addfld( 'FLASHENGY',     horiz_only,    'I',   '   ', 'lighting flash rate' )        ! flash frequency in grid box per minute (PPP)
    call addfld( 'CLDHGT',      horiz_only,    'I',      'KM', 'cloud top height' )           ! cloud top height
    call addfld( 'DCHGZONE',      horiz_only,    'I',    'KM', 'depth of discharge zone' )           ! depth of discharge zone
    call addfld( 'CGIC',   horiz_only,    'I',        'RATIO', 'ratio of cloud-ground/intracloud discharges' )        ! ratio of cloud-ground/intracloud discharges

  end subroutine lightning_inti

  subroutine lightning_no_prod( state, pbuf2d,  cam_in )
    !----------------------------------------------------------------------
    !	... set no production from lightning
    !----------------------------------------------------------------------
    use physics_types,    only : physics_state
    
    use physics_buffer,   only : pbuf_get_index, physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use physconst,        only : rga
    use phys_grid,        only : get_rlat_all_p, get_lat_all_p, get_lon_all_p, get_wght_all_p
    use cam_history,      only : outfld
    use camsrfexch,       only : cam_in_t
    use shr_reprosum_mod, only : shr_reprosum_calc
    use mo_constants,  only : rearth, d2r
    implicit none

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    type(physics_state), intent(in) :: state(begchunk:endchunk) ! physics state
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(cam_in_t), intent(in) :: cam_in(begchunk:endchunk) ! physics state

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    real(r8), parameter    :: land   = 1._r8
    real(r8), parameter    :: secpyr = 365._r8 * 8.64e4_r8

    integer :: i, c
    integer :: cldtind             ! level index for cloud top
    integer :: cldbind             ! level index for cloud base > 273k
    integer :: surf_type
    integer :: file                ! file index
    integer :: k, kk, zlow_ind, zhigh_ind, itype
    real(r8) :: glob_flashfreq     ! global flash frequency [s-1]
    real(r8) :: glob_noprod        ! global rate of no production [as tgn/yr]
    real(r8) :: frac_sum           ! work variable
    real(r8) :: zlow
    real(r8) :: zhigh
    real(r8) :: zlow_scal
    real(r8) :: zhigh_scal
    real(r8) :: fraction
    real(r8) :: dchgz
    real(r8) :: dchgzone(pcols,begchunk:endchunk)           ! depth of discharge zone [km]
    real(r8) :: cldhgt(pcols,begchunk:endchunk)             ! cloud top height [km]
    real(r8) :: cgic(pcols,begchunk:endchunk)               ! cloud-ground/intracloud discharge ratio
    real(r8) :: flash_energy(pcols,begchunk:endchunk)       ! energy of flashes per second
    real(r8) :: prod_no_col(pcols,begchunk:endchunk)        ! global no production rate for diagnostics
    real(r8) :: wrk, wrk1, wrk2(1)
    integer  :: ncol                                    ! columns per chunk
    integer  :: lchnk                                   ! columns per chunk
    real(r8),pointer :: cldtop(:)                       ! cloud top level index
    real(r8),pointer :: cldbot(:)                       ! cloud bottom level index
    real(r8) :: zmid(pcols,pver)                        ! geopot height above surface at midpoints (m)
    real(r8) :: zint(pcols,pver+1,begchunk:endchunk)    ! geopot height above surface at interfaces (m)
    real(r8) :: zsurf(pcols)                            ! geopot height above surface at interfaces (m)
    real(r8) :: rlats(pcols,begchunk:endchunk)          ! column latitudes in chunks
    real(r8) :: wght(pcols)

    !----------------------------------------------------------------------
    ! 	... parameters to determine cg/ic ratio [price and rind, 1993]
    !----------------------------------------------------------------------
    real(r8), parameter  :: ca = .021_r8
    real(r8), parameter  :: cb = -.648_r8
    real(r8), parameter  :: cc = 7.49_r8
    real(r8), parameter  :: cd = -36.54_r8
    real(r8), parameter  :: ce = 64.09_r8
    real(r8), parameter  :: t0 = 273._r8
    real(r8), parameter  :: m2km  = 1.e-3_r8
    real(r8), parameter  :: km2cm = 1.e5_r8
    real(r8), parameter  :: lat25 = 25._r8*d2r      ! 25 degrees latitude in radians
    integer  :: cldtop_ndx, cldbot_ndx
    integer  :: istat
    real(r8) :: flash_freq_land, flash_freq_ocn

    if (.not.has_no_lightning_prod) return

    !----------------------------------------------------------------------
    !	... initialization
    !----------------------------------------------------------------------

    flash_freq(:,:)       = 0._r8
    prod_no(:,:,:)        = 0._r8
    prod_no_col(:,:)      = 0._r8
    cldhgt(:,:)           = 0._r8
    dchgzone(:,:)         = 0._r8
    cgic(:,:)             = 0._r8
    flash_energy(:,:)     = 0._r8
    glob_prod_no_col(:,:) = 0._r8

    cldtop_ndx = pbuf_get_index('CLDTOP')
    cldbot_ndx = pbuf_get_index('CLDBOT')

    !--------------------------------------------------------------------------------
    !	... estimate flash frequency and resulting no emissions
    !           [price, penner, prather, 1997 (jgr)]
    !    lightning only occurs in convective clouds with a discharge zone, i.e.
    !    an altitude range where liquid water, ice crystals, and graupel coexist.
    !    we test this by examining the temperature at the cloud base.
    !    it is assumed that only one thunderstorm occurs per grid box, and its
    !    flash frequency is determined by the maximum cloud top height (not the
    !    depth of the discharge zone). this is somewhat speculative but yields
    !    reasonable results.
    !
    !       the cg/ic ratio is determined by an empirical formula from price and
    !    rind [1993]. the average energy of a cg flash is estimated as 6.7e9 j,
    !    and the average energy of a ic flash is assumed to be 1/10 of that value.
    !       the no production rate is assumed proportional to the discharge energy
    !    with 1e17 n atoms per j. the total number of n atoms is then distributed
    !    over the complete column of grid boxes.
    !--------------------------------------------------------------------------------
    Chunk_loop : do c = begchunk,endchunk
       ncol  = state(c)%ncol
       lchnk = state(c)%lchnk
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), cldtop_ndx, cldtop )
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), cldbot_ndx, cldbot )
       zsurf(:ncol) = state(c)%phis(:ncol)*rga
       call get_rlat_all_p( c, ncol, rlats(1,c) )
       call get_wght_all_p(c, ncol, wght)

       do k = 1,pver
          zmid(:ncol,k)   = state(c)%zm(:ncol,k) + zsurf(:ncol)
          zint(:ncol,k,c) = state(c)%zi(:ncol,k) + zsurf(:ncol)
       end do
       zint(:ncol,pver+1,c) = state(c)%zi(:ncol,pver+1) + zsurf(:ncol)

       col_loop : do i = 1,ncol
          !--------------------------------------------------------------------------------
          ! 	... find cloud top and bottom level above 273k
          !--------------------------------------------------------------------------------
          cldtind = nint( cldtop(i) )
          cldbind = nint( cldbot(i) )
          do
             if( cldbind <= cldtind .or. state(c)%t(i,cldbind) < t0 ) then
                exit
             end if
             cldbind = cldbind - 1
          end do
          cloud_layer : if( cldtind < pver .and. cldtind > 0 .and. cldtind < cldbind ) then
             !--------------------------------------------------------------------------------
             !       ... compute cloud top height and depth of charging zone
             !--------------------------------------------------------------------------------
             cldhgt(i,c)   = m2km*max( 0._r8,zint(i,cldtind,c) )
             dchgz = cldhgt(i,c) - m2km*zmid(i,cldbind)
             dchgzone(i,c) = dchgz
             !--------------------------------------------------------------------------------
             !       ... compute flash frequency for given cloud top height
             !           (flashes storm^-1 min^-1)
             !--------------------------------------------------------------------------------
             flash_freq_land = 3.44e-5_r8 * cldhgt(i,c)**4.9_r8 
             flash_freq_ocn  = 6.40e-4_r8 * cldhgt(i,c)**1.7_r8
             flash_freq(i,c) = cam_in(c)%landfrac(i)*flash_freq_land + &
                               cam_in(c)%ocnfrac(i) *flash_freq_ocn

             !--------------------------------------------------------------------------------
             !       ... compute cg/ic ratio
             !           cgic = proportion of cg flashes (=pg from ppp paper)
             !--------------------------------------------------------------------------------
             cgic(i,c) = 1._r8/((((ca*dchgz + cb)*dchgz + cc) *dchgz + cd)*dchgz + ce)
             if( dchgz < 5.5_r8 ) then
                cgic(i,c) = 0._r8
             else if( dchgz > 14._r8 ) then
                cgic(i,c) = .02_r8
             end if
             !--------------------------------------------------------------------------------
             !       ... compute flash energy (cg*6.7e9 + ic*6.7e8)
             !           and convert to total energy per second
             !           set ic = cg
             !--------------------------------------------------------------------------------
             flash_energy(i,c) = 6.7e9_r8 * flash_freq(i,c)/60._r8
             !--------------------------------------------------------------------------------
             !       ... LKE Aug 23, 2005: scale production to account for different grid
             !           box sizes. This requires a reduction in the overall fudge factor 
             !           (e.g., from 1.2 to 0.5) 
             !--------------------------------------------------------------------------------
             flash_energy(i,c) =  flash_energy(i,c) * wght(i) * geo_factor
             !--------------------------------------------------------------------------------
             ! 	... compute number of n atoms produced per second
             !           and convert to n atoms per second per cm2 and apply fudge factor
             !--------------------------------------------------------------------------------
             prod_no_col(i,c) = 1.e17_r8*flash_energy(i,c)/(1.e4_r8*rearth*rearth*wght(i)) * factor

             !--------------------------------------------------------------------------------
             ! 	... compute global no production rate in tgn/yr:
             !           tgn per second: * 14.00674 * 1.65979e-24 * 1.e-12
             !             nb: 1.65979e-24 = 1/avo
             !           tgn per year: * secpyr
             !--------------------------------------------------------------------------------
             glob_prod_no_col(i,c) = 1.e17_r8*flash_energy(i,c) &
                  * 14.00674_r8 * 1.65979e-24_r8 * 1.e-12_r8 * secpyr * factor

          end if cloud_layer
       end do Col_loop
    end do Chunk_loop
    !--------------------------------------------------------------------------------
    ! 	... Accumulate global total, convert to flashes per second
    ! 	... Accumulate global NO production rate
    !--------------------------------------------------------------------------------
    kk = pcols*(endchunk-begchunk+1)
    call shr_reprosum_calc( flash_freq, wrk2,kk,kk,1, commid=mpicom)
    glob_flashfreq=wrk2(1)/60._r8
    call shr_reprosum_calc( glob_prod_no_col, wrk2,kk,kk,1, commid=mpicom)
    glob_noprod = wrk2(1)
    if( masterproc ) then
       write(iulog,*) ' '
       write(iulog,'(''Global flash freq (/s), lightning NOx (TgN/y) = '',2f10.4)') &
            glob_flashfreq, glob_noprod
    end if

    if( glob_noprod > 0._r8 ) then
       !--------------------------------------------------------------------------------
       !	... Distribute production up to cloud top [Pickering et al., 1998 (JGR)]
       !--------------------------------------------------------------------------------
       do c = begchunk,endchunk
          ncol  = state(c)%ncol
          lchnk = state(c)%lchnk
          call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), cldtop_ndx, cldtop )
          do i = 1,ncol
             cldtind = nint( cldtop(i) )
             if( prod_no_col(i,c) > 0._r8 ) then
                if( cldhgt(i,c) > 0._r8 ) then
                   if( abs( rlats(i,c) ) > lat25 ) then
                      itype = 1                                                    ! midlatitude continental
                   else if( nint( cam_in(c)%landfrac(i) ) == land ) then
                      itype = 3                                                    ! tropical continental
                   else
                      itype = 2                                                    ! topical marine
                   end if
                   frac_sum = 0._r8
                   do k = cldtind,pver
                      zlow       = zint(i,k+1,c) * m2km                            ! lower interface height (km)
                      zlow_scal  = zlow * 16._r8/cldhgt(i,c)                       ! scale to 16 km convection height
                      zlow_ind   = max( 1,INT(zlow_scal)+1 )                       ! lowest vdist index to include in layer
                      zhigh      = zint(i,k,c) * m2km                              ! upper interface height (km)
                      zhigh_scal = zhigh * 16._r8/cldhgt(i,c)                      ! height (km) scaled to 16km convection height
                      zhigh_ind  = max( 1,MIN( 16,INT(zhigh_scal)+1 ) )            ! highest vdist index to include in layer
                      do kk = zlow_ind,zhigh_ind
                         wrk  = kk
                         wrk1 = kk - 1
                         fraction = min( zhigh_scal,wrk ) &                         ! fraction of vdist in this model layer
                              - max( zlow_scal,wrk1 )
                         fraction = max( 0._r8, min( 1._r8,fraction ) )
                         frac_sum = frac_sum + fraction*vdist(kk,itype)
                         prod_no(i,k,c) = prod_no(i,k,c) &                         ! sum the fraction of column NOx in layer k
                              + fraction*vdist(kk,itype)*.01_r8
                      end do
                      prod_no(i,k,c) = prod_no_col(i,c) * prod_no(i,k,c) &         ! multiply fraction by column amount
                           / (km2cm*(zhigh - zlow))                    ! and convert to atom N cm^-3 s^-1
                   end do
                end if
             end if
          end do
       end do
    end if

    !--------------------------------------------------------------------------------
    !       ... output lightning no production to history file
    !--------------------------------------------------------------------------------
    do c = begchunk,endchunk
       lchnk = state(c)%lchnk
       call outfld( 'LNO_PROD',     prod_no(:,:,c),        pcols, lchnk )
       call outfld( 'LNO_COL_PROD', glob_prod_no_col(:,c), pcols, lchnk )
       call outfld( 'FLASHFRQ',     flash_freq(:,c),       pcols, lchnk )
       call outfld( 'FLASHENGY',    flash_energy(:,c),     pcols, lchnk )
       call outfld( 'CLDHGT',       cldhgt(:,c),           pcols, lchnk )
       call outfld( 'DCHGZONE',     dchgzone(:,c),         pcols, lchnk )
       call outfld( 'CGIC',         cgic(:,c),             pcols, lchnk )
    enddo

  end subroutine lightning_no_prod

end module mo_lightning
