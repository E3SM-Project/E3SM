
module mo_setext

  use cam_logfile, only: iulog

  private
  public :: setext_inti, setext, has_ions

  save

  integer :: co_ndx, no_ndx, synoz_ndx, xno_ndx
  integer :: op_ndx, o2p_ndx, np_ndx, n2p_ndx, n2d_ndx, n_ndx, e_ndx, oh_ndx
  logical :: has_ions = .false.

contains

  subroutine setext_inti
    !--------------------------------------------------------
    !	... Initialize the external forcing module
    !--------------------------------------------------------

    use mo_chem_utls, only : get_extfrc_ndx
    use ppgrid,       only : pver
    use cam_history,  only : addfld, phys_decomp
    use spmd_utils,   only : masterproc

    implicit none

    co_ndx    = get_extfrc_ndx( 'CO' )
    no_ndx    = get_extfrc_ndx( 'NO' )
    synoz_ndx = get_extfrc_ndx( 'SYNOZ' )
    xno_ndx   = get_extfrc_ndx( 'XNO' )

    op_ndx   = get_extfrc_ndx( 'Op' )
    o2p_ndx  = get_extfrc_ndx( 'O2p' )
    np_ndx   = get_extfrc_ndx( 'Np' )
    n2p_ndx  = get_extfrc_ndx( 'N2p' )
    n2d_ndx  = get_extfrc_ndx( 'N2D' )
    n_ndx    = get_extfrc_ndx( 'N' )
    e_ndx    = get_extfrc_ndx( 'e' )
    oh_ndx   = get_extfrc_ndx( 'OH' )

    has_ions = op_ndx > 0 .and. o2p_ndx > 0 .and. np_ndx > 0 .and. n2p_ndx > 0 .and. e_ndx > 0

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'setext_inti: diagnostics: co_ndx, no_ndx, synoz_ndx, xno_ndx'
       write(iulog,'(10i5)') co_ndx, no_ndx, synoz_ndx, xno_ndx
    endif

    call addfld( 'NO_Lightning','molec/cm3/s', pver, 'A', 'lightning NO source', phys_decomp )
    call addfld( 'NO_Aircraft', 'molec/cm3/s', pver, 'A', 'aircraft NO source',  phys_decomp )
    call addfld( 'CO_Aircraft', 'molec/cm3/s', pver, 'A', 'aircraft CO source',  phys_decomp )

    call addfld( 'N4S_SPE', 'molec/cm3/s', pver, 'I', 'solar proton event N(4S) source',phys_decomp )
    call addfld( 'N2D_SPE', 'molec/cm3/s', pver, 'I', 'solar proton event N(2S) source',phys_decomp )
    call addfld( 'OH_SPE',  'molec/cm3/s', pver, 'I', 'solar proton event HOx source',  phys_decomp )

    if ( has_ions ) then
       call addfld( 'P_Op', '/s ', pver, 'I', 'production o+', phys_decomp )
       call addfld( 'P_O2p', '/s ', pver, 'I', 'production o2+', phys_decomp )
       call addfld( 'P_N2p', '/s ', pver, 'I', 'production n2+', phys_decomp )
       call addfld( 'P_Np', '/s ', pver, 'I', 'production n+', phys_decomp )
       call addfld( 'P_IONS', '/s ', pver, 'I', 'total ion production', phys_decomp )
    endif

  end subroutine setext_inti

  subroutine setext( extfrc, zint_abs, zint_rel, cldtop, &
       zmid, lchnk, tfld, o2mmr, ommr, &
       pmid, mbar, rlats, calday, ncol, rlons, pbuf )
    !--------------------------------------------------------
    !     ... for this latitude slice:
    !         - form the production from datasets
    !         - form the nox (xnox) production from lighing
    !         - form the nox (xnox) production from airplanes
    !         - form the co production from airplanes
    !--------------------------------------------------------

    use cam_history,  only : outfld
    use shr_kind_mod, only : r8 => shr_kind_r8
    use ppgrid,       only : pver, pcols
    use mo_airplane,  only : airpl_set
    use chem_mods,    only : extcnt
    use mo_lightning, only : prod_no

    use mo_extfrc,    only : extfrc_set
    use chem_mods,    only : extcnt
    use tracer_srcs,  only : num_tracer_srcs, tracer_src_flds, get_srcs_data
    use mo_chem_utls, only : get_extfrc_ndx
    use mo_synoz,     only : po3

    use mo_aurora,      only : aurora
    use mo_solarproton, only : spe_prod 
    use physics_buffer, only : physics_buffer_desc

    implicit none

    !--------------------------------------------------------
    !     ... dummy arguments
    !--------------------------------------------------------
    !--------------------------------------------------------
    integer,  intent(in)  ::   lchnk                       ! chunk id
    integer,  intent(in)  ::   ncol                        ! columns in chunk
    real(r8), intent(in)  ::   zint_abs(ncol,pver+1)           ! interface geopot height ( km )
    real(r8), intent(in)  ::   zint_rel(ncol,pver+1)           ! interface geopot height ( km )
    real(r8), intent(in)  ::   cldtop(ncol)                ! cloud top index
    real(r8), intent(out) ::   extfrc(ncol,pver,extcnt)    ! the "extraneous" forcing

    real(r8), intent(in)  ::   calday                      ! calendar day of year
    real(r8), intent(in)  ::   rlats(ncol)                 ! column latitudes (radians)
    real(r8), intent(in)  ::   rlons(ncol)                 ! column longitudes (radians)
    real(r8), intent(in)  ::   zmid(ncol,pver)             ! midpoint geopot height ( km )
    real(r8), intent(in)  ::   pmid(pcols,pver)            ! midpoint pressure (Pa)
    real(r8), intent(in)  ::   tfld(pcols,pver)            ! midpoint temperature (K)
    real(r8), intent(in)  ::   mbar(ncol,pver)             ! mean molecular mass (g/mole)
    real(r8), intent(in)  ::   o2mmr(ncol,pver)            ! o2 concentration (kg/kg)
    real(r8), intent(in)  ::   ommr(ncol,pver)             ! o concentration (kg/kg)

    type(physics_buffer_desc),pointer :: pbuf(:)

    !--------------------------------------------------------
    !     ... local variables
    !--------------------------------------------------------
    integer :: i, k, k1, kk, cldind, ii, jj, nlev
    real(r8) :: srcs_offline( ncol, pver )
    integer :: ndx

    real(r8), dimension(ncol,pver) :: no_lgt

    real(r8)    :: spe_nox(ncol,pver)        ! Solar Proton Event NO production
    real(r8)    :: spe_hox(ncol,pver)        ! Solar Proton Event HOx production

    extfrc(:,:,:) = 0._r8

    no_lgt(:,:) = 0._r8

    !--------------------------------------------------------
    !     ... set frcing from datasets
    !--------------------------------------------------------
    call extfrc_set( lchnk, zint_rel, extfrc, ncol )
    
    !--------------------------------------------------------
    !     ... set nox production from lighting
    !         note: from ground to cloud top production is c shaped
    !--------------------------------------------------------
    if ( no_ndx > 0 ) then
       do i = 1,ncol
          cldind = nint( cldtop(i) )
          if( cldind < pver .and. cldind > 0  ) then
             extfrc(i,cldind:pver,no_ndx) = extfrc(i,cldind:pver,no_ndx) &
                  + prod_no(i,cldind:pver,lchnk)
             no_lgt(i,cldind:pver) = prod_no(i,cldind:pver,lchnk)
          end if
       end do
    endif
    if ( xno_ndx > 0 ) then
       do i = 1,ncol
          cldind = nint( cldtop(i) )
          if( cldind < pver .and. cldind > 0  ) then
             extfrc(i,cldind:pver,xno_ndx) = extfrc(i,cldind:pver,xno_ndx) &
                  + prod_no(i,cldind:pver,lchnk)
          end if
       end do
    endif

    call outfld( 'NO_Lightning', no_lgt(:ncol,:), ncol, lchnk )

    call airpl_set( lchnk, ncol, no_ndx, co_ndx, xno_ndx, cldtop, zint_abs, extfrc)

    !---------------------------------------------------------------------
    ! 	... synoz production
    !---------------------------------------------------------------------
    if( synoz_ndx > 0 ) then
       do k = 1,pver
          extfrc(:ncol,k,synoz_ndx) = extfrc(:ncol,k,synoz_ndx) + po3(:ncol,k,lchnk)
       end do
    end if

    do i = 1,num_tracer_srcs

       ndx =  get_extfrc_ndx( tracer_src_flds(i) )
       call get_srcs_data( tracer_src_flds(i), srcs_offline,  ncol, lchnk, pbuf )
       do k = 1,pver
          extfrc(:ncol,k,ndx) = extfrc(:ncol,k,ndx) + srcs_offline(:ncol,k)
       enddo

    enddo


    if ( has_ions ) then
       !---------------------------------------------------------------------
       !     ... set ion auroral production
       !---------------------------------------------------------------------

       call aurora( tfld, o2mmr, ommr, mbar, rlats, &
            extfrc(:,:,o2p_ndx), extfrc(:,:,op_ndx), extfrc(:,:,n2p_ndx), extfrc(:,:,np_ndx), pmid, &
            lchnk, calday, ncol, rlons, pbuf )
       !---------------------------------------------------------------------
       !     ... set n(2d) and n(4s) auroral production
       !         Stan Solomon HAO
       !         include production of N by secondary auroral hot electrons (e_s*):
       !         this is not a "real" reaction; instead, the production is parameterized in terms
       !         of the production rate of N2+ by primary electrons, QN2P (which is in the model),
       !         as follows:
       !---------------------------------------------------------------------
       do k = 1,pver
          extfrc(:,k,n2d_ndx) = 1.57_r8*.6_r8*extfrc(:,k,n2p_ndx)
          extfrc(:,k,n_ndx) = 1.57_r8*.4_r8*extfrc(:,k,n2p_ndx)
       end do
       !---------------------------------------------------------------------
       !     ... set electron auroral production
       !---------------------------------------------------------------------
       do k = 1,pver
          extfrc(:,k,e_ndx) = extfrc(:,k,op_ndx) + extfrc(:,k,o2p_ndx) &
               + extfrc(:,k,np_ndx) + extfrc(:,k,n2p_ndx)
       end do
       !---------------------------------------------------------------------
       !     ... set SPE NOx and HOx production
       !     Jackman et al., JGR, 2005 
       !     production of 1.25 Nitrogen atoms/ion pair with branching ratios 
       !     of 0.55 N(4S) and 0.7 N(2D).
       !---------------------------------------------------------------------
       call spe_prod( spe_nox, spe_hox, pmid, zmid, lchnk, ncol)

       call outfld( 'N2D_SPE', 0.7_r8*spe_nox,  ncol, lchnk ) ! N(2D) produciton (molec/cm3/s)
       call outfld( 'N4S_SPE',0.55_r8*spe_nox,  ncol, lchnk ) ! N(4S) produciton (molec/cm3/s)
       call outfld( 'OH_SPE' ,        spe_hox,  ncol, lchnk ) ! HOX produciton (molec/cm3/s)

       extfrc(:ncol,:pver,n2d_ndx) = extfrc(:ncol,:pver,n2d_ndx) +  0.7_r8*spe_nox(:ncol,:pver)
       extfrc(:ncol,:pver,  n_ndx) = extfrc(:ncol,:pver,  n_ndx) + 0.55_r8*spe_nox(:ncol,:pver)
       extfrc(:ncol,:pver, oh_ndx) = extfrc(:ncol,:pver, oh_ndx)         + spe_hox(:ncol,:pver)

       call outfld( 'P_Op',  extfrc(:,:,op_ndx), ncol, lchnk )
       call outfld( 'P_O2p', extfrc(:,:,o2p_ndx), ncol, lchnk )
       call outfld( 'P_Np',  extfrc(:,:,np_ndx), ncol, lchnk )
       call outfld( 'P_N2p', extfrc(:,:,n2p_ndx), ncol, lchnk )
       call outfld( 'P_IONS',extfrc(:,:,n2d_ndx), ncol, lchnk )

    endif


  end subroutine setext

end module mo_setext
