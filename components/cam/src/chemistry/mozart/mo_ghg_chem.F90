!-----------------------------------------------------------------
! Manages the reaction rates of the green house gas species.
! This is used with the reduced ghg chemical mechanism.
!
! Created by: Francis Vitt -- 20 Aug 2008
!-----------------------------------------------------------------
module mo_ghg_chem

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use boundarydata,   only : boundarydata_type, boundarydata_init, boundarydata_update
  use physics_types,  only : physics_state
  use cam_abortutils,     only : endrun
  use ppgrid,         only : pcols, pver, begchunk, endchunk

  implicit none

  private
  save

  public :: ghg_chem_set_rates
  public :: ghg_chem_set_flbc
  public :: ghg_chem_init
  public :: ghg_chem_timestep_init
  public :: ghg_chem_final

  integer, parameter :: ncnst=4                      ! number of constituents
  type(boundarydata_type) :: chemdata
  character(len=6), dimension(ncnst), parameter :: nc_names = & ! constituent names
       (/'TN2O  ', 'TCH4  ', 'TCFC11', 'TCFC12'/)

  integer :: n2o_rxt, ch4_rxt, cfc11_rxt, cfc12_rxt, lyman_alpha_rxt
  integer :: n2o_ndx, ch4_ndx, cfc11_ndx, cfc12_ndx
  integer :: ghg_ndx(ncnst)
  character(len=6) :: ghg_bnd_names(ncnst)

  logical :: lyman_alpha = .false.
  type(boundarydata_type) :: h2orate_data
  character(len=4), parameter :: h2orate_name = 'jh2o'

contains

!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine ghg_chem_init(phys_state, bndtvg, h2orates)
    use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
    use cam_history,  only : addfld

    implicit none

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    character(len=*),    intent(in) :: bndtvg ! pathname for greenhouse gas loss rate
    character(len=*),optional,    intent(in) :: h2orates ! lyman-alpha h2o loss rates

    integer    ::  ids(8)
    integer :: m,mm

    n2o_rxt = get_rxt_ndx( 'n2o_loss' )
    ch4_rxt = get_rxt_ndx( 'ch4_loss' )
    cfc11_rxt = get_rxt_ndx( 'cfc11_loss' )
    cfc12_rxt = get_rxt_ndx( 'cfc12_loss' )
    lyman_alpha_rxt = get_rxt_ndx( 'lyman_alpha' )

    n2o_ndx = get_spc_ndx('N2O')
    ch4_ndx = get_spc_ndx('CH4')
    cfc11_ndx = get_spc_ndx('CFC11')
    cfc12_ndx = get_spc_ndx('CFC12')

    ids(1) = n2o_rxt   
    ids(2) = ch4_rxt   
    ids(3) = cfc11_rxt 
    ids(4) = cfc12_rxt 
    ids(5) = n2o_ndx   
    ids(6) = ch4_ndx   
    ids(7) = cfc11_ndx 
    ids(8) = cfc12_ndx 

    if( any( ids < 1 ) ) then
       call endrun('need to configure with ghg chemistry mechanism')
    endif

    call boundarydata_init(bndtvg,phys_state,nc_names,ncnst,chemdata)

    if ( present( h2orates ) ) then
       if ( len_trim( h2orates ) > 0 .and. lyman_alpha_rxt > 0 ) then
          lyman_alpha = .true.
          call boundarydata_init(h2orates,phys_state,(/h2orate_name/),1,h2orate_data)
       endif
    endif

    call addfld( 'GHG_CFC11_R', (/ 'lev' /), 'I', '1/sec', 'prescribed cfc11 loss rate for ghg chem' )
    call addfld( 'GHG_CFC12_R', (/ 'lev' /), 'I', '1/sec', 'prescribed cfc12 loss rate for ghg chem' )
    call addfld( 'GHG_N2O_R', (/ 'lev' /), 'I',   '1/sec', 'prescribed n2o loss rate for ghg chem' )
    call addfld( 'GHG_CH4_R', (/ 'lev' /), 'I',   '1/sec', 'prescribed ch4 loss rate for ghg chem' )
    call addfld( 'GHG_H2O_R', (/ 'lev' /), 'I',   '1/sec', 'prescribed h2o loss rate for ghg chem' )

    ghg_ndx =       (/  n2o_ndx,  ch4_ndx, cfc11_ndx, cfc12_ndx /)
    ghg_bnd_names = (/ 'N2OVMR', 'CH4VMR',  'F11VMR',  'F12VMR' /)

  end subroutine ghg_chem_init

!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine ghg_chem_timestep_init(phys_state)
    implicit none

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    call boundarydata_update(phys_state,chemdata)
    if (lyman_alpha) then
       call boundarydata_update(phys_state,h2orate_data)
    endif
  end subroutine ghg_chem_timestep_init

!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine ghg_chem_set_rates( rxn_rates, latmapback, zen_angle, ncol, lchnk )
    use chem_mods,    only : rxntot
    use chem_mods,    only : gas_pcnst
    use cam_history,  only : outfld
    use mo_constants, only : pi

    implicit none

    integer,  intent(in)    :: ncol                           ! number columns in chunk
    real(r8), intent(inout) :: rxn_rates(ncol,pver,rxntot)   ! ghg loss rates
    integer,  intent(in)    :: latmapback(pcols)
    real(r8), intent(in)    :: zen_angle(ncol)
    integer,  intent(in)    :: lchnk                          ! chunk index

    integer :: i,k
    real(r8), parameter :: half_pi = pi/2._r8 

    do k=1,pver-2
       do i=1,ncol
          rxn_rates(i,k,n2o_rxt)   = chemdata%datainst(latmapback(i),k,lchnk,1)
          rxn_rates(i,k,ch4_rxt)   = chemdata%datainst(latmapback(i),k,lchnk,2)
          rxn_rates(i,k,cfc11_rxt) = chemdata%datainst(latmapback(i),k,lchnk,3)
          rxn_rates(i,k,cfc12_rxt) = chemdata%datainst(latmapback(i),k,lchnk,4)
       enddo
    enddo

    rxn_rates(:ncol,pver-1:pver,n2o_rxt)   = 0._r8
    rxn_rates(:ncol,pver-1:pver,ch4_rxt)   = 0._r8
    rxn_rates(:ncol,pver-1:pver,cfc11_rxt) = 0._r8
    rxn_rates(:ncol,pver-1:pver,cfc12_rxt) = 0._r8

    call outfld( 'GHG_CFC11_R', rxn_rates(:ncol,:,cfc11_rxt), ncol, lchnk )
    call outfld( 'GHG_CFC12_R', rxn_rates(:ncol,:,cfc12_rxt), ncol, lchnk )
    call outfld( 'GHG_N2O_R',   rxn_rates(:ncol,:,n2o_rxt),   ncol, lchnk )
    call outfld( 'GHG_CH4_R',   rxn_rates(:ncol,:,n2o_rxt),   ncol, lchnk )
    
    if (lyman_alpha_rxt > 0) then
       rxn_rates(:ncol,:,lyman_alpha_rxt) = 0._r8 
    endif

    if (lyman_alpha) then
       do i=1,ncol
          if (zen_angle(i) < half_pi) then
             rxn_rates(i,:,lyman_alpha_rxt) = h2orate_data%datainst(latmapback(i),:,lchnk,1)
          endif
       enddo

       call outfld( 'GHG_H2O_R', rxn_rates(:ncol,:,lyman_alpha_rxt), ncol, lchnk )
    endif

  endsubroutine ghg_chem_set_rates

!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine ghg_chem_set_flbc( vmr, ncol )
    use chem_surfvals, only : chem_surfvals_get
    use chem_mods,     only : gas_pcnst
    use mo_flbc,       only : has_flbc
    implicit none

    integer,  intent(in)    :: ncol                           ! number columns in chunk
    real(r8), intent(inout) ::  vmr(ncol,pver,gas_pcnst)              ! xported species (vmr)
    integer  :: i,ndx

    do i = 1,ncnst
       ndx = ghg_ndx(i)
       if ( has_flbc(ndx)) then
          vmr(:ncol, pver-1, ndx) = vmr(:ncol, pver, ndx)
       else
          vmr(:ncol, pver-1:pver, ndx) = chem_surfvals_get(ghg_bnd_names(i))
       endif
    enddo

  endsubroutine ghg_chem_set_flbc

!-----------------------------------------------------------------
!-----------------------------------------------------------------
  subroutine ghg_chem_final
    implicit none

    deallocate(chemdata%fields)
    deallocate(chemdata%datainst)
    deallocate(chemdata%cdates)
    deallocate(chemdata%lat)
    deallocate(chemdata%zi)

    if (lyman_alpha) then

       deallocate(h2orate_data%fields)
       deallocate(h2orate_data%datainst)
       deallocate(h2orate_data%cdates)
       deallocate(h2orate_data%lat)
       deallocate(h2orate_data%zi)

    endif

  end subroutine ghg_chem_final

end module mo_ghg_chem
