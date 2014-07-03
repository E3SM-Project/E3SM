module mo_setsoa

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use abortutils,   only : endrun

  implicit none
  private
  public :: soa_inti, setsoa, has_soa
  public :: soa_register

  save

  integer, parameter :: NRX = 10                  ! number of SOA forming reactions
  integer, parameter :: NPR = 2                   ! number of products (always 2!)
  integer, target  :: spc_ndx(22)
  integer, pointer :: soam_ndx, soai_ndx, soab_ndx, soat_ndx, soax_ndx
  integer, pointer :: sogm_ndx, sogi_ndx, sogb_ndx, sogt_ndx, sogx_ndx
  integer, pointer :: oc1_ndx, oc2_ndx, c10h16_ndx, isop_ndx, o3_ndx, oh_ndx
  integer, pointer :: no3_ndx, no_ndx, ho2_ndx, tolo2_ndx, beno2_ndx, xylo2_ndx

  integer :: rxn_soa(nrx),react_ndx(NRX,NPR)
  real(r8), dimension(NRX,NPR) :: alpha         ! mass-based stoichiometric coefficients
  real(r8), dimension(NRX,NPR) :: k_om          ! equilibrium gas-particule partition
  real(r8), dimension(NRX)     :: T1, delH      ! Clausium Clayperson parameters (K, J/mol)
  real(r8), dimension(nrx,npr) :: fracsog_init ! mass fraction of each SOA class from each reaction
  real(r8), dimension(nrx,npr) :: fracsoa_init  ! mass fraction of each SOG class from each reaction

  integer :: fracsog_ndx = -1
  integer :: fracsoa_ndx = -1

  integer, pointer :: soa_ndx
  integer, pointer ::  bigalk_ndx, toluene_ndx

  real(r8), dimension(6)   :: bulk_yield               ! total yield of condensable compound (ug/m3/ppm)
  real(r8), dimension(6)   :: fraction                 ! fraction of VOC used in reaction

  real(r8), parameter :: avogadro = 6.023e23_r8        ! Avogadro number
  real(r8), parameter :: Rgas     = 8.314_r8           ! gas constant (J/K/mol)
  real(r8), parameter :: OMscale  = 2.1_r8             ! scaling factor for OM:OC [Turpin and Lim, 2001]
  logical :: has_soa = .false.
  logical :: has_soa_equil = .false.

contains

!===============================================================================
!===============================================================================
subroutine soa_inti(pbuf2d)
  use physics_buffer, only : physics_buffer_desc
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  if ( has_soa_equil) then
     call soa_inti_equil(pbuf2d)
  else
     call soa_inti_old()
  endif
endsubroutine soa_inti
!===============================================================================
!===============================================================================
subroutine setsoa(dt,reaction_rates,tfld,vmr,xhnm,ncol,lchnk,pbuf)
  use physics_buffer, only : physics_buffer_desc
  use chem_mods,      only : gas_pcnst, rxntot
  !-----------------------------------------------------------------------
  !      ... dummy arguments
  !-----------------------------------------------------------------------
  integer, intent(in)      :: ncol                              ! number columns in chunkx  
  integer, intent(in)      :: lchnk                             ! chunk index
  real(r8), intent(in)     :: dt                                ! time step
  real(r8), intent(in)     :: reaction_rates(ncol,pver,rxntot)  ! reaction rates
  real(r8), intent(inout)  :: vmr(ncol,pver,gas_pcnst)          ! xported species ( vmr )
  real(r8), intent(in)     :: tfld(pcols,pver) &                ! temperature (K)
       ,xhnm(ncol,pver)                   ! total atms density (molec/cm**3)
  type(physics_buffer_desc), pointer :: pbuf(:)

  if ( has_soa_equil) then
     call setsoa_equil(dt,reaction_rates,tfld,vmr,xhnm,ncol,lchnk,pbuf)
  else
     call setsoa_old(dt,reaction_rates,tfld,vmr,xhnm,ncol,lchnk)
  endif

end subroutine setsoa
!===============================================================================
!===============================================================================
subroutine soa_inti_old

  use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
  use cam_history,  only : addfld, phys_decomp
  use ppgrid,       only : pver
  use spmd_utils,   only : masterproc

  implicit none

!-----------------------------------------------------------------------      
! 	... check if this is an aerosol simulation
!-----------------------------------------------------------------------      
  if( .not. has_soa ) then 
    return
  end if

!-----------------------------------------------------------------------      
! 	... set reaction indicies
!-----------------------------------------------------------------------      
  rxn_soa(1) = get_rxt_ndx( 'soa1' )
  rxn_soa(2) = get_rxt_ndx( 'soa2' )
  rxn_soa(3) = get_rxt_ndx( 'soa3' )
  rxn_soa(4) = get_rxt_ndx( 'soa4' )
  rxn_soa(5) = get_rxt_ndx( 'soa4' )
  rxn_soa(6) = get_rxt_ndx( 'soa5' )
  if( all( rxn_soa(:) < 1 ) ) then
     has_soa = .false.
     return
  else
     if (masterproc) then
        write(iulog,*) '-----------------------------------------'
        write(iulog,*) 'mozart will do soa aerosols'
        write(iulog,*) '-----------------------------------------'
     endif
  end if

!
! define reactants
!
  react_ndx(1,1) = c10h16_ndx
  react_ndx(1,2) = o3_ndx
  react_ndx(2,1) = c10h16_ndx
  react_ndx(2,2) = oh_ndx
  react_ndx(3,1) = c10h16_ndx
  react_ndx(3,2) = no3_ndx
  react_ndx(4,1) = toluene_ndx
  react_ndx(4,2) = oh_ndx
  react_ndx(5,1) = toluene_ndx
  react_ndx(5,2) = oh_ndx
  react_ndx(6,1) = bigalk_ndx
  react_ndx(6,2) = oh_ndx

  if ( masterproc ) then
     write(iulog,*)'soa_inti ',c10h16_ndx, o3_ndx, oh_ndx, no3_ndx, bigalk_ndx, toluene_ndx
     write(iulog,*)'soa_inti ',soa_ndx, oc1_ndx, oc2_ndx
     write(iulog,*)'soa_inti ',react_ndx
  endif
!
! define partitioning coefficients for each reaction
! bulk yields are from Seinfeld and Pandis (1998)
!
! c10h16 + o3 (from Chung and Seinfeld, JGR, 107, 2002)
!
  alpha(1,1)    = 0.067_r8
  alpha(1,2)    = 0.354_r8
  k_om (1,1)    = 0.184_r8
  k_om (1,2)    = 0.0043_r8
  fraction(1)   = 1._r8
  bulk_yield(1) = 762._r8
!
! c10h16 + oh (from Chung and Seinfeld, JGR, 107, 2002)
!
  alpha(2,1)    = 0.067_r8
  alpha(2,2)    = 0.354_r8
  k_om (2,1)    = 0.184_r8
  k_om (2,2)    = 0.0043_r8
  fraction(2)   = 1._r8
  bulk_yield(2) = 762._r8
! 
! c10h16 + no3 (from Chung and Seinfeld, JGR, 107, 2002)
!
  alpha(3,1)    = 1.000_r8
  alpha(3,2)    = 0.000_r8
  k_om (3,1)    = 0.0163_r8
  k_om (3,2)    = 0.0000_r8
  fraction(3)   = 1._r8
  bulk_yield(3) = 762._r8
!
! toluene + oh : toluene (from Odum et al., Environ. Sci. Technol., 1892, 1997)
!
  alpha(4,1)    = 0.038_r8
  alpha(4,2)    = 0.167_r8
  k_om (4,1)    = 0.042_r8
  k_om (4,2)    = 0.0014_r8
  fraction(4)   = 0.7_r8
  bulk_yield(4) = 424._r8
!
! toluene + oh : m-xylene (from Cocker et al., Atmos. Environ., 6079, 2001)
!
  alpha(5,1)    = 0.120_r8
  alpha(5,2)    = 0.019_r8
  k_om (5,1)    = 0.060_r8
  k_om (5,2)    = 0.010_r8
  fraction(5)   = 0.2_r8
  bulk_yield(5) = 419._r8
!
! bigalk + oh : only for alkanes >= heptane (assume low-yield aromatics as in Lack et al.)
!               (from Odum et al., Environ. Sci. Technol., 1892, 1997)
!
  alpha(6,1)    = 0.071_r8
  alpha(6,2)    = 0.138_r8
  k_om (6,1)    = 0.053_r8
  k_om (6,2)    = 0.0019_r8
  fraction(6)   = 0.1_r8
  bulk_yield(6) = 200._r8
!
  call addfld( 'SOA_PROD', 'kg/kg/s', pver, 'A', 'production of SOA', phys_decomp )

  return
end subroutine soa_inti_old
!
subroutine setsoa_old(dt,reaction_rates,tfld,vmr,xhnm,ncol,lchnk)
!
! secondary organic aerosol for mozart v2.5
!
! based on Lack et al., JGR, 109, D03203, 2004
!
! rewritten by Jean-Francois Lamarque for updated chemical
! mechanism (March 2004)
!
! adapted to CAM (May 2004)
!
  use ppgrid,       only : pcols, pver
  use chem_mods,    only : adv_mass, gas_pcnst, rxntot
  use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
  use cam_history,  only : outfld
  use abortutils,   only : endrun
!
  implicit none
!
!-----------------------------------------------------------------------
!      ... dummy arguments
!-----------------------------------------------------------------------
  integer, intent(in)  :: ncol                              ! number columns in chunkx  
  integer, intent(in)  :: lchnk                             ! chunk index
  real(r8), intent(in)     :: dt                                ! time step
  real(r8), intent(in)     :: reaction_rates(ncol,pver,rxntot)  ! reaction rates
  real(r8), intent(inout)  :: vmr(ncol,pver,gas_pcnst)          ! xported species ( vmr )
  real(r8), intent(in)     :: tfld(pcols,pver) &                ! temperature (K)
                             ,xhnm(ncol,pver)                   ! total atms density (mol/cm**3)

!-----------------------------------------------------------------------
!      ... local variables
!-----------------------------------------------------------------------
  integer :: i,k,n
  real(r8) :: m_0
  real(r8) :: mw_soa,yield,prod,soa_mass
  real(r8) :: soa_prod(ncol,pver)
!
! find molecular weight of SOA
!
  mw_soa = adv_mass(soa_ndx)
!
  do k=1,pver
    do i=1,ncol
!
! calculate initial mass of organic aerosols from OC1 and OC2 
! and convert to ug/m3
!
      m_0 = (vmr(i,k,oc1_ndx)+vmr(i,k,oc2_ndx)) * xhnm(i,k) * adv_mass(oc1_ndx)/avogadro * 1.e12_r8
!
! switch based on a minimum value of m_0.  The bulk approach is
! used to initiate the process
!
      if ( m_0 <= 0.2_r8 ) then
!
! bulk theory
!
        soa_mass = 0._r8
        do n=1,6
!
          if ( rxn_soa(n) <= 0 ) cycle
!
          yield = bulk_yield(n)
!
! define chemical production from gas-phase chemistry
!
          prod  = reaction_rates(i,k,rxn_soa(n)) * fraction(n) &
                * vmr(i,k,react_ndx(n,1)) * vmr(i,k,react_ndx(n,2)) * dt
!
! convert from mixing ratio to ppm
!
          prod = 1e6_r8 * prod
!
! collect into total SOA mass
!
          soa_mass = soa_mass + yield * prod
!
        end do
!
      else
!
! partitioning theory
!
        soa_mass = 0._r8
        do n=1,6
!
          if ( rxn_soa(n) <= 0 ) cycle
!
! define yield from available m_0
!
          yield = soa_yield(m_0,alpha(n,1:2),k_om(n,1:2))
!
! define chemical production from gas-phase chemistry
!
          prod  = reaction_rates(i,k,rxn_soa(n)) * fraction(n) &
                * vmr(i,k,react_ndx(n,1)) * vmr(i,k,react_ndx(n,2)) * dt
!
! convert from mixing ratio to mass (ug/m3)       
!
          prod = prod * xhnm(i,k) * mw_soa/avogadro * 1.e12_r8
!
! collect into total SOA mass
!
          soa_mass = soa_mass + yield * prod
!
        end do
!
      endif
!
! convert from ug/m3 to mixing ratio and update vmr
!
      vmr(i,k,soa_ndx) = vmr(i,k,soa_ndx) + soa_mass * 1.e-12_r8 * avogadro/(mw_soa*xhnm(i,k))
      if ( vmr(i,k,soa_ndx) > 1.e0_r8 ) then
        write(iulog,*)i,k,soa_mass,m_0
        call endrun('soa_yield: vmr(i,k,soa_ndx) > 1.e0_r8')
      endif
!
      soa_prod(i,k) = soa_mass*1.e-12_r8*avogadro/(28.966_r8*xhnm(i,k)*dt)
    end do
  end do
!
  call outfld('SOA_PROD',soa_prod(:ncol,:),ncol, lchnk)
  return
end subroutine setsoa_old
!
real(r8) function soa_yield(m_0,alpha,k)
!
  implicit none
!
  real(r8) :: m_0
  real(r8), dimension(2) :: alpha,k
!
!!$  soa_yield = m_0 * ( ((alpha(1)*k(1))/(1._r8+(alpha(1)*k(1)*m_0))) &
!!$                    + ((alpha(2)*k(2))/(1._r8+(alpha(2)*k(2)*m_0))) )
  soa_yield = m_0 * ( ((alpha(1)*k(1))/(1._r8+k(1)*m_0)) &
                  +   ((alpha(2)*k(2))/(1._r8+k(2)*m_0)) ) 
!
  return
end function soa_yield
!

  !===============================================================================
  !===============================================================================
  subroutine soa_register
    use physics_buffer, only : pbuf_add_field, dtype_r8
    use mo_chem_utls,   only : get_spc_ndx

    !-----------------------------------------------------------------------      
    !       ... set species indices
    !-----------------------------------------------------------------------      

    oc1_ndx     => spc_ndx(1)
    oc2_ndx     => spc_ndx(2)
    soam_ndx    => spc_ndx(3)
    soai_ndx    => spc_ndx(4)
    soat_ndx    => spc_ndx(5)
    soab_ndx    => spc_ndx(6)
    soax_ndx    => spc_ndx(7)
    c10h16_ndx  => spc_ndx(8)
    isop_ndx    => spc_ndx(9)
    tolo2_ndx   => spc_ndx(10)
    beno2_ndx   => spc_ndx(11)
    xylo2_ndx   => spc_ndx(12)
    ho2_ndx     => spc_ndx(13)
    no_ndx      => spc_ndx(14)
    o3_ndx      => spc_ndx(15)
    oh_ndx      => spc_ndx(16)
    no3_ndx     => spc_ndx(17)
    sogm_ndx    => spc_ndx(18)
    sogi_ndx    => spc_ndx(19)
    sogt_ndx    => spc_ndx(20)
    sogb_ndx    => spc_ndx(21)
    sogx_ndx    => spc_ndx(22)

    oc1_ndx      = get_spc_ndx( 'OC1' )
    oc2_ndx      = get_spc_ndx( 'OC2' )
    c10h16_ndx   = get_spc_ndx( 'C10H16')
    isop_ndx     = get_spc_ndx( 'ISOP' )
    tolo2_ndx    = get_spc_ndx( 'TOLO2' )
    beno2_ndx    = get_spc_ndx( 'BENO2' )
    xylo2_ndx    = get_spc_ndx( 'XYLO2' )
    ho2_ndx      = get_spc_ndx( 'HO2' )
    no_ndx       = get_spc_ndx( 'NO' )
    o3_ndx       = get_spc_ndx( 'OX' )
    if( o3_ndx < 1 ) then
       o3_ndx =  get_spc_ndx( 'O3' )
    end if
    oh_ndx       = get_spc_ndx( 'OH' )
    no3_ndx      = get_spc_ndx( 'NO3' )

    soam_ndx     = get_spc_ndx( 'SOAM' )
    soai_ndx     = get_spc_ndx( 'SOAI' )
    soat_ndx     = get_spc_ndx( 'SOAT' )
    soab_ndx     = get_spc_ndx( 'SOAB' )
    soax_ndx     = get_spc_ndx( 'SOAX' )

    sogm_ndx     = get_spc_ndx( 'SOGM' )
    sogi_ndx     = get_spc_ndx( 'SOGI' )
    sogt_ndx     = get_spc_ndx( 'SOGT' )
    sogb_ndx     = get_spc_ndx( 'SOGB' )
    sogx_ndx     = get_spc_ndx( 'SOGX' )

    has_soa_equil = all( spc_ndx(1:22) > 0 )
    has_soa = has_soa_equil

    if ( has_soa_equil ) then
       ! fracsog and fracsoa are added to pbuffer for persistence across restarts
       call pbuf_add_field( 'FRACSOG' ,'global',dtype_r8,(/pcols,pver,nrx,npr/), fracsog_ndx)
       call pbuf_add_field( 'FRACSOA' ,'global',dtype_r8,(/pcols,pver,nrx,npr/), fracsoa_ndx)
    else
       ! reassign these ndx pointers...
       soa_ndx     => spc_ndx(1)
       oc1_ndx     => spc_ndx(2)
       oc2_ndx     => spc_ndx(3)
       c10h16_ndx  => spc_ndx(4)
       o3_ndx      => spc_ndx(5)
       oh_ndx      => spc_ndx(6)
       no3_ndx     => spc_ndx(7)
       bigalk_ndx  => spc_ndx(8)
       toluene_ndx => spc_ndx(9)

       soa_ndx      = get_spc_ndx( 'SOA' )
       oc1_ndx      = get_spc_ndx( 'OC1' )
       oc2_ndx      = get_spc_ndx( 'OC2' )
       c10h16_ndx   = get_spc_ndx( 'C10H16')
       o3_ndx       = get_spc_ndx( 'OX' )
       if( o3_ndx < 1 ) then
          o3_ndx =  get_spc_ndx( 'O3' )
       end if
       oh_ndx       = get_spc_ndx( 'OH' )
       no3_ndx      = get_spc_ndx( 'NO3' )

       bigalk_ndx   = get_spc_ndx( 'BIGALK' )
       toluene_ndx  = get_spc_ndx( 'TOLUENE' )

       has_soa = all( spc_ndx(1:7) > 0 )
    endif

  end subroutine soa_register

  !===============================================================================
  !===============================================================================
  subroutine soa_inti_equil(pbuf2d)

    use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
    use cam_history,  only : addfld, phys_decomp
    use ppgrid,       only : pver
    use spmd_utils,   only : masterproc
    use cam_control_mod, only: nsrest
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_get_field

    implicit none

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_ptr(:)

    integer :: astat
    integer :: i,k,j, c

    real(r8), pointer :: fracsog(:,:,:,:)  ! mass fraction of each SOA class from each reaction
    real(r8), pointer :: fracsoa(:,:,:,:)  ! mass fraction of each SOG class from each reaction

    !-----------------------------------------------------------------------      
    !       ... check if this is an aerosol simulation
    !-----------------------------------------------------------------------      
    if( .not. has_soa ) then 
       return
    end if

    !-----------------------------------------------------------------------      
    !       ... set reaction indicies
    !-----------------------------------------------------------------------      
    rxn_soa(1) = get_rxt_ndx( 'C10H16_O3' )
    rxn_soa(2) = get_rxt_ndx( 'C10H16_OH' )
    rxn_soa(3) = get_rxt_ndx( 'C10H16_NO3' )
    rxn_soa(4) = get_rxt_ndx( 'ISOP_OH' )
    rxn_soa(5) = get_rxt_ndx( 'TOLO2_HO2' )
    rxn_soa(6) = get_rxt_ndx( 'ox_p12' )
    rxn_soa(7) = get_rxt_ndx( 'BENO2_HO2' )
    rxn_soa(8) = get_rxt_ndx( 'BENO2_NO' )
    rxn_soa(9) = get_rxt_ndx( 'XYLO2_HO2' )
    rxn_soa(10) = get_rxt_ndx( 'XYLO2_NO' )
    if( all( rxn_soa(:) < 1 ) ) then
       has_soa = .false.
       return
    else
       if (masterproc) then
          write(*,*) '-----------------------------------------'
          write(*,*) 'mozart will do soa aerosols'
          write(*,*) '-----------------------------------------'
       endif
    end if

    !
    ! define reactants
    !
    react_ndx(1,1) = c10h16_ndx
    react_ndx(1,2) = o3_ndx
    react_ndx(2,1) = c10h16_ndx
    react_ndx(2,2) = oh_ndx
    react_ndx(3,1) = c10h16_ndx
    react_ndx(3,2) = no3_ndx
    react_ndx(4,1) = isop_ndx
    react_ndx(4,2) = oh_ndx
    react_ndx(5,1) = tolo2_ndx
    react_ndx(5,2) = ho2_ndx
    react_ndx(6,1) = tolo2_ndx
    react_ndx(6,2) = no_ndx
    react_ndx(7,1) = beno2_ndx
    react_ndx(7,2) = ho2_ndx
    react_ndx(8,1) = beno2_ndx
    react_ndx(8,2) = no_ndx
    react_ndx(9,1) = xylo2_ndx
    react_ndx(9,2) = ho2_ndx
    react_ndx(10,1) = xylo2_ndx
    react_ndx(10,2) = no_ndx

    if ( masterproc ) then
       print *,'soa_inti ',c10h16_ndx, isop_ndx, tolo2_ndx, beno2_ndx, xylo2_ndx,o3_ndx, oh_ndx, no3_ndx
       print *,'soa_inti ',soam_ndx, soai_ndx, soab_ndx, soat_ndx, soax_ndx, oc1_ndx, oc2_ndx
       print *,'soa_inti ',sogm_ndx, sogi_ndx, sogb_ndx, sogt_ndx, sogx_ndx
       print *,'soa_inti ',react_ndx
    endif
    !
    ! define partitioning coefficients for each reaction
    ! bulk yields are from Seinfeld and Pandis (1998)
    !
    ! c10h16 + o3 (from Chung and Seinfeld, JGR, 107, 2002)
    !
    alpha(1,1)    = 0.067_r8
    alpha(1,2)    = 0.354_r8
    k_om(1,1)     = 0.184_r8
    k_om(1,2)     = 0.0043_r8
    T1(1)         = 310._r8
    delH(1)       = 42.e3_r8

    ! c10h16 + oh (from Chung and Seinfeld, JGR, 107, 2002)
    !
    alpha(2,1)    = 0.067_r8
    alpha(2,2)    = 0.354_r8
    k_om(2,1)     = 0.184_r8
    k_om(2,2)     = 0.0043_r8
    T1(2)         = 310._r8
    delH(2)       = 42.e3_r8
    ! 
    ! c10h16 + no3 (from Chung and Seinfeld, JGR, 107, 2002)
    !
    alpha(3,1)    = 1.000_r8
    alpha(3,2)    = 0.000_r8
    k_om(3,1)     = 0.0163_r8
    k_om(3,2)     = 0.0000_r8
    T1(3)         = 310._r8
    delH(3)       = 42.e3_r8
    !
    !! isop + oh (from Henze and Seinfeld, GRL, 2006): low NOx
    !!
    !  alpha(4,1)    = 0.232_r8
    !  alpha(4,2)    = 0.0288_r8
    !  k_om(4,1)     = 0.00862_r8
    !  k_om(4,2)     = 1.62_r8
    !  T1(4)         = 295._r8
    !  delH(4)       = 42.e3_r8
    !!
    ! isop + oh (from Henze and Seinfeld, GRL, 2006): high NOx
    !
    alpha(4,1)    = 0.264_r8
    alpha(4,2)    = 0.0173_r8
    k_om(4,1)     = 0.00115_r8
    k_om(4,2)     = 1.52_r8
    T1(4)         = 295._r8
    delH(4)       = 42.e3_r8
    !
    ! toluene + oh (pers comm Seinfeld and Henze): low NOx (TOLO2 + HO2) 
    !
    alpha(5,1)    = 0.2349_r8
    alpha(5,2)    = 0.0_r8
    k_om(5,1)     = 1000.0_r8
    k_om(5,2)     = 0.0_r8
    T1(5)         = 295._r8
    delH(5)       = 42.e3_r8
    !
    ! toluene + oh (pers comm Seinfeld and Henze): high NOx (TOLO2 + NO) 
    !
    alpha(6,1)    = 0.0378_r8
    alpha(6,2)    = 0.0737_r8
    k_om(6,1)     = 0.4300_r8
    k_om(6,2)     = 0.0470_r8
    T1(6)         = 295._r8
    delH(6)       = 42.e3_r8
    !
    ! benzene + oh (pers comm Seinfeld and Henze): low NOx (BENO2 + HO2) 
    !
    alpha(7,1)    = 0.2272_r8
    alpha(7,2)    = 0.0_r8
    k_om (7,1)    = 1000.0_r8
    k_om (7,2)    = 0.0_r8
    T1(7)         = 295._r8
    delH(7)       = 42.e3_r8
    !
    ! benzene + oh (pers comm Seinfeld and Henze): high NOx (BENO2 + NO) 
    !
    alpha(8,1)    = 0.0442_r8
    alpha(8,2)    = 0.5454_r8
    k_om(8,1)     = 3.3150_r8
    k_om(8,2)     = 0.0090_r8
    T1(8)         = 295._r8
    delH(8)       = 42.e3_r8
    !
    ! xylene + oh (pers comm Seinfeld and Henze): low NOx (XYLO2 + HO2) 
    !
    alpha(9,1)    = 0.2052_r8
    alpha(9,2)    = 0.0_r8
    k_om(9,1)     = 1000.0_r8
    k_om(9,2)     = 0.0_r8
    T1(9)         = 295._r8
    delH(9)       = 42.e3_r8
    !
    ! xylene + oh (pers comm Seinfeld and Henze): high NOx (XYLO2 + NO) 
    !
    alpha(10,1)    = 0.0212_r8
    alpha(10,2)    = 0.0615_r8
    k_om(10,1)     = 0.7610_r8
    k_om(10,2)     = 0.0290_r8
    T1(10)         = 295._r8
    delH(10)       = 42.e3_r8
    !
    call addfld( 'SOAM_PROD', 'molec/molec/s', pver, 'A', 'production of SOAM', phys_decomp )
    call addfld( 'SOAI_PROD', 'molec/molec/s', pver, 'A', 'production of SOAI', phys_decomp )
    call addfld( 'SOAT_PROD', 'molec/molec/s', pver, 'A', 'production of SOAT', phys_decomp )
    call addfld( 'SOAB_PROD', 'molec/molec/s', pver, 'A', 'production of SOAB', phys_decomp )
    call addfld( 'SOAX_PROD', 'molec/molec/s', pver, 'A', 'production of SOAX', phys_decomp )

    call addfld( 'SOAM_dens', 'ug/m3', pver, 'A', 'density of SOAM', phys_decomp )
    call addfld( 'SOAI_dens', 'ug/m3', pver, 'A', 'density of SOAI', phys_decomp )
    call addfld( 'SOAT_dens', 'ug/m3', pver, 'A', 'density of SOAT', phys_decomp )
    call addfld( 'SOAB_dens', 'ug/m3', pver, 'A', 'density of SOAB', phys_decomp )
    call addfld( 'SOAX_dens', 'ug/m3', pver, 'A', 'density of SOAX', phys_decomp )

    !
    !initialize fracsoa for first timestep and store values for future
    fracsoa_init(1:3,:)=0.2_r8
    fracsoa_init(3,2)=0._r8
    fracsoa_init(4,:)=0.5_r8
    fracsoa_init(5:10,:)=0.33_r8
    fracsoa_init(5,2)=0._r8
    fracsoa_init(7,2)=0._r8
    fracsoa_init(9,2)=0._r8

    !initialize fracsog for first timestep and store values for future
    fracsog_init(1:3,:)=0.2_r8
    fracsog_init(3,2)=0._r8
    fracsog_init(4,:)=0.5_r8
    fracsog_init(5:10,:)=0.33_r8
    fracsog_init(5,2)=0._r8
    fracsog_init(7,2)=0._r8
    fracsog_init(9,2)=0._r8
    
    if (nsrest == 0) then ! initial run
       do c=begchunk, endchunk
          pbuf_ptr=>pbuf_get_chunk(pbuf2d, c)
          call pbuf_get_field(pbuf_ptr, fracsoa_ndx, fracsoa, start=(/1,1,1,1/),kount=(/pcols,pver,nrx,npr/))
          call pbuf_get_field(pbuf_ptr, fracsog_ndx, fracsog )
          do i = 1,pcols
             do k = 1,pver   
                fracsoa( i,k, :,: ) = fracsoa_init(:,:)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine soa_inti_equil
  !===============================================================================
  ! clh (08/01/06): added T dependence of gas-aerosol phase partitioning
  !                 added isoprene as SOA precursor [Henze and Seinfeld, 2006]
  !                 added anthropogenic SOA according to [Henze et al., 2007]
  !                 split SOA into classes (SOAM=monoterpenes,SOAI=isoprene,SOAB=benzene,SOAT=toluene
  !                                         SOAX=xylene)
  !                 fixed error in yield calculation
  !                 added scaling for OC to OM in pre-existing aerosol mass
  !                 modified yield calculation to allow for re-evaporation
  ! Note: do not need to subtract formed aerosol product from reactants, because gas-phase
  !       products do not account for entire mass
  !===============================================================================
  subroutine setsoa_equil(dt,reaction_rates,tfld,vmr,xhnm,ncol,lchnk,pbuf)
    !
    ! updated SOA mechanism 
    ! based on Chung and Seinfeld, JGR, 2002
    !
    ! implemented in CAM by Colette Heald (summer 2007)
    !
    use ppgrid,       only : pcols, pver
    use chem_mods,    only : adv_mass, gas_pcnst, rxntot
    use cam_history,  only : outfld
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field
    !
    implicit none
    !
    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)      :: ncol                              ! number columns in chunkx  
    integer, intent(in)      :: lchnk                             ! chunk index
    real(r8), intent(in)     :: dt                                ! time step
    real(r8), intent(in)     :: reaction_rates(ncol,pver,rxntot)  ! reaction rates
    real(r8), intent(inout)  :: vmr(ncol,pver,gas_pcnst)          ! xported species ( vmr )
    real(r8), intent(in)     :: tfld(pcols,pver) &                ! temperature (K)
         ,xhnm(ncol,pver)                   ! total atms density (molec/cm**3)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------------
    !      ... local variables
    !-----------------------------------------------------------------------
    integer :: i,k,n,p,iter
    real(r8) :: T_fac
    real(r8) :: delHC, sumorg, poa, mnew, numer, maxM, minM, tol, mw_soa,m_air
    real(r8) :: k_om_T(NRX,NPR)                   ! equilibrium gas-particule partition (T dependent)
    real(r8) :: prod(NRX,NPR)                     ! oxidized produced (alpha_i*delHC)
    real(r8) :: sog0(NRX,NPR)                     ! pre-existing SOG
    real(r8) :: soa0(NRX,NPR)                     ! pre-exisiting SOA
    real(r8) :: orggas(NRX,NPR)                   ! intermediate (SOG0+prod of ox)
    real(r8) :: sog(NRX,NPR)                      ! final SOG
    real(r8) :: soa(NRX,NPR)                      ! final SOA
    real(r8), dimension(pcols,pver) :: soam_mass,soai_mass,soat_mass,soab_mass,soax_mass
    real(r8), dimension(pcols,pver) :: sogm_mass,sogi_mass,sogt_mass,sogb_mass,sogx_mass
    real(r8) :: soam_prod(ncol,pver), soai_prod(ncol,pver),soat_prod(ncol,pver),soab_prod(ncol,pver),soax_prod(ncol,pver)
    
    real(r8), pointer, dimension(:,:,:,:) :: fracsog  ! mass fraction of each SOA class from each reaction
    real(r8), pointer, dimension(:,:,:,:) :: fracsoa  ! mass fraction of each SOG class from each reaction
    
    soam_mass(:,:)=0._r8
    soai_mass(:,:)=0._r8
    soat_mass(:,:)=0._r8
    soab_mass(:,:)=0._r8
    soax_mass(:,:)=0._r8
    sogm_mass(:,:)=0._r8
    sogi_mass(:,:)=0._r8
    sogt_mass(:,:)=0._r8
    sogb_mass(:,:)=0._r8
    sogx_mass(:,:)=0._r8

    call pbuf_get_field(pbuf, fracsoa_ndx, fracsoa, start=(/1,1,1,1/),kount=(/ncol,pver,nrx,npr/))
    call pbuf_get_field(pbuf, fracsog_ndx, fracsog, start=(/1,1,1,1/),kount=(/ncol,pver,nrx,npr/))

    do i=1,ncol
       do k=1,pver
          !
          ! INIALIZATION AND LUMPING
          !
          !       calculate mass concentration of air (ug/m3)
          m_air = xhnm(i,k)*28.966_r8/avogadro*1.e12_r8
          !
          !       calculate initial mass of POA from OC1 and OC2 (in ug/m3)
          poa = (vmr(i,k,oc1_ndx)+vmr(i,k,oc2_ndx)) * OMscale * xhnm(i,k) * adv_mass(oc1_ndx)/avogadro * 1.e12_r8 
          !
          !       specify pre-existing SOG/SOA for each class (in ug/m3)


          sog0(1,1)=vmr(i,k,sogm_ndx) * xhnm(i,k) * adv_mass(sogm_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 1,1)
          sog0(1,2)=vmr(i,k,sogm_ndx) * xhnm(i,k) * adv_mass(sogm_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 1,2)
          sog0(2,1)=vmr(i,k,sogm_ndx) * xhnm(i,k) * adv_mass(sogm_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 2,1)
          sog0(2,2)=vmr(i,k,sogm_ndx) * xhnm(i,k) * adv_mass(sogm_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 2,2)
          sog0(3,1)=vmr(i,k,sogm_ndx) * xhnm(i,k) * adv_mass(sogm_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 3,1)
          sog0(3,2)=0._r8
          sog0(4,1)=vmr(i,k,sogi_ndx) * xhnm(i,k) * adv_mass(sogi_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 4,1)
          sog0(4,2)=vmr(i,k,sogi_ndx) * xhnm(i,k) * adv_mass(sogi_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 4,2)
          sog0(5,1)=vmr(i,k,sogt_ndx) * xhnm(i,k) * adv_mass(sogt_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 5,1)
          sog0(5,2)=0._r8
          sog0(6,1)=vmr(i,k,sogt_ndx) * xhnm(i,k) * adv_mass(sogt_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 6,1)
          sog0(6,2)=vmr(i,k,sogt_ndx) * xhnm(i,k) * adv_mass(sogt_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 6,2)
          sog0(7,1)=vmr(i,k,sogb_ndx) * xhnm(i,k) * adv_mass(sogb_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 7,1)
          sog0(7,2)=0._r8
          sog0(8,1)=vmr(i,k,sogb_ndx) * xhnm(i,k) * adv_mass(sogb_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 8,1)
          sog0(8,2)=vmr(i,k,sogb_ndx) * xhnm(i,k) * adv_mass(sogb_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 8,2)
          sog0(9,1)=vmr(i,k,sogx_ndx) * xhnm(i,k) * adv_mass(sogx_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 9,1)
          sog0(9,2)=0._r8
          sog0(10,1)=vmr(i,k,sogx_ndx) * xhnm(i,k) * adv_mass(sogx_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 10,1)
          sog0(10,2)=vmr(i,k,sogx_ndx) * xhnm(i,k) * adv_mass(sogx_ndx)/avogadro * 1.e12_r8 * fracsog(i,k, 10,2)
          !
          soa0(1,1)=vmr(i,k,soam_ndx) * xhnm(i,k) * adv_mass(soam_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 1,1)
          soa0(1,2)=vmr(i,k,soam_ndx) * xhnm(i,k) * adv_mass(soam_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 1,2)
          soa0(2,1)=vmr(i,k,soam_ndx) * xhnm(i,k) * adv_mass(soam_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 2,1)
          soa0(2,2)=vmr(i,k,soam_ndx) * xhnm(i,k) * adv_mass(soam_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 2,2)
          soa0(3,1)=vmr(i,k,soam_ndx) * xhnm(i,k) * adv_mass(soam_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 3,1)
          soa0(3,2)=0._r8
          soa0(4,1)=vmr(i,k,soai_ndx) * xhnm(i,k) * adv_mass(soai_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 4,1)
          soa0(4,2)=vmr(i,k,soai_ndx) * xhnm(i,k) * adv_mass(soai_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 4,2)
          soa0(5,1)=vmr(i,k,soat_ndx) * xhnm(i,k) * adv_mass(soat_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 5,1)
          soa0(5,2)=0._r8
          soa0(6,1)=vmr(i,k,soat_ndx) * xhnm(i,k) * adv_mass(soat_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 6,1)
          soa0(6,2)=vmr(i,k,soat_ndx) * xhnm(i,k) * adv_mass(soat_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 6,2)
          soa0(7,1)=vmr(i,k,soab_ndx) * xhnm(i,k) * adv_mass(soab_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 7,1)
          soa0(7,2)=0._r8
          soa0(8,1)=vmr(i,k,soab_ndx) * xhnm(i,k) * adv_mass(soab_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 8,1)
          soa0(8,2)=vmr(i,k,soab_ndx) * xhnm(i,k) * adv_mass(soab_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 8,2)
          soa0(9,1)=vmr(i,k,soax_ndx) * xhnm(i,k) * adv_mass(soax_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 9,1)
          soa0(9,2)=0._r8
          soa0(10,1)=vmr(i,k,soax_ndx) * xhnm(i,k) * adv_mass(soax_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 10,1)
          soa0(10,2)=vmr(i,k,soax_ndx) * xhnm(i,k) * adv_mass(soax_ndx)/avogadro * 1.e12_r8 * fracsoa(i,k, 10,2)
          !
          !
          !--------------------------------------------------------------
          ! CHEMISTRY
          !
          sumorg=0._r8
          !
          do n=1,NRX
             !
             !         temperature dependence of paritioning via Clausiaus Clayperon [C&S, 2002]: clh (02/24/06)
             T_fac = tfld(i,k)/T1(n) * exp(delH(n)/Rgas*(1/tfld(i,k) - 1/T1(n)))
             k_om_T(n,:)=k_om(n,:)*T_fac
             !
             !         find molecular weight of SOA
             if ( (n >=1) .and. (n < 4) ) then
                mw_soa = adv_mass(soam_ndx)
             else if (n == 4) then
                mw_soa = adv_mass(soai_ndx)
             else if ( (n > 4) .and. (n < 7) ) then
                mw_soa = adv_mass(soat_ndx)
             else if ( (n >= 7) .and. (n < 9) ) then
                mw_soa = adv_mass(soab_ndx)
             else if ( (n >= 9) .and. (n < 11) ) then
                mw_soa = adv_mass(soax_ndx)
             end if
             !
             !         define chemical production from gas-phase chemistry (and convert to ug/m3)
             delHC  = reaction_rates(i,k,rxn_soa(n)) &
                  * vmr(i,k,react_ndx(n,1)) * vmr(i,k,react_ndx(n,2))* dt &
                  * xhnm(i,k) * mw_soa/avogadro * 1.e12_r8
             !
             !         specify the new total of gas-phase products (before re-partitioning)
             !         and total up all SOA in gas and aerosol phase
             do p=1,nPR
                orggas(n,p) = sog0(n,p)+alpha(n,p)*delHC
                sumorg=sumorg+orggas(n,p)+soa0(n,p)
                prod(n,p)=alpha(n,p)*delHC
             enddo
             !
          enddo
          !
          !--------------------------------------------------------------
          !       check to see if no organics no partitioning!
          !             (set here to previous timestep concentration+oxprod to preserve mass)
          if (sumorg < 1.e-10_r8) then
             do n=1,NRX
                do p=1,NPR
                   sog(n,p)=sog0(n,p)
                   soa(n,p)=soa0(n,p)+prod(n,p)
                enddo
             enddo
          else
             ! PARTITION (Need to iteratively solve for MNEW, set tolerances to 0.1 ng/m3 or 1% of MNEW)
             !
             numer=0._r8
             maxM=0._r8
             !           If POA is essentially zero, equations simplify
             if (POA < 1.e-10_r8) then
                do n=1,NRX
                   do p=1,NPR
                      numer = numer + k_om_T(n,p)* (orggas(n,p)+soa0(n,p))
                   enddo
                enddo
                !              if numerator is less than 1 then MNEW must be zero
                if (numer <= 1._r8) then 
                   mnew=0._r8
                   iter=0
                else
                   minM = 1.e-40_r8
                   maxM = poa + sumorg
                   tol = 1.e-4_r8
                   mnew = zeroin(minM,maxM,tol,poa,soa0,orggas,k_om_T,sumorg,iter)
                end if
                !
                !           if additional organic mass is less than 1% of POA, or difference between the
                !           two is less than 0.1 ng/m3 (the tolerance) then POA is partitioning mass
             else if ( (sumorg < 0.01_r8*poa) .or. (abs(sumorg-poa) < 1.e-4_r8) ) then
                mnew = poa
                iter = 0
                !           otherwise solve for MNEW iteratively on interval [poa, poa+sumorg]
             else
                maxM = poa+sumorg
                minM = poa
                tol = 1.e-4_r8
                mnew = zeroin(minM,maxM,tol,poa,soa0,orggas,k_om_T,sumorg,iter) 
             end if
             !
             !
             !           Now equilibrium partitioning with new MNEW
             !           If no MNEW then all the SOA evaporates to gas-phase
             if (mnew > 0._r8) then 
                do n=1,NRX
                   do p=1,NPR
                      soa(n,p)=k_om_T(n,p)*mnew*(orggas(n,p)+soa0(n,p))/(1._r8+k_om_T(n,p)*mnew)
                      if (k_om_T(n,p) /= 0._r8) then
                         sog(n,p)=soa(n,p)/(k_om_T(n,p)*mnew)
                      else
                         sog(n,p)=0._r8
                      end if
                   enddo
                enddo

             else 
                do n=1,NRX
                   do p=1,NPR
                      sog(n,p)=orggas(n,p)+soa0(n,p)
                      soa(n,p)=1.e-20_r8
                   enddo
                enddo
             end if
             !
          end if
          !
          !--------------------------------------------------------------
          ! LUMP INTO ARRAYS
          do n=1,NRX
             do p=1,NPR
                if ( (n >=1) .and. (n < 4) ) then
                   soam_mass(i,k) = soam_mass(i,k) + soa(n,p)
                   sogm_mass(i,k) = sogm_mass(i,k) + sog(n,p)
                else if (n == 4) then
                   soai_mass(i,k) = soai_mass(i,k) + soa(n,p)
                   sogi_mass(i,k) = sogi_mass(i,k) + sog(n,p)
                else if ( (n > 4) .and. (n < 7) ) then
                   soat_mass(i,k) = soat_mass(i,k) + soa(n,p)
                   sogt_mass(i,k) = sogt_mass(i,k) + sog(n,p)
                else if ( (n >= 7) .and. (n < 9) ) then
                   soab_mass(i,k) = soab_mass(i,k) + soa(n,p)
                   sogb_mass(i,k) = sogb_mass(i,k) + sog(n,p)
                else if ( (n >= 9) .and. (n < 11) ) then
                   soax_mass(i,k) = soax_mass(i,k) + soa(n,p)
                   sogx_mass(i,k) = sogx_mass(i,k) + sog(n,p)
                end if
             enddo
          enddo
          !
          !       Save mass fraction of each SOA rxn to SOA class
          !       (but if sumorg essentially zero revert to init fracs OR
          !        if mnew = 0 (all evap) then fracsoa revert to old, calculate fracsog only)
          if (sumorg < 1.e-10_r8) then
             do n=1,NRX
                do p=1,NPR
                   fracsoa(i,k,n,p)=fracsoa_init(n,p)
                   fracsog(i,k,n,p)=fracsog_init(n,p)
                enddo
             enddo
          else
             if (mnew < 1.e-20_r8) then
                do n=1,NRX
                   do p=1,NPR
                      fracsoa(i,k,n,p)=fracsoa_init(n,p)
                      if ( (n >=1) .and. (n < 4) ) then
                         if (sogm_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogm_mass(i,k)
                         end if
                      else if (n == 4) then
                         if (sogi_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogi_mass(i,k)
                         end if
                      else if ( (n > 4) .and. (n < 7) ) then
                         if (sogt_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogt_mass(i,k)
                         end if
                      else if ( (n >= 7) .and. (n < 9) ) then
                         if (sogb_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogb_mass(i,k)
                         end if
                      else if ( (n >= 9) .and. (n < 11) ) then
                         if (sogx_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogx_mass(i,k)
                         end if
                      end if
                      if ( (p==2) .and. (n==3 .or. n==5 .or. n==7 .or. n==9) ) then
                         fracsog(i,k,n,p)=0._r8
                      end if
                   enddo
                enddo
             else
                do n=1,NRX
                   do p=1,NPR
                      if ( (n >=1) .and. (n < 4) ) then
                         if (soam_mass(i,k) == 0._r8) then
                            fracsoa(i,k,n,p)=fracsoa_init(n,p)
                         else
                            fracsoa(i,k,n,p)=soa(n,p)/soam_mass(i,k)
                         end if
                         if (sogm_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogm_mass(i,k)
                         end if
                      else if (n == 4) then
                         if (soai_mass(i,k) == 0._r8) then
                            fracsoa(i,k,n,p)=fracsoa_init(n,p)
                         else
                            fracsoa(i,k,n,p)=soa(n,p)/soai_mass(i,k)
                         end if
                         if (sogi_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogi_mass(i,k)
                         end if
                      else if ( (n > 4) .and. (n < 7) ) then
                         if (soat_mass(i,k) == 0._r8) then
                            fracsoa(i,k,n,p)=fracsoa_init(n,p)
                         else
                            fracsoa(i,k,n,p)=soa(n,p)/soat_mass(i,k)
                         end if
                         if (sogt_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogt_mass(i,k)
                         end if
                      else if ( (n >= 7) .and. (n < 9) ) then
                         if (soab_mass(i,k) == 0._r8) then
                            fracsoa(i,k,n,p)=fracsoa_init(n,p)
                         else
                            fracsoa(i,k,n,p)=soa(n,p)/soab_mass(i,k)
                         end if
                         if (sogb_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogb_mass(i,k)
                         end if
                      else if ( (n >= 9) .and. (n < 11) ) then
                         if (soax_mass(i,k) == 0._r8) then
                            fracsoa(i,k,n,p)=fracsoa_init(n,p)
                         else
                            fracsoa(i,k,n,p)=soa(n,p)/soax_mass(i,k)
                         end if
                         if (sogx_mass(i,k) == 0._r8) then
                            fracsog(i,k,n,p)=fracsog_init(n,p)
                         else
                            fracsog(i,k,n,p)=sog(n,p)/sogx_mass(i,k)
                         end if
                      end if
                      if ( (p==2) .and. (n==3 .or. n==5 .or. n==7 .or. n==9) ) then
                         fracsoa(i,k,n,p)=0._r8
                         fracsog(i,k,n,p)=0._r8
                      end if
                   enddo
                enddo
             end if
          end if
          !
          !--------------------------------------------------------------
          !
          ! calculate NET production in kg/kg/s (subtract initial mass) 
          !
          soam_prod(i,k) = ( soam_mass(i,k)*1.e-12_r8*avogadro/(adv_mass(soam_ndx)*xhnm(i,k)) - &
                           vmr(i,k,soam_ndx) )/dt
          soai_prod(i,k) = ( soai_mass(i,k)*1.e-12_r8*avogadro/(adv_mass(soai_ndx)*xhnm(i,k)) - &
                           vmr(i,k,soai_ndx) )/dt
          soat_prod(i,k) = ( soat_mass(i,k)*1.e-12_r8*avogadro/(adv_mass(soat_ndx)*xhnm(i,k)) - &
                           vmr(i,k,soat_ndx) )/dt
          soab_prod(i,k) = ( soab_mass(i,k)*1.e-12_r8*avogadro/(adv_mass(soab_ndx)*xhnm(i,k)) - &
                           vmr(i,k,soab_ndx) )/dt
          soax_prod(i,k) = ( soax_mass(i,k)*1.e-12_r8*avogadro/(adv_mass(soax_ndx)*xhnm(i,k)) - &
                           vmr(i,k,soax_ndx) )/dt
          !
          ! convert from ug/m3 to mixing ratio and update vmr
          !
          vmr(i,k,soam_ndx) = soam_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(soam_ndx)*xhnm(i,k))
          vmr(i,k,soai_ndx) = soai_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(soai_ndx)*xhnm(i,k))
          vmr(i,k,soat_ndx) = soat_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(soat_ndx)*xhnm(i,k))
          vmr(i,k,soab_ndx) = soab_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(soab_ndx)*xhnm(i,k))
          vmr(i,k,soax_ndx) = soax_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(soax_ndx)*xhnm(i,k))
          vmr(i,k,sogm_ndx) = sogm_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(sogm_ndx)*xhnm(i,k))
          vmr(i,k,sogi_ndx) = sogi_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(sogi_ndx)*xhnm(i,k))
          vmr(i,k,sogt_ndx) = sogt_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(sogt_ndx)*xhnm(i,k))
          vmr(i,k,sogb_ndx) = sogb_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(sogb_ndx)*xhnm(i,k))
          vmr(i,k,sogx_ndx) = sogx_mass(i,k) * 1.e-12_r8 * avogadro/(adv_mass(sogx_ndx)*xhnm(i,k))
       enddo
    enddo
    !
    call outfld('SOAM_PROD',soam_prod(:ncol,:),ncol, lchnk)
    call outfld('SOAI_PROD',soai_prod(:ncol,:),ncol, lchnk)
    call outfld('SOAT_PROD',soat_prod(:ncol,:),ncol, lchnk)
    call outfld('SOAB_PROD',soab_prod(:ncol,:),ncol, lchnk)
    call outfld('SOAX_PROD',soax_prod(:ncol,:),ncol, lchnk)

    call outfld('SOAM_dens',soam_mass(:ncol,:), ncol, lchnk)
    call outfld('SOAI_dens',soai_mass(:ncol,:), ncol, lchnk)
    call outfld('SOAT_dens',soat_mass(:ncol,:), ncol, lchnk)
    call outfld('SOAB_dens',soab_mass(:ncol,:), ncol, lchnk)
    call outfld('SOAX_dens',soax_mass(:ncol,:), ncol, lchnk)

    !
    !
    return
  end subroutine setsoa_equil

  !===============================================================================
  !===============================================================================
  real(r8) function zeroin(x1,x2,tol,poa,aer,gas,k,totorg,iter)
    ! function to iteratively solve function using bilinear method
    !
    implicit none
    !
    integer             :: iter
    real(r8),intent(in) :: x1, x2                      ! min/max of interval
    real(r8),intent(in) :: tol                         ! tolerance (interval of uncertainty)
    real(r8),intent(in) :: poa,totorg                         
    real(r8),intent(in) :: aer(NRX,NPR), gas(NRX,NPR)  ! aerosol and gas phase concentrations
    real(r8),intent(in) :: k(NRX,NPR)                  ! partitioning coeff
    ! local vars
    real(r8)            :: xa,xb,xm,fa,fb,fm
    !
    xa=x1
    xb=x2
    xm=0._r8
    fa=soa_function(xa,poa,aer,gas,k)
    fb=soa_function(xb,poa,aer,gas,k)
    !
    ! check that functions have opposite signs
    if (fa >= 0._r8) then
       if (fb >=0._r8) then 
          write(iulog,*) 'ABORT IN ZEROIN: SAME SIGN ON FUNCTION',poa,totorg,x1,x2,fa,fb,aer,gas,k
          call endrun
       end if
    else
       if (fb <=0._r8) then
          write(iulog,*) 'ABORT IN ZEROIN: SAME SIGN ON FUNCTION',poa,totorg,x1,x2,fa,fb,aer,gas,k
          call endrun
       end if
    end if
    !
    iter=0
    do while ((abs(xa-xb) > 2._r8*tol) .and. (abs(xa-xb) > 0.01_r8*xm) )
       xm=(xa+xb)/2
       fm=soa_function(xm,poa,aer,gas,k)
       if (fa >=0._r8) then
          if (fm >=0._r8) then 
             xa=xm
             fa=fm
          else
             xb=xm
             fb=fm
          end if
       else
          if (fm < 0._r8) then
             xa=xm
             fa=fm
          else
             xb=xm
             fb=fm
          end if
       end if
       iter=iter+1        
    enddo
    !
    zeroin = (xa+xb)/2
    !
    return
  end function zeroin
  !===============================================================================
  !===============================================================================
  real(r8) function soa_function(m0,poa,aer,gas,k)
    ! function which calculates SOAeqn (trying to minimize to zero)
    !
    implicit none
    !
    real(r8),intent(in) :: m0,poa                         
    real(r8),intent(in) :: aer(NRX,NPR), gas(NRX,NPR)  ! aerosol and gas phase concentrations
    real(r8),intent(in) :: k(NRX,NPR)                  ! partitioning coeff
    ! local vars
    integer             :: n,p
    real(r8)            :: value
    !
    value=0._r8
    !
    do n=1,NRX
       do p=1,NPR
          value = value + k(n,p)*(gas(n,p)+aer(n,p))/(1._r8 + k(n,p)*m0)
       enddo
    enddo
    soa_function = value + (poa/m0) - 1._r8
    !
    !  write(iulog,*) 'clh soa_function out: ',m0,soa_function
    return
  end function soa_function
  !===============================================================================

end module mo_setsoa
