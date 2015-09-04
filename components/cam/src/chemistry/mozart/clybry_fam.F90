!-----------------------------------------------------------------------
!
! Manages the adjustment of ClOy and BrOy family components in response
! to conservation issues resulting from advection.
!
! Created by: Francis Vitt
! Date: 21 May 2008
! Modified by Stacy Walters
! Date: 13 August 2008
!-----------------------------------------------------------------------

module clybry_fam

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use ppgrid,        only : pcols, pver
  use chem_mods,     only : gas_pcnst, adv_mass
  use constituents,  only : pcnst
  use short_lived_species,only: set_short_lived_species,get_short_lived_species

  implicit none

  save

  private
  public :: clybry_fam_set
  public :: clybry_fam_adj
  public :: clybry_fam_init

  integer :: id_cly,id_bry

  integer :: id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2
  integer :: id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr

  logical :: has_clybry

contains

  !------------------------------------------
  !------------------------------------------
  subroutine clybry_fam_init

    use mo_chem_utls, only : get_spc_ndx
    implicit none

    integer :: ids(16)

    id_cly = get_spc_ndx('CLY')
    id_bry = get_spc_ndx('BRY')

    id_cl = get_spc_ndx('CL')
    id_clo = get_spc_ndx('CLO')
    id_hocl = get_spc_ndx('HOCL')
    id_cl2 = get_spc_ndx('CL2')
    id_cl2o2 = get_spc_ndx('CL2O2')
    id_oclo = get_spc_ndx('OCLO')
    id_hcl = get_spc_ndx('HCL')
    id_clono2 = get_spc_ndx('CLONO2')

    id_br = get_spc_ndx('BR')
    id_bro = get_spc_ndx('BRO')
    id_hbr = get_spc_ndx('HBR')
    id_brono2 = get_spc_ndx('BRONO2')
    id_brcl = get_spc_ndx('BRCL')
    id_hobr = get_spc_ndx('HOBR')

    ids = (/ id_cly,id_bry, &
             id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2, &
             id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr /)

    has_clybry = all( ids(:) > 0 )

  endsubroutine clybry_fam_init

!--------------------------------------------------------------
! set the ClOy and BrOy mass mixing ratios
!  - this is call before advection
!--------------------------------------------------------------
  subroutine clybry_fam_set( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : get_nstep
    use physics_buffer, only : physics_buffer_desc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8) :: wrk(ncol,pver,2)
    real(r8) :: mmr(pcols,pver,gas_pcnst)
    integer  :: n, m

    if (.not. has_clybry) return

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          mmr(:ncol,:,m) = q(:ncol,:, n)
       endif
    enddo
    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

!--------------------------------------------------------------
!       ... form updated chlorine, bromine atom mass mixing ratios
!--------------------------------------------------------------
    wrk(:,:,1) = cloy( mmr, pcols, ncol )
    wrk(:,:,2) = broy( mmr, pcols, ncol )

    mmr(:ncol,:,id_cly) = wrk(:,:,1)
    mmr(:ncol,:,id_bry) = wrk(:,:,2)

    call set_short_lived_species( mmr, lchnk, ncol, pbuf )
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          q(:ncol,:, n) = mmr(:ncol,:,m)
       endif
    enddo

  end subroutine clybry_fam_set

!--------------------------------------------------------------
! adjust the ClOy and BrOy individual family members 
!  - this is call after advection
!--------------------------------------------------------------
  subroutine clybry_fam_adj( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : is_first_step
    use physics_buffer, only : physics_buffer_desc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: factor(ncol,pver)
    real(r8) :: wrk(ncol,pver)
    real(r8) :: mmr(pcols,pver,gas_pcnst)

    integer  :: n, m

    if (.not. has_clybry) return

!--------------------------------------------------------------
!       ... CLY,BRY are not adjusted until the end of the first timestep
!--------------------------------------------------------------
    if (is_first_step()) return

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          mmr(:ncol,:,m) = q(:ncol,:, n)
       endif
    enddo
    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

!--------------------------------------------------------------
!       ... form updated chlorine atom mass mixing ratio
!--------------------------------------------------------------
    wrk(:,:) = cloy( mmr, pcols, ncol )

    factor(:ncol,:) = mmr(:ncol,:,id_cly) / wrk(:ncol,:)
!--------------------------------------------------------------
!       ... adjust "group" members
!--------------------------------------------------------------
    mmr(:ncol,:,id_cl)     = factor(:ncol,:)*mmr(:ncol,:,id_cl)
    mmr(:ncol,:,id_clo)    = factor(:ncol,:)*mmr(:ncol,:,id_clo)
    mmr(:ncol,:,id_hocl)   = factor(:ncol,:)*mmr(:ncol,:,id_hocl)
    mmr(:ncol,:,id_cl2)    = factor(:ncol,:)*mmr(:ncol,:,id_cl2)
    mmr(:ncol,:,id_cl2o2)  = factor(:ncol,:)*mmr(:ncol,:,id_cl2o2)
    mmr(:ncol,:,id_oclo)   = factor(:ncol,:)*mmr(:ncol,:,id_oclo)
    mmr(:ncol,:,id_hcl)    = factor(:ncol,:)*mmr(:ncol,:,id_hcl)
    mmr(:ncol,:,id_clono2) = factor(:ncol,:)*mmr(:ncol,:,id_clono2)

!--------------------------------------------------------------
!        ... form updated bromine atom mass mixing ratio
!--------------------------------------------------------------
    wrk(:,:) = broy( mmr, pcols, ncol )

    factor(:ncol,:) = mmr(:ncol,:,id_bry) / wrk(:ncol,:)
!--------------------------------------------------------------
!       ... adjust "group" members
!--------------------------------------------------------------
    mmr(:ncol,:,id_br)     = factor(:ncol,:)*mmr(:ncol,:,id_br)
    mmr(:ncol,:,id_bro)    = factor(:ncol,:)*mmr(:ncol,:,id_bro)
    mmr(:ncol,:,id_hbr)    = factor(:ncol,:)*mmr(:ncol,:,id_hbr)
    mmr(:ncol,:,id_brono2) = factor(:ncol,:)*mmr(:ncol,:,id_brono2)
    mmr(:ncol,:,id_brcl)   = factor(:ncol,:)*mmr(:ncol,:,id_brcl)
    mmr(:ncol,:,id_hobr)   = factor(:ncol,:)*mmr(:ncol,:,id_hobr)

    call set_short_lived_species( mmr, lchnk, ncol, pbuf )
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          q(:ncol,:, n) = mmr(:ncol,:,m)
       endif
    enddo

  end subroutine clybry_fam_adj

!--------------------------------------------------------------
! private methods
!--------------------------------------------------------------

!--------------------------------------------------------------
! compute the mass mixing retio of ClOy
!--------------------------------------------------------------
  function cloy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: cloy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_cl)           /adv_mass(id_cl) &
              + q(:ncol,k,id_clo)          /adv_mass(id_clo) &
              + q(:ncol,k,id_hocl)         /adv_mass(id_hocl) &
              + 2._r8*( q(:ncol,k,id_cl2)  /adv_mass(id_cl2) &
                      + q(:ncol,k,id_cl2o2)/adv_mass(id_cl2o2) ) &
              + q(:ncol,k,id_oclo)         /adv_mass(id_oclo) &
              + q(:ncol,k,id_hcl)          /adv_mass(id_hcl) &
              + q(:ncol,k,id_clono2)       /adv_mass(id_clono2) 
       cloy(:,k) = adv_mass(id_cl) * wrk(:)
    end do

  end function cloy

!--------------------------------------------------------------
! compute the mass mixing retio of BrOy
!--------------------------------------------------------------
  function broy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: broy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_br)    /adv_mass(id_br) &
              + q(:ncol,k,id_bro)   /adv_mass(id_bro) &
              + q(:ncol,k,id_hbr)   /adv_mass(id_hbr) &
              + q(:ncol,k,id_brono2)/adv_mass(id_brono2) &
              + q(:ncol,k,id_brcl)  /adv_mass(id_brcl) &
              + q(:ncol,k,id_hobr)  /adv_mass(id_hobr)
       broy(:,k) = adv_mass(id_br) * wrk(:)
    end do

  end function broy

end module clybry_fam
