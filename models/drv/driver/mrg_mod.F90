module mrg_mod

  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use mct_mod
  use seq_cdata_mod
  use seq_comm_mct
  use seq_infodata_mod
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! TODO - write summary of naming convention here as well  
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: mrg_x2a_run_mct
  public :: mrg_x2i_run_mct
  public :: mrg_x2l_run_mct
  public :: mrg_x2r_run_mct
  public :: mrg_x2o_run_mct
  public :: mrg_x2g_run_mct
  public :: mrg_x2s_run_mct
  public :: mrg_x2w_run_mct

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: getfld 

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

  subroutine mrg_x2a_run_mct( cdata_a, l2x_a, o2x_a, xao_a, i2x_a, fractions_a, x2a_a )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_a
    type(mct_aVect), intent(in)     :: l2x_a
    type(mct_aVect), intent(in)     :: o2x_a
    type(mct_aVect), intent(in)     :: xao_a
    type(mct_aVect), intent(in)     :: i2x_a
    type(mct_aVect), intent(in)     :: fractions_a
    type(mct_aVect), intent(inout)  :: x2a_a
    !----------------------------------------------------------------------- 
    !
    ! Local workspace
    !
    real(r8) :: fracl, fraci, fraco
    integer  :: n,ka,ki,kl,ko,kx,kof,kif,klf
    integer  :: lsize       
    integer  :: index_x2a_Sf_lfrac
    integer  :: index_x2a_Sf_ifrac
    integer  :: index_x2a_Sf_ofrac
    character(CL) :: field_atm   ! string converted to char
    character(CL) :: field_lnd   ! string converted to char
    character(CL) :: field_ice   ! string converted to char
    character(CL) :: field_xao   ! string converted to char
    character(CL) :: field_ocn   ! string converted to char
    character(CL) :: itemc_atm   ! string converted to char
    character(CL) :: itemc_lnd   ! string converted to char
    character(CL) :: itemc_ice   ! string converted to char
    character(CL) :: itemc_xao   ! string converted to char
    character(CL) :: itemc_ocn   ! string converted to char
    logical :: iamroot  
    logical :: first_time = .true.
    logical, pointer, save :: lmerge(:),imerge(:),xmerge(:),omerge(:)
    integer, pointer, save :: lindx(:), iindx(:), oindx(:),xindx(:)
    integer, save          :: naflds, klflds,niflds,noflds,nxflds
    !-----------------------------------------------------------------------
    !
    call seq_comm_setptrs(CPLID, iamroot=iamroot)

    if (first_time) then
          
       naflds = mct_aVect_nRattr(x2a_a)
       klflds = mct_aVect_nRattr(l2x_a)
       niflds = mct_aVect_nRattr(i2x_a)
       noflds = mct_aVect_nRattr(o2x_a)
       nxflds = mct_aVect_nRattr(xao_a)

       allocate(lindx(naflds), lmerge(naflds))
       allocate(iindx(naflds), imerge(naflds))
       allocate(xindx(naflds), xmerge(naflds))
       allocate(oindx(naflds), omerge(naflds))

       lindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       oindx(:) = 0
       lmerge(:)  = .true.
       imerge(:)  = .true.
       xmerge(:)  = .true.
       omerge(:)  = .true.

       ! Field naming rules
       ! Only atm states that are Sx_... will be merged
       ! Only fluxes that are F??x_... will be merged 
       ! All fluxes will be multiplied by corresponding component fraction

       do ka = 1,naflds
          call getfld(ka, x2a_a, field_atm, itemc_atm)
          if (field_atm(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through 
          end if
          if (field_atm(1:1) == 'S' .and. field_atm(2:2) /= 'x') then
             cycle ! any state fields that are not Sx_ will just be copied
          end if

          do kl = 1,klflds
             call getfld(kl, l2x_a, field_lnd, itemc_lnd)
             if (trim(itemc_atm) == trim(itemc_lnd)) then
                if ((trim(field_atm) == trim(field_lnd))) then
                   if (field_lnd(1:1) == 'F') lmerge(ka) = .false.
                end if
                lindx(ka) = kl
                exit 
             end if
          end do
          do ki = 1,niflds
             call getfld(ki, i2x_a, field_ice, itemc_ice)
             if (field_ice(1:1) == 'F' .and. field_ice(2:4) == 'ioi') then
                cycle ! ignore all fluxes that are ice/ocn fluxes
             end if
             if (trim(itemc_atm) == trim(itemc_ice)) then
                if ((trim(field_atm) == trim(field_ice))) then
                   if (field_ice(1:1) == 'F') imerge(ka) = .false.
                end if
                iindx(ka) = ki
                exit 
             end if
          end do
          do kx = 1,nxflds
             call getfld(kx, xao_a, field_xao, itemc_xao)
             if (trim(itemc_atm) == trim(itemc_xao)) then
                if ((trim(field_atm) == trim(field_xao))) then
                   if (field_xao(1:1) == 'F') xmerge(ka) = .false.
                end if
                xindx(ka) = kx
                exit 
             end if
          end do
          do ko = 1,noflds
             call getfld(ko, o2x_a, field_ocn, itemc_ocn)
             if (trim(itemc_atm) == trim(itemc_ocn)) then
                if ((trim(field_atm) == trim(field_ocn))) then
                   if (field_ocn(1:1) == 'F') omerge(ka) = .false.
                end if
                oindx(ka) = ko
                exit 
             end if
          end do
          if (lindx(ka) == 0) itemc_lnd = 'unset'
          if (iindx(ka) == 0) itemc_ice = 'unset'
          if (xindx(ka) == 0) itemc_xao = 'unset'
          if (oindx(ka) == 0) itemc_ocn = 'unset'

          if (iamroot) then
             write(logunit,10)trim(itemc_atm),trim(itemc_lnd),&
                  trim(itemc_ice),trim(itemc_xao),trim(itemc_ocn)
10           format(' ',' atm field: ',a15,', lnd merge: ',a15, &
             ', ice merge: ',a15,', xao merge: ',a15,', ocn merge: ',a15)
             write(logunit, *)'field_atm,lmerge, imerge, xmerge, omerge= ',&
                  trim(field_atm),lmerge(ka),imerge(ka),xmerge(ka),omerge(ka)
         end if
       end do
       first_time = .false.
    end if

    ! Zero attribute vector

    call mct_avect_zero(x2a_a)

    ! Update surface fractions

    kif=mct_aVect_indexRA(fractions_a,"ifrac")
    klf=mct_aVect_indexRA(fractions_a,"lfrac")
    kof=mct_aVect_indexRA(fractions_a,"ofrac")
    lsize = mct_avect_lsize(x2a_a)

    index_x2a_Sf_lfrac = mct_aVect_indexRA(x2a_a,'Sf_lfrac')
    index_x2a_Sf_ifrac = mct_aVect_indexRA(x2a_a,'Sf_ifrac')
    index_x2a_Sf_ofrac = mct_aVect_indexRA(x2a_a,'Sf_ofrac')
    do n = 1,lsize
       x2a_a%rAttr(index_x2a_Sf_lfrac,n) = fractions_a%Rattr(klf,n)
       x2a_a%rAttr(index_x2a_Sf_ifrac,n) = fractions_a%Rattr(kif,n)
       x2a_a%rAttr(index_x2a_Sf_ofrac,n) = fractions_a%Rattr(kof,n)
    end do

    ! Copy attributes that do not need to be merged
    ! These are assumed to have the same name in 
    ! (o2x_a and x2a_a) and in (l2x_a and x2a_a), etc.

    call mct_aVect_copy(aVin=l2x_a, aVout=x2a_a, vector=mct_usevector)
    call mct_aVect_copy(aVin=o2x_a, aVout=x2a_a, vector=mct_usevector)
    call mct_aVect_copy(aVin=i2x_a, aVout=x2a_a, vector=mct_usevector) 
    call mct_aVect_copy(aVin=xao_a, aVout=x2a_a, vector=mct_usevector)

    ! If flux to atm is coming only from the ocean (based on field being in o2x_a) - 
    ! -- then scale by both ocean and ice fraction
    ! If flux to atm is coming only from the land or ice or coupler
    ! -- then do scale by fraction above
    
    do ka = 1,naflds
       do n = 1,lsize
          fracl = fractions_a%Rattr(klf,n)
          fraci = fractions_a%Rattr(kif,n)
          fraco = fractions_a%Rattr(kof,n)
          if (lindx(ka) > 0 .and. fracl > 0._r8) then
             if (lmerge(ka)) then 
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + l2x_a%rAttr(lindx(ka),n) * fracl
             else
                x2a_a%rAttr(ka,n) = l2x_a%rAttr(lindx(ka),n) * fracl
             end if
          end if
          if (iindx(ka) > 0 .and. fraci > 0._r8) then
             if (imerge(ka)) then 
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + i2x_a%rAttr(iindx(ka),n) * fraci
             else
                x2a_a%rAttr(ka,n) = i2x_a%rAttr(iindx(ka),n) * fraci
             end if
          end if
          if (xindx(ka) > 0 .and. fraco > 0._r8) then
             if (xmerge(ka)) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + xao_a%rAttr(xindx(ka),n) * fraco
             else
                x2a_a%rAttr(ka,n) = xao_a%rAttr(xindx(ka),n) * fraco
             end if
          end if
          if (oindx(ka) > 0) then
             if (omerge(ka) .and. fraco > 0._r8) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco
             end if
             if (.not. omerge(ka)) then
                !--- NOTE: This IS using the ocean fields and ice fraction !! ---
                x2a_a%rAttr(ka,n) = o2x_a%rAttr(oindx(ka),n) * fraci
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco 
             end if
          end if
       end do
    end do

  end subroutine mrg_x2a_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2i_run_mct( cdata_i, a2x_i, o2x_i, x2i_i )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_i
    type(mct_aVect),intent(in) :: a2x_i
    type(mct_aVect),intent(in) :: o2x_i
    type(mct_aVect),intent(inout):: x2i_i
    !
    ! Local variables
    !
    integer :: i
    real(r8):: flux_epbalfact
    character(len=cl) :: flux_epbal
    type(seq_infodata_type),pointer :: infodata
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_x2i_Faxa_rain
    integer, save :: index_x2i_Faxa_snow
    logical, save :: first_time = .true.
    !----------------------------------------------------------------------- 

    if (first_time) then
       index_a2x_Faxa_snowc = mct_aVect_indexRA(a2x_i,'Faxa_snowc')
       index_a2x_Faxa_snowl = mct_aVect_indexRA(a2x_i,'Faxa_snowl')
       index_a2x_Faxa_rainc = mct_aVect_indexRA(a2x_i,'Faxa_rainc')
       index_a2x_Faxa_rainl = mct_aVect_indexRA(a2x_i,'Faxa_rainl')
       index_x2i_Faxa_rain  = mct_aVect_indexRA(x2i_i,'Faxa_rain' )
       index_x2i_Faxa_snow  = mct_aVect_indexRA(x2i_i,'Faxa_snow' )
       first_time = .false.
    end if

    ! Apply correction to precipitation of requested driver namelist
    call seq_cdata_setptrs(cdata_i,infodata=infodata)
    call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)

    call mct_aVect_copy(aVin=o2x_i, aVout=x2i_i, vector=mct_usevector)
    call mct_aVect_copy(aVin=a2x_i, aVout=x2i_i, vector=mct_usevector)

    ! Merge total snow and precip for ice input
    ! Scale total precip and runoff by flux_epbalfact 

    do i = 1,mct_aVect_lsize(x2i_i)
       x2i_i%rAttr(index_x2i_Faxa_rain,i) = a2x_i%rAttr(index_a2x_Faxa_rainc,i) + &
	                                    a2x_i%rAttr(index_a2x_Faxa_rainl,i)
       x2i_i%rAttr(index_x2i_Faxa_snow,i) = a2x_i%rAttr(index_a2x_Faxa_snowc,i) + &
	                                    a2x_i%rAttr(index_a2x_Faxa_snowl,i) 

       x2i_i%rAttr(index_x2i_Faxa_rain,i) = x2i_i%rAttr(index_x2i_Faxa_rain,i) * flux_epbalfact
       x2i_i%rAttr(index_x2i_Faxa_snow,i) = x2i_i%rAttr(index_x2i_Faxa_snow,i) * flux_epbalfact
    end do

  end subroutine mrg_x2i_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2r_run_mct( cdata_r, l2x_r, fractions_r, x2r_r)

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_r
    type(mct_aVect),intent(in) :: l2x_r
    type(mct_aVect),intent(in) :: fractions_r
    type(mct_aVect),intent(inout):: x2r_r
    !
    ! Local variables
    !
    integer :: i
    type(seq_infodata_type),pointer :: infodata
    integer, save :: index_l2x_Flrl_rofliq
    integer, save :: index_l2x_Flrl_rofice
    integer, save :: index_x2r_Flrl_rofliq
    integer, save :: index_x2r_Flrl_rofice
    integer, save :: index_lfrac
    logical, save :: first_time = .true.
    real(r8) :: lfrac
    !----------------------------------------------------------------------- 

    if (first_time) then
       index_l2x_Flrl_rofliq = mct_aVect_indexRA(l2x_r,'Flrl_rofliq' )
       index_l2x_Flrl_rofice = mct_aVect_indexRA(l2x_r,'Flrl_rofice' )
       index_x2r_Flrl_rofliq = mct_aVect_indexRA(x2r_r,'Flrl_rofliq' )
       index_x2r_Flrl_rofice = mct_aVect_indexRA(x2r_r,'Flrl_rofice' )
       index_lfrac = mct_aVect_indexRA(fractions_r,"lfrac")
       first_time = .false.
    end if

    ! Merge land rof and ice forcing for rof input

    do i = 1,mct_aVect_lsize(x2r_r)
       lfrac = fractions_r%rAttr(index_lfrac,i)
       x2r_r%rAttr(index_x2r_Flrl_rofliq,i) = l2x_r%rAttr(index_l2x_Flrl_rofliq,i) * lfrac
       x2r_r%rAttr(index_x2r_Flrl_rofice,i) = l2x_r%rAttr(index_l2x_Flrl_rofice,i) * lfrac
    end do

  end subroutine mrg_x2r_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2l_run_mct( cdata_l, a2x_l, r2l_l, x2l_l )

    !----------------------------------------------------------------------- 
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_l
    type(mct_aVect), intent(in)     :: a2x_l  ! input
    type(mct_aVect), intent(in)     :: r2l_l  ! input
    type(mct_aVect), intent(inout)  :: x2l_l  ! output
    !----------------------------------------------------------------------- 

    ! Create input land state directly from atm and runoff outputs 
    call mct_aVect_copy(aVin=a2x_l, aVout=x2l_l, vector=mct_usevector)
    call mct_aVect_copy(aVin=r2l_l, aVout=x2l_l, vector=mct_usevector)

  end subroutine mrg_x2l_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2o_run_mct( cdata_o, a2x_o, i2x_o, w2x_o, xao_o, fractions_o, x2o_o )

    !----------------------------------------------------------------------- 
    ! Arguments
    type(seq_cdata), intent(in)    :: cdata_o
    type(mct_aVect), intent(in)    :: a2x_o
    type(mct_aVect), intent(in)    :: i2x_o
    type(mct_aVect), intent(in)    :: w2x_o
    type(mct_aVect), intent(in)    :: xao_o
    type(mct_aVect), intent(in)    :: fractions_o
    type(mct_aVect), intent(inout) :: x2o_o
    !
    ! Local variables
    !
    integer  :: n,ka,ki,ko,kir,kor
    integer  :: lsize
    real(r8) :: ifrac,ifracr
    real(r8) :: afrac,afracr
    real(r8) :: flux_epbalfact
    real(r8) :: frac_sum
    real(r8) :: avsdr, anidr, avsdf, anidf   ! albedos
    real(r8) :: fswabsv, fswabsi             ! sw
    integer  :: noflds,naflds,niflds,nxflds
    integer  :: kof,kaf,kif,kxf
    character(len=cl) :: flux_epbal
    character(CL) :: field_ocn   ! string converted to char
    character(CL) :: field_atm   ! string converted to char
    character(CL) :: field_ice   ! string converted to char
    character(CL) :: field_xao   ! string converted to char
    character(CL) :: itemc_ocn   ! string converted to char
    character(CL) :: itemc_atm   ! string converted to char
    character(CL) :: itemc_ice   ! string converted to char
    character(CL) :: itemc_xao   ! string converted to char
    logical :: iamroot  
    type(seq_infodata_type),pointer :: infodata
    integer, save :: index_a2x_Faxa_swvdr
    integer, save :: index_a2x_Faxa_swvdf
    integer, save :: index_a2x_Faxa_swndr
    integer, save :: index_a2x_Faxa_swndf
    integer, save :: index_i2x_Fioi_swpen
    integer, save :: index_xao_So_avsdr
    integer, save :: index_xao_So_anidr
    integer, save :: index_xao_So_avsdf
    integer, save :: index_xao_So_anidf
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_x2o_Foxx_swnet
    integer, save :: index_x2o_Faxa_snow 
    integer, save :: index_x2o_Faxa_rain 
    integer, save :: index_x2o_Faxa_prec  
    logical, save, pointer :: amerge(:),imerge(:),xmerge(:)
    integer, save, pointer :: aindx(:), iindx(:), oindx(:), xindx(:)
    logical, save :: first_time = .true.
    character(*),parameter :: subName = '(mrg_x2o_run_mct) '
    !----------------------------------------------------------------------- 

    call seq_comm_setptrs(CPLID, iamroot=iamroot)

    noflds = mct_aVect_nRattr(x2o_o)
    naflds = mct_aVect_nRattr(a2x_o)
    niflds = mct_aVect_nRattr(i2x_o)
    nxflds = mct_aVect_nRattr(xao_o)

    if (first_time) then
       index_a2x_Faxa_swvdr     = mct_aVect_indexRA(a2x_o,'Faxa_swvdr')
       index_a2x_Faxa_swvdf     = mct_aVect_indexRA(a2x_o,'Faxa_swvdf')
       index_a2x_Faxa_swndr     = mct_aVect_indexRA(a2x_o,'Faxa_swndr')
       index_a2x_Faxa_swndf     = mct_aVect_indexRA(a2x_o,'Faxa_swndf')
       index_i2x_Fioi_swpen     = mct_aVect_indexRA(i2x_o,'Fioi_swpen') 
       index_xao_So_avsdr       = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_anidr       = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_avsdf       = mct_aVect_indexRA(xao_o,'So_avsdf')
       index_xao_So_anidf       = mct_aVect_indexRA(xao_o,'So_anidf')
       index_x2o_Foxx_swnet     = mct_aVect_indexRA(x2o_o,'Foxx_swnet')

       index_a2x_Faxa_snowc     = mct_aVect_indexRA(a2x_o,'Faxa_snowc')
       index_a2x_Faxa_snowl     = mct_aVect_indexRA(a2x_o,'Faxa_snowl')
       index_a2x_Faxa_rainc     = mct_aVect_indexRA(a2x_o,'Faxa_rainc')
       index_a2x_Faxa_rainl     = mct_aVect_indexRA(a2x_o,'Faxa_rainl')
       index_x2o_Faxa_snow      = mct_aVect_indexRA(x2o_o,'Faxa_snow')
       index_x2o_Faxa_rain      = mct_aVect_indexRA(x2o_o,'Faxa_rain')
       index_x2o_Faxa_prec      = mct_aVect_indexRA(x2o_o,'Faxa_prec') 

       ! Compute all other quantities based on standardized naming convention (see below)
       ! Only ocn field states that have the name-prefix Sx_ will be merged
       ! Only field names have the same name-suffix (after the "_") will be merged
       !    (e.g. Si_fldname, Sa_fldname => merged to => Sx_fldname)
       ! All fluxes will be scaled by the corresponding afrac or ifrac 
       !   EXCEPT for 
       !    -- Faxa_snnet, Faxa_snow, Faxa_rain, Faxa_prec (derived)
       !    -- Forr_* (treated in ccsm_comp_mod)
       ! All i2x_o fluxes that have the name-suffix "Faii" (atm/ice fluxes) will be ignored 
       ! - only ice fluxes that are Fioi_... will be used in the ocean merges

       allocate(aindx(noflds), amerge(noflds))
       allocate(iindx(noflds), imerge(noflds))
       allocate(xindx(noflds), xmerge(noflds))
       aindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       amerge(:) = .true.
       imerge(:) = .true.
       xmerge(:) = .true.

       do kof = 1,noflds
          call getfld(kof, x2o_o, field_ocn, itemc_ocn)
          if (field_ocn(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through 
          end if
          if (field_ocn(1:1) == 'S' .and. field_ocn(2:2) /= 'x') then
             cycle ! ignore all ocn states that do not have a Sx_ prefix 
          end if
          if (trim(field_ocn) == 'Foxx_swnet'.or. &
              trim(field_ocn) == 'Faxa_snow' .or. &
              trim(field_ocn) == 'Faxa_rain' .or. &
              trim(field_ocn) == 'Faxa_prec') then
             cycle ! ignore swnet, snow, rain, prec - treated explicitly above
          end if
          if (trim(field_ocn(1:5)) == 'Forr_') then
             cycle ! ignore runoff fields from land - treated in coupler
          end if

          do kaf = 1,naflds
             call getfld(kaf, a2x_o, field_atm, itemc_atm)
             if (trim(itemc_ocn) == trim(itemc_atm)) then
                if ((trim(field_ocn) == trim(field_atm))) then
                   if (field_atm(1:1) == 'F') amerge(kof) = .false.
                end if
                aindx(kof) = kaf
                exit
             end if
          end do
          do kif = 1,niflds
             call getfld(kif, i2x_o, field_ice, itemc_ice)
             if (field_ice(1:1) == 'F' .and. field_ice(2:4) == 'aii') then
                cycle ! ignore all i2x_o fluxes that are ice/atm fluxes
             end if
             if (trim(itemc_ocn) == trim(itemc_ice)) then
                if ((trim(field_ocn) == trim(field_ice))) then
                   if (field_ice(1:1) == 'F') imerge(kof) = .false.
                end if
                iindx(kof) = kif
                exit
             end if
          end do
          do kxf = 1,nxflds
             call getfld(kxf, xao_o, field_xao, itemc_xao)
             if (trim(itemc_ocn) == trim(itemc_xao)) then
                if ((trim(field_ocn) == trim(field_xao))) then
                   if (field_xao(1:1) == 'F') xmerge(kof) = .false.
                end if
                xindx(kof) = kxf
                exit
             end if
          end do
          if (aindx(kof) == 0) itemc_atm = 'unset'
          if (iindx(kof) == 0) itemc_ice = 'unset'
          if (xindx(kof) == 0) itemc_xao = 'unset'

          if (iamroot) then
             write(logunit,10)trim(itemc_ocn),&
                  trim(itemc_xao),trim(itemc_ice),trim(itemc_atm)
10           format(' ',' ocn field: ',a15,', xao merge: ',a15, &
                  ', ice merge: ',a15,', atm merge: ',a15)
             write(logunit, *)'field_ocn,kof,imerge,amerge,xmerge= ',&
                  trim(field_ocn),kof,imerge(kof),xmerge(kof),amerge(kof) 
         end if
       end do

       first_time = .false.
    end if
    
    call seq_cdata_setptrs(cdata_o, infodata=infodata)
    call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)

    call mct_aVect_zero(x2o_o)

    call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=w2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=mct_usevector)

    ! Compute input ocn state (note that this only applies to non-land portion of gridcell)

    ki  = mct_aVect_indexRa(fractions_o,"ifrac",perrWith=subName)
    ko  = mct_aVect_indexRa(fractions_o,"ofrac",perrWith=subName)
    kir = mct_aVect_indexRa(fractions_o,"ifrad",perrWith=subName)
    kor = mct_aVect_indexRa(fractions_o,"ofrad",perrWith=subName)
    lsize = mct_aVect_lsize(x2o_o)
    do n = 1,lsize

       ifrac = fractions_o%rAttr(ki,n)
       afrac = fractions_o%rAttr(ko,n)
       frac_sum = ifrac + afrac
       if ((frac_sum) /= 0._r8) then
          ifrac = ifrac / (frac_sum)
          afrac = afrac / (frac_sum)
       endif

       ifracr = fractions_o%rAttr(kir,n)
       afracr = fractions_o%rAttr(kor,n)
       frac_sum = ifracr + afracr
       if ((frac_sum) /= 0._r8) then
          ifracr = ifracr / (frac_sum)
          afracr = afracr / (frac_sum)
       endif

       ! Derived: compute net short-wave
       avsdr = xao_o%rAttr(index_xao_So_avsdr,n)  
       anidr = xao_o%rAttr(index_xao_So_anidr,n)  
       avsdf = xao_o%rAttr(index_xao_So_avsdf,n)  
       anidf = xao_o%rAttr(index_xao_So_anidf,n)  
       fswabsv  =  a2x_o%rAttr(index_a2x_Faxa_swvdr,n) * (1.0_R8 - avsdr) &
                 + a2x_o%rAttr(index_a2x_Faxa_swvdf,n) * (1.0_R8 - avsdf)
       fswabsi  =  a2x_o%rAttr(index_a2x_Faxa_swndr,n) * (1.0_R8 - anidr) &
                 + a2x_o%rAttr(index_a2x_Faxa_swndf,n) * (1.0_R8 - anidf)
       x2o_o%rAttr(index_x2o_Foxx_swnet,n) = (fswabsv + fswabsi)                 * afracr + &
                                             i2x_o%rAttr(index_i2x_Fioi_swpen,n) * ifrac

       ! Derived: compute total precipitation - scale total precip 
       ! Note that runoff is scaled by flux_epbalfact in ccsm_comp_mod
       x2o_o%rAttr(index_x2o_Faxa_snow ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc,n) * afrac + &
                                             a2x_o%rAttr(index_a2x_Faxa_snowl,n) * afrac 
       x2o_o%rAttr(index_x2o_Faxa_rain ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc,n) * afrac + &
                                             a2x_o%rAttr(index_a2x_Faxa_rainl,n) * afrac

       x2o_o%rAttr(index_x2o_Faxa_snow ,n) = x2o_o%rAttr(index_x2o_Faxa_snow ,n) * flux_epbalfact
       x2o_o%rAttr(index_x2o_Faxa_rain ,n) = x2o_o%rAttr(index_x2o_Faxa_rain ,n) * flux_epbalfact

       x2o_o%rAttr(index_x2o_Faxa_prec ,n) = x2o_o%rAttr(index_x2o_Faxa_rain ,n) + &
                                             x2o_o%rAttr(index_x2o_Faxa_snow ,n)
    end do

    do kof = 1,noflds
       do n = 1,lsize
          ifrac = fractions_o%rAttr(ki,n)
          afrac = fractions_o%rAttr(ko,n)
          frac_sum = ifrac + afrac
          if ((frac_sum) /= 0._r8) then
             ifrac = ifrac / (frac_sum)
             afrac = afrac / (frac_sum)
          endif
          if (iindx(kof) > 0) then
             if (imerge(kof)) then 
                x2o_o%rAttr(kof,n) = x2o_o%rAttr(kof,n) + i2x_o%rAttr(iindx(kof),n) * ifrac
             else
                x2o_o%rAttr(kof,n) = i2x_o%rAttr(iindx(kof),n) * ifrac
             end if
          end if
          if (aindx(kof) > 0) then
             if (amerge(kof)) then
                x2o_o%rAttr(kof,n) = x2o_o%rAttr(kof,n) + a2x_o%rAttr(aindx(kof),n) * afrac
             else
                x2o_o%rAttr(kof,n) = a2x_o%rAttr(aindx(kof),n) * afrac
             end if
          end if
          if (xindx(kof) > 0) then
             if (xmerge(kof)) then
                x2o_o%rAttr(kof,n) = x2o_o%rAttr(kof,n) + xao_o%rAttr(xindx(kof),n) * afrac
             else
                x2o_o%rAttr(kof,n) = xao_o%rAttr(xindx(kof),n) * afrac
             end if
          end if
       end do
    end do
       
  end subroutine mrg_x2o_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2g_run_mct( cdata_g, s2x_g, x2g_g )

    !----------------------------------------------------------------------- 
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_g
    type(mct_aVect), intent(inout)  :: s2x_g  ! input
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !----------------------------------------------------------------------- 

    ! Create input glc state directly from land snow output state
    call mct_aVect_copy(aVin=s2x_g, aVout=x2g_g, vector=mct_usevector)

  end subroutine mrg_x2g_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2s_run_mct( cdata_s, g2x_s, x2s_s )

    !----------------------------------------------------------------------- 
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_s
    type(mct_aVect), intent(inout)  :: g2x_s  ! input
    type(mct_aVect), intent(inout)  :: x2s_s  ! output
    !----------------------------------------------------------------------- 

    ! Create input land state directly from glc output state
    call mct_aVect_copy(aVin=g2x_s, aVout=x2s_s, vector=mct_usevector)

  end subroutine mrg_x2s_run_mct

!--------------------------------------------------------------------------

  subroutine mrg_x2w_run_mct( cdata_w, a2x_w, o2x_w, i2x_w, frac_w, x2w_w)

    !-----------------------------------------------------------------------
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_w
    type(mct_aVect), intent(inout)  :: a2x_w  ! input
    type(mct_aVect), intent(inout)  :: o2x_w  ! input
    type(mct_aVect), intent(inout)  :: i2x_w  ! input
    type(mct_aVect), intent(inout)  :: frac_w ! input
    type(mct_aVect), intent(inout)  :: x2w_w  ! output
    !-----------------------------------------------------------------------

    ! Create input wave state directly from atm, ocn, ice output state

    call mct_aVect_copy(aVin=a2x_w, aVout=x2w_w, vector=mct_usevector)
    call mct_aVect_copy(aVin=o2x_w, aVout=x2w_w, vector=mct_usevector)
    call mct_aVect_copy(aVin=i2x_w, aVout=x2w_w, vector=mct_usevector)

  end subroutine mrg_x2w_run_mct

!--------------------------------------------------------------------------

  subroutine getfld(n, av, field, suffix)
    integer         , intent(in)    :: n
    type(mct_aVect) , intent(in)    :: av 
    character(len=*), intent(out)   :: field
    character(len=*), intent(out)   :: suffix

    type(mct_string) :: mstring     ! mct char type

    call mct_aVect_getRList(mstring,n,av)
    field  = mct_string_toChar(mstring)
    suffix = trim(field(scan(field,'_'):))
    call mct_string_clean(mstring)

    if (field(1:1) /= 'S' .and. field(1:1) /= 'F' .and. field(1:2) /= 'PF') then
       write(6,*)'field attribute',trim(field),' must start with S or F or PF' 
       call shr_sys_abort()
    end if
  end subroutine getfld

end module mrg_mod

