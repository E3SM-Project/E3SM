module crmclouds_camaerosols
#if (defined MODAL_AERO) 
!---------------------------------------------------------------------------------------------------
! Purpose: Provide methods for MMF to allow CRM clouds to interact with GCM aerosols
!---------------------------------------------------------------------------------------------------
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use cam_abortutils,     only: endrun
  use ppgrid

  implicit none
  private

  public :: crmclouds_mixnuc_tend 

contains 
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine crmclouds_mixnuc_tend (state, ptend, dtime, cflx, pblht, pbuf, &
                                  wwqui_cen, wwqui_cloudy_cen, wwqui_bnd, &
                                  wwqui_cloudy_bnd,  species_class )
  !-----------------------------------------------------------------------------
  ! Purpose: calculate aerosol tendency from droplet activation and mixing
  !          Adopted from mmicro_pcond in cldwat2m.F90
  !-----------------------------------------------------------------------------
  use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,   only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
  use physconst,        only: gravit, rair, karman, spec_class_gas
  use constituents,     only: cnst_get_ind, pcnst
  use time_manager,     only: is_first_step
  use ndrop,            only: dropmixnuc
  use phys_control,     only: phys_getopts
  use modal_aero_data,  only: ntot_amode, nspec_amode, numptr_amode, lmassptr_amode
  use rad_constituents, only: rad_cnst_get_info
  !-----------------------------------------------------------------------------
  ! Input 
  type(physics_state), intent(in)               :: state             ! state variables
  type(physics_buffer_desc), pointer            :: pbuf(:)           ! physics buffer
  real(r8), intent(in), dimension(pcols)        :: pblht             ! PBL height (m)
  real(r8), intent(in)                          :: dtime             ! timestep
  real(r8), intent(in), dimension(pcols,pcnst)  :: cflx              ! constituent flux from surface
  real(r8), intent(in), dimension(pcols,pver)   :: wwqui_cen         ! vert velocity variance in quiescent class        (m2/s2)
  real(r8), intent(in), dimension(pcols,pver)   :: wwqui_cloudy_cen  ! vert velocity variance in quiescent+cloudy class (m2/s2)
  real(r8), intent(in), dimension(pcols,pver+1) :: wwqui_bnd         ! vert velocity variance in quiescent class        (m2/s2)
  real(r8), intent(in), dimension(pcols,pver+1) :: wwqui_cloudy_bnd  ! vert velocity variance in quiescent+cloudy class (m2/s2)
  integer,  intent(in)                          :: species_class(:)  ! 
  ! Output
  type(physics_ptend), intent(out)              :: ptend             ! output tendencies
  !-----------------------------------------------------------------------------
  ! Local variables
  integer :: i, k, m, k1, k2
  integer :: itim
  integer :: ixcldliq, ixcldice, ixnumliq
  integer :: l, lnum, lnumcw, lmass, lmasscw
  integer :: lchnk, ncol
  integer :: nmodes

  logical :: lq(pcnst)
  logical :: use_ECPP

  real(r8), parameter :: qsmall = 1.e-18_r8

  real(r8), dimension(pcols,pver ) :: nc        ! droplet number concentration (#/kg)
  real(r8), dimension(pcols,pver ) :: nctend    ! change in droplet number concentration
  real(r8), dimension(pcols,pver ) :: omega     ! grid-averaaged vertical velocity 
  real(r8), dimension(pcols,pver ) :: qc        ! liquid water content (kg/kg)
  real(r8), dimension(pcols,pver ) :: qi        ! ice water content (kg/kg) 
  real(r8), dimension(pcols,pver ) :: lcldn     ! new cloud fraction
  real(r8), dimension(pcols,pver ) :: lcldo     ! old cloud fraction
  real(r8), dimension(pcols,pver ) :: wsub      ! subgrid vertical velocity
  real(r8), dimension(pcols,pverp) :: ekd_crm   ! diffusivity
  real(r8), dimension(pcols,pverp) :: kkvh_crm  ! eddy diffusivity
  real(r8), dimension(pcols,pver ) :: zs        ! inverse of distance between levels (meter)
  real(r8), dimension(pcols,pver ) :: dz        ! layer depth (m)
  real(r8), dimension(pcols,pver ) :: cs        ! air density
  real(r8), dimension(pcols,pverp) :: lc        ! mixing length (m)
  real(r8), dimension(pcols,pverp) :: zheight   ! height at lay interface (m)
  real(r8), dimension(pcols,pverp) :: alc       ! asymptotic length scale (m)
  real(r8), dimension(pcols,pver ) :: tendnd    ! tendency of cloud droplet number concentrations (not used in the MMF) 
  
  real(r8), allocatable :: factnum(:,:,:)       ! activation fraction for aerosol number
  real(r8) :: qcld

  ! Variables in physics buffer
  real(r8), pointer, dimension(:,:) :: cldn       ! cloud fractin at the current time step
  real(r8), pointer, dimension(:,:) :: cldo       ! cloud fraction at the previous time step
  real(r8), pointer, dimension(:,:) :: acldy_cen  ! liquid cloud fraction at the previous time step from ECPP
  real(r8), pointer, dimension(:,:) :: kkvh       ! vertical diffusivity
  real(r8), pointer, dimension(:,:) :: tke        ! turbulence kenetic energy 
  real(r8), pointer, dimension(:,:) :: tk_crm     ! m2/s
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call phys_getopts(use_ECPP_out = use_ECPP)

  lchnk = state%lchnk
  ncol  = state%ncol

  call rad_cnst_get_info(0, nmodes=nmodes)
  allocate(factnum(pcols,pver,nmodes))

  !-----------------------------------------------------------------------------
  ! Initialize ptend
  !-----------------------------------------------------------------------------
  lq(:)=.false.
  do m=1,ntot_amode

    lnum = numptr_amode(m)
    if (lnum>0) lq(lnum) = .true.

    do l=1,nspec_amode(m)
      lmass = lmassptr_amode(l,m)
      lq(lmass) = .true.
    end do ! l

  end do ! m
 
  call physics_ptend_init(ptend,state%psetcols,'crmclouds_mixnuc', lq=lq)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  ! In the MMF, turbulent mixing for tracer species are disabled in tphysac,
  ! so the turbulent for gas species mixing are added here
  do m=1, pcnst
    if(species_class(m).eq.spec_class_gas) then
      if (use_ECPP) then
        ptend%lq(m) = .false.
      else
        ptend%lq(m) = .true.
      end if
    end if
  end do ! m

  itim = pbuf_old_tim_idx ()
  call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, pbuf_get_index('CLDO'), cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, pbuf_get_index('ACLDY_CEN'), acldy_cen)
  call pbuf_get_field(pbuf, pbuf_get_index('kvh'), kkvh)
  call pbuf_get_field(pbuf, pbuf_get_index('tke'), tke)
  call pbuf_get_field(pbuf, pbuf_get_index('TK_CRM'), tk_crm)

  if (is_first_step()) then
    kkvh(:,:)= 0.0_r8
    tke(:,:) = 0.0_r8
  endif

  do i=1, ncol
    do k=1, pver-1
    zs(i,k) = 1._r8/(state%zm(i,k)-state%zm(i,k+1))
    end do
    zs(i,pver) = zs(i,pver-1)

    ! calculate height at layer interface (simple calculation)
    zheight(i,pverp) = 0.0
    do k=pver, 1, -1
      zheight(i,k) = zheight(i,k+1) + state%pdel(i,k)/state%pmid(i,k)*(rair*state%t(i,k)/gravit)
    end do

    ! calculate mixing length
    ! from Holtslag and Boville, 1993, J. Climate. 
    do k=1, pverp
      if(zheight(i,k).le.pblht(i)) then
        alc(i,k) = 300.
      else
        alc(i,k) = 30+270*exp(1.-zheight(i,k)/pblht(i))
      endif
      lc(i,k) = alc(i,k)*karman*zheight(i,k)/(alc(i,k)+karman*zheight(i,k))
    enddo 
  end do

  kkvh_crm = 0._r8
  do i=1, ncol
    do k=1, pver
      ! wsub(i,k) = sqrt(tke(i,k)/3.)  ! tke from CRM is located in the middle of model level
      ! should be tke or tke/3. in cldwat2m.F90, it is tke.wsub seems too large from this approach. 

      ! wsub(i,k) = tk_crm(i,k) * zs(i,k)
      ! wsub(i,k) = min(wsub(i,k),10._r8)

      ! from vertical variance in the quiescent class, which excldues 
      ! the contribution from strong updraft and downdraft. 
      ! wsub(i,k) = sqrt(wwqui_cen(i,k))        ! use variance in quiescent area
      wsub(i,k) = sqrt(wwqui_cloudy_cen(i,k))  ! use variance in cloudy quiescent area

      ! from tke in CAM
      ! wsub(i,k) = sqrt(0.5_r8*(tke(i,k)+tke(i,k+1)))

      wsub(i,k) = min(wsub(i,k), 10._r8) 
      wsub(i,k) = max(0.20_r8, wsub(i,k))
    end do   ! k

    do k=1, pver+1

      k1=min(k, pver)
      k2=max(k-1, 1)
      ! calculate ekd_crm from wsub in the cloudy quiescent class (following a part of ndrop.F90)
      ! ekd_crm(i,k) = wsub(i,k) / zs(i,k)
      ! ekd_crm(i,k) = min(10.0_r8, max(0.20_r8, sqrt(wwqui_cloudy_bnd(i,k))))*2.0 / (zs(i,k1)+zs(i,k2))  ! use wsub in layer boundary large ekd at free troposphere. 
      ekd_crm(i,k) = min(10.0_r8, max(0.20_r8, sqrt(wwqui_cloudy_bnd(i,k))))* lc(i,k) 
      kkvh_crm(i,k) = ekd_crm(i,k)

      ! kkvh_crm is from kvh in CAM
      ! kkvh_crm(i,k) = kkvh(i,k)

      ! set kkvh to kkvh_crm so this will be used in dropmixnuc in the mmf call
      kkvh(i,k) = kkvh_crm(i,k)

    end do   ! k
  end do ! i

  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)
  call cnst_get_ind('NUMLIQ', ixnumliq)

  qc(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
  qi(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
  nc(:ncol,:pver) = state%q(:ncol,:pver,ixnumliq)

  do k=1,pver
   do i=1,ncol
      qcld=qc(i,k)+qi(i,k)
      if(qcld.gt.qsmall)then
        if (use_ECPP) then
          ! When ECPP is called, activation associated with cloud fraction 
          ! change is treated in ECPP, so set two cloud fractio be the same here. 
          ! But ECPP still did not treat activation associated with turbulent 
          ! scale motion, and is done in dropmixnuc
          ! lcldn(i,k) = cldo(i,k)*qc(i,k)/qcld
          ! lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
          lcldn(i,k) = acldy_cen(i,k)
          lcldo(i,k) = acldy_cen(i,k)
        else
          lcldn(i,k) = cldn(i,k)*qc(i,k)/qcld
          lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
        end if
      else
        lcldn(i,k) = 0._r8
        lcldo(i,k) = 0._r8
      endif
    end do ! i 
  end do ! k

  omega(:ncol,:) = state%omega(:ncol,:) ! should we set omega to be zero ??

  call dropmixnuc(state, ptend, dtime, pbuf, wsub, lcldn, lcldo, tendnd,factnum, species_class, .true. )
  deallocate(factnum)

  ! this part is moved into tphysbc after aerosol stuff
  ! cldo(:ncol,:) = cldn(:ncol,:)

end subroutine crmclouds_mixnuc_tend
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
#endif /* MODAL_AERO */

end module crmclouds_camaerosols
