module crmclouds_camaerosols
#if (defined MODAL_AERO) 
!---------------------------------------------------------------------------------------------
! Purpose: 
! 
!  Provides the necessary subroutines to use cloud fields from the CRM model to drive the 
!  aerosol-related subroutines in CAM. Several taskes:
!     i) to fill the physics buffers with those diagnosed from the CRM clouds.  
!    ii) to provide the interface for some physics prcoesses, such as droplet activaiton, 
!         and convetive transport. 
!
!  An alternative (and better?) approach is to use the ECPP (explicit-cloud parameterized-pollutant). 
!  This will be done later.
!
!  Revision history: 
!  July, 27, 2009: Minghuai Wang
! 
!-------------------------------------------------------------------------------------------- 
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use cam_abortutils,      only: endrun

   implicit none
   private

   public :: crmclouds_mixnuc_tend 
   public :: crmclouds_diag 
   public :: crmclouds_convect_tend

!======================================================================================================
contains 
!======================================================================================================

!------------------------------------------------------------------------------------------------------
subroutine crmclouds_mixnuc_tend (state, ptend, dtime, cflx, pblht, pbuf,   &
                   wwqui_cen, wwqui_cloudy_cen, wwqui_bnd, wwqui_cloudy_bnd,  species_class )
!-----------------------------------------------------------------------------------------------------
!
! Purpose: to calculate aerosol tendency from dropelt activation and mixing. 
!          Adopted from mmicro_pcond in cldwat2m.F90
!
!------------------------------------------------------------------------------------------------------
  use physics_types,    only: physics_state, physics_ptend, physics_tend, physics_ptend_init
  use physics_buffer,   only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
  use physconst,        only: gravit, rair, karman, spec_class_gas
  use constituents,     only: cnst_get_ind, pcnst
  use time_manager,     only: is_first_step
  use cam_history,      only: outfld
  use ndrop,            only: dropmixnuc
  use modal_aero_data
  use rad_constituents, only: rad_cnst_get_info

  !!! Input 
  type(physics_state), intent(in)    :: state               ! state variables
  type(physics_buffer_desc), pointer :: pbuf(:)             ! physics buffer
  real(r8), intent(in) :: pblht(pcols)                      ! PBL height (meter)
  real(r8), intent(in)  :: dtime                            ! timestep
  real(r8), intent(in) :: cflx(pcols,pcnst)                 ! constituent flux from surface
  real(r8), intent(in) :: wwqui_cen(pcols, pver)            ! vertical velocity variance in quiescent class (m2/s2)
  real(r8), intent(in) :: wwqui_cloudy_cen(pcols, pver)     ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
  real(r8), intent(in) :: wwqui_bnd(pcols, pver+1)          ! vertical velocity variance in quiescent class (m2/s2)
  real(r8), intent(in) :: wwqui_cloudy_bnd(pcols, pver+1)   ! vertical velocity variance in quiescent, and cloudy class (m2/s2)
  integer,  intent(in) :: species_class(:)

  !!! output
  type(physics_ptend), intent(out) :: ptend   ! package tendencies

  !!! Local variables
  integer i,k,m, k1, k2
  integer ifld, itim
  integer ixcldliq, ixcldice, ixnumliq
  integer l,lnum,lnumcw,lmass,lmasscw

  integer :: lchnk    ! chunk identifier
  integer :: ncol     ! number of atmospheric columns
  integer :: nmodes
 
  
  real(r8) :: nc(pcols, pver)             ! droplet number concentration (#/kg)
  real(r8) :: nctend(pcols, pver)         ! change in droplet number concentration
  real(r8) :: omega(pcols, pver)          ! grid-averaaged vertical velocity 
  real(r8) :: qc(pcols, pver)             ! liquid water content (kg/kg)
  real(r8) :: qi(pcols, pver)             ! ice water content (kg/kg) 
  real(r8) :: lcldn(pcols, pver)
  real(r8) :: lcldo(pcols, pver) 

  real(r8) :: wsub(pcols, pver)           ! subgrid vertical velocity
  real(r8) :: ekd_crm(pcols, pverp)       ! diffusivity
  real(r8) :: kkvh_crm(pcols, pverp)      ! eddy diffusivity
  real(r8) :: zs(pcols, pver)             ! inverse of distance between levels (meter)
  real(r8) :: dz(pcols, pver)             ! layer depth (m)
  real(r8) :: cs(pcols, pver)             ! air density
  real(r8) :: lc(pcols, pverp)            ! mixing length (m)
  real(r8) :: zheight(pcols, pverp)       ! height at lay interface (m)
  
  real(r8) :: alc(pcols, pverp)           ! asymptotic length scale (m)
  real(r8) :: tendnd(pcols, pver)         ! tendency of cloud droplet number concentrations (not used in the MMF) 

  real(r8),allocatable :: factnum(:,:,:)  ! activation fraction for aerosol number

  real(r8) :: qcld, qsmall

  logical :: dommf=.true.                 ! value insignificant, if present, means that dropmixnuc is called the mmf part. 

  !!! Variables in the physics buffer:
  real(r8), pointer, dimension(:,:) :: cldn       ! cloud fractin at the current time step
  real(r8), pointer, dimension(:,:) :: cldo       ! cloud fraction at the previous time step
  real(r8), pointer, dimension(:,:) :: acldy_cen  ! liquid cloud fraction at the previous time step from ECPP
  real(r8), pointer, dimension(:,:) ::  kkvh      ! vertical diffusivity
  real(r8), pointer, dimension(:,:) :: tke        ! turbulence kenetic energy 
  real(r8), pointer, dimension(:,:) :: tk_crm     ! m2/s
  logical :: lq(pcnst)

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  lchnk = state%lchnk
  ncol  = state%ncol

  qsmall = 1.e-18_r8

  call rad_cnst_get_info(0, nmodes=nmodes)
  allocate(factnum(pcols,pver,nmodes))

  !----------------------------------------------------------------------------
  ! Initialize ptend
  !----------------------------------------------------------------------------
  lq(:)=.false.
  do m=1,ntot_amode

    lnum = numptr_amode(m)
    if (lnum>0)then
       lq(lnum) = .true.
    end if

    do l=1,nspec_amode(m)
      lmass = lmassptr_amode(l,m)
      lq(lmass) = .true.
    end do ! l

  end do ! m
 
  call physics_ptend_init(ptend,state%psetcols,'crmclouds_mixnuc', lq=lq)

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

   ! In the MMF model, turbulent mixing for tracer species are turned off in tphysac.
   ! So the turbulent for gas species mixing are added here.
   do m=1, pcnst
      if(species_class(m).eq.spec_class_gas) then
#ifdef ECPP
        ptend%lq(m) = .false.
#else
        ptend%lq(m) = .true.
#endif
      end if
   end do

  itim = pbuf_old_tim_idx ()
  ifld = pbuf_get_index ('CLD')
  call pbuf_get_field(pbuf, ifld, cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  ifld = pbuf_get_index ('CLDO')
  call pbuf_get_field(pbuf, ifld, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  ifld = pbuf_get_index ('ACLDY_CEN')
  call pbuf_get_field(pbuf, ifld, acldy_cen)
  ifld = pbuf_get_index('kvh')
  call pbuf_get_field(pbuf, ifld, kkvh)

  ifld=pbuf_get_index('tke')
  call pbuf_get_field(pbuf, ifld, tke)

  ifld = pbuf_get_index('TK_CRM')
  call pbuf_get_field(pbuf, ifld, tk_crm)


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
!
    do k=1, pverp
      if(zheight(i,k).le.pblht(i)) then
        alc(i,k) = 300.
      else
        alc(i,k) = 30+270*exp(1.-zheight(i,k)/pblht(i))
      endif
      lc(i,k) = alc(i,k)*karman*zheight(i,k)/(alc(i,k)+karman*zheight(i,k))
    enddo 
  end do

  call outfld('LENGC', lc, pcols, lchnk)

  kkvh_crm = 0._r8
  do i=1, ncol
    do k=1, pver
!      wsub(i,k) = sqrt(tke(i,k)/3.)  ! tke from CRM is located in the middle of model level
                                  ! should be tke or tke/3. 
                                  ! in cldwat2m.F90, it is tke.
                                  ! wsub seems too large from this approach. 

!      wsub(i,k) = tk_crm(i,k) * zs(i,k)
!      wsub(i,k) = min(wsub(i,k),10._r8)

! from vertical variance in the quiescent class, which excldues 
! the contribution from strong updraft and downdraft. 
!       wsub(i,k) = sqrt(wwqui_cen(i,k))        ! use variance in quiescent area
       wsub(i,k) = sqrt(wwqui_cloudy_cen(i,k))  ! use variance in cloudy quiescent area

! from tke in CAM
!       wsub(i,k) = sqrt(0.5_r8*(tke(i,k)+tke(i,k+1)))

       wsub(i,k) = min(wsub(i,k), 10._r8) 
       wsub(i,k) = max(0.20_r8, wsub(i,k))
    end do   ! end k

    do k=1, pver+1

      k1=min(k, pver)
      k2=max(k-1, 1)
!
! calculate ekd_crm from wsub in the cloudy quiescent class (following a part of ndrop.F90)
!     ekd_crm(i,k) = wsub(i,k) / zs(i,k)
!      ekd_crm(i,k) = min(10.0_r8, max(0.20_r8, sqrt(wwqui_cloudy_bnd(i,k))))*2.0 / (zs(i,k1)+zs(i,k2))  ! use wsub in layer boundary
                                                                                                         ! large ekd at free troposphere. 
      ekd_crm(i,k) = min(10.0_r8, max(0.20_r8, sqrt(wwqui_cloudy_bnd(i,k))))* lc(i,k) 
      kkvh_crm(i,k) = ekd_crm(i,k)

! kkvh_crm is from kvh in CAM
!      kkvh_crm(i,k) = kkvh(i,k)

! set kkvh to kkvh_crm so this will be used in dropmixnuc in the mmf call
       kkvh(i,k) = kkvh_crm(i,k)

    end do   !end k

  end do

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

#ifdef ECPP
!
!      When ECPP is called, activation associated with cloud fraction change is treated in ECPP.
!      so set two cloud fractio be the same here. 
!      But ECPP still did not treat activation associated with turbulent scale motion, and is
!      done in dropmixnuc
!         lcldn(i,k)=cldo(i,k)*qc(i,k)/qcld
!         lcldo(i,k)=cldo(i,k)*qc(i,k)/qcld
         lcldn(i,k)=acldy_cen(i,k)
         lcldo(i,k)=acldy_cen(i,k)
#else
         lcldn(i,k)=cldn(i,k)*qc(i,k)/qcld
         lcldo(i,k)=cldo(i,k)*qc(i,k)/qcld
#endif
      else
         lcldn(i,k)=0._r8
         lcldo(i,k)=0._r8
      endif
    enddo
  enddo

! should we set omega to be zero ??
  omega(:ncol, :) = state%omega(:ncol, :)

  call dropmixnuc(state, ptend, dtime, pbuf, wsub, lcldn, lcldo, tendnd,factnum, species_class,dommf )
  deallocate(factnum)


! this part is moved into tphysbc after aerosol stuffs. 
!
!  cldo(:ncol,:)=cldn(:ncol,:)

end subroutine crmclouds_mixnuc_tend
!======================================================================================================

!------------------------------------------------------------------------------------------------------
subroutine crmclouds_convect_tend(state,  ptend,  ztodt,  pbuf)
!-----------------------------------------------------------------
!
! Purpose: to do convective transport of tracer species using the cloud fields from CRM and using the 
!          subroutine of convtran. 
!
! Minghuai Wang, July, 2009: adopted from zm_conv_tend_2 
!
!------------------------------------------------------------------------------------------------------
   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
   use constituents,  only: pcnst, cnst_get_ind
   use zm_conv,       only: convtran
   use error_messages, only: alloc_err  

! Arguments
! Input variables:
   type(physics_state), intent(in ) :: state          ! Physics state variables
   real(r8), intent(in) :: ztodt

! Output variables:
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer

! Local variables
   integer :: i, lchnk, istat
   integer :: ncol
   integer :: nstep
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
   real(r8), dimension(pcols,pver) :: dpdry
   real(r8), dimension(pcols,pver) :: dp   ! layer thickness in mbs (between upper/lower interface).
   real(r8), dimension(pcols) :: dsubcld   !  wg layer thickness in mbs between lcl and maxi.

! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble


   real(r8), pointer, dimension(:,:) :: mu  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: eu  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: du  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: md  !(pcols,pver,begchunk:endchunk)
   real(r8), pointer, dimension(:,:) :: ed  !(pcols,pver,begchunk:endchunk)

   real(r8), pointer, dimension(:) :: jtr8   !(pcols,begchunk:endchunk)
        ! wg top  level index of deep cumulus convection.
   real(r8), pointer, dimension(:) :: maxgr8 !(pcols,begchunk:endchunk)
        ! wg gathered values of maxi.
   real(r8), pointer, dimension(:) :: ideepr8 !(pcols,begchunk:endchunk)               
        ! w holds position of gathered points vs longitude index

   integer :: jt(pcols)
   integer :: maxg(pcols)
   integer :: ideep(pcols) 
   integer  :: lengath !(begchunk:endchunk)

!==Guangxing Lin
   logical :: lq(pcnst)
!
! Initialize
!
 ! call physics_ptend_init(ptend)
   lq(:) = .true.
   lq(1)        = .false.
   lq(ixcldice) = .false.
   lq(ixcldliq) = .false.

   call physics_ptend_init(ptend,state%psetcols,'convtran2',lq=lq)
!==Guangxing Lin    

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols,pver,pcnst/) )

   ifld = pbuf_get_index('MU_CRM')
   call pbuf_get_field(pbuf, ifld, mu)
   ifld = pbuf_get_index('MD_CRM')
   call pbuf_get_field(pbuf, ifld, md)
   ifld = pbuf_get_index('DU_CRM')
   call pbuf_get_field(pbuf, ifld, du)
   ifld = pbuf_get_index('EU_CRM')
   call pbuf_get_field(pbuf, ifld, eu)
   ifld = pbuf_get_index('ED_CRM')
   call pbuf_get_field(pbuf, ifld, ed)
   ifld = pbuf_get_index('JT_CRM')
   call pbuf_get_field(pbuf, ifld, jtr8)
   ifld = pbuf_get_index('MX_CRM')
   call pbuf_get_field(pbuf, ifld, maxgr8)
   ifld = pbuf_get_index('IDEEP_CRM')
   call pbuf_get_field(pbuf, ifld, ideepr8)


! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk
  ncol = state%ncol

   nstep = get_nstep()

!
!     Convective transport of all trace species except cloud liquid 
!     and cloud ice done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)


   ptend%name  = 'convtran2'
   ptend%lq(:) = .true.
   ptend%lq(1)        = .false.
   ptend%lq(ixcldice) = .false.
   ptend%lq(ixcldliq) = .false.

   
! Is this ok to get the index???
   jt = int(jtr8+0.5_r8)
   maxg = int(maxgr8+0.5_r8)
   ideep = int(ideepr8+0.5_r8)

! calculate lengath from ideep
   lengath = 0
   do i=1, ncol
    if(ideep(i).ge.1) then
      lengath = lengath + 1
    endif
   end do

!
! initialize dpdry for call to convtran 
! it is used for tracers of dry smixing ratio type
!
   dpdry = 0._r8
   do i = 1,lengath
     dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
     dp(i,:) = state%pdel(ideep(i),:)/100._r8
  end do

! dsubdld is not used in convtran, and is set to be zero. 
   dsubcld = 0._r8


!   call t_startf ('crmclouds_convtran')
   call convtran (lchnk,                                        &
                  ptend%lq,state%q, pcnst,  mu(:,:), md(:,:),   &
                  du(:,:), eu(:,:), ed(:,:), dp(:,:), dsubcld(:),  &
                  jt(:),maxg(:),ideep(:), 1, lengath,  &
                  nstep,   fracis,  ptend%q, dpdry  )
!   call t_stopf ('crm_clouds_convtran')

end subroutine crmclouds_convect_tend
!=====================================================================================================

!------------------------------------------------------------------------------------------------------
subroutine crmclouds_diag

end subroutine crmclouds_diag
!======================================================================================================

#endif

end module crmclouds_camaerosols
