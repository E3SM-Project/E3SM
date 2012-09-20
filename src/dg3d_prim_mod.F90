#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
  
  !=======================================================================================================!
  !	DG3D Dynamical core development:								!
  ! 		Creating By:  	R. D. Nair (SCD/01/2006)						!
  ! 		Modifying By:	H.-W. Choi (SCD/01/2006)						!
  ! 		New addtions:	RDN & JG   (CISL/11/2011)						!
  ! 	 (u,v) momentum form:	RDN & JG   (CISL/02/2012)						!
  !=======================================================================================================!
  module dg3d_prim_mod
    !=========================================!
    use kinds, only : real_kind
    !=========================================!
    use physical_constants, only : rearth , g
    !=========================================!
    use dimensions_mod, only : np, nlev
    !=========================================!
    use element_mod, only : element_t
    !=========================================!
    use edge_mod, only : EdgeBuffer_t, edgevpack,edgerotate,edgevunpack,edgeDGVpack,edgeDGVunpack
    !=========================================!
    use filter_mod, only : filter_t, filter_P
    !=========================================!
    use hybrid_mod, only : hybrid_t    
    !=========================================! 
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    !=========================================!
    use quadrature_mod, only : quadrature_t, gauss, gausslobatto
    !=========================================!
    use derivative_mod, only : derivative_t, gradient_wk, vorticity
    !=========================================!
    use time_mod, only : timelevel_t, smooth
    !=========================================!
    use bndry_mod, only : bndry_exchangev
    !=========================================!   
    use control_mod, only : nu, remap_type, filter_freq, filter_counter, topology, test_case, remapfreq, statefreq
    !=========================================!
    use dg3d_dynamics_mod, only:  dg3d_rhs_terms, pres_grad_term, horizontal_diff, implicit_diff, &
        dg3d_diff_grads, dg3d_diff_grads_uv, dg3d_diff_flux ,cov_vorticity , dg3d_uvform_rhs,  &
        gradient_p3d, diffusion_theta , diffusion_hypr 
    !=========================================!
    use dg3d_vertical_mod, only : lagrangian_surfvars, etalevel_temp, pt2temp, temp2pt
    !=========================================!
    use dg3d_remap_mod, only :  linear_remap, parabolic_remap, energy_remap, monotonic_remap
    !=========================================!
    use dg3d_tests_mod, only : jw_baroclinic, heldsuarez_initial,heldsuarez_uv_forcing, &
         heldsuarez_th_correction,heldsuarez_pt_forcing 
    !=========================================!
    use dg3d_errors_mod, only : slice850, slice_at_prlevel
    !=========================================!
    use dg_core_mod, only: co2contra,contra2co,sphere2contra,contra2sphere,height2phi,phi2height,  &
         psi2height,dp2pt,pt2dp, mono_filter
    !=======================================================================================================!
    use parallel_mod, only: parallel_t, syncmp
    !=======================================================================================================! 
    !=======================================================================================================!
    real (kind=real_kind), dimension(:,:), pointer :: fcor,rmv,mv,metdet,rmetdetp
    real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
    real (kind=real_kind), dimension(np,np):: mmx,mmxi    
    real (kind=real_kind), dimension(np,np,2,nlev) :: tmp  ! temporary variable
    !=======================================================================================================!
    public ::  dg3d_advance 
    public ::  dg3d_ssprk2 
    public ::  dg3d_uv_ssprk3, dg3d_uv_ssprk2
    !public ::  dg3d_rk3 
    !=======================================================================================================!
  contains
    !=======================================================================================================!
    !  	DG3D_ADVANCE											!
    !=======================================================================================================!
    subroutine dg3d_advance(elem, edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
      !=======================================================================================================!
      implicit none
      type (element_t)     , intent(inout) :: elem(:)
      type (EdgeBuffer_t)  , intent(in) :: edge3
      type (derivative_t)  , intent(in) :: deriv
      type (filter_t)                   :: flt
      type (hybrid_t)      , intent(in) :: hybrid
      real (kind=real_kind), intent(in) :: dt
      real (kind=real_kind), intent(in) :: pmean
      type (TimeLevel_t)   , intent(in) :: tl
      integer              , intent(in) :: nets
      integer              , intent(in) :: nete
      !=======================================================================================================!
      ! 	Local
      !=======================================================================================================!   
      integer:: i,j,k,ie
      integer:: nm1,n0,np1,nstep
      !=======================================================================================================!
      !=======================================================================================================!
      nm1   = tl%nm1
      n0    = tl%n0
      np1   = tl%np1
      nstep = tl%nstep
      !=======================================================================================================!
      !  Primitive DG Initialization										!
      !=======================================================================================================!
      !=======================================================================================================!
      if (nstep == 0) then
         if (hybrid%par%masterproc) then
            if (test_case=='jw_bcl')   then 
               print *,'DG 3D JW_BCL Data Initialization with Q'
            elseif (test_case=='heldsuarez') then
               print *,'DG 3D Held Suarez Data Initialization'   
            endif
         endif
         !=======================================================================================================!
         ! 	Data Initialization,  mass and wind fields for Baroclinic model 				!
         !=======================================================================================================!    
         if (topology == 'cube') then   
            do ie=nets,nete 
               if (test_case=='jw_bcl') then    
                  call jw_baroclinic(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,          &
                       elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),        &
                       elem(ie)%state%pt3d,elem(ie)%state%dp3d, elem(ie)%state%qt3d)
               elseif (test_case=='heldsuarez') then  
                  call heldsuarez_initial(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,     &
                       elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),   &
                       elem(ie)%state%pt3d,elem(ie)%state%dp3d)
               endif
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,            &
                    elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,              &
                    elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               elem(ie)%state%pr3d_ref(:,:,:)= elem(ie)%state%pr3d(:,:,:)  
               ! Set initial value for Q based on qt3d
               ! I need to pass qt3d by lagrangian_surfvars later
               elem(ie)%state%Q(:,:,:) = elem(ie)%state%qt3d(:,:,:)
               do k=1,nlev    
                  elem(ie)%state%uv(:,:,:,k)  = contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
                  elem(ie)%state%uv0(:,:,:,k) = elem(ie)%state%uv(:,:,:,k)
               enddo
            enddo
         endif
         !=======================================================================================================!
      endif
      !=======================================================================================================!
      !	Time integration with RK3 option 
      !=======================================================================================================!
      if (nstep >= 1) then
         !=======================================================================================================!
         !=======================================================================================================!
         if (hybrid%par%masterproc .and. mod(nstep,statefreq)==0 .and. statefreq>0) then
            print *,'Time-Stepping at timestep:',nstep
         endif
         !call dg3d_rk3(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
         !call dg3d_ssprk2(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
         !call dg3d_uv_ssprk2(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
          call dg3d_uv_ssprk3(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
         
         !=======================================================================================================!
      endif

      !=======================================================================================================!
    end subroutine dg3d_advance

    !=======================================================================================================!
    !      (u,v) test-zone       
    !=======================================================================================================!


    !=======================================================================================================!
    subroutine dg3d_uv_ssprk3(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
      !=======================================================================================================!
      implicit none
      !=======================================================================================================!
      type (Element_t)     , intent(inout), target :: elem(:)
      type (EdgeBuffer_t)  , intent(in) :: edge3
      type (derivative_t)  , intent(in) :: deriv
      type (filter_t)                   :: flt
      type (hybrid_t)      , intent(in) :: hybrid
      type (quadrature_t)               :: gll
      real (kind=real_kind), intent(in) :: dt
      real (kind=real_kind), intent(in) :: pmean
      type (TimeLevel_t)   , intent(in) :: tl
      integer              , intent(in) :: nets
      integer              , intent(in) :: nete
      integer                              :: ig, ntime
      !=======================================================================================================!
      ! Local
      !=======================================================================================================!
      integer, parameter                                                :: neqn=5   !!(# Eqns to be solved)
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dp_rk,dp_zero 
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: ht_rk,ht_zero
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: pt_rk,pt_zero, qt_rk, qt_zero 
      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: uv_rk,uv_zero
      real (kind=real_kind), dimension(np,np,2,nlev)                    :: hs_uv_force 
      real (kind=real_kind), dimension(np,np,nlev)                      :: hs_th_force, tforce, rtlnp
      real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: dpbuf,ptbuf,htbuf, qtbuf 
      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf 
      real (kind=real_kind), dimension(np,np,neqn,nlev,nets:nete)       :: dg3d_rhs, rhs_force
      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: dubuf, dvbuf, pgbuf
      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: dif_uv, difgr1, difgr2, prgr
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dif_th
      real (kind=real_kind), dimension(np,np,2)                         :: gr1,gr2 , lapuv 
      real (kind=real_kind), dimension(np,np)                           :: difu, difv, htop
      !=======================================================================================================!
      real (kind=real_kind) :: zero,one,two,three,four
      real (kind=real_kind) :: lenscale,grv,hdt
      integer    :: i,j,k,ie
      integer    :: kptr
      integer    :: nm1,n0,np1,nstep
      integer    :: eqn
      !=======================================================================================================!
      nm1   = tl%nm1
      n0    = tl%n0
      np1   = tl%np1
      nstep = tl%nstep
      !=======================================================================================================!
      zero  = 0.0D0
      one   = 1.0D0
      two   = 2.0D0
      three = 3.0D0
      four  = 4.0D0
      !=======================================================================================================!
      hdt = dt
      grv = g
      lenscale=rearth
      !=======================================================================================================!
      gll= gausslobatto(np)
      !=======================================================================================================!
      if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
         filter_counter = filter_counter + 1
         uvbuf= zero
         dpbuf= zero
         ptbuf= zero
         qtbuf= zero
         htbuf= zero
         uv_zero= zero
         dp_zero = zero
         pt_zero = zero
         qt_zero = zero
         uv_rk= zero
         dp_rk = zero
         pt_rk = zero
         qt_rk = zero
         dg3d_rhs= zero
         !=======================================================================================================!
         !  Starting Time Integration                                                                            !
         !=======================================================================================================!
         if (nstep >= 1) Then
            !=======================================================================================================!
            !  RK-3 Time Marching                                                                                   !
            !=======================================================================================================!
            do ie=nets, nete
               !=======================================================================================================!
               !  Data Initialization
               !=======================================================================================================!
               if (nstep == 1 ) then
                  if (topology == 'cube' .and. test_case=='jw_bcl') then
                     call jw_baroclinic(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,   &
                          elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),    & 
                          elem(ie)%state%pt3d,elem(ie)%state%dp3d,elem(ie)%state%qt3d)
                  elseif (test_case=='heldsuarez') then
                     call heldsuarez_initial(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,&
                          elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),&
                          elem(ie)%state%pt3d,elem(ie)%state%dp3d)

                          elem(ie)%state%qt3d(:,:,:)  = 0.0D0       !No moisture for HS
                  endif

                     call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,  &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,       &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
                     elem(ie)%state%pr3d_ref(:,:,:)= elem(ie)%state%pr3d(:,:,:)
            !=======================================================================================================!
                  do k=1,nlev
                     elem(ie)%state%uv(:,:,:,k)  = contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
                     elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
                  enddo
               endif
            !=======================================================================================================!
            !       Filter
            !=======================================================================================================!
               do k=1,nlev
                  call filter_P(elem(ie)%state%uv(:,:,1,k),flt)
                  call filter_P(elem(ie)%state%uv(:,:,2,k),flt)
                  call filter_P(elem(ie)%state%dp3d(:,:,k),flt)
                  call filter_P(elem(ie)%state%pt3d(:,:,k),flt)
                  call filter_P(elem(ie)%state%qt3d(:,:,k),flt)
               enddo

               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= elem(ie)%state%uv(:,:,:,k)
               enddo

            !=======================================================================================================!
              ! RHS forcing (psuedo physics) 

              if (test_case=='jw_bcl') then
                           rhs_force(:,:,:,:,ie) = 0.0D0 

               elseif (test_case == 'heldsuarez') then

                  hs_uv_force(:,:,:,:) = heldsuarez_uv_forcing(elem(ie)%state%pr3d,elem(ie)%state%couv)

                  tforce(:,:,:)= heldsuarez_pt_forcing(hdt,elem(ie)%spherep,elem(ie)%state%pr3d, &
                                       elem(ie)%state%T(:,:,:,n0))  

                  do k=1,nlev
                     do j=1,np
                        do i=1,np
                           rhs_force(i,j,1,k,ie) = hs_uv_force(i,j,1,k) 
                           rhs_force(i,j,2,k,ie) = hs_uv_force(i,j,2,k) 
                           rhs_force(i,j,3,k,ie) = 0.0D0
                           rhs_force(i,j,4,k,ie) = tforce(i,j,k) 
                        end do
                     end do
                  end do

             endif

         !=======================================================================================================!
         !       Temporary Storage for fiels from previous time step
         !=======================================================================================================!
               do k=1,nlev
                  uv_zero(:,:,:,k,ie)= elem(ie)%state%uv(:,:,:,k)
                  dp_zero(:,:,k,ie)= elem(ie)%state%dp3d(:,:,k) * elem(ie)%metdet(:,:)
                  pt_zero(:,:,k,ie)= elem(ie)%state%pt3d(:,:,k) * dp_zero(:,:,k,ie)
                  qt_zero(:,:,k,ie)= elem(ie)%state%qt3d(:,:,k) * dp_zero(:,:,k,ie)
               enddo

               !=======================================================================================================!
               !       DG: Packing u,v  velocity & height fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)
               !=======================================================================================================!
            enddo

            call bndry_exchangeV(hybrid,edge3)
       
         if (nu .ne. 0.0D0 )  then      !LDG diffusion activated 
            do ie=nets,nete
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)

             do k=1,nlev
               call dg3d_diff_grads_uv(elem(ie),deriv,uvbuf(:,:,:,k,ie),elem(ie)%state%v(:,:,:,k,n0),gr1,gr2)
                  difgr1(:,:,:,k,ie) = gr1(:,:,:)
                  difgr2(:,:,:,k,ie) = gr2(:,:,:)
               end do
            end do

            do ie = nets, nete
               kptr=6*nlev
               call edgeDGVpack(edge3,reshape(difgr1(:,:,:,:,ie),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=8*nlev
               call edgeDGVpack(edge3,reshape(difgr2(:,:,:,:,ie),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
            end do

            call bndry_exchangeV(hybrid,edge3)

            do ie=nets,nete
               kptr=6*nlev
               call edgeDGVunpack(edge3, dubuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr=8*nlev
               call edgeDGVunpack(edge3, dvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

               do k=1,nlev
                  call dg3d_diff_flux(elem(ie),deriv,dubuf(0,0,1,k,ie),difgr1(1,1,1,k,ie),difu)
                  call dg3d_diff_flux(elem(ie),deriv,dvbuf(0,0,1,k,ie),difgr2(1,1,1,k,ie),difv)

                    ! Converting to [u,v] = A^-T [u_1,u_2] 
                 do j=1,np
                 do i=1,np
                  dif_uv(i,j,1,k,ie) = difu(i,j) *elem(ie)%Dinv(1,1,i,j) + difv(i,j) * elem(ie)%Dinv(2,1,i,j)
                  dif_uv(i,j,2,k,ie) = difu(i,j) *elem(ie)%Dinv(1,2,i,j) + difv(i,j) * elem(ie)%Dinv(2,2,i,j)
                 enddo
                 enddo

               dif_th(:,:,k,ie) = diffusion_theta(elem(ie),deriv,elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k))
               end do
           end do
         
          else             !No diffusion 
              dif_uv = 0.0D0 
              dif_th = 0.0D0 
          
            do ie=nets,nete
               !=======================================================================================================!
               ! SSPDG-RK1: Unpack the edges for uvcomp and height  (nair)
               !=======================================================================================================!
                  kptr=0*nlev
                  call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
                  kptr=2*nlev
                  call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
                  kptr=3*nlev
                  call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
                  kptr=4*nlev
                  call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
                  kptr=5*nlev
                  call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
            enddo
          endif       ! LDG diffusion 

            do ie=nets,nete
               !=======================================================================================================!
               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_uvform_rhs(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),          &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%uv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               !! SSP-RK3 stage-1 computation 
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        uv_rk(i,j,1,k,ie)= uv_zero(i,j,1,k,ie) + hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))
                        uv_rk(i,j,2,k,ie)= uv_zero(i,j,2,k,ie) + hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))
                         pt_rk(i,j,k,ie) = pt_zero(i,j,k,ie) + hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))
                         dp_rk(i,j,k,ie) = dp_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,3,k,ie)
                         qt_rk(i,j,k,ie) = qt_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,5,k,ie) 
                     end do
                  end do
               enddo
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            do ie= nets, nete
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%uv(i,j,1,k)= uv_rk(i,j,1,k,ie)
                        elem(ie)%state%uv(i,j,2,k)= uv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,      &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,     &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               !       Converting covarinat to contravariant and Psi to Height
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= elem(ie)%state%uv(:,:,:,k)
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
               !SSPDG-RK2: Packing  contravariant vectors,  ht-field
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)

            enddo

            call bndry_exchangeV(hybrid,edge3)

            !=======================================================================================================!
            do ie = nets, nete
               !=======================================================================================================!
               !       DG-RK3: Unpack the edges for  flux vectors and scalar fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !=======================================================================================================!
               
               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_uvform_rhs(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),           &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%uv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               !  SSP-RK3 Stage-2                                                                                             !
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        uv_rk(i,j,1,k,ie)=(3.0D0*uv_zero(i,j,1,k,ie) + uv_rk(i,j,1,k,ie)+  &
                                               hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))) *0.25D0
                        uv_rk(i,j,2,k,ie)=(3.0D0*uv_zero(i,j,2,k,ie) + uv_rk(i,j,2,k,ie)+  &
                                               hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))) *0.25D0
                        pt_rk(i,j,k,ie) = (3.0D0*pt_zero(i,j,k,ie) + pt_rk(i,j,k,ie)+ &
                                               hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))) *0.25D0
                        dp_rk(i,j,k,ie) = (3.0D0*dp_zero(i,j,k,ie) + dp_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,3,k,ie)) *0.25D0
                        qt_rk(i,j,k,ie) = (3.0D0*qt_zero(i,j,k,ie) + qt_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,5,k,ie)) *0.25D0
                     end do
                  end do
               end do
            
           enddo


            do ie= nets, nete
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%uv(i,j,1,k)= uv_rk(i,j,1,k,ie)
                        elem(ie)%state%uv(i,j,2,k)= uv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,      &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,     &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= elem(ie)%state%uv(:,:,:,k)
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
               ! Packing  u,v  vectors,  ht-field
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)

            enddo

            call bndry_exchangeV(hybrid,edge3)

            !=======================================================================================================!
            do ie = nets, nete
               !=======================================================================================================!
               !       DG-RK3: Unpack the edges for  flux vectors and scalar fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !=======================================================================================================!

               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_uvform_rhs(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),          &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%uv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               !! SSP-RK3 stage-3 computations 
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        uv_rk(i,j,1,k,ie)=(uv_zero(i,j,1,k,ie) + 2.0D0*(uv_rk(i,j,1,k,ie) +  &
                                            hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))) )/ 3.0D0 
                        uv_rk(i,j,2,k,ie)=(uv_zero(i,j,2,k,ie) + 2.0D0*(uv_rk(i,j,2,k,ie) +  &
                                            hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))) )/ 3.0D0 
                        pt_rk(i,j,k,ie) = (pt_zero(i,j,k,ie) + 2.0D0*(pt_rk(i,j,k,ie)+ &
                                            hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))) )/ 3.0D0 
                        dp_rk(i,j,k,ie) = (dp_zero(i,j,k,ie) + 2.0D0*(dp_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,3,k,ie)) )/ 3.0D0 
                        qt_rk(i,j,k,ie) = (qt_zero(i,j,k,ie) + 2.0D0*(qt_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,5,k,ie)) )/ 3.0D0 
                     end do
                  end do
               end do
            
           enddo
            !=======================================================================================================!
            ! Extract new evolved variables 
            !=======================================================================================================!
            do ie= nets, nete
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%uv(i,j,1,k)= uv_rk(i,j,1,k,ie)
                        elem(ie)%state%uv(i,j,2,k)= uv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               !       Remapping Zone (ssp-rk2)
               !=======================================================================================================!
               if ((nstep > 1).and.(mod(nstep,remapfreq) == 0)) then
                  !=======================================================================================================!
                  if (remap_type == "linear") then
                     call linear_remap(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,                     &
                          elem(ie)%state%uv,elem(ie)%state%dp3d)
                  else if (remap_type == "parabolic") then
                     call parabolic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                          elem(ie)%state%uv,elem(ie)%state%dp3d,elem(ie)%state%qt3d)
                  else if (remap_type == "monotonic") then
                     call monotonic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                          elem(ie)%state%uv,elem(ie)%state%dp3d)
                  else if (remap_type == "energy") then
                     call energy_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,	&
                          elem(ie)%state%uv,elem(ie)%state%dp3d)
                  endif
                  !=======================================================================================================!
               endif

               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,              &
                                        elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,                &
                                        elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               !       Converting covarinat to contravariant and Psi to Height
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= sphere2contra(elem(ie)%state%uv(:,:,:,k),elem(ie)%Dinv(:,:,:,:))
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !   Ending Time Step nstep >= 1                                                                         !
            !=======================================================================================================!
         ENDIF
      endif
      !=======================================================================================================!

      ! Remap the potential temperature to real temperature for this time step
      do ie= nets, nete
         elem(ie)%state%T(:,:,:,n0) = pt2temp(elem(ie)%state%pr3d,elem(ie)%state%pt3d)
         ! The moist is just the qt3d field which may be converted later to Q by something else.
         ! if not, this is just a waste of memory.
         elem(ie)%state%Q(:,:,:) = elem(ie)%state%qt3d(:,:,:)
      end do

    end subroutine dg3d_uv_ssprk3  
!+++++++++++++++++++++++++++++++++++++++++++++++ (u,v) formulation end +++++++++++++++
!========= =============================================================!
    subroutine dg3d_uv_ssprk2(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
      !=======================================================================================================!
      !=======================================================================================================!
      implicit none
      !=======================================================================================================!
      type (Element_t)     , intent(inout), target :: elem(:)
      type (EdgeBuffer_t)  , intent(in) :: edge3
      type (derivative_t)  , intent(in) :: deriv
      type (filter_t)                   :: flt
      type (hybrid_t)      , intent(in) :: hybrid
      type (quadrature_t)               :: gll
      real (kind=real_kind), intent(in) :: dt
      real (kind=real_kind), intent(in) :: pmean
      type (TimeLevel_t)   , intent(in) :: tl
      integer              , intent(in) :: nets
      integer              , intent(in) :: nete
      integer                              :: ig, ntime
      !=======================================================================================================!
      ! Local
      !=======================================================================================================!
      !=======================================================================================================!
      integer, parameter                                                :: neqn=5   !!(# Eqns to be solved)
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dp_rk,dp_zero 
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: ht_rk,ht_zero
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: pt_rk,pt_zero, qt_rk, qt_zero 
      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: cuv_rk,cuv_zero
      real (kind=real_kind), dimension(np,np,2,nlev)                    :: hs_uv_force 
      real (kind=real_kind), dimension(np,np,nlev)                      :: hs_th_force, tforce, rtlnp
      real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: dpbuf,ptbuf,htbuf, qtbuf 
      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf 
      real (kind=real_kind), dimension(np,np,neqn,nlev,nets:nete)       :: dg3d_rhs, rhs_force
      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: dubuf, dvbuf, pgbuf
      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: dif_uv, difgr1, difgr2, prgr
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dif_th
      real (kind=real_kind), dimension(np,np,2)                         :: gr1,gr2 , lapuv 
      real (kind=real_kind), dimension(np,np)                           :: difu, difv, htop
      !=======================================================================================================!
      real (kind=real_kind) :: zero,one,two,three,four
      real (kind=real_kind) :: lenscale,grv,hdt
      integer    :: i,j,k,ie
      integer    :: kptr
      integer    :: nm1,n0,np1,nstep
      integer    :: eqn
      !=======================================================================================================!
      !=======================================================================================================!
      nm1   = tl%nm1
      n0    = tl%n0
      np1   = tl%np1
      nstep = tl%nstep
      !=======================================================================================================!
      zero  = 0.0D0
      one   = 1.0D0
      two   = 2.0D0
      three = 3.0D0
      four  = 4.0D0
      !=======================================================================================================!
      hdt = dt
      grv = g
      lenscale=rearth
      !=======================================================================================================!
      gll= gausslobatto(np)
      !=======================================================================================================!
      if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
         filter_counter = filter_counter + 1
         uvbuf= zero
         dpbuf= zero
         ptbuf= zero
         qtbuf= zero
         htbuf= zero
         cuv_zero= zero
         dp_zero = zero
         pt_zero = zero
         qt_zero = zero
         cuv_rk= zero
         dp_rk = zero
         pt_rk = zero
         qt_rk = zero
         dg3d_rhs= zero
         !=======================================================================================================!
         !=======================================================================================================!
         !  Starting Time Integration                                                                            !
         !=======================================================================================================!
         if (nstep >= 1) Then
            !=======================================================================================================!
            !  RK-3 Time Marching                                                                                   !
            !=======================================================================================================!
            do ie=nets, nete
               !=======================================================================================================!
               !  Data Initialization
               !=======================================================================================================!
               if (nstep == 1 ) then
                  if (topology == 'cube' .and. test_case=='jw_bcl') then
                     call jw_baroclinic(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,   &
                          elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),    & 
                          elem(ie)%state%pt3d,elem(ie)%state%dp3d,elem(ie)%state%qt3d)

                     call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,  &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,       &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
                     elem(ie)%state%pr3d_ref(:,:,:)= elem(ie)%state%pr3d(:,:,:)
                  endif
                  !=======================================================================================================!
                  !  DG: Compute  covariant components for the first time
                  !=======================================================================================================!
                  do k=1,nlev
                     elem(ie)%state%uv(:,:,:,k)  = contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
                     !elem(ie)%state%uv0(:,:,:,k) = elem(ie)%state%uv(:,:,:,k)
                     !elem(ie)%state%couv(:,:,:,k) = elem(ie)%state%uv(:,:,:,k)
                     !elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
                     elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
                  enddo
               endif
               !=======================================================================================================!
               !       Filter
               !=======================================================================================================!
               do k=1,nlev
                  call filter_P(elem(ie)%state%uv(:,:,1,k),flt)
                  call filter_P(elem(ie)%state%uv(:,:,2,k),flt)
                  call filter_P(elem(ie)%state%dp3d(:,:,k),flt)
                  call filter_P(elem(ie)%state%pt3d(:,:,k),flt)
                  call filter_P(elem(ie)%state%qt3d(:,:,k),flt)
               enddo
               !=======================================================================================================!
               !       Converting covarinat to contravariant (second timestep onwards)
               !=======================================================================================================!
               do k=1,nlev
                  !lem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
                  elem(ie)%state%v(:,:,:,k,n0)= elem(ie)%state%uv(:,:,:,k)
               enddo
               !=======================================================================================================!
               !       Temporary Storage for fiels from previous time step
               !=======================================================================================================!
               do k=1,nlev
                  cuv_zero(:,:,:,k,ie)= elem(ie)%state%uv(:,:,:,k)
                  dp_zero(:,:,k,ie)= elem(ie)%state%dp3d(:,:,k) * elem(ie)%metdet(:,:)
                  pt_zero(:,:,k,ie)= elem(ie)%state%pt3d(:,:,k) * dp_zero(:,:,k,ie)
                  qt_zero(:,:,k,ie)= elem(ie)%state%qt3d(:,:,k) * dp_zero(:,:,k,ie)
                  rhs_force(:,:,:,k,ie)= 0.0D0
               enddo

               !=======================================================================================================!
               !       DG: Packing Contravariant velocity & height fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            call bndry_exchangeV(hybrid,edge3)
            !=======================================================================================================!
            !=======================================================================================================!
            !   LDG-Step1
            !=======================================================================================================!
            do ie=nets,nete
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)

             do k=1,nlev
               call dg3d_diff_grads_uv(elem(ie),deriv,uvbuf(:,:,:,k,ie),elem(ie)%state%v(:,:,:,k,n0),gr1,gr2)
                  difgr1(:,:,:,k,ie) = gr1(:,:,:)
                  difgr2(:,:,:,k,ie) = gr2(:,:,:)
               end do
            end do

            do ie = nets, nete
               kptr=6*nlev
               call edgeDGVpack(edge3,reshape(difgr1(:,:,:,:,ie),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=8*nlev
               call edgeDGVpack(edge3,reshape(difgr2(:,:,:,:,ie),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
            end do

            call bndry_exchangeV(hybrid,edge3)

            do ie=nets,nete
               kptr=6*nlev
               call edgeDGVunpack(edge3, dubuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr=8*nlev
               call edgeDGVunpack(edge3, dvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

               do k=1,nlev
                  call dg3d_diff_flux(elem(ie),deriv,dubuf(0,0,1,k,ie),difgr1(1,1,1,k,ie),difu)
                  call dg3d_diff_flux(elem(ie),deriv,dvbuf(0,0,1,k,ie),difgr2(1,1,1,k,ie),difv)
                  !dif_uv(:,:,1,k,ie) = difu(:,:)
                  !dif_uv(:,:,2,k,ie) = difv(:,:)

             do j=1,np
             do i=1,np
                 dif_uv(i,j,1,k,ie) = difu(i,j) *elem(ie)%Dinv(1,1,i,j) + difv(i,j) * elem(ie)%Dinv(2,1,i,j)
                 dif_uv(i,j,2,k,ie) = difu(i,j) *elem(ie)%Dinv(1,2,i,j) + difv(i,j) * elem(ie)%Dinv(2,2,i,j)
             enddo
             enddo

               dif_th(:,:,k,ie) = diffusion_theta(elem(ie),deriv,elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k))
               end do
            end do
            !=======================================================================================================!
            !     End of LDG stages
            !=======================================================================================================!
            !=======================================================================================================!
            !  RK-1 Stage                                                                                           !
            !=======================================================================================================!
         
          
            do ie=nets,nete
               !=======================================================================================================!
               ! SSPDG-RK1: Unpack the edges for uvcomp and height  (nair)
               !=======================================================================================================!
               !  kptr=0*nlev
               !  call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               !  kptr=2*nlev
               !  call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !  kptr=3*nlev
               !  call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !  kptr=4*nlev
               !  call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !  kptr=5*nlev
               !  call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)

               !=======================================================================================================!
               tforce(:,:,:) = 0.0D0
               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               !=======================================================================================================!
               rhs_force(:,:,:,:,ie) = 0.0D0
             ! dif_uv(:,:,:,:,ie) = 0.0D0
             ! dif_th(:,:,:,ie) = 0.0D0
               !=======================================================================================================!
               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_uvform_rhs(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),          &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%uv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        cuv_rk(i,j,1,k,ie)= cuv_zero(i,j,1,k,ie) + hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))
                        cuv_rk(i,j,2,k,ie)= cuv_zero(i,j,2,k,ie) + hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))
                        dp_rk(i,j,k,ie)   = dp_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,3,k,ie)
                        pt_rk(i,j,k,ie)   = pt_zero(i,j,k,ie) + hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))
                        qt_rk(i,j,k,ie)   = qt_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,5,k,ie) 
                     end do
                  end do
               enddo
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !=======================================================================================================!
            do ie= nets, nete
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%uv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
                        elem(ie)%state%uv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,      &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,     &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               !       Converting covarinat to contravariant and Psi to Height
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= elem(ie)%state%uv(:,:,:,k)
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
               !SSPDG-RK2: Packing  contravariant vectors,  ht-field
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)

            enddo
            !=======================================================================================================!
            call bndry_exchangeV(hybrid,edge3)
            !=======================================================================================================!
            !  SSP RK-2 STAGE                                                                                               !
            !=======================================================================================================!
            do ie = nets, nete
               !=======================================================================================================!
               !       DG-RK3: Unpack the edges for  flux vectors and scalar fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !=======================================================================================================!
               tforce(:,:,:) = 0.0D0
               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               rhs_force(:,:,:,:,ie) = 0.0D0
               !dif_uv(:,:,:,:,ie) = 0.0D0
               !dif_th(:,:,:,ie) = 0.0D0
               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_uvform_rhs(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),           &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%uv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        cuv_rk(i,j,1,k,ie)=(cuv_zero(i,j,1,k,ie) + cuv_rk(i,j,1,k,ie)+  &
                                            hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))) *0.5D0
                        cuv_rk(i,j,2,k,ie)=(cuv_zero(i,j,2,k,ie) + cuv_rk(i,j,2,k,ie)+  &
                                            hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))) *0.5D0
                        dp_rk(i,j,k,ie) = (dp_zero(i,j,k,ie) + dp_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,3,k,ie)) *0.5D0
                        pt_rk(i,j,k,ie) = (pt_zero(i,j,k,ie) + pt_rk(i,j,k,ie)+ &
                                            hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))) *0.5D0
                        qt_rk(i,j,k,ie) = (qt_zero(i,j,k,ie) + qt_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,5,k,ie)) *0.5D0
                     end do
                  end do
               end do
            
           enddo
            !=======================================================================================================!
            ! Extract new evolved variables 
            !=======================================================================================================!
            do ie= nets, nete
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%uv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
                        elem(ie)%state%uv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               !       Remapping Zone (ssp-rk2)
               !=======================================================================================================!
               if ((nstep > 1).and.(mod(nstep,remapfreq) == 0)) then
                  !=======================================================================================================!
                  if (remap_type == "linear") then
                     call linear_remap(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,                     &
                          elem(ie)%state%uv,elem(ie)%state%dp3d)
                  else if (remap_type == "parabolic") then
                     call parabolic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                          elem(ie)%state%uv,elem(ie)%state%dp3d,elem(ie)%state%qt3d)
                  else if (remap_type == "monotonic") then
                     call monotonic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                          elem(ie)%state%uv,elem(ie)%state%dp3d)
                  else if (remap_type == "energy") then
                     call energy_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,	&
                          elem(ie)%state%uv,elem(ie)%state%dp3d)
                  endif
                  !=======================================================================================================!
               endif
               !=======================================================================================================!
               !=======================================================================================================!
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,              &
                                        elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,                &
                                        elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               !       Converting covarinat to contravariant and Psi to Height
               !=======================================================================================================!
               do k=1,nlev
                  !lem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
                  elem(ie)%state%v(:,:,:,k,n0)= sphere2contra(elem(ie)%state%uv(:,:,:,k),elem(ie)%Dinv(:,:,:,:))
                  !elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
               !=======================================================================================================!
               !       Temprature field: from 1 to nlev                                                                !
               !       Surface Pressure: at nlev+1                                                                     !
               !=======================================================================================================!
               do k=1,nlev
                  !elem(ie)%state%uv(:,:,:,k)= contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
                  !elem(ie)%state%zeta(:,:,k)= vorticity(elem(ie)%state%couv(:,:,:,k),deriv)/elem(ie)%metdet(:,:)
                  !elem(ie)%state%zeta(:,:,k)= cov_vorticity(deriv,elem(ie)%state%couv(:,:,:,k))/elem(ie)%metdet(:,:)
               end do
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !   Ending Time Step nstep >= 1                                                                         !
            !=======================================================================================================!
         ENDIF
         !=======================================================================================================!
      endif
      !=======================================================================================================!
      !=======================================================================================================!

      ! Remap the potential temperature to real temperature for this time step
      do ie= nets, nete
         elem(ie)%state%T(:,:,:,n0) = pt2temp(elem(ie)%state%pr3d,elem(ie)%state%pt3d)
         ! The moist is just the qt3d field which may be converted later to Q by something else.
         ! if not, this is just a waste of memory.
         elem(ie)%state%Q(:,:,:) = elem(ie)%state%qt3d(:,:,:)
      end do

    end subroutine dg3d_uv_ssprk2
!+++++++++++++++++++++++++++++++++++++++++++++++ (u,v) formulation end +++++++++++++++

!=======================================================================================================!
!========= covariant (classic) formulation =============================================================!
    subroutine dg3d_ssprk2(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
      !=======================================================================================================!
      !=======================================================================================================!
      implicit none
      !=======================================================================================================!
      type (Element_t)     , intent(inout), target :: elem(:)
      type (EdgeBuffer_t)  , intent(in) :: edge3
      type (derivative_t)  , intent(in) :: deriv
      type (filter_t)                   :: flt
      type (hybrid_t)      , intent(in) :: hybrid
      type (quadrature_t)               :: gll
      real (kind=real_kind), intent(in) :: dt
      real (kind=real_kind), intent(in) :: pmean
      type (TimeLevel_t)   , intent(in) :: tl
      integer              , intent(in) :: nets
      integer              , intent(in) :: nete
      integer                              :: ig, ntime
      !=======================================================================================================!
      ! Local
      !=======================================================================================================!
      !=======================================================================================================!
      integer, parameter                                                :: neqn=5   !!(# Eqns to be solved)
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dp_rk,dp_zero 
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: ht_rk,ht_zero
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: pt_rk,pt_zero, qt_rk, qt_zero 
      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: cuv_rk,cuv_zero
      real (kind=real_kind), dimension(np,np,2,nlev)                    :: hs_uv_force 
      real (kind=real_kind), dimension(np,np,nlev)                      :: hs_th_force, tforce, rtlnp
      real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: dpbuf,ptbuf,htbuf, qtbuf 
      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf 
      real (kind=real_kind), dimension(np,np,neqn,nlev,nets:nete)       :: dg3d_rhs, rhs_force
      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: dubuf, dvbuf, pgbuf
      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: dif_uv, difgr1, difgr2, prgr
      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dif_th
      real (kind=real_kind), dimension(np,np,2)                         :: gr1,gr2 , lapuv 
      real (kind=real_kind), dimension(np,np)                           :: difu, difv, htop
      !=======================================================================================================!
      real (kind=real_kind) :: zero,one,two,three,four
      real (kind=real_kind) :: lenscale,grv,hdt
      integer    :: i,j,k,ie
      integer    :: kptr
      integer    :: nm1,n0,np1,nstep
      integer    :: eqn
      !=======================================================================================================!
      !=======================================================================================================!
      nm1   = tl%nm1
      n0    = tl%n0
      np1   = tl%np1
      nstep = tl%nstep
      !=======================================================================================================!
      zero  = 0.0D0
      one   = 1.0D0
      two   = 2.0D0
      three = 3.0D0
      four  = 4.0D0
      !=======================================================================================================!
      hdt = dt
      grv = g
      lenscale=rearth
      !=======================================================================================================!
      gll= gausslobatto(np)
      !=======================================================================================================!
      if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
         filter_counter = filter_counter + 1
         uvbuf= zero
         dpbuf= zero
         ptbuf= zero
         qtbuf= zero
         htbuf= zero
         cuv_zero= zero
         dp_zero = zero
         pt_zero = zero
         qt_zero = zero
         cuv_rk= zero
         dp_rk = zero
         pt_rk = zero
         qt_rk = zero
         dg3d_rhs= zero
         !=======================================================================================================!
         !=======================================================================================================!
         !  Starting Time Integration                                                                            !
         !=======================================================================================================!
         if (nstep >= 1) Then
            !=======================================================================================================!
            !  RK-3 Time Marching                                                                                   !
            !=======================================================================================================!
            do ie=nets, nete
               !=======================================================================================================!
               !  Data Initialization
               !=======================================================================================================!
               if (nstep == 1 ) then
                  if (topology == 'cube' .and. test_case=='jw_bcl') then
                     call jw_baroclinic(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,   &
                          elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),    & 
                          elem(ie)%state%pt3d,elem(ie)%state%dp3d,elem(ie)%state%qt3d)
                   ! call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,          &
                   !      elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,            &
                   !      elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
                     call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,  &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,       &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
                     elem(ie)%state%pr3d_ref(:,:,:)= elem(ie)%state%pr3d(:,:,:)
                  endif
                  !=======================================================================================================!
                  !  DG: Compute  covariant components for the first time
                  !=======================================================================================================!
                  do k=1,nlev
                     elem(ie)%state%uv(:,:,:,k)  = contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
                     elem(ie)%state%uv0(:,:,:,k) = elem(ie)%state%uv(:,:,:,k)
                     elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
                     elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
                  enddo
               endif
               !=======================================================================================================!
               !       Filter
               !=======================================================================================================!
               do k=1,nlev
                  call filter_P(elem(ie)%state%couv(:,:,1,k),flt)
                  call filter_P(elem(ie)%state%couv(:,:,2,k),flt)
                  call filter_P(elem(ie)%state%dp3d(:,:,k),flt)
                  call filter_P(elem(ie)%state%pt3d(:,:,k),flt)
                  call filter_P(elem(ie)%state%qt3d(:,:,k),flt)
               enddo
               !=======================================================================================================!
               !       Converting covarinat to contravariant (second timestep onwards)
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
               enddo
               !=======================================================================================================!
               !       Temporary Storage for fiels from previous time step
               !=======================================================================================================!
               do k=1,nlev
                  cuv_zero(:,:,:,k,ie)= elem(ie)%state%couv(:,:,:,k)
                  dp_zero(:,:,k,ie)= elem(ie)%state%dp3d(:,:,k) * elem(ie)%metdet(:,:)
                  pt_zero(:,:,k,ie)= elem(ie)%state%pt3d(:,:,k) * dp_zero(:,:,k,ie)
                  qt_zero(:,:,k,ie)= elem(ie)%state%qt3d(:,:,k) * dp_zero(:,:,k,ie)
                  rhs_force(:,:,:,k,ie)= 0.0D0
               enddo

               !=======================================================================================================!
               !       DG: Packing Contravariant velocity & height fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)
               !=======================================================================================================!
               !       DG: Rotating velocity (contra)
               !=======================================================================================================!
               kptr=0*nlev
               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !       Insert communications here: for shared memory, just a single
            !       sync is required
            !=======================================================================================================!
            call bndry_exchangeV(hybrid,edge3)
            !=======================================================================================================!
            !=======================================================================================================!
            !   LDG-Step1
            !=======================================================================================================!
            do ie=nets,nete

               !      !============================================================
               !      ! Unpack the edges for uvcomp and height  (nair)
               !      !============================================================
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)

               kptr=3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)

               do k=1,nlev

                  call dg3d_diff_grads(elem(ie),deriv,uvbuf(:,:,:,k,ie),elem(ie)%state%v(:,:,:,k,n0), &
                       elem(ie)%state%couv(:,:,:,k), gr1,gr2)


                  difgr1(:,:,:,k,ie) = gr1(:,:,:)
                  difgr2(:,:,:,k,ie) = gr2(:,:,:)
               end do
            end do
            do ie = nets, nete
               kptr=6*nlev
               call edgeDGVpack(edge3,reshape(difgr1(:,:,:,:,ie),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=6*nlev
               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
               kptr=8*nlev
               call edgeDGVpack(edge3,reshape(difgr2(:,:,:,:,ie),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=8*nlev
               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
            end do

            call bndry_exchangeV(hybrid,edge3)

            do ie=nets,nete
               kptr=6*nlev
               call edgeDGVunpack(edge3, dubuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr=8*nlev
               call edgeDGVunpack(edge3, dvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)


               do k=1,nlev
                  call dg3d_diff_flux(elem(ie),deriv,dubuf(0,0,1,k,ie),difgr1(1,1,1,k,ie),difu)
                  call dg3d_diff_flux(elem(ie),deriv,dvbuf(0,0,1,k,ie),difgr2(1,1,1,k,ie),difv)
                  dif_uv(:,:,1,k,ie) = difu(:,:)
                  dif_uv(:,:,2,k,ie) = difv(:,:)

                  !lapuv(:,:,1) = difu(:,:) 
                  !lapuv(:,:,2) = difv(:,:) 
                  !dif_uv(:,:,:,k,ie) = diffusion_hypr(elem(ie),deriv,lapuv)

                  dif_th(:,:,k,ie) = diffusion_theta(elem(ie),deriv,elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k))

               end do
            end do
            !=======================================================================================================!
            !     End of LDG stages
            !=======================================================================================================!

            !  RK-1 Stage                                                                                           !
            !=======================================================================================================!
            !=======================================================================================================!
            !=======================================================================================================!
            do ie=nets,nete
               !=======================================================================================================!
               ! SSPDG-RK1: Unpack the edges for uvcomp and height  (nair)
               !=======================================================================================================!
               !  kptr=0*nlev
               !  call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               !  kptr=2*nlev
               !  call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !  kptr=3*nlev
               !  call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !  kptr=4*nlev
               !  call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !=======================================================================================================!
               tforce(:,:,:) = 0.0D0
               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               !=======================================================================================================!
               rhs_force(:,:,:,:,ie) = 0.0D0
               !     dif_uv(:,:,:,:,ie) = 0.0D0
               !=======================================================================================================!
               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_rhs_terms(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),           &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%couv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        cuv_rk(i,j,1,k,ie)= cuv_zero(i,j,1,k,ie) + hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))
                        cuv_rk(i,j,2,k,ie)= cuv_zero(i,j,2,k,ie) + hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))
                        dp_rk(i,j,k,ie)   = dp_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,3,k,ie)
                        pt_rk(i,j,k,ie)   = pt_zero(i,j,k,ie) + hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))
                        qt_rk(i,j,k,ie)   = qt_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,5,k,ie) 
                     end do
                  end do
               enddo
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !=======================================================================================================!
            do ie= nets, nete
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%couv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
                        elem(ie)%state%couv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,      &
                                           elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,     &
                                           elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               !       Converting covarinat to contravariant and Psi to Height
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
               !SSPDG-RK2: Packing  contravariant vectors,  ht-field
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
               kptr=2*nlev
               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
               kptr=3*nlev
               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
               kptr=4*nlev
               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
               kptr=5*nlev
               call edgeDGVpack(edge3,elem(ie)%state%qt3d,nlev,kptr,elem(ie)%desc)
               !=======================================================================================================!
               !       DG: Rotating velocity (contra)
               !=======================================================================================================!
               kptr=0*nlev
               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !=======================================================================================================!
            !       Insert communications here: for shared memory, just a single
            !       sync is required
            !=======================================================================================================!
            call bndry_exchangeV(hybrid,edge3)
            !=======================================================================================================!
            !=======================================================================================================!
            !  SSP RK-2 STAGE                                                                                               !
            !=======================================================================================================!
            !=======================================================================================================!
            !=======================================================================================================!
            do ie = nets, nete
               !=======================================================================================================!
               !       DG-RK3: Unpack the edges for  flux vectors and scalar fields
               !=======================================================================================================!
               kptr=0*nlev
               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
               kptr =2*nlev
               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =3*nlev
               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =4*nlev
               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               kptr =5*nlev
               call edgeDGVunpack(edge3, qtbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
               !=======================================================================================================!
               tforce(:,:,:) = 0.0D0
               call gradient_p3d(elem(ie),deriv,elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0),elem(ie)%state%pgrads)

               rhs_force(:,:,:,:,ie) = 0.0D0
               !   dif_uv(:,:,:,:,ie) = 0.0D0
               !=======================================================================================================!
               !       DG3D Model RHS                                                                                  !
               !=======================================================================================================!
               htop(:,:)=elem(ie)%state%ht(:,:,1)
               do k=1,nlev
                  call dg3d_rhs_terms(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),           &
                         dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),qtbuf(:,:,k,ie),                                &
                         elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%couv(:,:,:,k),elem(ie)%fcor,        &
                         elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
                         elem(ie)%state%qt3d(:,:,k),htop,elem(ie)%state%pgrads(:,:,:,k),                 & 
                         rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
               enddo
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        cuv_rk(i,j,1,k,ie)=(cuv_zero(i,j,1,k,ie) + cuv_rk(i,j,1,k,ie)+  &
                                            hdt * (dg3d_rhs(i,j,1,k,ie) + dif_uv(i,j,1,k,ie))) *0.5D0
                        cuv_rk(i,j,2,k,ie)=(cuv_zero(i,j,2,k,ie) + cuv_rk(i,j,2,k,ie)+  &
                                            hdt * (dg3d_rhs(i,j,2,k,ie) + dif_uv(i,j,2,k,ie))) *0.5D0
                        dp_rk(i,j,k,ie) = (dp_zero(i,j,k,ie) + dp_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,3,k,ie)) *0.5D0
                        pt_rk(i,j,k,ie) = (pt_zero(i,j,k,ie) + pt_rk(i,j,k,ie)+ &
                                            hdt * (dg3d_rhs(i,j,4,k,ie) + dif_th(i,j,k,ie))) *0.5D0
                        qt_rk(i,j,k,ie) = (qt_zero(i,j,k,ie) + qt_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,5,k,ie)) *0.5D0
                     end do
                  end do
               end do
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            do ie= nets, nete
               !=======================================================================================================!
               do k=1,nlev
                  do j=1,np
                     do i=1,np
                        elem(ie)%state%couv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
                        elem(ie)%state%couv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                        elem(ie)%state%qt3d(i,j,k)  = qt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
                     end do
                  end do
               end do

               !=======================================================================================================!
               !       Remapping Zone (ssp-rk2)
               !=======================================================================================================!
               if ((nstep > 1).and.(mod(nstep,remapfreq) == 0)) then
                  !=======================================================================================================!
                  if (remap_type == "linear") then
                     call linear_remap(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,                     &
                          elem(ie)%state%couv,elem(ie)%state%dp3d)
                  else if (remap_type == "parabolic") then
                     call parabolic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                          elem(ie)%state%couv,elem(ie)%state%dp3d,elem(ie)%state%qt3d)
                   ! call linear_remap(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,                     &
                   !       elem(ie)%state%couv,elem(ie)%state%dp3d)
                   ! call monotonic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                   !      elem(ie)%state%couv,elem(ie)%state%dp3d)
                   ! call energy_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,	&
                   !      elem(ie)%state%couv,elem(ie)%state%dp3d)
                  else if (remap_type == "monotonic") then
                     call monotonic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,     &
                          elem(ie)%state%couv,elem(ie)%state%dp3d)
                  else if (remap_type == "energy") then
                     call energy_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,	&
                          elem(ie)%state%couv,elem(ie)%state%dp3d)
                  endif
                  !=======================================================================================================!
               endif
               !=======================================================================================================!
               !=======================================================================================================!
               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,              &
                                        elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,                &
                                        elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
               !=======================================================================================================!
               !       Converting covarinat to contravariant and Psi to Height
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))
                  !elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
               enddo
               !=======================================================================================================!
               !=======================================================================================================!
               !       Temprature field: from 1 to nlev                                                                !
               !       Surface Pressure: at nlev+1                                                                     !
               !=======================================================================================================!
               do k=1,nlev
                  elem(ie)%state%uv(:,:,:,k)= contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
                  !elem(ie)%state%zeta(:,:,k)= vorticity(elem(ie)%state%couv(:,:,:,k),deriv)/elem(ie)%metdet(:,:)
                  elem(ie)%state%zeta(:,:,k)= cov_vorticity(deriv,elem(ie)%state%couv(:,:,:,k))/elem(ie)%metdet(:,:)
               end do
               !=======================================================================================================!
            enddo
            !=======================================================================================================!
            !   Ending Time Step nstep >= 1                                                                         !
            !=======================================================================================================!
         ENDIF
         !=======================================================================================================!
      endif
      !=======================================================================================================!
      !=======================================================================================================!

      ! Remap the potential temperature to real temperature for this time step
      do ie= nets, nete
         elem(ie)%state%T(:,:,:,n0) = pt2temp(elem(ie)%state%pr3d,elem(ie)%state%pt3d)
         ! The moist is just the qt3d field which may be converted later to Q by something else.
         ! if not, this is just a waste of memory.
         elem(ie)%state%Q(:,:,:) = elem(ie)%state%qt3d(:,:,:)
      end do

    end subroutine dg3d_ssprk2
    !+++++++++++++++++++++++++
    !=======================================================================================================!
    !=======================================================================================================!
    !	DG3D_RK3											!
    !=======================================================================================================!
!    subroutine dg3d_rk3(elem,edge3,deriv,flt,hybrid,dt,pmean,tl,nets,nete)
!      !=======================================================================================================!
!      !=======================================================================================================!
!      implicit none
!      !=======================================================================================================!
!      type (Element_t)     , intent(inout), target :: elem(:)
!      type (EdgeBuffer_t)  , intent(in) :: edge3
!      type (derivative_t)  , intent(in) :: deriv
!      type (filter_t)                   :: flt
!      type (hybrid_t)      , intent(in) :: hybrid
!      type (quadrature_t)               :: gll
!      real (kind=real_kind), intent(in) :: dt
!      real (kind=real_kind), intent(in) :: pmean
!      type (TimeLevel_t)   , intent(in) :: tl
!      integer              , intent(in) :: nets
!      integer              , intent(in) :: nete
!      integer                              :: ig, ntime
!      !=======================================================================================================!
!      ! Local
!      !=======================================================================================================! 
!      !=======================================================================================================!  
!      integer, parameter                                                :: neqn=4   !!(# Eqns to be solved)
!      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: dp_rk,dp_zero
!      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: ht_rk,ht_zero, old_t3d
!      real (kind=real_kind), dimension(np,np,nlev,nets:nete)            :: pt_rk,pt_zero
!      real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)          :: cuv_rk,cuv_zero
!      real (kind=real_kind), dimension(np,np,2,nlev)                    :: hs_uv_force
!      real (kind=real_kind), dimension(np,np,nlev)                      :: hs_th_force,tforce
!      real (kind=real_kind), dimension(0:np+1,0:np+1,nlev,nets:nete)    :: dpbuf,ptbuf,htbuf
!      real (kind=real_kind), dimension(0:np+1,0:np+1,2,nlev,nets:nete)  :: uvbuf 
!      real (kind=real_kind), dimension(np,np,neqn,nlev,nets:nete)       :: dg3d_rhs, rhs_force
!      real (kind=real_kind), dimension(np,np,2)                         :: diffuv              
!      real (kind=real_kind), dimension(np,np)                           :: diff_pott , htop 
!      !real (kind=real_kind), dimension(np,np,5,nlev,nets:nete):: dg3d_rhs
!      !=======================================================================================================!
!      real (kind=real_kind) :: zero,one,two,three,four
!      real (kind=real_kind) :: lenscale,grv,hdt 
!      integer    :: i,j,k,ie
!      integer    :: kptr
!      integer    :: nm1,n0,np1,nstep
!      integer    :: eqn
!      !=======================================================================================================!
!      !=======================================================================================================!
!      nm1   = tl%nm1
!      n0    = tl%n0
!      np1   = tl%np1
!      nstep = tl%nstep
!      !=======================================================================================================!   
!      zero  = 0.0D0
!      one   = 1.0D0
!      two   = 2.0D0
!      three = 3.0D0
!      four  = 4.0D0
!      !=======================================================================================================!
!      hdt = dt
!      grv = g  
!      lenscale=rearth
!      !=======================================================================================================!
!      !=======================================================================================================!
!      if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0) then
!         filter_counter = filter_counter + 1
!         uvbuf= zero
!         dpbuf= zero
!         ptbuf= zero
!         htbuf= zero
!         cuv_zero= zero
!         dp_zero = zero
!         pt_zero = zero   
!         cuv_rk= zero
!         dp_rk = zero
!         pt_rk = zero
!         dg3d_rhs= zero
!         !=======================================================================================================!
!         !=======================================================================================================!
!         !  Starting Time Integration										!
!         !=======================================================================================================!
!         IF (nstep >= 1) Then
!            !=======================================================================================================!
!            !  RK-3 Time Marching											!
!            !=======================================================================================================!     
!            do ie=nets, nete
!               !=======================================================================================================!
!               !  Data Initialization
!               !=======================================================================================================!  
!               if (nstep == 1 ) then 
!                  if (topology == 'cube' ) then 
!                     if (test_case=='jw_bcl') then    
!                        call jw_baroclinic(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,            &
!                             elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),          &
!                             elem(ie)%state%pt3d,elem(ie)%state%dp3d, elem(ie)%state%qt3d)
!                     elseif (test_case=='heldsuarez') then  
!                        call heldsuarez_initial(ie,elem(ie)%spherep,elem(ie)%Dinv,elem(ie)%fcor,elem(ie)%state%sgp,       &
!                             elem(ie)%state%ptop,elem(ie)%state%tbar,elem(ie)%state%v(:,:,:,:,n0),     &
!                             elem(ie)%state%pt3d,elem(ie)%state%dp3d)
!                     endif
!                     call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,		&
!                          elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,		&
!                          elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))
!                     elem(ie)%state%pr3d_ref(:,:,:)= elem(ie)%state%pr3d(:,:,:)
!                  endif
!                  !=======================================================================================================!
!                  !  DG: Compute  covariant components for the first time 
!                  !=======================================================================================================!    
!                  do k=1,nlev    
!                     elem(ie)%state%uv(:,:,:,k)  = contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
!                     elem(ie)%state%uv0(:,:,:,k) = elem(ie)%state%uv(:,:,:,k)
!                     elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))
!                     elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv)
!                  enddo
!               endif
!               !=======================================================================================================!
!               !	Filter
!               !=======================================================================================================! 
!               do k=1,nlev      
!                  call filter_P(elem(ie)%state%couv(:,:,1,k),flt)
!                  call filter_P(elem(ie)%state%couv(:,:,2,k),flt)
!                  call filter_P(elem(ie)%state%dp3d(:,:,k),flt)
!                  !if (mod(nstep,6) == 0 ) then 
!                  call filter_P(elem(ie)%state%pt3d(:,:,k),flt)
!                  !endif 
!               enddo
!               !=======================================================================================================!
!               !	Converting covarinat to contravariant (second timestep onwards)
!               !=======================================================================================================!
!               do k=1,nlev
!                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))  
!               enddo
!
!               !=============================================================
!               rhs_force(:,:,:,:,ie)= 0.0D0  
!
!               if (test_case=='jw_bcl') then    
!
!                  diffuv(:,:,:) = 0.0D0 
!
!                  do k=1,nlev   
!                     diffuv(:,:,:) =  horizontal_diff(k,deriv,elem(ie)%metdet,elem(ie)%metinv,elem(ie)%state%couv(:,:,:,k))  
!
!                     do j=1,np
!                        do i=1,np
!                           rhs_force(i,j,1,k,ie) = diffuv(i,j,1)
!                           rhs_force(i,j,2,k,ie) = diffuv(i,j,2)
!                           rhs_force(i,j,3,k,ie) = 0.0D0  
!                           rhs_force(i,j,4,k,ie) = 0.0D0
!                        end do
!                     end do
!                  end do
!               endif
!               !=======================================================================================================!
!               !	Forcing Vector for Held-Suarez
!               !=======================================================================================================!   
!               if (test_case=='heldsuarez') then    
!                  diffuv(:,:,:) = 0.0D0 
!                  hs_uv_force(:,:,:,:) = heldsuarez_uv_forcing(elem(ie)%state%pr3d,elem(ie)%state%couv)  
!
!                  do k=1,nlev   
!                     diffuv(:,:,:) =  horizontal_diff(k,deriv,elem(ie)%metdet,elem(ie)%metinv,elem(ie)%state%couv(:,:,:,k))  
!                     old_t3d(:,:,k,ie) = elem(ie)%state%T(:,:,k,n0)
!                     ! diff_pott(:,:) = elem(ie)%state%pt3d(:,:,k)
!                     ! elem(ie)%state%pt3d(:,:,k) = implicit_diff(hdt,deriv,elem(ie)%metdet,elem(ie)%metinv,diff_pott(:,:))  
!
!                     do j=1,np
!                        do i=1,np
!                           rhs_force(i,j,1,k,ie) = hs_uv_force(i,j,1,k) + diffuv(i,j,1)
!                           rhs_force(i,j,2,k,ie) = hs_uv_force(i,j,2,k) + diffuv(i,j,2)
!                           rhs_force(i,j,3,k,ie) = 0.0D0  
!                           rhs_force(i,j,4,k,ie) = 0.0D0
!                        end do
!                     end do
!                  end do
!               endif
!               !=======================================================================================================!
!               !=======================================================================================================!
!               !	Temporary Storage for fiels from previous time step
!               !=======================================================================================================!
!               do k=1,nlev      
!                  cuv_zero(:,:,:,k,ie)= elem(ie)%state%couv(:,:,:,k)
!                  dp_zero(:,:,k,ie)= elem(ie)%state%dp3d(:,:,k) * elem(ie)%metdet(:,:)
!                  pt_zero(:,:,k,ie)= elem(ie)%state%pt3d(:,:,k) * dp_zero(:,:,k,ie)
!               enddo
!               !=======================================================================================================!
!               !	DG: Packing Contravariant velocity & height fields
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
!               kptr=2*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
!               kptr=3*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
!               kptr=4*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
!               !=======================================================================================================!
!               ! 	DG: Rotating velocity (contra)
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
!               !=======================================================================================================!
!            enddo
!            !=======================================================================================================!
!            ! 	Insert communications here: for shared memory, just a single
!            ! 	sync is required
!            !=======================================================================================================!
!            call bndry_exchangeV(hybrid,edge3)
!            !=======================================================================================================!
!            !  RK-1 Stage												!
!            !=======================================================================================================!     
!            !=======================================================================================================!
!            do ie=nets,nete
!               !=======================================================================================================!
!               !	DG-RK1: Unpack the edges for uvcomp and height  (nair)
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
!               kptr=2*nlev
!               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               kptr=3*nlev
!               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               kptr=4*nlev
!               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               !=======================================================================================================! 
!               tforce(:,:,:) = 0.0D0  
!               do k=1,nlev
!                  call pres_grad_term(elem(ie),k,hdt,deriv,elem(ie)%state%peta(:,:,k),elem(ie)%state%T(:,:,k,n0),			&
!                       tforce(:,:,k),elem(ie)%state%pgrads(:,:,:,k))
!               enddo
!               !=======================================================================================================! 
!               !	DG3D Model RHS											!
!               !=======================================================================================================!
!               htop(:,:)=elem(ie)%state%ht(:,:,1)
!               do k=1,nlev
!                  call dg3d_rhs_terms(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),		&
!                       elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%couv(:,:,:,k),elem(ie)%fcor,		&
!                       elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k),		&
!                       htop,	&
!                       elem(ie)%state%pgrads(:,:,:,k),rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
!               enddo
!               !=======================================================================================================! 
!               do k=1,nlev
!                  do j=1,np
!                     do i=1,np
!                        cuv_rk(i,j,1,k,ie)= cuv_zero(i,j,1,k,ie) + hdt * dg3d_rhs(i,j,1,k,ie) 
!                        cuv_rk(i,j,2,k,ie)= cuv_zero(i,j,2,k,ie) + hdt * dg3d_rhs(i,j,2,k,ie)
!                        dp_rk(i,j,k,ie)   = dp_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,3,k,ie) 
!                        pt_rk(i,j,k,ie)   = pt_zero(i,j,k,ie) + hdt * dg3d_rhs(i,j,4,k,ie) 
!                     end do
!                  end do
!               enddo
!               !=======================================================================================================! 
!            enddo
!            !=======================================================================================================!
!            do ie= nets, nete
!               !=======================================================================================================!  
!               do k=1,nlev  
!                  do j=1,np
!                     do i=1,np
!                        elem(ie)%state%couv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
!                        elem(ie)%state%couv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
!                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
!                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
!                     end do
!                  end do
!               end do
!               !=======================================================================================================!
!               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,		&
!                    elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,			&
!                    elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0)) 
!               !=======================================================================================================!
!               !	Converting covarinat to contravariant and Psi to Height
!               !=======================================================================================================!
!               do k=1,nlev
!                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))       
!                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv) 
!               enddo
!               !=======================================================================================================!	      
!               !	DG-RK2: Packing  contravariant vectors,  ht-field 
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
!               kptr=2*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
!               kptr=3*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
!               kptr=4*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
!               !=======================================================================================================!
!               !	DG: Rotating vectors (nair)
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
!            end do
!            !=======================================================================================================!
!            ! 	Insert communications here: for shared memory, just a single
!            ! 	sync is required
!            !=======================================================================================================!
!            call bndry_exchangeV(hybrid,edge3)
!            !=======================================================================================================!
!            !  RK-2 STAGE												!
!            !=======================================================================================================! 
!            !=======================================================================================================!
!            do ie = nets, nete
!               !=======================================================================================================!
!               ! 	Unpack the edges for  vectors and scalar fields  (nair)
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
!               kptr=2*nlev
!               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               kptr=3*nlev
!               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               kptr=4*nlev
!               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               !=======================================================================================================!   
!               tforce(:,:,:) = 0.0D0  
!               do k=1,nlev
!                  call pres_grad_term(elem(ie),k,hdt,deriv,elem(ie)%state%peta(:,:,k),elem(ie)%state%T(:,:,k,n0),			&
!                       tforce(:,:,k),elem(ie)%state%pgrads(:,:,:,k))
!               enddo
!               !=======================================================================================================! 
!               !	DG3D Model RHS											!
!               !=======================================================================================================! 
!               htop(:,:)=elem(ie)%state%ht(:,:,1)
!               do k=1,nlev
!                  call dg3d_rhs_terms(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),&
!                       elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%couv(:,:,:,k),elem(ie)%fcor, &
!                       elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
!                       htop, &
!                       elem(ie)%state%pgrads(:,:,:,k),rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
!               enddo
!               !=======================================================================================================! 
!               do k=1,nlev   
!                  do j=1,np
!                     do i=1,np
!                        cuv_rk(i,j,1,k,ie)=(three*cuv_zero(i,j,1,k,ie) + cuv_rk(i,j,1,k,ie)+ hdt * dg3d_rhs(i,j,1,k,ie) )/four
!                        cuv_rk(i,j,2,k,ie)=(three*cuv_zero(i,j,2,k,ie) + cuv_rk(i,j,2,k,ie)+ hdt * dg3d_rhs(i,j,2,k,ie) )/four
!                        dp_rk(i,j,k,ie) =  (three*dp_zero(i,j,k,ie) + dp_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,3,k,ie) )/four
!                        pt_rk(i,j,k,ie) =  (three*pt_zero(i,j,k,ie) + pt_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,4,k,ie) )/four
!                     end do
!                  end do
!               enddo
!               !=======================================================================================================!
!            enddo
!            !=======================================================================================================!
!            do ie= nets, nete
!               !=======================================================================================================!  
!               do k=1,nlev  
!                  do j=1,np
!                     do i=1,np
!                        elem(ie)%state%couv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
!                        elem(ie)%state%couv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
!                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
!                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
!                     end do
!                  end do
!               end do
!               !=======================================================================================================!
!               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,		&
!                    elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,		&
!                    elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0)) 
!               !=======================================================================================================!
!               !	Converting covarinat to contravariant and Psi to Height
!               !=======================================================================================================!
!               do k=1,nlev
!                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))       
!                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv) 
!               enddo
!               !=======================================================================================================!
!               !	DG-Rk3: Packing  flux vectors,  psi-field  (nair)
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgeDGVpack(edge3,reshape(elem(ie)%state%v(:,:,:,:,n0),(/np,np,2*nlev/)),2*nlev,kptr,elem(ie)%desc)
!               kptr=2*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%ht,nlev,kptr,elem(ie)%desc)
!               kptr=3*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%dp3d,nlev,kptr,elem(ie)%desc)
!               kptr=4*nlev
!               call edgeDGVpack(edge3,elem(ie)%state%pt3d,nlev,kptr,elem(ie)%desc)
!               !=======================================================================================================!
!               !	DG: Rotating contravariant  vectors (nair)
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgerotate(edge3,2*nlev,kptr,elem(ie)%desc)
!            end do
!            !=======================================================================================================!
!            ! 	Insert communications here: for shared memory, just a single
!            ! 	sync is required
!            !=======================================================================================================!
!            call bndry_exchangeV(hybrid,edge3)
!            !=======================================================================================================!
!            do ie = nets, nete
!               !=======================================================================================================!
!               !	DG-RK3: Unpack the edges for  flux vectors and scalar fields
!               !=======================================================================================================!
!               kptr=0*nlev
!               call edgeDGVunpack(edge3, uvbuf(0,0,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
!               kptr =2*nlev
!               call edgeDGVunpack(edge3, htbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               kptr =3*nlev
!               call edgeDGVunpack(edge3, dpbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               kptr =4*nlev
!               call edgeDGVunpack(edge3, ptbuf(0,0,1,ie), nlev, kptr, elem(ie)%desc)
!               !=======================================================================================================!
!               tforce(:,:,:) = 0.0D0  
!               do k=1,nlev
!                  call pres_grad_term(elem(ie),k,hdt,deriv,elem(ie)%state%peta(:,:,k),elem(ie)%state%T(:,:,k,n0),	&
!                       tforce(:,:,k),elem(ie)%state%pgrads(:,:,:,k))
!               enddo
!               !=======================================================================================================! 
!               !	DG3D Model RHS											!
!               !=======================================================================================================! 
!               htop(:,:)=elem(ie)%state%ht(:,:,1)
!               do k=1,nlev
!                  call dg3d_rhs_terms(elem(ie),k,neqn,deriv,uvbuf(:,:,:,k,ie),htbuf(:,:,k,ie),dpbuf(:,:,k,ie),ptbuf(:,:,k,ie),	&
!                       elem(ie)%state%v(:,:,:,k,n0),elem(ie)%state%couv(:,:,:,k),elem(ie)%fcor, &
!                       elem(ie)%state%ht(:,:,k),elem(ie)%state%dp3d(:,:,k),elem(ie)%state%pt3d(:,:,k), &
!                       htop,	&
!                       elem(ie)%state%pgrads(:,:,:,k),rhs_force(:,:,:,k,ie),dg3d_rhs(:,:,:,k,ie))
!               enddo
!               !=======================================================================================================! 
!               do k=1,nlev   
!                  do j=1,np
!                     do i=1,np
!                        cuv_rk(i,j,1,k,ie)=(cuv_zero(i,j,1,k,ie) + two*(cuv_rk(i,j,1,k,ie)+ hdt * dg3d_rhs(i,j,1,k,ie)) )/three 
!                        cuv_rk(i,j,2,k,ie)=(cuv_zero(i,j,2,k,ie) + two*(cuv_rk(i,j,2,k,ie)+ hdt * dg3d_rhs(i,j,2,k,ie)) )/three 
!                        dp_rk(i,j,k,ie) = (dp_zero(i,j,k,ie) + two*(dp_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,3,k,ie)) )/three 
!                        pt_rk(i,j,k,ie) = (pt_zero(i,j,k,ie) + two*(pt_rk(i,j,k,ie)+ hdt * dg3d_rhs(i,j,4,k,ie)) )/three  
!                     end do
!                  end do
!               end do
!               !=======================================================================================================!
!            enddo
!            !=======================================================================================================!
!            do ie= nets, nete
!               !=======================================================================================================!  
!               do k=1,nlev  
!                  do j=1,np
!                     do i=1,np
!                        elem(ie)%state%couv(i,j,1,k)= cuv_rk(i,j,1,k,ie)
!                        elem(ie)%state%couv(i,j,2,k)= cuv_rk(i,j,2,k,ie)
!                        elem(ie)%state%dp3d(i,j,k)  = dp_rk(i,j,k,ie) / elem(ie)%metdet(i,j)
!                        elem(ie)%state%pt3d(i,j,k)  = pt_rk(i,j,k,ie) / dp_rk(i,j,k,ie)
!                     end do
!                  end do
!               end do
!               !=======================================================================================================!
!               !	Remapping Zone
!               !=======================================================================================================!
!               if ((nstep > 1).and.(mod(nstep,remapfreq) == 0)) then
!                  !=======================================================================================================!   
!                  if (remap_type == "linear") then
!                     call linear_remap(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,			&
!                          elem(ie)%state%couv,elem(ie)%state%dp3d)
!                  else if (remap_type == "parabolic") then
!                     call parabolic_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,	&
!                          elem(ie)%state%couv,elem(ie)%state%dp3d)
!                  else if (remap_type == "energy") then
!                     call energy_remap(elem(ie)%metinv,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,	&
!                          elem(ie)%state%couv,elem(ie)%state%dp3d)
!                  endif
!                  !=======================================================================================================! 			 
!               endif
!               !=======================================================================================================!
!               !=======================================================================================================!
!               call lagrangian_surfvars(ie,elem(ie)%state%sgp,elem(ie)%state%ptop,elem(ie)%state%pt3d,		&
!                    elem(ie)%state%dp3d,elem(ie)%state%pr3d,elem(ie)%state%gp3d,		&
!                    elem(ie)%state%psi,elem(ie)%state%peta,elem(ie)%state%T(:,:,:,n0))        
!               !=======================================================================================================!
!               !	Converting covarinat to contravariant and Psi to Height
!               !=======================================================================================================!
!               do k=1,nlev
!                  elem(ie)%state%v(:,:,:,k,n0)= co2contra(elem(ie)%state%couv(:,:,:,k),elem(ie)%metinv(:,:,:,:))       
!                  !elem(ie)%state%couv(:,:,:,k)= contra2co(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%met(:,:,:,:))  
!                  elem(ie)%state%ht(:,:,k)    = psi2height(elem(ie)%state%psi(:,:,k),grv) 
!               enddo
!               !=======================================================================================================!        
!               !	Temprature field: from 1 to nlev								!
!               !	Surface Pressure: at nlev+1									!
!               !=======================================================================================================!
!               do k=1,nlev      
!                  elem(ie)%state%uv(:,:,:,k)= contra2sphere(elem(ie)%state%v(:,:,:,k,n0),elem(ie)%D(:,:,:,:))
!                  elem(ie)%state%zeta(:,:,k)= vorticity(elem(ie)%state%couv(:,:,:,k),deriv)/(rearth*elem(ie)%metdet(:,:))
!               end do
!               elem(ie)%state%T(:,:,:,n0) = pt2temp(elem(ie)%state%pr3d,elem(ie)%state%pt3d)  
!               !Jose Garcia: For now we simply assign one field to the other until we complete the advection.
!               !elem(ie)%state%Q(:,:,:,n0) = qt2moist(elem(ie)%state%pr3d,elem(ie)%state%qt3d)  
!               elem(ie)%state%Q(:,:,:,n0) = elem(ie)%state%qt3d  
!               !=======================================================================================================!
!               !	Implicit Temperature Forcing for Held-Suarez
!               !=======================================================================================================!
!               if (test_case=='heldsuarez') then 
!                  elem(ie)%state%pt3d(:,:,:)= heldsuarez_th_correction(hdt,elem(ie)%state%pt3d,elem(ie)%spherep, &
!                       elem(ie)%state%pr3d,elem(ie)%state%T(:,:,:,n0))  
!                  elem(ie)%state%T(:,:,:,n0) = pt2temp(elem(ie)%state%pr3d,elem(ie)%state%pt3d)  
!               endif
!              !=======================================================================================================!
!            enddo
!            !=======================================================================================================!
!            !   Ending Time Step nstep >= 1										!
!            !=======================================================================================================!
!         ENDIF
!         !=======================================================================================================!
!      endif
!      !=======================================================================================================!
!      !=======================================================================================================!
!    end subroutine dg3d_rk3
     !=======================================================================================================!
     !=======================================================================================================!
   end module dg3d_prim_mod
