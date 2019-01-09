module fvm_consistent_se_cslam
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: nc, nhe, nlev, ntrac, np, nhr, nhc, ngpc, ns, nht, lbc,ubc
  use dimensions_mod,         only: irecons_tracer
  use dimensions_mod,         only: kmin_jet,kmax_jet
  use cam_abortutils,         only: endrun

  use time_mod,               only: timelevel_t
  use element_mod,            only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod,             only: hybrid_t
  use perf_mod,               only: t_startf, t_stopf 

  implicit none
  private
  save

  real (kind=r8), dimension(ngpc), private :: gsweights, gspts
  real (kind=r8),parameter       , private :: eps=1.0e-14_r8
  public :: run_consistent_se_cslam
contains
  !
  !**************************************************************************************
  !
  ! Consistent CSLAM-SE algorithm documented in
  !
  ! Lauritzen et al. (2017): CAM-SE-CSLAM: Consistent finite-volume transport with
  !                          spectral-element dynamics. Mon. Wea. Rev.
  !
  !
  !**************************************************************************************
  !
  subroutine run_consistent_se_cslam(elem,fvm,hybrid,dt_fvm,tl,nets,nete,hvcoord)
    ! ---------------------------------------------------------------------------------
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm
    use fvm_mod               , only:  fill_halo_fvm, ghostBufQnhc, ghostBufQ1, ghostBufFlux
    use fvm_reconstruction_mod, only: reconstruction
    use fvm_analytic_mod      , only: gauss_points
    use derivative_mod        , only: subcell_integration
    use edge_mod              , only: ghostpack, ghostunpack
    use bndry_mod             , only: ghost_exchange
    use hybvcoord_mod         , only: hvcoord_t
    use constituents          , only: qmin
    implicit none
    type (element_t)      , intent(inout) :: elem(:)
    type (fvm_struct)     , intent(inout) :: fvm(:)
    type (hybrid_t)       , intent(in)    :: hybrid   ! distributed parallel structure (shared)
    type (TimeLevel_t)    , intent(in)    :: tl              ! time level struct
    type (hvcoord_t)      , intent(in)    :: hvcoord
    integer               , intent(in)    :: nets  ! starting thread element number (private)
    integer               , intent(in)    :: nete  ! ending thread element number   (private)
    real (kind=r8)        , intent(in)    :: dt_fvm


    !high-order air density reconstruction
    real (kind=r8) :: ctracer(irecons_tracer,1-nhe:nc+nhe,1-nhe:nc+nhe,ntrac)
    real (kind=r8) :: inv_dp_area(nc,nc), ps_se(nc,nc)

    logical :: llimiter(ntrac)
    integer :: i,j,k,ie,itr,ntmp

    llimiter = .true.

    call gauss_points(ngpc,gsweights,gspts) !set gauss points/weights
    gspts = 0.5_r8*(gspts+1.0_r8) !shift location so in [0:1] instead of [-1:1]

    call t_startf('fvm:before_Qnhc')
    do ie=nets,nete
       do k=1,nlev
          elem(ie)%sub_elem_mass_flux(:,:,:,k) = dt_fvm*elem(ie)%sub_elem_mass_flux(:,:,:,k)*fvm(ie)%dp_ref_inverse(k)
          fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)   =         fvm(ie)%dp_fvm (1:nc,1:nc,k,n0_fvm)*fvm(ie)%dp_ref_inverse(k)
       end do
       call ghostpack(ghostbufQnhc,fvm(ie)%dp_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,1:nlev,n0_fvm),nlev,      0,ie)
       call ghostpack(ghostbufQnhc,fvm(ie)%c(1-nhc:nc+nhc,1-nhc:nc+nhc,1:nlev,1:ntrac,n0_fvm)   ,nlev*ntrac,nlev,ie)
    end do
    call t_stopf('fvm:before_Qnhc')

    call t_startf('fvm:ghost_exchange:Qnhc')
    call ghost_exchange(hybrid,ghostbufQnhc)
    call t_stopf('fvm:ghost_exchange:Qnhc')

    call t_startf('fvm:orthogonal_swept_areas')
    do ie=nets,nete
      fvm(ie)%se_flux    (1:nc,1:nc,:,:) = elem(ie)%sub_elem_mass_flux(:,:,:,:)
      call ghostunpack(ghostbufQnhc, fvm(ie)%dp_fvm(1-nhc:nc+nhc,1-nhc:nc+nhc,1:nlev,n0_fvm),nlev     ,0,ie)
      call ghostunpack(ghostbufQnhc, fvm(ie)%c(1-nhc:nc+nhc,1-nhc:nc+nhc,1:nlev,1:ntrac,n0_fvm),   nlev*ntrac,nlev,ie) 
      call compute_displacements_for_swept_areas (fvm(ie),fvm(ie)%dp_fvm(1-nhe:nc+nhe,1-nhe:nc+nhe,:,n0_fvm),1)
      call ghostpack(ghostBufFlux, fvm(ie)%se_flux(:,:,:,:),4*nlev,0,ie)
    end do
    call ghost_exchange(hybrid,ghostBufFlux)
    do ie=nets,nete
      call ghostunpack(ghostBufFlux, fvm(ie)%se_flux(:,:,:,:),4*nlev,0,ie)
      call ghost_flux_unpack(fvm(ie))
    enddo

    call t_stopf('fvm:orthogonal_swept_areas')

    do ie=nets,nete
      fvm(ie)%c(:,:,:,:,np1_fvm)    = 0.0_r8!to avoid problems when uninitialized variables are set to NaN
      fvm(ie)%dp_fvm(:,:,:,np1_fvm) = 0.0_r8!to avoid problems when uninitialized variables are set to NaN      
      do k=1,nlev
        call t_startf('fvm:tracers_reconstruct')
        call reconstruction(fvm(ie)%c(1-nhc:nc+nhc,1-nhc:nc+nhc,k,1:ntrac,n0_fvm),&
             ctracer(:,:,:,:),irecons_tracer,llimiter,ntrac,&
             nc,nhe,nhr,nhc,nht,ns,nhr+(nhe-1),&
             fvm(ie)%jx_min,fvm(ie)%jx_max,fvm(ie)%jy_min,fvm(ie)%jy_max,&
             fvm(ie)%cubeboundary,fvm(ie)%halo_interp_weight,fvm(ie)%ibase,&
             fvm(ie)%spherecentroid(:,1-nhe:nc+nhe,1-nhe:nc+nhe),&
             fvm(ie)%recons_metrics,fvm(ie)%recons_metrics_integral,&
             fvm(ie)%rot_matrix,fvm(ie)%centroid_stretch,&
             fvm(ie)%vertex_recons_weights,fvm(ie)%vtx_cart&
             )
        call t_stopf('fvm:tracers_reconstruct')
        call t_startf('fvm:swept_flux')
        call swept_flux(elem(ie),fvm(ie),k,ctracer)
        call t_stopf('fvm:swept_flux')
      end do
    end do
     !
     !***************************************
     !
     ! Large Courant number increment
     !
     !***************************************
     !
     ! In the jet region the effective Courant number
     ! in the cslam trajectory algorithm can be > 1
     ! (by up to 20%)
     !
     ! We limit the trajectories to < 1 but in this step
     ! we do a piecewise constant update for the
     ! amount of mass for which the Courant number is >1
     !
     !
     call t_startf('fvm:fill_halo_fvm:large_Courant')
     call fill_halo_fvm(ghostbufQ1,elem,fvm,hybrid,nets,nete,np1_fvm,1,kmin_jet,kmax_jet)
     call t_stopf('fvm:fill_halo_fvm:large_Courant')
     call t_startf('fvm:large_Courant_number_increment')
     do ie=nets,nete
       do k=kmin_jet,kmax_jet !1,nlev
          call large_courant_number_increment(fvm(ie),k)
        end do
      end do
     call t_stopf('fvm:large_Courant_number_increment')

     call t_startf('fvm:end_of_reconstruct_subroutine')
     do ie=nets,nete
       !
       ! convert to mixing ratio
       !
       do k=1,nlev         
         do j=1,nc
           do i=1,nc
             inv_dp_area(i,j) = 1.0_r8/fvm(ie)%dp_fvm(i,j,k,np1_fvm)
           end do
         end do

         do itr=1,ntrac
           do j=1,nc
             do i=1,nc
               ! convert to mixing ratio
               fvm(ie)%c(i,j,k,itr,np1_fvm) = fvm(ie)%c(i,j,k,itr,np1_fvm)*inv_dp_area(i,j)
               ! remove round-off undershoots
               fvm(ie)%c(i,j,k,itr,np1_fvm) = MAX(fvm(ie)%c(i,j,k,itr,np1_fvm),qmin(itr))
             end do
           end do
         end do
         !
         ! convert to dp and scale back dp
         !
         fvm(ie)%dp_fvm(1:nc,1:nc,k,np1_fvm) = fvm(ie)%dp_fvm(1:nc,1:nc,k,np1_fvm)*fvm(ie)%dp_ref(k)*fvm(ie)%inv_area_sphere
       end do
       !
       ! to avoid accumulation of truncation error overwrite CSLAM surface pressure with SE
       ! surface pressure
       !
       call subcell_integration(elem(ie)%state%psdry(:,:), np, nc, elem(ie)%metdet,ps_se)
       fvm(ie)%psc = ps_se*fvm(ie)%inv_se_area_sphere
!       do j=1,nc
!         do i=1,nc
!           fvm(ie)%psc(i,j) = sum(fvm(ie)%dp_fvm(i,j,:,np1_fvm)) +  hvcoord%hyai(1)*hvcoord%ps0
!         end do
!       end do
     end do
     call t_stopf('fvm:end_of_reconstruct_subroutine')
     !
     ! advance fvm time-levels
     !
     ntmp     = np1_fvm
     np1_fvm  = n0_fvm
     n0_fvm   = ntmp
  end subroutine run_consistent_se_cslam

  subroutine swept_flux(elem,fvm,ilev,ctracer)
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm
    use fvm_analytic_mod      , only: get_high_order_weights_over_areas
    use dimensions_mod, only : kmin_jet,kmax_jet
    implicit none
    type (element_t) , intent(in)   :: elem
    type (fvm_struct), intent(inout):: fvm
    integer          , intent(in) :: ilev
    real (kind=r8), intent(inout) :: ctracer(irecons_tracer,1-nhe:nc+nhe,1-nhe:nc+nhe,ntrac)

    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=r8)    , dimension(0:7       , imin:imax,imin:imax,num_sides) :: displ
    integer (kind=r8) , dimension(1:2,11    , imin:imax,imin:imax,num_sides) :: base_vec
    real (kind=r8)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides) :: base_vtx
    integer                  , dimension(2,num_area, imin:imax,imin:imax,num_sides) :: idx
    real (kind=r8)    , dimension(imin:imax,imin:imax,num_sides)             :: mass_flux_se
    real (kind=r8)    , dimension(irecons_tracer,num_area) :: weights
    real (kind=r8)                     :: gamma
    integer :: i,j,iside,iarea,iw

    integer, parameter :: num_seg_max=5
    REAL(KIND=r8), dimension(2,num_seg_max,num_area) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)               :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8) :: x_start, dgam_vec
    REAL(KIND=r8) :: gamma_max, displ_first_guess

    REAL(KIND=r8) :: flux,flux_tracer(ntrac)

    REAL(KIND=r8), dimension(num_area) :: dp_area

    logical :: tl1,tl2,tr1,tr2

    integer, dimension(4), parameter :: imin_side = (/1   ,0   ,1   ,1   /)
    integer, dimension(4), parameter :: imax_side = (/nc  ,nc  ,nc  ,nc+1/)
    integer, dimension(4), parameter :: jmin_side = (/1   ,1   ,0   ,1   /)
    integer, dimension(4), parameter :: jmax_side = (/nc+1,nc  ,nc  ,nc  /)

    integer :: iseg, iseg_tmp,flowcase,ii,jj,itr

    call define_swept_areas(fvm,ilev,displ,base_vec,base_vtx,idx)

    mass_flux_se(1:nc,1:nc,1:4)  = -elem%sub_elem_mass_flux(1:nc,1:nc,1:4,ilev)
    mass_flux_se(0   ,1:nc,2  )  =  elem%sub_elem_mass_flux(1   ,1:nc,4  ,ilev)
    mass_flux_se(nc+1,1:nc,4  )  =  elem%sub_elem_mass_flux(nc  ,1:nc,2  ,ilev)
    mass_flux_se(1:nc,0   ,3  )  =  elem%sub_elem_mass_flux(1:nc,1   ,1  ,ilev)
    mass_flux_se(1:nc,nc+1,1  )  =  elem%sub_elem_mass_flux(1:nc,nc  ,3  ,ilev)
    !
    ! prepare for air/tracer update
    !
    fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm) = fvm%dp_fvm(1:nc,1:nc,ilev,n0_fvm)*fvm%area_sphere
    do itr=1,ntrac
      fvm%c(1:nc,1:nc,ilev,itr,np1_fvm) = fvm%c(1:nc,1:nc,ilev,itr,n0_fvm)*fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm)
      do iw=1,irecons_tracer
        ctracer(iw,1-nhe:nc+nhe,1-nhe:nc+nhe,itr)=ctracer(iw,1-nhe:nc+nhe,1-nhe:nc+nhe,itr)*&
             fvm%dp_fvm(1-nhe:nc+nhe,1-nhe:nc+nhe,ilev,n0_fvm)
      end do
    end do

    do iside=1,4
      do j=jmin_side(iside),jmax_side(iside)
        do i=imin_side(iside),imax_side(iside)
           !DO NOT USE MASS_FLUX_SE AS THRESHOLD - THRESHOLD CONDITION MUST BE CONSISTENT WITH 
           !THE ONE USED IN DEFINE_SWEPT_AREAS
!          if (mass_flux_se(i,j,iside)>eps) then 
          if (fvm%se_flux(i,j,iside,ilev)>eps) then
            !
            !        ||             ||
            !  tl1   ||             || tr1
            !        ||             ||
            !  =============================
            !        ||             ||
            !  tl2   ||             || tr2
            !        ||             ||
            !
            tl1 = displ(3,i,j,iside)<0.0_r8.and.displ(6,i,j,iside).ge.0.0_r8 !departure point in tl1 quadrant
            tl2 = displ(6,i,j,iside)<0.0_r8.and.displ(7,i,j,iside)   >0.0_r8 !departure point in tl2 quadrant
            tr1 = displ(2,i,j,iside)<0.0_r8.and.displ(4,i,j,iside).ge.0.0_r8 !departure point in tr1 quadrant
            tr2 = displ(4,i,j,iside)<0.0_r8.and.displ(5,i,j,iside)   >0.0_r8 !departure point in tr2 quadrant

            !
            ! pathological cases
            !
            !        |  ||           ||                      ||           ||
            !        |  ||-----------||                      ||-----------||
            !        |  ||           ||                      ||           ||
            !  ================================     =================================
            !           ||           ||                   |  ||           ||
            !  ---------||           ||             ------|--||           ||
            !           ||           ||                   |  ||           ||
            !
            !                tl1=tl1.or.tl2
            !                tr1=tr1.or.tr2
            !                tl1=displ(3,i,j,iside)<0.0_r8.and..not.(tl1.and.tl2)
            !                tr1=displ(2,i,j,iside)<0.0_r8.and..not.(tr1.and.tr2)

            num_seg=-1; num_seg_static=-1 !initialization
            if (.not.tl1.and..not.tl2.and..not.tr1.and..not.tr2) then
              flowcase=0
              !
              !        ||             ||                 ||             ||                ||             ||
              !        ||  *       *  ||                 ||  *----------*                 |*----------*  ||
              !        || /         \ ||                 || /           ||                ||           \ ||
              !        ||/           \||                 ||/            ||                ||            \||
              !  =============================     =============================     =============================
              !        ||             ||                 ||             ||                ||             ||
              !
              !
              call define_area3_center (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                   num_seg_static,x_start, dgam_vec,fvm%se_flux(i,j,iside,ilev),displ_first_guess)

              gamma=1.0_r8!fvm%se_flux(i,j,iside,ilev)
              gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
            else
              if (tl1.and.tr1) then
                flowcase=1
                !
                !
                !  tl1   ||             || tr1             ||             ||                ||             ||
                !     *--||-------------||--*           *--||-------------||                ||-------------||--*
                !      \ ||             || /             \ ||             ||\              /||             || /
                !       \||             ||/               \||             || \            / ||             ||/
                !  =============================     =========================*===     ==*==========================
                !        ||             ||                 ||             ||                ||             ||
                !
                call define_area2           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static,&
                     num_seg, num_seg_static,x_start, dgam_vec,displ_first_guess)
                call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static,&
                     num_seg, num_seg_static,x_start, dgam_vec)
                call define_area4           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static,&
                     num_seg, num_seg_static,x_start, dgam_vec)
                gamma=1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tl1.and..not.tr1.and..not.tr2) then
                flowcase=2
                !
                !        ||             ||                 ||             ||                ||             ||
                !     *--||----------*  ||                /||----------*  ||             *--||-------------*
                !      \ ||           \ ||               / ||           \ ||              \ ||             ||
                !       \||            \||              /  ||            \||               \||             ||
                !  =============================     ==*==========================     =============================
                !        ||             ||                 ||             ||                ||             ||
                !
                call define_area2     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
                     x_start, dgam_vec,displ_first_guess)
                call define_area3_left(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
                     x_start, dgam_vec)
                gamma=1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tr1.and..not.tl1.and..not.tl2) then !displ(3).ge.0.0_r8) then
                flowcase=3
                !
                !        ||  *----------||--*              ||  *----------||\                *-------------||--*
                !        || /           || /               || /           || \              ||             || /
                !        ||/            ||/                ||/            ||  \             ||             ||/
                !  =============================     ==========================*==     =============================
                !        ||             ||                 ||             ||                ||             ||
                !        ||             ||                 ||             ||                ||             ||
                !        ||             ||                 ||             ||                ||             ||
                !
                call define_area3_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
                     num_seg_static, x_start, dgam_vec)
                call define_area4      (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
                     num_seg_static, x_start, dgam_vec,displ_first_guess)
                gamma=1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tl2.and..not.tr1.and..not.tr2) then !displ(2).ge.0.0_r8) then
                flowcase=4
                !
                !        ||----------*  ||                 ||-------------*
                !       /||           \ ||                /||             ||
                !      / ||            \||               / ||             ||
                !  ===/=========================     ===/=========================
                !     | /||             ||              | /||             ||
                !     |/ ||             ||              |/ ||             ||
                !     *  ||             ||              *  ||             ||
                !
                call define_area1_area2(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area3_left (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,&
                     x_start, dgam_vec,displ_first_guess)
                gamma = 1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tr2.and..not.tl1.and..not.tl2) then !displ(3).ge.0.0_r8) then
                flowcase=5
                !                case(5)
                !
                !
                !        ||  *-----2----||
                !        || /1         3||\
                !        ||/      4     || \
                !  =============================
                !        ||             ||\ |
                !        ||             || \|
                !        ||             ||  *
                !
                call define_area3_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area4_area5(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec,displ_first_guess)
                gamma=1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tl2.and.tr1.and..not.tr2) then
                flowcase=6
                !                case(6)
                !
                !
                !        ||-------------||--*
                !       /||             || /
                !      / ||             ||/
                !  ===/=========================
                !     | /||             ||
                !     |/ ||             ||
                !     *  ||             ||
                !
                !
                call define_area1_area2     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area4           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec,displ_first_guess)

                gamma=1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tr2.and.tl1.and..not.tl2) then
                flowcase=7
                !                case(7)
                !
                !
                !     *--||-------------||
                !      \ ||             ||\
                !       \||             || \
                !  =============================
                !        ||             ||\ |
                !        ||             || \|
                !        ||             ||  *
                !
                !
                call define_area2           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec,displ_first_guess)
                call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area4_area5     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                gamma =  1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else if (tl2.and.tr2) then
                flowcase=8
                !                case(8)
                !
                !
                !        ||-------------||
                !       /||             ||\
                !      / ||             || \
                !  =============================
                !     | /||             ||\ |
                !     |/ ||             || \|
                !     *  ||             ||  *
                !
                !
                !
                !
                !
                call define_area1_area2     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec)
                call define_area4_area5     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                     num_seg_static,x_start, dgam_vec,displ_first_guess)
                gamma =  1.0_r8
                gamma_max = fvm%displ_max(i,j,iside)/displ_first_guess
              else
                call endrun('ERROR - unknown flow case')
              end if
            end if
            !
            ! iterate to get flux area
            !
            call t_startf('fvm:swept_area:get_gamma')
            do iarea=1,num_area
              dp_area(iarea) = fvm%dp_fvm(idx(1,iarea,i,j,iside),idx(2,iarea,i,j,iside),ilev,n0_fvm)
            end do
            call get_flux_segments_area_iterate(x,x_static,dx_static,dx,x_start,dgam_vec,num_seg,num_seg_static,&
                 num_seg_max,num_area,dp_area,flowcase,gamma,mass_flux_se(i,j,iside),0.0_r8,gamma_max)
            call t_stopf('fvm:swept_area:get_gamma')
            !
            ! pack segments for high-order weights computation
            !
            do iarea=1,num_area
              do iseg=1,num_seg_static(iarea)
                iseg_tmp=num_seg(iarea)+iseg
                x (:,iseg_tmp,iarea)  = x_static (:,iseg,iarea)
                dx(:,iseg_tmp,iarea)  = dx_static(:,iseg,iarea)
              end do
              num_seg(iarea)=num_seg(iarea)+MAX(0,num_seg_static(iarea))
            end do
            !
            ! compute higher-order weights
            !
            call t_startf('fvm:swept_area:get_high_order_w')
            call get_high_order_weights_over_areas(x,dx,num_seg,num_seg_max,num_area,weights,ngpc,gsweights, gspts,irecons_tracer)
            call t_stopf('fvm:swept_area:get_high_order_w')
            !
            !**************************************************
            !
            ! remap air and tracers
            !
            !**************************************************
            !
            call t_startf('fvm:swept_area:remap')
            flux=0.0_r8; flux_tracer=0.0_r8
            do iarea=1,num_area
              if (num_seg(iarea)>0) then
                ii=idx(1,iarea,i,j,iside); jj=idx(2,iarea,i,j,iside)
                flux=flux+weights(1,iarea)*fvm%dp_fvm(ii,jj,ilev,n0_fvm)
                do itr=1,ntrac
                  do iw=1,irecons_tracer
                    flux_tracer(itr) = flux_tracer(itr)+weights(iw,iarea)*ctracer(iw,ii,jj,itr)
                  end do
                end do
              end if
            end do
            fvm%se_flux(i,j,iside,ilev) = mass_flux_se(i,j,iside)-flux
            if (fvm%se_flux(i,j,iside,ilev)>1.0E-13_r8.and.(ilev<kmin_jet.or.ilev>kmax_jet)) then
              write(*,*) "CN excess flux outside of pre-scribed jet region"
              write(*,*) "Increase jet region with kmin_jet and kmax_jet ",&
                   ilev,fvm%se_flux(i,j,iside,ilev),mass_flux_se(i,j,iside),flux,flowcase,&
                   kmin_jet,kmax_jet
            end if

            fvm%dp_fvm(i  ,j  ,ilev        ,np1_fvm) = fvm%dp_fvm(i  ,j  ,ilev        ,np1_fvm)-flux
            fvm%     c(i  ,j  ,ilev,1:ntrac,np1_fvm) = fvm%     c(i  ,j  ,ilev,1:ntrac,np1_fvm)-flux_tracer(1:ntrac)
            !
            ! update flux in nearest neighbor cells
            !
            if (iside==1) then
              fvm%dp_fvm(i,j-1,ilev        ,np1_fvm) = fvm%dp_fvm(i,j-1,ilev        ,np1_fvm)+flux
              fvm%     c(i,j-1,ilev,1:ntrac,np1_fvm) = fvm%     c(i,j-1,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            if (iside==2) then
              fvm%dp_fvm(i+1,j,ilev        ,np1_fvm) = fvm%dp_fvm(i+1,j,ilev        ,np1_fvm)+flux
              fvm%     c(i+1,j,ilev,1:ntrac,np1_fvm) = fvm%     c(i+1,j,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            if (iside==3) then
              fvm%dp_fvm(i,j+1,ilev        ,np1_fvm) = fvm%dp_fvm(i,j+1,ilev        ,np1_fvm)+flux
              fvm%     c(i,j+1,ilev,1:ntrac,np1_fvm) = fvm%     c(i,j+1,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            if (iside==4) then
              fvm%dp_fvm(i-1,j,ilev        ,np1_fvm) = fvm%dp_fvm(i-1,j,ilev        ,np1_fvm)+flux
              fvm%     c(i-1,j,ilev,1:ntrac,np1_fvm) = fvm%     c(i-1,j,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            call t_stopf('fvm:swept_area:remap')
          end if
        end do
      end do
    end do    
  end subroutine swept_flux


  subroutine large_courant_number_increment(fvm,ilev)
    use fvm_control_volume_mod, only: np1_fvm
    implicit none
    type (fvm_struct), intent(inout):: fvm
    integer          , intent(in) :: ilev

    integer, parameter :: num_sides=4, imin= 0, imax=nc+1

    integer, dimension(4), parameter :: imin_side = (/1   ,0   ,1   ,1   /)
    integer, dimension(4), parameter :: imax_side = (/nc  ,nc  ,nc  ,nc+1/)
    integer, dimension(4), parameter :: jmin_side = (/1   ,1   ,0   ,1   /)
    integer, dimension(4), parameter :: jmax_side = (/nc+1,nc  ,nc  ,nc  /)

    integer :: i,j,iside,itr
    real (kind=r8)    :: flux,flux_tracer(ntrac)
    real (kind=r8), dimension(0:nc+1,0:nc+1)      :: inv_dp_area
    real (kind=r8), dimension(0:nc+1,0:nc+1,ntrac):: c_tmp

    inv_dp_area=1.0_r8/fvm%dp_fvm(0:nc+1,0:nc+1,ilev,np1_fvm)
    c_tmp      = fvm%c(0:nc+1,0:nc+1,ilev,1:ntrac,np1_fvm)
    do iside=1,4
      do j=jmin_side(iside),jmax_side(iside)
        do i=imin_side(iside),imax_side(iside)
          if (fvm%se_flux(i,j,iside,ilev)>eps) then
            flux = fvm%se_flux(i,j,iside,ilev)
            do itr=1,ntrac
              flux_tracer(itr) = fvm%se_flux(i,j,iside,ilev)*c_tmp(i,j,itr)*inv_dp_area(i,j)
            end do
            fvm%dp_fvm(i  ,j  ,ilev        ,np1_fvm) = fvm%dp_fvm(i  ,j  ,ilev        ,np1_fvm)-flux
            fvm%     c(i  ,j  ,ilev,1:ntrac,np1_fvm) = fvm%     c(i  ,j  ,ilev,1:ntrac,np1_fvm)-flux_tracer(1:ntrac)
            !
            ! update flux in nearest neighbor cells
            !
            if (iside==1) then
              fvm%dp_fvm(i,j-1,ilev        ,np1_fvm) = fvm%dp_fvm(i,j-1,ilev        ,np1_fvm)+flux
              fvm%     c(i,j-1,ilev,1:ntrac,np1_fvm) = fvm%     c(i,j-1,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            if (iside==2) then
              fvm%dp_fvm(i+1,j,ilev        ,np1_fvm) = fvm%dp_fvm(i+1,j,ilev        ,np1_fvm)+flux
              fvm%     c(i+1,j,ilev,1:ntrac,np1_fvm) = fvm%     c(i+1,j,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            if (iside==3) then
              fvm%dp_fvm(i,j+1,ilev        ,np1_fvm) = fvm%dp_fvm(i,j+1,ilev        ,np1_fvm)+flux
              fvm%     c(i,j+1,ilev,1:ntrac,np1_fvm) = fvm%     c(i,j+1,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
            if (iside==4) then
              fvm%dp_fvm(i-1,j,ilev        ,np1_fvm) = fvm%dp_fvm(i-1,j,ilev        ,np1_fvm)+flux
              fvm%     c(i-1,j,ilev,1:ntrac,np1_fvm) = fvm%     c(i-1,j,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
            end if
          end if
        end do
      end do
    end do
  end subroutine large_courant_number_increment

  subroutine ghost_flux_unpack(fvm)
    use control_mod, only : neast, nwest, seast, swest
    implicit none
    type (fvm_struct), intent(inout) :: fvm

    integer :: i,j,k,ishft
    !
    ! rotate coordinates if needed
    !
    if (fvm%cubeboundary.NE.0) then
       do k=1,nlev
          do j=1-nhe,nc+nhe
             do i=1-nhe,nc+nhe
                ishft = NINT(fvm%flux_orient(2,i,j))
                fvm%se_flux(i,j,1:4,k) = cshift(fvm%se_flux(i,j,1:4,k),shift=ishft)
             end do
          end do
       end do
       !
       ! non-existent cells in physical space - necessary?
       !
       if (fvm%cubeboundary==nwest) then
          fvm%se_flux(1-nhe:0,nc+1 :nc+nhe,:,:) = 0.0_r8
       else if (fvm%cubeboundary==swest) then
          fvm%se_flux(1-nhe:0,1-nhe:0     ,:,:) = 0.0_r8
       else if (fvm%cubeboundary==neast) then
          fvm%se_flux(nc+1 :nc+nhe,nc+1 :nc+nhe,:,:) = 0.0_r8
       else if (fvm%cubeboundary==seast) then
          fvm%se_flux(nc+1 :nc+nhe,1-nhe:0,:,:) = 0.0_r8
       end if
    end if
  end subroutine ghost_flux_unpack

  subroutine compute_displacements_for_swept_areas(fvm,cair,irecons)
    implicit none
    type (fvm_struct), intent(inout)     :: fvm
    integer, intent(in) :: irecons
    real (kind=r8)                :: cair(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons,nlev) !high-order air density reconstruction
    !
    !   flux iside 1                     flux iside 3                    flux iside 2       flux iside 4
    !
    !   |          |                     |  ---1--> |                    |    --2-->|       |--1-->    |
    !  -4----------3-   /\              -4----------3-                  -4----------3-     -4----------3-   ||
    !   |          |   /||\              |\\\\\\\\\\|    ||              |   |\\\\\\|       |\\\\\\|   |
    !   |  --2-->  |    || dv(1)         |\\\\\\\\\\|    ||              |   |\\\\\\|       |\\\\\\|   |
    !   |----------|    ||               |----------|    || dv(3)        |   |\\\\\\|       |\\\\\\|   |
    !   |\\\\\\\\\\|    ||               | <--2---  |   \||/             |   |\\\\\\|       |\\\\\\|   |
    !   |\\\\\\\\\\|    ||               |          |    \/              |   |\\\\\\|       |\\\\\\|   |
    !  -1----------2-                   -1----------2-                  -1----------2-     -1----------2-
    !   |  <--1--  |                     |          |                    |    <--1--|       |<--2--
    !
    !                                                                     /                          \
    !   line-integral                                                    <==========         =========>
    !   from vertex 2                                                     \  dv(2)              dv(4)/
    !   to 1
    !
    !   Note vertical
    !   lines have
    !   zero line-
    !   integral!
    !
    integer               :: i,j,k,iside,ix
    integer, parameter :: num_area=1, num_seg_max=2
    REAL(KIND=r8), dimension(2,num_seg_max,num_area,4,nc,nc) :: x_static, dx_static
    REAL(KIND=r8), dimension(2,num_seg_max,num_area,4,nc,nc) :: x, dx
    REAL(KIND=r8), dimension(2,num_seg_max,num_area)         :: x_tmp, dx_tmp
    integer             , dimension(              num_area,4      ) :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8,                   4,nc,nc) :: x_start, dgam_vec
    REAL(KIND=r8), dimension(num_area) :: dp_area
    integer, dimension(4) :: flowcase
    REAL(KIND=r8)  :: gamma(4), flux_se

    num_seg_static(1,1) =  1; num_seg(1,1) = 1; flowcase(1) = -1
    num_seg_static(1,2) =  0; num_seg(1,2) = 2; flowcase(2) = -2
    num_seg_static(1,3) =  1; num_seg(1,3) = 1; flowcase(3) = -1
    num_seg_static(1,4) =  0; num_seg(1,4) = 2; flowcase(4) = -4

    do j=1,nc
       do i=1,nc
          do ix=1,2
             iside=1;
             x_static (ix,1,1,iside,i,j) = fvm%vtx_cart(2,ix,i,j)
             dx_static(ix,1,1,iside,i,j) = fvm%vtx_cart(1,ix,i,j)-fvm%vtx_cart(2,ix,i,j)
             x_start  (ix,1,  iside,i,j) = fvm%vtx_cart(1,ix,i,j)
             x_start  (ix,2,  iside,i,j) = fvm%vtx_cart(2,ix,i,j)
             dgam_vec (ix,1,  iside,i,j) = fvm%vtx_cart(4,ix,i,j)-fvm%vtx_cart(1,ix,i,j)
             !
             ! compute first guess
             !
             gamma(iside)                       = 0.5_r8
             x        (ix,1,1,iside,i,j) = x_start(ix,1,iside,i,j)+gamma(iside)*dgam_vec(ix,1,iside,i,j)
             dx       (ix,1,1,iside,i,j) = -dx_static(ix,1,1,iside,i,j)
             !
             ! side 2
             !
             iside=2;
             x_start  (ix,1,  iside,i,j) = fvm%vtx_cart(2,ix,i,j)
             x_start  (ix,2,  iside,i,j) = fvm%vtx_cart(3,ix,i,j)
             dgam_vec (ix,1,  iside,i,j) = fvm%vtx_cart(1,ix,i,j)-fvm%vtx_cart(2,ix,i,j)
             x        (ix,1,1,iside,i,j) = x_start(ix,1,iside,i,j)
             !
             ! compute first guess - gamma=1
             !
             gamma(iside)                       = 0.5_r8
             dx       (ix,1,1,iside,i,j) =  gamma(iside)*dgam_vec (ix,1,  iside,i,j)
             x        (ix,2,1,iside,i,j) =  x_start(ix,2,iside,i,j)+gamma(iside)*dgam_vec(ix,1,iside,i,j)
             dx       (ix,2,1,iside,i,j) = -gamma(iside)*dgam_vec (ix,1,  iside,i,j)
             !
             ! side 3
             !
             iside=3;
             x_static (ix,1,1,iside,i,j) = fvm%vtx_cart(4,ix,i,j)
             dx_static(ix,1,1,iside,i,j) = fvm%vtx_cart(3,ix,i,j)-fvm%vtx_cart(4,ix,i,j)
             x_start  (ix,1,  iside,i,j) = fvm%vtx_cart(3,ix,i,j)
             x_start  (ix,2,  iside,i,j) = fvm%vtx_cart(4,ix,i,j)
             dgam_vec (ix,1,  iside,i,j) = fvm%vtx_cart(2,ix,i,j)-fvm%vtx_cart(3,ix,i,j)
             !
             ! compute first guess - gamma(iside)=1
             !
             gamma(iside)                       = 0.5_r8
             x        (ix,1,1,iside,i,j) = x_start(ix,1,iside,i,j)+gamma(iside)*dgam_vec(ix,1,iside,i,j)
             dx       (ix,1,1,iside,i,j) = -dx_static(ix,1,1,iside,i,j)
             !
             ! side 4
             !
             iside=4;
             x_start  (ix,1,  iside,i,j) = fvm%vtx_cart(1,ix,i,j)
             x_start  (ix,2,  iside,i,j) = fvm%vtx_cart(4,ix,i,j)
             dgam_vec (ix,1,  iside,i,j) = fvm%vtx_cart(2,ix,i,j)-fvm%vtx_cart(1,ix,i,j)
             x        (ix,2,1,iside,i,j) = x_start(ix,2,iside,i,j)
             !
             ! compute first guess - gamma(iside)=1
             !
             gamma(iside)                       = 0.5_r8
             dx       (ix,2,1,iside,i,j) =  gamma(iside)*dgam_vec (ix,1,  iside,i,j)
             x        (ix,1,1,iside,i,j) =  x_start(ix,1,iside,i,j)+gamma(iside)*dgam_vec(ix,1,iside,i,j)
             dx       (ix,1,1,iside,i,j) = -gamma(iside)*dgam_vec (ix,1,  iside,i,j)
          end do
       end do
    end do

    do k=1,nlev
      do j=1,nc
        do i=1,nc
          dp_area = cair(i,j,1,k)
          do iside=1,4
            flux_se = -fvm%se_flux(i,j,iside,k)
            if (flux_se>eps) then
              gamma(iside)=0.5_r8
              !
              ! this copying is necessary since get_flux_segments_area_iterate change x and dx
              !
              x_tmp (:,1:num_seg(1,iside),:)=x (:,1:num_seg(1,iside),:,iside,i,j)
              dx_tmp(:,1:num_seg(1,iside),:)=dx(:,1:num_seg(1,iside),:,iside,i,j)
              call get_flux_segments_area_iterate(&
                   x_tmp(:,:,:),x_static(:,:,:,iside,i,j),dx_static(:,:,:,iside,i,j),dx_tmp(:,:,:),&
                   x_start(:,:,iside,i,j),dgam_vec(:,:,iside,i,j),num_seg(:,iside),num_seg_static(:,iside),&
                   num_seg_max,num_area,dp_area,flowcase(iside),gamma(iside),flux_se,0.0_r8,1.0_r8)
              fvm%se_flux(i,j,iside,k) = ABS(SUM(gamma(iside)*dgam_vec(:,1,iside,i,j)))
              if (gamma(iside)>1_r8) then
                gamma(iside)=1.0_r8-eps
              end if              
            else
              fvm%se_flux(i,j,iside,k) = 0.0_r8
            end if
          enddo
        end do
      end do
    end do
  end subroutine compute_displacements_for_swept_areas



  subroutine get_flux_segments_area_iterate(x,x_static,dx_static,dx,x_start,dgam_vec,num_seg,num_seg_static,&
       num_seg_max,num_area,c,flow_case,gamma,flux,gamma_min,gamma_max)
    implicit none
    integer                                                , intent(in)    :: num_area, num_seg_max
    REAL(KIND=r8), dimension(2,num_seg_max,num_area), intent(in)    :: x_static, dx_static
    REAL(KIND=r8), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx
    integer             , dimension(num_area              ), intent(in) :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8)                   , intent(in) :: x_start, dgam_vec
    REAL(KIND=r8)                                   , intent(inout) :: gamma
    REAL(KIND=r8)                                   , intent(in) :: flux,gamma_min,gamma_max
    integer                                                , intent(in) :: flow_case

    real (kind=r8), dimension(num_area)             , intent(in) :: c

    real (kind=r8)                                :: flux_static
    real (kind=r8)                                :: weight_area(num_area), xtmp(2), xtmp2(2)
    real (kind=r8)                                :: gamma1, gamma2, gamma3, dgamma, f1, f2

    real (kind=r8), dimension(  ngpc  ) :: xq,yq
    real (kind=r8), dimension(  ngpc,1) :: F !linear

    real (kind=r8) :: xq2,xq2i, rho, rhoi, yrh, w_static(num_area)

    integer :: iseg,iarea,iter,ipt
    integer, parameter :: iter_max=20
    logical :: lexit_after_one_more_iteration

    lexit_after_one_more_iteration = .false.
    !
    ! compute static line-integrals (not necessary to recompute them for every iteration)
    !
    flux_static = 0.0_r8
    w_static    = 0.0_r8
    weight_area = 0.0_r8
    do iarea=1,num_area
       do iseg=1,num_seg_static(iarea)

!rck vector directive needed here
!DIR$ SIMD
          do ipt=1,ngpc
             xq(ipt) = x_static(1,iseg,iarea)+dx_static(1,iseg,iarea)*gspts(ipt)! create quadrature point locations
             yq(ipt) = x_static(2,iseg,iarea)+dx_static(2,iseg,iarea)*gspts(ipt)
             F(ipt,1) = yq(ipt)/(SQRT(1.0_r8+xq(ipt)*xq(ipt) + yq(ipt)*yq(ipt))*(1.0_r8+xq(ipt)*xq(ipt)))! potential ! potential
          enddo
          weight_area(iarea) = weight_area(iarea)+sum(gsweights(:)*F(:,1))*0.5_r8*dx_static(1,iseg,iarea) !integral
       end do
       w_static(iarea)= weight_area(iarea)
       flux_static = flux_static+weight_area(iarea)*c(iarea)      !add to swept flux
    end do
    !
    ! initilization
    !
    gamma1=0.0_r8; f1=-flux   ! zero flux guess 1
    !
    ! compute flux integrals of first guess passed to subroutine
    !
    gamma2=gamma
    f2 = 0.0_r8
    weight_area=w_static
    do iarea=1,num_area
       do iseg=1,num_seg(iarea)
!rck vector directive needed here
!DIR$ SIMD
          do ipt=1,ngpc
             xq(ipt)  = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(ipt)! create quadrature point locations
             yq(ipt)  = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(ipt)
             xq2      =  xq(ipt)*xq(ipt)
             xq2i     =  1.0_r8/(1.0_r8+xq2)
             rho      =  SQRT(1.0_r8+xq2+yq(ipt)*yq(ipt))
             rhoi     =  1.0_r8/rho
             yrh      =  yq(ipt)*rhoi
             F(ipt,1) =  yrh*xq2i
          enddo
          weight_area(iarea) = weight_area(iarea)+sum(gsweights(:)*F(:,1))*0.5_r8*dx(1,iseg,iarea)! integral
       end do
       f2 = f2+weight_area(iarea)*c(iarea)
    end do
    f2 = f2-flux !integral error
    iter=0
    if (abs(f2-f1)<eps) then
      !
      ! in case the first guess is converged
      !
      return
    end if
    

    dgamma=(gamma2-gamma1)*f2/(f2-f1);
    gamma3 = gamma2-dgamma;                    ! Newton "guess" for gamma
    gamma1 = gamma2; f1 = f2; gamma2 = gamma3; ! prepare for iteration
    do iter=1,iter_max
       !
       ! update vertex location: flow_case dependent to avoid many zero operations
       !
       select case(flow_case)
       case(-4)
          iarea=1
          dx       (:,2,1) =  gamma3*dgam_vec (:,1)
          x        (:,1,1) =  x_start(:,1)+gamma3*dgam_vec(:,1)
          dx       (:,1,1) = -gamma3*dgam_vec (:,1)

       case(-2)
          iarea=1
          dx       (:,1,iarea) =  gamma3*dgam_vec (:,1)
          x        (:,2,iarea) =  x_start(:,2)+gamma3*dgam_vec(:,1)
          dx       (:,2,iarea) = -gamma3*dgam_vec (:,1)
       case(-1)
          !
          ! to compute first-guess perpendicular displacements for iside=1
          !
          iarea=1          
          x        (:,1,iarea) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx       (:,1,iarea) = -dx_static(:,1,iarea)
          x        (:,2,iarea) = x_start(:,2)+gamma3*dgam_vec(:,1)
          dx       (:,2,iarea) = x_start(:,2)-x(:,2,iarea)
       case(0)
          iarea=3
          xtmp = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx       (:,1,iarea) = xtmp(:  )-x(:,1,iarea)           !dynamic - line 2
          x        (:,2,iarea) = xtmp(:  )                        !dynamic - line 3
          dx       (:,2,iarea) = x_static(:,2,iarea)-x(:,2,iarea) !dynamic - line 3
       case(1)
          iarea=2
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 2
          x   (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x    (:,2,iarea) = xtmp (:)         !dynamic
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(2)
          iarea=2
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 2
          x   (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3

          iarea=3
          xtmp(:        ) = x_start(:,4)+gamma3*dgam_vec(:,4)!
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 1
          x   (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(3)
          iarea         = 3
          xtmp    (:  ) = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,2,iarea)-xtmp(:) !dynamic - line 2

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(4)
          iarea           = 1
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x (:,2,iarea) = xtmp(:)                      !dynamic
          dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic

          iarea         = 2
          xtmp    (:  ) = x_start(:,2)+gamma3*dgam_vec(:,2)
          xtmp2   (:  ) = x_start(:,3)+gamma3*dgam_vec(:,3)

          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic

          x (:,2,iarea) = xtmp (:)          !dynamic
          dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic

          x (:,3,iarea) = xtmp2(:)                !dynamic
          dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic

          iarea            = 3
          xtmp (:        ) = x_start(:,4)+gamma3*dgam_vec(:,4)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic - line 1
          x    (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx   (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(5)
          iarea                = 3
          xtmp    (:  )        = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,2,iarea)-xtmp(:) !dynamic - line 2

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          xtmp2   (:  ) = x_start(:,7)+gamma3*dgam_vec(:,7)

          dx(:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x (:,2,iarea) = xtmp(:)          !dynamic -line 2
          dx       (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic - line 2
          x        (:,3,iarea) = xtmp2(:)               !dynamic  -line 1
          dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic - line 1

          iarea             = 5
          xtmp  (:  )       = x_start(:,8)+gamma3*dgam_vec(:,8)

          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x        (:,2,iarea) = xtmp(:)                     !dynamic -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(6)
          iarea = 1
          xtmp(:  ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x (:,2,iarea) = xtmp(:)                      !dynamic
          dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic

          iarea         = 2
          xtmp    (:  ) = x_start(:,2)+gamma3*dgam_vec(:,2)
          xtmp2   (:  ) = x_start(:,3)+gamma3*dgam_vec(:,3)

          dx(:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic
          x (:,2,iarea) = xtmp (:)          !dynamic
          dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic
          x (:,3,iarea) = xtmp2(:)                !dynamic
          dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x    (:,2,iarea) = xtmp (:)         !dynamic
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(7)
          iarea=2
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 2
          x   (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x    (:,2,iarea) = xtmp (:)         !dynamic
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea      = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          xtmp2   (:  ) = x_start(:,7)+gamma3*dgam_vec(:,7)

          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic
          x        (:,2,iarea) = xtmp(:)              !dynamic
          dx       (:,2,iarea) = xtmp2(:)-xtmp(:)     !dynamic
          x        (:,3,iarea) = xtmp2(:)               !dynamic
          dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic

          iarea      = 5
          xtmp (:  ) = x_start(:,8)+gamma3*dgam_vec(:,8)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x    (:,2,iarea) = xtmp(:)                     !dynamic -line 2
          dx   (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(8)
          iarea = 1
          xtmp(:  ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x (:,2,iarea) = xtmp(:)                      !dynamic
          dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic

          iarea         = 2
          xtmp    (:  ) = x_start(:,2)+gamma3*dgam_vec(:,2)
          xtmp2   (:  ) = x_start(:,3)+gamma3*dgam_vec(:,3)

          dx(:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic
          x (:,2,iarea) = xtmp (:)          !dynamic
          dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic
          x (:,3,iarea) = xtmp2(:)                !dynamic
          dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x    (:,2,iarea) = xtmp (:)         !dynamic
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea      = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          xtmp2   (:  ) = x_start(:,7)+gamma3*dgam_vec(:,7)

          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic
          x        (:,2,iarea) = xtmp(:)              !dynamic
          dx       (:,2,iarea) = xtmp2(:)-xtmp(:)     !dynamic
          x        (:,3,iarea) = xtmp2(:)               !dynamic
          dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic

          iarea      = 5
          xtmp (:  ) = x_start(:,8)+gamma3*dgam_vec(:,8)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x    (:,2,iarea) = xtmp(:)                     !dynamic -line 2
          dx   (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case default
          call endrun('flow case not defined in get_flux_segments_area_iterate')
       end select
       !
       ! compute flux integral
       !
       f2 = 0.0_r8
       weight_area=w_static
       do iarea=1,num_area
         do iseg=1,num_seg(iarea)
!rck vector directive needed here
!DIR$ SIMD
           do ipt=1,ngpc

             xq(ipt) = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(ipt)! create quadrature point locations
             yq(ipt) = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(ipt)

             xq2      =  xq(ipt)*xq(ipt)
             xq2i     =  1.0_r8/(1.0_r8+xq2)
             rho      =  SQRT(1.0_r8+xq2+yq(ipt)*yq(ipt))
             rhoi     =  1.0_r8/rho
             yrh      =  yq(ipt)*rhoi
             F(ipt,1) =  yrh*xq2i
           end do
           weight_area(iarea) = weight_area(iarea)+sum(gsweights(:)*F(:,1))*0.5_r8*dx(1,iseg,iarea)! integral
         end do
         f2 = f2+weight_area(iarea)*c(iarea)
       end do
       f2 = f2-flux !integral error

       !
       ! uncommented logic leads to noise in PS in FKESSLER at element boundary
       !
       !       if (ABS(f2)<eps.or.ABS((gamma2-gamma1)*f2)<eps.or.lexit_after_one_more_iteration) then
       !
       if (ABS(f2)<eps.or.lexit_after_one_more_iteration) then
         gamma=gamma3
         if (gamma>gamma_max) then
           lexit_after_one_more_iteration=.true.
           gamma=gamma_max
           gamma3=gamma_max
         else
           exit
         end if
       else
          !
          ! Newton increment
          !
         if (abs(f2-f1)<eps) then
           !
           ! if entering here abs(f2)>eps and abs(f1)>eps but abs(f2-f1)<eps
           !
           dgamma=-0.5_r8*(gamma2-gamma1)
           lexit_after_one_more_iteration=.true.
         else
           dgamma=(gamma2-gamma1)*f2/(f2-f1)
         endif
         if (ABS(dgamma)>eps) then
           gamma3 = gamma2-dgamma;
         else
           !
           ! dgamma set to minimum displacement to avoid f2-f1=0
           !
           gamma3=gamma2-SIGN(1.0_r8,dgamma)*eps
           write(*,*) "WARNING: setting gamma to min",gamma3,iter
         end if
         gamma3=MAX(gamma3,gamma_min)
         !
         ! prepare for next iteration
         !
         gamma1 = gamma2; f1 = f2; gamma2 = gamma3;
       endif
     end do
     if (iter>iter_max) write(*,*) "WARNING: iteration not converged",&
          ABS(f2),flux,gamma1,gamma2,gamma3
  end subroutine get_flux_segments_area_iterate

  subroutine define_swept_areas(fvm,ilev,displ,base_vec,base_vtx,idx)
    use control_mod, only : neast, nwest, seast, swest
    implicit none
    type (fvm_struct), intent(inout) :: fvm
    integer          , intent(in)    :: ilev


    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=r8)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(out) :: displ
    integer (kind=r8) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(out) :: base_vec
    real (kind=r8)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(out) :: base_vtx
    integer                  , dimension(2,num_area, imin:imax,imin:imax,num_sides), intent(out) :: idx

    real (kind=r8) :: flux_sum     (0:nc+1,0:nc+1,2)
    integer               :: degenerate   (1:nc+1,1:nc+1  )
    integer               :: circular_flow(1:nc+1,1:nc+1  )
    integer               :: illcond      (1:nc+1,1:nc+1)
    integer               :: ib,i,j,sgn, iside, iarea

    !
    ! set where reconstruction function is as a function of area and side
    !
    integer, dimension(num_area*4), parameter :: idx_shift_tmp = (/-1,-1, 0, 1, 1,&  !iside=1
                                                                    1, 0, 0, 0, 1,&  !iside=2
                                                                    1, 1, 0,-1,-1,&  !iside=3
                                                                   -1, 0, 0, 0,-1/)  !iside=4

    integer, dimension(num_area*4), parameter :: idy_shift_tmp = (/-1, 0, 0, 0,-1,&  !iside=1
                                                                   -1,-1, 0, 1, 1,&  !iside=2
                                                                    1, 0, 0, 0, 1,&  !iside=3
                                                                    1, 1, 0,-1,-1/)  !iside=4

    integer, dimension(num_area,4), parameter :: idx_shift = RESHAPE(idx_shift_tmp,(/num_area,4/))
    integer, dimension(num_area,4), parameter :: idy_shift = RESHAPE(idy_shift_tmp,(/num_area,4/))

    integer, dimension(4), parameter :: iside_m1 = (/4,1,2,3/)
    integer, dimension(4), parameter :: iside_p1 = (/2,3,4,1/)
    integer, dimension(4), parameter :: iside_p2 = (/3,4,1,2/)
    integer, dimension(4), parameter :: iside_p3 = (/4,1,2,3/)

    integer, dimension(4), parameter :: imin_side = (/1   ,0   ,1   ,1   /)
    integer, dimension(4), parameter :: imax_side = (/nc  ,nc  ,nc  ,nc+1/)
    integer, dimension(4), parameter :: jmin_side = (/1   ,1   ,0   ,1   /)
    integer, dimension(4), parameter :: jmax_side = (/nc+1,nc  ,nc  ,nc  /)



    integer :: iur,jur,ilr,jlr,iul,jul,ill,jll

    ib = fvm%cubeboundary
    flux_sum(0:nc+1,1:nc+1,1) = fvm%se_flux(0:nc+1,0:nc  ,3,ilev)-fvm%se_flux(0:nc+1,1:nc+1,1,ilev)
    flux_sum(1:nc+1,0:nc+1,2) = fvm%se_flux(0:nc  ,0:nc+1,2,ilev)-fvm%se_flux(1:nc+1,0:nc+1,4,ilev)

    !
    ! Degenerate case ("two departure points")
    !
    !           ||  |                        || no change in this situation ||  no change in this situation
    !           ||  |                        ||                             ||
    !           ||--------                   ||----------                   ||----------
    !           ||  |                        ||                             ||
    ! =======================      =======================         =====================
    !       |   ||                       |   ||                             ||
    !  -----|---||                 ------|---||                    ---------||
    !       |   ||                       |   ||                             ||
    !       |   ||                       |   ||                             ||
    !
    !
    where (flux_sum(0:nc,1:nc+1,1)*flux_sum(1:nc+1,1:nc+1,1)<0.0_r8.and.flux_sum(1:nc+1,0:nc,2)*flux_sum(1:nc+1,1:nc+1,2)<0.0_r8)
       degenerate(:,:) = 0
    elsewhere
       degenerate(:,:) = 1
    end where

    if (ib>0) then
       if (ib==swest) degenerate(1   ,1   ) = 1
       if (ib==nwest) degenerate(1   ,nc+1) = 1
       if (ib==neast) degenerate(nc+1,nc+1) = 1
       if (ib==seast) degenerate(nc+1,1   ) = 1
    end if

    do j=1,nc+1
       do i=1,nc+1
          do sgn=-1,1,2
             if (&
                  sgn*flux_sum(i-1,j,1)<0.0_r8.and.sgn*flux_sum(i,j-1,2)>0.0_r8.and.&
                  sgn*flux_sum(i  ,j,1)>0.0_r8.and.sgn*flux_sum(i,j  ,2)<0.0_r8) then
                circular_flow(i,j) = 0
             else
                circular_flow(i,j) = 1
             end if
          end do
       end do
    end do
    !
    ! wrap around corners
    !
    if (ib==nwest) then
       flux_sum(0,nc+1,1) = fvm%se_flux(0,nc,3,ilev)-fvm%se_flux(1,nc+1,4,ilev)
       flux_sum(1,nc+1,2) = fvm%se_flux(0,nc,3,ilev)-fvm%se_flux(1,nc+1,4,ilev)

       i=1;j=nc+1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i,j-1,2)>0.0_r8.and.&
               sgn*flux_sum(i  ,j,1)>0.0_r8.and.sgn*flux_sum(i,j  ,2)<0.0_r8) then
             circular_flow(i,j) = 0
          end if
       end do
    else if (ib==swest) then
       flux_sum(0,1,1) = fvm%se_flux(1,0,4,ilev)-fvm%se_flux(0,1,1,ilev)
       flux_sum(1,0,2) = fvm%se_flux(0,1,1,ilev)-fvm%se_flux(1,0,4,ilev)
       i=1;j=1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i-1,j,1)<0.0_r8.and.&
               sgn*flux_sum(i  ,j,1)>0.0_r8.and.sgn*flux_sum(i,j  ,2)<0.0_r8) then
             circular_flow(i,j) = 0
          end if
       end do
    else if (ib==neast) then
       flux_sum(nc+1,nc+1,1) = fvm%se_flux(nc+1,nc,3,ilev)-fvm%se_flux(nc,nc+1,2,ilev)
       flux_sum(nc+1,nc+1,2) = fvm%se_flux(nc,nc+1,2,ilev)-fvm%se_flux(nc+1,nc,3,ilev)
       i=nc+1;j=nc+1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i-1,j,1)<0.0_r8.and.sgn*flux_sum(i,j-1,2)>0.0_r8.and.&
               sgn*flux_sum(i,j  ,2)<0.0_r8) then
             circular_flow(i,j) = 0
          end if
       end do
    else if (ib==seast) then
       flux_sum(nc+1,1   ,1) = fvm%se_flux(nc,0,2,ilev)-fvm%se_flux(nc+1,1,1,ilev)
       flux_sum(nc+1,0   ,2) = fvm%se_flux(nc,0,2,ilev)-fvm%se_flux(nc+1,1,1,ilev)
       i=nc+1;j=1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i-1,j,1)<0.0_r8.and.sgn*flux_sum(i,j-1,2)>0.0_r8.and.&
               sgn*flux_sum(i,j  ,2)<0.0_r8) then
             circular_flow(i,j) = 0
          end if
       end do
    end if
    illcond = circular_flow*degenerate
    !
    !
    !
    !
    do iside=1,4
       do j=jmin_side(iside),jmax_side(iside)
          do i=imin_side(iside),imax_side(iside)
             if (fvm%se_flux(i,j,iside,ilev)>eps) then
                iur = i+idx_shift(4,iside); jur = j+idy_shift(4,iside) !(i,j) index of upper right quadrant
                ilr = i+idx_shift(5,iside); jlr = j+idy_shift(5,iside) !(i,j) index of lower left  quadrant
                iul = i+idx_shift(2,iside); jul = j+idy_shift(2,iside) !(i,j) index of upper right quadrant
                ill = i+idx_shift(1,iside); jll = j+idy_shift(1,iside) !(i,j) index of lower left  quadrant

                !iside=1
                if (iside==1) then
                displ(0,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i,j)     !center left
                displ(1,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i+1,j)   !center right
                displ(2,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j)   !c2
                displ(3,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j)   !c3
                displ(4,i,j,iside) = -flux_sum   (i+1,j  ,1)*illcond(i+1,j)   !r1
                displ(5,i,j,iside) = -flux_sum   (i+1,j-1,2)*illcond(i+1,j)   !r2
                displ(6,i,j,iside) = -flux_sum   (i-1,j  ,1)*illcond(i  ,j)   !l1
                displ(7,i,j,iside) =  flux_sum   (i  ,j-1,2)*illcond(i  ,j)   !l2

                end if
                if (iside==2) then
                !iside=2
                displ(0,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j  )     !center left
                displ(1,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j+1)   !center right
                displ(2,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i+1,j+1)   !c2
                displ(3,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i+1,j  )   !c3
                displ(4,i,j,iside) =  flux_sum   (i+1,j+1,2)*illcond(i+1,j+1)   !r1
                displ(5,i,j,iside) = -flux_sum   (i+1,j+1,1)*illcond(i+1,j+1)   !r2
                displ(6,i,j,iside) =  flux_sum   (i+1,j-1,2)*illcond(i+1,j)   !l1
                displ(7,i,j,iside) =  flux_sum   (i+1,j  ,1)*illcond(i+1,j)   !l2
                end if
                !iside=3
                if (iside==3) then
                displ(0,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i+1,j+1)     !center left
                displ(1,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i  ,j+1)   !center right
                displ(2,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j+1)   !c2
                displ(3,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j+1)   !c3
                displ(4,i,j,iside) =  flux_sum   (i-1,j+1,1)*illcond(i  ,j+1)   !r1
                displ(5,i,j,iside) =  flux_sum   (i  ,j+1,2)*illcond(i  ,j+1)   !r2
                displ(6,i,j,iside) =  flux_sum   (i+1,j+1,1)*illcond(i+1,j+1)   !l1
                displ(7,i,j,iside) = -flux_sum   (i+1,j+1,2)*illcond(i+1,j+1)   !l2
                end if
                if (iside==4) then
                !iside=4
                displ(0,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j+1)     !center left
                displ(1,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j  )   !center right
                displ(2,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i  ,j  )   !c2
                displ(3,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i  ,j+1)   !c3
                displ(4,i,j,iside) = -flux_sum   (i  ,j-1,2)*illcond(i  ,j  )   !r1
                displ(5,i,j,iside) =  flux_sum   (i-1,j  ,1)*illcond(i  ,j  )   !r2
                displ(6,i,j,iside) = -flux_sum   (i  ,j+1,2)*illcond(i  ,j+1)   !l1
                displ(7,i,j,iside) = -flux_sum   (i-1,j+1,1)*illcond(i  ,j+1)   !l2
                end if

                base_vtx(:,1,i,j,iside) = fvm%vtx_cart(iside,:,i  ,j            )       !vertex center left
                base_vtx(:,2,i,j,iside) = fvm%vtx_cart(iside_p1(iside),:,i  ,j  )       !vertex center right
                base_vtx(:,3,i,j,iside) = fvm%vtx_cart(iside,:,iur,jur          )       !vertex upper right
                base_vtx(:,4,i,j,iside) = fvm%vtx_cart(iside_p3(iside),:,ilr,jlr)       !vertex lower right
                base_vtx(:,5,i,j,iside) = fvm%vtx_cart(iside_p1(iside),:,iul,jul)       !vertex upper left
                base_vtx(:,6,i,j,iside) = fvm%vtx_cart(iside_p2(iside),:,ill,jll)       !vertex lower left

                base_vec(:, 1,i,j,iside) = fvm%flux_vec    (:,i  ,j  ,iside          )      !vector center
                base_vec(:, 2,i,j,iside) = fvm%flux_vec    (:,i  ,j  ,iside_p1(iside))      !vector center right
                base_vec(:, 3,i,j,iside) = fvm%flux_vec    (:,i  ,j  ,iside_p3(iside))      !vector center left
                base_vec(:, 4,i,j,iside) = fvm%flux_vec    (:,iur,jur,iside          )      !vector upper right 1
                base_vec(:, 5,i,j,iside) = fvm%flux_vec    (:,iur,jur,iside_p3(iside))      !vector upper right 2
                base_vec(:, 6,i,j,iside) = fvm%flux_vec    (:,ilr,jlr,iside_p3(iside))      !vector lower right 1
                base_vec(:, 7,i,j,iside) = fvm%flux_vec    (:,ilr,jlr,iside_p2(iside))      !vector lower right 2
                base_vec(:, 8,i,j,iside) = fvm%flux_vec    (:,iul,jul,iside          )      !vector upper left 1
                base_vec(:, 9,i,j,iside) = fvm%flux_vec    (:,iul,jul,iside_p1(iside))      !vector upper left 2
                base_vec(:,10,i,j,iside) = fvm%flux_vec    (:,ill,jll,iside_p1(iside))      !vector lower left 1
                base_vec(:,11,i,j,iside) = fvm%flux_vec    (:,ill,jll,iside_p2(iside))      !vector lower left 2

                do iarea=1,5
                   idx(1,iarea,i,j,iside) = i+idx_shift(iarea,iside)
                   idx(2,iarea,i,j,iside) = j+idy_shift(iarea,iside)
                end do
             else
                displ(:,i,j,iside) = 9D99!for debugging
             end if
          end do
       end do
    end do
    !
    ! wrap around corners here
    !

  end subroutine define_swept_areas


  !
  ! Notation conventions used in define_area subroutines
  !
  !
  !
  !   ^    ||--->   ^   <---||    ^
  !  /|\   || 3    /|\    2 ||   /|\
  !   | 6  ||     1 |       ||    | 4
  !   |    ||       |       ||    |
  ! =================================
  !        ||               ||
  !        ||               ||
  !      7 ||               || 5
  !    <---||               ||--->
  !

  subroutine define_area1_area2(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
       x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=r8)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=r8), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=r8)    , dimension(2,3) :: xdep !departure points
    real (kind=r8)                     :: gamma
    integer :: iarea


    REAL(KIND=r8) :: xtmp(2),xtmp2(2)
    !
    !
    !        ||-----        ||
    !       /||             ||
    !      / ||             ||
    !  ===X=========================
    !     | /||             ||
    !     |/ ||             ||
    !     *  ||             ||
    !
    !
    ! crossing X
    if (SUM(ABS(base_vec(:,9,i,j,iside))).NE.0) then
       gamma = displ(0,i,j,iside)*displ(7,i,j,iside)/(displ(0,i,j,iside)-displ(6,i,j,iside))
!       gamma = MAX(MIN(gamma,displ(7,i,j,iside),-displ(3,i,j,iside)),0.0_r8)!MWR manuscript
       gamma = MAX(MIN(gamma,displ(7,i,j,iside),-0.25_r8*displ(3,i,j,iside)),0.0_r8)!dbgxxx
    else
       !
       ! corner case
       !
       gamma=displ(0,i,j,iside)
    end if


    xdep    (:,1) = base_vtx(:, 6,i,j,iside)+displ(7,i,j,iside)*base_vec(:,10,i,j,iside)-displ(6,i,j,iside)*base_vec(:,11,i,j,iside)
    x_start (:,1) = base_vtx(:, 6,i,j,iside)
    dgam_vec(:,1) = base_vec(:,10,i,j,iside)*gamma

    xdep(:,2) = base_vtx(:,2,i,j,iside)+displ(1,i,j,iside)*base_vec(:, 1,i,j,iside)+displ(2,i,j,iside)*base_vec(:, 2,i,j,iside)

    iarea                  = 1
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1

    x_static (:,1,iarea) = base_vtx(:,6,i,j,iside)       !static
    dx_static(:,1,iarea) = xdep(:,1)-x_static(:,1,iarea) !static

    xtmp(:        ) = x_start(:,1)+dgam_vec(:,1)
    x   (:,1,iarea) = xdep(:,1)                  !static
    dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic

    x (:,2,iarea) = xtmp(:)                      !dynamic
    dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic
    !
    !
    !
    iarea                  = 2
    num_seg       (iarea)  = 3

    x_start (:,2) = base_vtx(:,5,i,j,iside)
    dgam_vec(:,2) = base_vec(:,9,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,2)+dgam_vec(:,2)

    x_start (:,3) = base_vtx(:,5,i,j,iside)
    dgam_vec(:,3) = base_vec(:,8,i,j,iside)*displ(0,i,j,iside)
    xtmp2   (:  ) = x_start(:,3)+dgam_vec(:,3)

    x   (:,1,iarea) = base_vtx(:,5,i,j,iside) !static
    dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic

    x (:,2,iarea) = xtmp (:)          !dynamic
    dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic

    x (:,3,iarea) = xtmp2(:)                !dynamic
    dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic
  end subroutine define_area1_area2


  subroutine define_area2(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,x_start, dgam_vec,&
       displ_first_guess)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=r8)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=r8), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=r8)    , dimension(2,3) :: xdep !departure points
    real (kind=r8), optional, intent(out)        :: displ_first_guess
    real (kind=r8) :: gamma
    integer :: iarea


    REAL(KIND=r8) :: xtmp(2)
    ! *: xdep(:,1)
    ! x: xtmp
    !
    !      2 ||             ||
    !     *--x              ||
    !     1\3||1            ||
    !       \||             ||
    !  =============================
    !        ||             ||
    !
    !
    ! compute departure points (xdep(1) is left; xdep(3) is right and xdep(2) is midway
    !
    xdep(:,1) = base_vtx(:,5,i,j,iside)+&
         MAX(0.0_r8,displ(6,i,j,iside))*base_vec(:,8,i,j,iside)-displ(3,i,j,iside)*base_vec(:,9,i,j,iside)
    x_start (:,1) = base_vtx(:,5,i,j,iside)
    gamma         = displ(0,i,j,iside)
    dgam_vec(:,1) = base_vec(:,8,i,j,iside)*gamma
    if (present(displ_first_guess)) displ_first_guess = gamma

    iarea                  = 2
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1

    x_static (:,1,iarea) = base_vtx(:,5,i,j,iside)       !static  - line 1
    dx_static(:,1,iarea) = xdep(:,1)-x_static(:,1,iarea) !static  - line 1

    xtmp     (:        ) = x_start(:,1)+dgam_vec(:,1)
    x        (:,1,iarea) = xdep(:,1)                  !static  - line 2
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic - line 2

    x        (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3
  end subroutine define_area2


  subroutine define_area3_left(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, &
       num_seg, num_seg_static,x_start, dgam_vec,displ_first_guess)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=r8)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=r8), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8)                   , intent(inout):: x_start, dgam_vec
    real (kind=r8), optional, intent(out)        :: displ_first_guess

    real (kind=r8)    , dimension(2,3) :: xdep !departure points
    real (kind=r8)                     :: gamma
    integer :: iarea


    REAL(KIND=r8) :: xtmp(2)

    ! iarea = 3
    !-------------------------------------------------------------------------------------------
    !
    !          xtmp         xdep(2)
    !           |x-----2------*   ||
    !           ||             \  ||
    !           |1              3 ||
    !           ||               \||
    !        ===========4==============
    !
    !
    xdep(:,2) = base_vtx(:,2,i,j,iside)+displ(1,i,j,iside)*base_vec(:,1,i,j,iside)&
         +MAX(0.0_r8,displ(2,i,j,iside))*base_vec(:,2,i,j,iside)
    x_start (:,4) = base_vtx(:,1,i,j,iside)
    gamma         = displ(0,i,j,iside)
    dgam_vec(:,4) = base_vec(:,1,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,4)+dgam_vec(:,4)

    if (present(displ_first_guess)) displ_first_guess = gamma

    iarea                  = 3
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 2

    x_static (:,1,iarea) = xdep(:,2)                         !static  - line 3
    dx_static(:,1,iarea) = base_vtx(:,2,i,j,iside)-xdep(:,2) !static  - line 3

    x_static (:,2,iarea) = base_vtx(:,2,i,j,iside)                         !static  - line 4
    dx_static(:,2,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static  - line 4

    x        (:,1,iarea) = base_vtx(:,1,i,j,iside)    !static  - line 1
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic - line 1

    x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area3_left

  subroutine define_area3_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=r8)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=r8), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=r8), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=r8)    , dimension(2,3) :: xdep !departure points
    real (kind=r8)                     :: gamma
    integer :: iarea

    REAL(KIND=r8) :: xtmp(2)
    !
    !
    !        ||  *-----2----||\
    !        || /1         3|| \
    !        ||/      4     ||
    !  =============================
    !        ||             ||
    !        ||             ||
    !        ||             ||
    !
    xdep(:,1) = base_vtx(:,1,i,j,iside)+displ(0,i,j,iside)*base_vec(:,1,i,j,iside)&
         +MAX(0.0_r8,displ(3,i,j,iside))*base_vec(:,3,i,j,iside)
    x_start (:,5) = base_vtx(:,2,i,j,iside)
    gamma         = displ(1,i,j,iside)
    dgam_vec(:,5) = base_vec(:,1,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,5)+dgam_vec(:,5)

    iarea                  = 3
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 2

    x_static (:,1,iarea) = base_vtx(:,1,i,j,iside)           !static  - line 1
    dx_static(:,1,iarea) = xdep(:,1)-base_vtx(:,1,i,j,iside) !static  - line 1

    x_static (:,2,iarea) = base_vtx(:,2,i,j,iside)                         !static  - line 4
    dx_static(:,2,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static  - line 4

    x        (:,1,iarea) = xdep(:,1)            !static  - line 2
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic - line 2

    x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
    dx       (:,2,iarea) = x_static(:,2,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area3_right


  subroutine define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec)
    implicit none
    integer, parameter  :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    integer, parameter  :: num_seg_max=5
    integer,                                                                 intent(in)   :: i,j,iside
    real (kind=r8),    dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout):: displ
    integer (kind=r8), dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout):: base_vec
    real (kind=r8),    dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout):: base_vtx
    real(KIND=r8),     dimension(2,num_seg_max,num_area),                    intent(inout):: x, dx, x_static, dx_static
    integer,           dimension(num_area),                                  intent(inout):: num_seg, num_seg_static
    real(KIND=r8),     dimension(2,8),                                       intent(inout):: x_start, dgam_vec

    real (kind=r8)      :: gamma
    integer             :: iarea
    real(KIND=r8)       :: xtmp(2),xtmp2(2)
    !
    !        ||-------------||
    !       /||             ||\
    !        ||             ||
    !  =============================
    !        ||             ||
    !        ||             ||
    !        ||             ||
    !
    x_start (:,4) = base_vtx(:,1,i,j,iside)
    x_start (:,5) = base_vtx(:,2,i,j,iside)
    gamma         = displ(0,i,j,iside)
    dgam_vec(:,4) = base_vec(:,1,i,j,iside)*gamma
    dgam_vec(:,5) = base_vec(:,1,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,4)+dgam_vec(:,4)
    xtmp2   (:  ) = x_start(:,5)+dgam_vec(:,5)

    iarea                  = 3
    num_seg       (iarea)  = 3
    num_seg_static(iarea)  = 1

    x_static (:,1,iarea) = base_vtx(:,2,i,j,iside)                         !static
    dx_static(:,1,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static

    x        (:,1,iarea) = base_vtx(:,1,i,j,iside)    !static
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic

    x        (:,2,iarea) = xtmp (:)         !dynamic
    dx       (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic

    x        (:,3,iarea) = xtmp2(:)              !dynamic
    dx       (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic
  end subroutine define_area3_left_right

  subroutine define_area4_area5(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec,displ_first_guess)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    integer, parameter :: num_seg_max=5
    real (kind=r8),    dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8), dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8),    dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    real(KIND=r8),     dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer,           dimension(num_area),               intent(inout) :: num_seg, num_seg_static
    real(KIND=r8),     dimension(2,8),                    intent(inout) :: x_start, dgam_vec
    real(KIND=r8),     optional,                          intent(out)   :: displ_first_guess


    real (kind=r8)    , dimension(2,3) :: xdep !departure points
    real (kind=r8)                     :: gamma
    integer :: iarea

    real(KIND=r8) :: xtmp(2),xtmp2(2)
    !
    !        ||     --------||
    !        ||             ||\
    !        ||             || \
    !  =============================
    !        ||             ||\ |
    !        ||             || \|
    !        ||             ||  *
    !
    !
    ! iarea  = 4
    !
    iarea                  = 4
    num_seg       (iarea)  = 3

    if (SUM(ABS(base_vec(:,5,i,j,iside))).NE.0) then
       gamma = displ(1,i,j,iside)*displ(5,i,j,iside)/(displ(1,i,j,iside)-displ(4,i,j,iside))
!       gamma = MAX(MIN(gamma,displ(5,i,j,iside),-displ(2,i,j,iside)),0.0_r8)!MWR manuscript
       gamma = MAX(MIN(gamma,displ(5,i,j,iside),-0.25_r8*displ(2,i,j,iside)),0.0_r8)
    else
       !
       ! corner case
       !
       gamma = displ(1,i,j,iside)
    end if

    if (present(displ_first_guess)) displ_first_guess = displ(1,i,j,iside)

    x_start (:,6) = base_vtx(:,3,i,j,iside)
    dgam_vec(:,6) = base_vec(:,4,i,j,iside)*displ(1,i,j,iside)
    xtmp    (:  ) = x_start(:,6)+dgam_vec(:,6)
    x_start (:,7) = base_vtx(:,3,i,j,iside)
    dgam_vec(:,7) = base_vec(:,5,i,j,iside)*gamma
    xtmp2   (:  ) = x_start(:,7)+dgam_vec(:,7)

    x        (:,1,iarea) = base_vtx(:,3,i,j,iside)!static   -line 1
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1

    x        (:,2,iarea) = xtmp(:)          !dynamic -line 2
    dx       (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic - line 2

    x        (:,3,iarea) = xtmp2(:)               !static   -line 1
    dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic - line 1
    !
    !iarea = 5
    !
    xdep(:,1) = base_vtx(:,4,i,j,iside)+displ(5,i,j,iside)*base_vec(:,6,i,j,iside)&
         -displ(4,i,j,iside)*base_vec(:,7,i,j,iside)
    x_start (:,8) = base_vtx(:,4,i,j,iside)
    dgam_vec(:,8) = base_vec(:,6,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,8)+dgam_vec(:,8)

    iarea                  = 5
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1

    x        (:,1,iarea) = base_vtx(:,4,i,j,iside)!static   -line 1
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1

    x_static (:,1,iarea) = xdep(:,1)                        !static - line 1
    dx_static(:,1,iarea) = x(:,1,iarea)-x_static(:,1,iarea) !static - line 1

    x        (:,2,iarea) = xtmp(:)                     !dynamic -line 2
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area4_area5


  subroutine define_area4(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec,displ_first_guess)
    implicit none
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    integer, parameter :: num_seg_max=5
    integer,                                                                 intent(in)    :: i,j,iside
    real (kind=r8),    dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8), dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8),    dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx

    real(KIND=r8),     dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer,           dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    real(KIND=r8),     dimension(2,8)                   , intent(inout) :: x_start, dgam_vec
    real(KIND=r8), optional,                              intent(out)   :: displ_first_guess



    real (kind=r8), dimension(2,3) :: xdep !departure points
    real (kind=r8)                 :: gamma
    integer                        :: iarea
    real(KIND=r8)                  :: xtmp(2)

    iarea                  = 4
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1

    xdep(:,1) = base_vtx(:,3,i,j,iside)+MAX(0.0_r8,displ(4,i,j,iside))*base_vec(:,4,i,j,iside)&
         -displ(2,i,j,iside)*base_vec(:,5,i,j,iside)
    x_start (:,6) = base_vtx(:,3,i,j,iside)
    gamma         = displ(1,i,j,iside)
    dgam_vec(:,6) = base_vec(:,4,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,6)+dgam_vec(:,6)

    if (present(displ_first_guess)) displ_first_guess = gamma

    x_static (:,1,iarea) = xdep(:,1)                         !static
    dx_static(:,1,iarea) = base_vtx(:,3,i,j,iside)-xdep(:,1) !static

    x        (:,1,iarea) = base_vtx(:,3,i,j,iside) !static  - line 2
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2

    x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area4

  subroutine define_area3_center(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
       x_start, dgam_vec,se_flux_center,displ_first_guess)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    integer, parameter :: num_seg_max=5
    real (kind=r8),    dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=r8), dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=r8),    dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx

    real(KIND=r8),     dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer,           dimension(num_area),               intent(inout) :: num_seg, num_seg_static
    real(KIND=r8),     dimension(2,8),                    intent(inout) :: x_start, dgam_vec
    real(KIND=r8) ,                                       intent(in   ) :: se_flux_center
    real(KIND=r8),     optional,                          intent(out)   :: displ_first_guess

    real (kind=r8)    , dimension(2,3) :: xdep !departure points
    real (kind=r8)                     :: gamma
    integer :: iarea
    !
    !                 xdep(2)
    !                 ______X______
    !        ||      /             \      ||
    !        ||  *--/               \--*  ||
    !        || /xdep(1)         xdep(3)\ ||
    !        ||/                         \||
    !  ========================================
    !        ||                           ||
    !
    !
    ! compute departure points (xdep(1) is left; xdep(3) is right and xdep(2) is midway
    !

    xdep(:,1) = base_vtx(:,1,i,j,iside)+&
         displ(0,i,j,iside)*base_vec(:,1,i,j,iside)+displ(3,i,j,iside)*base_vec(:,3,i,j,iside)
    xdep(:,3) = base_vtx(:,2,i,j,iside)+&
         displ(1,i,j,iside)*base_vec(:,1,i,j,iside)+displ(2,i,j,iside)*base_vec(:,2,i,j,iside)
    xdep(:,2) = 0.5_r8*(xdep(:,1)+xdep(:,3))

    gamma= se_flux_center
    x_start(:,1) = ABS(base_vec(:,3,i,j,iside))*((xdep(:,2)-base_vtx(:,1,i,j,iside)))+&
         base_vtx(:,1,i,j,iside) !xdep(2) - midway between departure points projected to side 1

    dgam_vec(:,1) = gamma*base_vec(:,1,i,j,iside)

    if (present(displ_first_guess)) displ_first_guess = gamma

    xdep(:,2)     = x_start(:,1)+dgam_vec(:,1)
    iarea                  = 3
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 3

    !                 ______X______
    !        ||    2 /             \ 3    ||
    !        ||  *--/               \--*  ||
    !        || /                       \ ||
    !        ||/ 1          5           4\||
    !  ========================================
    !        ||                           ||
    !
    x_static (:,1,iarea) = base_vtx(:,1,i,j,iside)       !static  - line 1
    dx_static(:,1,iarea) = xdep(:,1)-x_static(:,1,iarea) !static  - line 1

    x        (:,1,iarea) = xdep(:,1)                     !static  - line 2
    dx       (:,1,iarea) = xdep(:,2)-x(:,1,iarea)        !dynamic - line 2

    x        (:,2,iarea) = xdep(:,2)                     !dynamic - line 3
    dx       (:,2,iarea) = xdep(:,3)-x(:,2,iarea)        !dynamic - line 3

    x_static (:,2,iarea) = xdep(:,3)                                  !static  - line 4
    dx_static(:,2,iarea) = base_vtx(:,2,i,j,iside)-x_static(:,2,iarea)!static  - line 4

    x_static (:,3,iarea) = base_vtx(:,2,i,j,iside)                         !static - line 5
    dx_static(:,3,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static - line 5

  end subroutine define_area3_center
end module fvm_consistent_se_cslam
