module flux_avg
!-------------------------------------------------------------------------------
! Purpose: Code for smoothing surface fluxes to aviod instabilities
!-------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: begchunk, endchunk, pcols
  
  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t    
  use phys_grid,        only: get_ncols_p
  use physics_buffer, only : pbuf_add_field, dtype_r8
  implicit none
  private
  save

  ! Public interfaces

  public :: flux_avg_register
  public :: flux_avg_init
  public :: flux_avg_run
  
  ! Private module data

  integer :: lhflx_idx      ! lhflx index in physics buffer
  integer :: shflx_idx      ! shflx index in physics buffer
  integer :: qflx_idx       ! qflx index in physics buffer
  integer :: taux_idx       ! taux index in physics buffer
  integer :: tauy_idx       ! tauy index in physics buffer
  integer :: lhflx_res_idx  ! lhflx_res index in physics buffer
  integer :: shflx_res_idx  ! shflx_res index in physics buffer
  integer :: qflx_res_idx   ! qflx_res index in physics buffer
  integer :: taux_res_idx   ! taux_res index in physics buffer
  integer :: tauy_res_idx   ! tauy_res index in physics buffer

!===============================================================================
contains
!===============================================================================

subroutine flux_avg_register()
   !----------------------------------------------------------------------------
   ! Purpose: Register the fluxes in the physics buffer
   !----------------------------------------------------------------------------

   ! Request physics buffer space for fields that persist across timesteps.
   call pbuf_add_field('LHFLX',    'global',dtype_r8,(/pcols,1/),lhflx_idx)
   call pbuf_add_field('SHFLX',    'global',dtype_r8,(/pcols,1/),shflx_idx)
   call pbuf_add_field('TAUX',     'global',dtype_r8,(/pcols,1/),taux_idx)
   call pbuf_add_field('TAUY',     'global',dtype_r8,(/pcols,1/),tauy_idx)
   call pbuf_add_field('QFLX',     'global',dtype_r8,(/pcols,1/),qflx_idx)
   call pbuf_add_field('LHFLX_RES','global',dtype_r8,(/pcols,1/),lhflx_res_idx)
   call pbuf_add_field('SHFLX_RES','global',dtype_r8,(/pcols,1/),shflx_res_idx)
   call pbuf_add_field('TAUX_RES', 'global',dtype_r8,(/pcols,1/),taux_res_idx)
   call pbuf_add_field('TAUY_RES', 'global',dtype_r8,(/pcols,1/),tauy_res_idx)
   call pbuf_add_field('QFLX_RES', 'global',dtype_r8,(/pcols,1/),qflx_res_idx)

end subroutine flux_avg_register

!===============================================================================

subroutine flux_avg_init(cam_in,  pbuf2d)
   use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk
   !----------------------------------------------------------------------------
   ! Initialize the surface fluxes in the physics buffer using the cam import state
   !----------------------------------------------------------------------------
   ! Input arguments
   type(cam_in_t),      intent(in)    :: cam_in(begchunk:endchunk)
   ! Local variables
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(physics_buffer_desc), pointer :: pbuf2d_chunk(:)
   integer :: lchnk
   integer :: ncol
   !---------------------------------------------------------------------------- 

   do lchnk = begchunk, endchunk
      ncol = get_ncols_p(lchnk)
      pbuf2d_chunk => pbuf_get_chunk(pbuf2d, lchnk)
      
      call pbuf_set_field(pbuf2d_chunk, lhflx_idx,  cam_in(lchnk)%lhf(:ncol))
      call pbuf_set_field(pbuf2d_chunk, shflx_idx,  cam_in(lchnk)%shf(:ncol))
      call pbuf_set_field(pbuf2d_chunk, qflx_idx,   cam_in(lchnk)%cflx(:ncol,1))
      call pbuf_set_field(pbuf2d_chunk, taux_idx,   cam_in(lchnk)%wsx(:ncol))
      call pbuf_set_field(pbuf2d_chunk, tauy_idx,   cam_in(lchnk)%wsy(:ncol))

      call pbuf_set_field(pbuf2d,       shflx_res_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, lhflx_res_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, qflx_res_idx,  0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, taux_res_idx,  0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, tauy_res_idx,  0.0_r8)
   end do

end subroutine flux_avg_init

!===============================================================================

subroutine flux_avg_run(state, cam_in,  pbuf, nstep, deltat, timescale)
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field
   !---------------------------------------------------------------------------- 
   ! Purpose: smooth fluxes from surface components
   !----------------------------------------------------------------------------
   ! Input arguments
   type(physics_state), intent(in)    :: state
   type(cam_in_t),      intent(inout) :: cam_in
   type(physics_buffer_desc), pointer :: pbuf(:)
   
   integer,             intent(in)    :: nstep
   real(r8),            intent(in)    :: deltat
   real(r8),            intent(in)    :: timescale

   ! Local variables
   integer :: lchnk                  ! chunk identifier
   integer :: ncol                   ! number of atmospheric columns

   ! physics buffer fields
   real(r8), pointer, dimension(:) :: lhflx   ! latent heat flux
   real(r8), pointer, dimension(:) :: shflx   ! sensible heat flux
   real(r8), pointer, dimension(:) :: qflx    ! water vapor heat flux
   real(r8), pointer, dimension(:) :: taux    ! x momentum flux
   real(r8), pointer, dimension(:) :: tauy    ! y momentum flux
   real(r8), pointer, dimension(:) :: lhflx_res   ! latent heat flux
   real(r8), pointer, dimension(:) :: shflx_res   ! sensible heat flux
   real(r8), pointer, dimension(:) :: qflx_res    ! water vapor heat flux
   real(r8), pointer, dimension(:) :: taux_res    ! x momentum flux
   real(r8), pointer, dimension(:) :: tauy_res    ! y momentum flux
   !----------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, lhflx_idx,     lhflx )
   call pbuf_get_field(pbuf, shflx_idx,     shflx )
   call pbuf_get_field(pbuf, qflx_idx,      qflx  )
   call pbuf_get_field(pbuf, taux_idx,      taux  )
   call pbuf_get_field(pbuf, tauy_idx,      tauy  )

   call pbuf_get_field(pbuf, lhflx_res_idx, lhflx_res )
   call pbuf_get_field(pbuf, shflx_res_idx, shflx_res )
   call pbuf_get_field(pbuf, qflx_res_idx,  qflx_res  )
   call pbuf_get_field(pbuf, taux_res_idx,  taux_res  )
   call pbuf_get_field(pbuf, tauy_res_idx,  tauy_res  )

   call smooth(cam_in%lhf,           lhflx,lhflx_res,nstep,deltat,timescale,ncol)
   call smooth(cam_in%shf,           shflx,shflx_res,nstep,deltat,timescale,ncol)
   call smooth(cam_in%wsx,           taux, taux_res, nstep,deltat,timescale,ncol)
   call smooth(cam_in%wsy,           tauy, tauy_res, nstep,deltat,timescale,ncol)
   call smooth(cam_in%cflx(:pcols,1),qflx, qflx_res, nstep,deltat,timescale,ncol)

end subroutine flux_avg_run

!===============================================================================

subroutine smooth(new, old, res, nstep, deltat, timescale, ncol)
   !----------------------------------------------------------------------------
   ! Purpose: smooth the input data in time
   !----------------------------------------------------------------------------
   ! Input arguments
   real(r8), intent(inout) :: new(pcols)
   real(r8), intent(inout) :: old(pcols)
   real(r8), intent(inout) :: res(pcols)
   real(r8), intent(in)    :: deltat
   real(r8), intent(in)    :: timescale
   integer,  intent(in)    :: nstep
   integer,  intent(in)    :: ncol
   ! Local variables
   real(r8) :: temp(pcols)
   real(r8) :: maxres
   integer i
   !----------------------------------------------------------------------------
   temp(1:ncol) = new(1:ncol)
   if (nstep > 0) then
      new(1:ncol) = 0.5_r8*(new(1:ncol)+old(1:ncol))
   else
      old(1:ncol) = new(1:ncol)
      res(1:ncol) = 0._r8
   endif

   ! first add the flux that the surface model wanted to provide less
   ! the flux the atmosphere will actually see to the residual
   res(1:ncol) = res(1:ncol) + temp(1:ncol)-new(1:ncol)

   ! now calculate the amount that we might increment the new flux
   ! to include some of the residual. If the residual is small 
   ! we will just add it all, but if it is large we will add it 
   ! at the rate required to put the residual back into the flux 
   ! over a period specified by "timescale"
   do i = 1,ncol
      maxres = max( abs(new(i)) , abs(old(i)) )
      if (abs(res(i)).lt.maxres*0.05_r8) then
         temp(i) = res(i)
         res(i) = 0._r8
      else
         temp(i) = res(i)*deltat/timescale
         res(i) = res(i)-temp(i)
      endif
   end do

   ! dont do conservative smoothing for first 12 hours
   if (nstep*deltat/86400._r8 < 0.5_r8) then
      temp = 0._r8
      res = 0._r8
   endif

   ! make the new flux the average of the sfc model and last timestep
   ! plus some of the residual
   new(1:ncol) = new(1:ncol) + temp(1:ncol)
   old(1:ncol) = new(1:ncol)

end subroutine smooth

!===============================================================================

end module flux_avg

