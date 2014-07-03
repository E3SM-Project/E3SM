!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module barotropic

!BOP
! !MODULE: barotropic
!
! !DESCRIPTION:
!  This module contains the routine for solving the barotropic 
!  equations.
!
! !REVISION HISTORY:
!  SVN:$Id: barotropic.F90 44198 2013-02-25 22:43:22Z mlevy@ucar.edu $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod 
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_SolversMod

   use kinds_mod, only: int_kind, i4, r8
   use blocks, only: nx_block, ny_block, block, get_block
!   use distribution, only: 
   use domain_size
   use domain, only: distrb_clinic, blocks_clinic, nblocks_clinic,  &
       POP_haloClinic
   use constants, only: field_type_vector, field_type_scalar,       &
       grav, c1, c0, field_loc_NEcorner, field_loc_center
   use prognostic, only: max_blocks_clinic, GRADPX, GRADPY, UBTROP, VBTROP, &
       PSURF, curtime, oldtime, newtime, PGUESS
   use operators, only: grad, div
   use grid, only: sfc_layer_type, sfc_layer_varthick, TAREA, REGION_MASK,    &
          KMT, FCOR, HU, CALCT, sfc_layer_rigid, sfc_layer_oldfree 
   use time_management, only: mix_pass, leapfrogts, impcor, c2dtu, theta,     &
          gamma, f_euler_ts, beta, c2dtp, dtp
   use global_reductions, only: global_sum
   use forcing_fields, only: ATM_PRESS, FW
   use forcing_ap, only: ap_data_type
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now

   use overflows
   use overflow_type

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_barotropic,  &
             barotropic_driver

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     tavg_SU,            &! tavg id for vertically-integrated U
     tavg_SV              ! tavg id for vertically-integrated V

!-----------------------------------------------------------------------
!
!  nullspace removal
!
!  CHECKER +/- checkerboard field for removing global
!          checkerboard nullspace from surface pressure field.
!          zero on land and marginal seas.
!  CONSTNT constant field for removing constant nullspace
!          from surface presure field.  zero on land and in marginal
!          seas, one in open ocean.
!
!  sum_check = sum(CHECKER)
!  sum_const = sum(CONSTNT)
!
!  acheck = sum(CHECKER*TAREA)/sum(CONSTNT*TAREA)
!
!  rcheck = acheck/(sum_const - acheck*sum_check)
!  rconst =      1/(sum_const - acheck*sum_check)
!
!-----------------------------------------------------------------------

   integer (i4), dimension (:,:,:), allocatable :: &  
      CHECKER,  &! checkerboard nullspace field
      CONSTNT    ! constant     nullspace field 

   real (r8) :: &
      rcheck, rconst  ! scalar constants for checkboard removal

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_barotropic
! !INTERFACE:

   subroutine init_barotropic

! !DESCRIPTION:
!  This routine initializes barotropic quantities - mostly tavg
!  diagnostics related to barotropic velocities.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,             &! loop indices
      iblock,            &! local block index
      sum_check,         &! global sum of checkboard field 
      sum_const           ! global sum of constant   field

   real (r8) :: & 
      acheck    ! sum(CHECKER*TAREA)/sum(CONSTNT*TAREA)

   real (r8), dimension(:,:,:), allocatable :: & 
      CHECK_AREA,    &! TAREA*CHEKER
      CONST_AREA      ! TAREA*CONSTNT

   type (block) :: &
      this_block   ! block information for this block

!-----------------------------------------------------------------------
!
!  define tavg diagnostics related to barotropic velocities
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_SU,'SU',2,                                   &
            long_name='Vertically Integrated Velocity in grid-x direction', &
                          units='centimeter^2/s', grid_loc='2221',          &
                          field_type = field_type_vector,                   &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_SV,'SV',2,                                  &
           long_name='Vertically Integrated Velocity in grid-y direction', &
                          units='centimeter^2/s', grid_loc='2221',         &
                          field_type = field_type_vector,                   &
                          coordinates='ULONG ULAT time')

!-----------------------------------------------------------------------
!
!  initialize nullspace removal fields
!
!-----------------------------------------------------------------------

   if (sfc_layer_type == sfc_layer_varthick) then

      allocate (   CHECKER(nx_block,ny_block,max_blocks_clinic), &
                   CONSTNT(nx_block,ny_block,max_blocks_clinic), &
                CHECK_AREA(nx_block,ny_block,max_blocks_clinic), &
                CONST_AREA(nx_block,ny_block,max_blocks_clinic))

      !$OMP PARALLEL DO PRIVATE(iblock, this_block, i, j, n)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)  

         do j = 1,ny_block
         do i = 1,nx_block

            n = this_block%i_glob(i) + abs(this_block%j_glob(j))
            CHECKER(i,j,iblock) = 2*mod(n,2) - 1
            CHECK_AREA(i,j,iblock) = CHECKER(i,j,iblock)* &
                                       TAREA(i,j,iblock)

         enddo
         enddo

         if (allocated(REGION_MASK)) then
            where(KMT(:,:,iblock) > 0 .and. &
                  REGION_MASK(:,:,iblock) > 0) 
               CONSTNT(:,:,iblock) = 1
               CONST_AREA(:,:,iblock) = TAREA(:,:,iblock)
            elsewhere
               CHECKER(:,:,iblock) = 0
               CONSTNT(:,:,iblock) = 0
               CHECK_AREA(:,:,iblock) = 0
               CONST_AREA(:,:,iblock) = 0
            endwhere
         else
            where(KMT(:,:,iblock) > 0) 
               CONSTNT(:,:,iblock) = 1
               CONST_AREA(:,:,iblock) = TAREA(:,:,iblock)
            elsewhere
               CHECKER(:,:,iblock) = 0
               CONSTNT(:,:,iblock) = 0
               CHECK_AREA(:,:,iblock) = 0
               CONST_AREA(:,:,iblock) = 0
            endwhere
         endif
      end do
      !$OMP END PARALLEL DO

      sum_check = global_sum(CHECKER,distrb_clinic,field_loc_center)
      sum_const = global_sum(CONSTNT,distrb_clinic,field_loc_center)

      acheck = global_sum(CHECK_AREA,distrb_clinic,field_loc_center) &
              /global_sum(CONST_AREA,distrb_clinic,field_loc_center)

      rcheck = acheck/(sum_const - acheck*sum_check)
      rconst =     c1/(sum_const - acheck*sum_check)

      deallocate(CHECK_AREA, CONST_AREA)

   endif ! varthick

!-----------------------------------------------------------------------
!EOC

   end subroutine init_barotropic

!***********************************************************************
!BOP
! !IROUTINE: barotropic_driver
! !INTERFACE:

 subroutine barotropic_driver(ZX,ZY,errorCode)

! !DESCRIPTION:
!  This routine solves the barotropic equations for the surface 
!  pressure and barotropic velocity field using the implicit 
!  free-surface formulation.
!
!  For leapfrog steps, the momentum equations for the auxiliary 
!  velocity $(u',v')$ are
!  \begin{eqnarray}
!     (u' - u^{n-1}) - 2 \Delta t \alpha f(v' - v^{n-1}) 
!         &=& 2 \Delta t [F_x - \gamma \nabla_x p - 
!                            (1-\gamma)\nabla_x p^{n-1}] \\ 
!     (v' - v^{n-1}) + 2 \Delta t \alpha f(u' - u^{n-1}) 
!         &=& 2 \Delta t [F_y - \gamma \nabla_y p - 
!                            (1-\gamma)\nabla_y p^{n-1}].
!  \end{eqnarray}
!  The elliptic equation for new pressure $p^{n+1}$ is
!  \begin{equation}
!     \nabla\cdot(H \nabla(p^{n+1}) - p^{n+1}/(\alpha 2 \Delta t^2 g)
!           = \nabla\cdot(H(u',v')/(\alpha 2\Delta t) + \nabla p^{n-1})
!              - p^n/(\alpha 2\Delta t^2 g - F_w/(\alpha 2\Delta t)
!  \end{equation}
!  New velocities $(U^{n+1},V^{n+1})$ are then constructed using
!  \begin{eqnarray}
!     U^{n+1} & = & u' - \alpha 2\Delta t \nabla_x(p^{n+1} - p^{n-1}) \\
!     V^{n+1} & = & v' - \alpha 2\Delta t \nabla_y(p^{n+1} - p^{n-1})
!  \end{eqnarray}
!
!  On the first pass of Matsuno steps, the auxiliary velocity is
!  \begin{eqnarray}
!         (u' - u^n) - \Delta t\theta f(v' - v^n)
!           &=& \Delta t (F_x - \nabla_x p^n) \\
!         (v' - v^n) + \Delta t\theta f(u' - u^n)
!           &=& \Delta t (F_y - \nabla_y p^n),
!  \end{eqnarray}
!  the elliptic equation for new pressure is
!  \begin{equation}
!         \nabla\cdot(H\nabla({p`}^{n+1}) - 
!         {p`}^{n+1}/(\theta\Delta t^2 g)
!           = \nabla\cdot(H[(u',v')/(\theta\Delta t) + \nabla p^n])
!              - p^n/(\theta*\Delta t^2 g) - F_w/(\theta\Delta t),
!  \end{equation}
!  and the velocities are constructed using
!  \begin{eqnarray}
!         U^{n+1} &=& u' - \theta\Delta t \nabla_x(p^{n+1} - p) \\
!         V^{n+1} &=& v' - \theta\Delta t \nabla_y(p^{n+1} - p)
!  \end{eqnarray}
!
!  On the second pass of Matsuno steps, the auxiliary velocity is
!  \begin{eqnarray}
!         (u' - u^n) - \Delta t \theta f(v' - v^n)
!           &=& \Delta t (F_x' - \theta\nabla_x {p'}^{n+1} - 
!                             (1-\theta)\nabla_x p^n) \\
!         (v' - v^n) + \Delta t \theta f(u' - u^n)
!           &=& \Delta t (F_y' - \theta\nabla_y {p'}^{n+1} - 
!                             (1-\theta)\nabla_y p^n),
!  \end{eqnarray}
!  the elliptic equation for new pressure is
!  \begin{equation}
!         \nabla\cdot(H\nabla p^{n+1}) - p^{n+1}/(\theta\Delta t^2 g)
!           = \nabla\cdot(H[(u',v')/(\theta\Delta t) + \nabla p^{n+1}])
!              - p^n/(\theta\Delta t^2 g) - F_w/(\theta\Delta t)
!  \end{equation}
!  and the velocities are constructed using
!  \begin{eqnarray}
!         U^{n+1} &=& u' - \theta\Delta t\nabla_x(p^{n+1} - {p'}^{n+1})\\
!         V^{n+1} &=& v' - \theta\Delta t\nabla_y(p^{n+1} - {p'}^{n+1})
!  \end{eqnarray}
!
!  The above equations are written for the case of implicit treatment
!  of the coriolis terms.  The parameters $\alpha$ and $\gamma$ are
!  used to vary the time-centering of the coriolis terms and surface 
!  pressure gradients, which enter the equations centered in time
!  as follows:
!  \begin{eqnarray}
!     \alpha Q^{n+1} + \gamma Q^n + (1-\alpha-\gamma)Q^{n-1} & &
!                   {\rm (leapfrog\ timesteps)} \\
!     \theta Q^{n+1} + (1-\theta)Q^n & & {\rm (matsuno\ timesteps)}
!  \end{eqnarray}
!  where Q is the coriolis term or surface-pressure gradient. 
!  The force terms $(F_x,F_y)$ contain r.h.s. coriolis terms which
!  vary depending on the type of timestep (see comments in clinic).
!  If the coriolis terms are treated explicitly, then they are 
!  simply evaluated at time (n), and appear only in the force terms.
!
!  In the 2nd pass of a matsuno timestep, the force terms
!  $(F_x',F_y')$ are constructed in baroclinic using the predicted
!  prognostic fields $({U'}^{n+1},{V'}^{n+1}),T^{n+1}$ from the first 
!  pass.
!
!  The auxiliary velocities $(u',v')$ are solutions of the momentum
!  equation using an earlier known pressure.  The elliptic equation
!  solves for the correction to this pressure, and its gradient is
!  added to $(u',v')$ to obtain the new velocites $(U^{n+1},V^{n+1})$.
!
!  The elliptic equation is multiplied by the T-cell area to make
!  the operator matrix used by the cg solver symmetric.
!  $F_w$ is the surface freshwater flux in cm/s for the variable
!  thickness surface layer option.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout) :: &
      ZX, ZY             ! vertical integrals of forcing

! !INPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode          ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
   integer (int_kind) :: iblock

   real (r8) ::          &
      xcheck              ! global sum of checkerboard

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      RHS,               &! r.h.s. of elliptic eqn times T-cell area
      UH,VH,             &! auxiliary velocities (see description above)
      PCHECK,            &! array for computing null space
      WORKX,WORKY         ! local temp space

   real (r8), dimension(nx_block,ny_block) :: &
      diagonalCorrection,     &! time dependent correction to operator
      WORK1,WORK2,WORK3,WORK4  ! local work space

   type (block) ::     &
      this_block        ! block information for current block

!-----------------------------------------------------------------------
!
!  calculate r.h.s. of barotropic momentum equations.
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock, this_block, diagonalCorrection, &
   !$OMP                     WORK1, WORK2, WORK3, WORK4)

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)  

      if (leapfrogts) then               ! leapfrog

         WORK3 = c2dtp*(ZX(:,:,iblock) -                          &
                              gamma *GRADPX(:,:,curtime,iblock) - & 
                        (c1 - gamma)*GRADPX(:,:,oldtime,iblock))
         WORK4 = c2dtp*(ZY(:,:,iblock) -                          &
                              gamma *GRADPY(:,:,curtime,iblock) - &
                        (c1 - gamma)*GRADPY(:,:,oldtime,iblock))

      elseif (mix_pass == 1 .or. f_euler_ts) then   ! matsuno 1st pass

 
         WORK3 = c2dtp*(ZX(:,:,iblock) - GRADPX(:,:,curtime,iblock))
         WORK4 = c2dtp*(ZY(:,:,iblock) - GRADPY(:,:,curtime,iblock))
 
      else ! (mix_pass == 2)                ! matsuno 2nd pass
 

         WORK3 = c2dtp*(ZX(:,:,iblock) -                          &
                              theta *GRADPX(:,:,newtime,iblock) - &
                        (c1 - theta)*GRADPX(:,:,curtime,iblock))
         WORK4 = c2dtp*(ZY(:,:,iblock) -                          &
                              theta *GRADPY(:,:,newtime,iblock) - &
                        (c1 - theta)*GRADPY(:,:,curtime,iblock))

      endif

!-----------------------------------------------------------------------
!
!     calculate negative gradient of surface atmospheric pressure
!     and add it to r.h.s. forcing
!
!-----------------------------------------------------------------------

      if (ap_data_type /= 'none') then

         call grad(1, WORK1, WORK2, ATM_PRESS(:,:,iblock), this_block)

         WORK3 = WORK3 - c2dtp*WORK1
         WORK4 = WORK4 - c2dtp*WORK2
      endif

!-----------------------------------------------------------------------
!
!     solve for auxiliary velocities ([Uh],[Vh]) 
!
!-----------------------------------------------------------------------

      if (impcor) then          ! implicit coriolis

         WORK1 = c2dtp*beta*FCOR(:,:,iblock)
         WORK2 = c1/(c1 + WORK1**2)
         UH(:,:,iblock) = WORK2*(WORK3 + WORK1*WORK4) + & 
                          UBTROP(:,:,oldtime,iblock)
         VH(:,:,iblock) = WORK2*(WORK4 - WORK1*WORK3) + &
                          VBTROP(:,:,oldtime,iblock)

      else                      ! explicit coriolis

         UH(:,:,iblock) = WORK3 + UBTROP(:,:,oldtime,iblock)
         VH(:,:,iblock) = WORK4 + VBTROP(:,:,oldtime,iblock)

      endif

!-----------------------------------------------------------------------
!
!     calculate r.h.s. of elliptic equation
!
!-----------------------------------------------------------------------

      if (leapfrogts) then               ! leapfrog
 
         WORK3 = HU(:,:,iblock)*(UH(:,:,iblock) + &
                                 beta*c2dtp*GRADPX(:,:,oldtime,iblock))
         WORK4 = HU(:,:,iblock)*(VH(:,:,iblock) + &
                                 beta*c2dtp*GRADPY(:,:,oldtime,iblock))

      elseif (mix_pass == 1 .or. f_euler_ts) then   ! matsuno 1st pass
 
         WORK3 = HU(:,:,iblock)*(UH(:,:,iblock) + &
                                 beta*c2dtp*GRADPX(:,:,curtime,iblock))
         WORK4 = HU(:,:,iblock)*(VH(:,:,iblock) + &
                                 beta*c2dtp*GRADPY(:,:,curtime,iblock))
 
      else ! (mix_pass == 2)                ! matsuno 2nd pass
 
         WORK3 = HU(:,:,iblock)*(UH(:,:,iblock) + &
                                 beta*c2dtp*GRADPX(:,:,newtime,iblock))
         WORK4 = HU(:,:,iblock)*(VH(:,:,iblock) + & 
                                 beta*c2dtp*GRADPY(:,:,newtime,iblock))
 
      endif

      if ( overflows_on .and. overflows_interactive ) then
         call ovf_brtrpc_renorm(WORK3,WORK4,iblock)
      endif

      !*** div returns T-cell area * divergence
      call div(1,RHS(:,:,iblock),WORK3,WORK4,this_block)
      RHS(:,:,iblock) = RHS(:,:,iblock)/(beta*c2dtp)

      if ( overflows_on .and. overflows_interactive ) then
         call ovf_rhs_brtrpc_continuity(RHS,iblock)
      endif

!-----------------------------------------------------------------------
!
!     add diagonal term to central coefficient in implicit free-surface 
!     formulation, and add correction to r.h.s.
!
!-----------------------------------------------------------------------

      select case (sfc_layer_type)

      case(sfc_layer_varthick)
         diagonalCorrection(:,:) =                                     & 
                        merge(TAREA(:,:,iblock)/(beta*c2dtp*dtp*grav), &
                              c0,CALCT(:,:,iblock))
         RHS(:,:,iblock) = RHS(:,:,iblock) -                           &
                   diagonalCorrection(:,:)*PSURF(:,:,curtime,iblock) - &
                   FW(:,:,iblock)*TAREA(:,:,iblock)/(beta*c2dtp)

      case(sfc_layer_rigid)
         diagonalCorrection(:,:) = 0.0_POP_r8

      case(sfc_layer_oldfree)
         diagonalCorrection(:,:) =                                     &
                        merge(TAREA(:,:,iblock)/(beta*c2dtp*dtp*grav), &
                              c0,CALCT(:,:,iblock))
         RHS(:,:,iblock) = RHS(:,:,iblock) -                           &
                     diagonalCorrection(:,:)*PSURF(:,:,curtime,iblock)

      end select

      call POP_SolversDiagonal(diagonalCorrection, iblock, errorCode)

!-----------------------------------------------------------------------
!
!     initial guess for solver
!     in matsuno 2nd pass, temporarily store press gradient from
!     1st pass
!
!-----------------------------------------------------------------------

      PSURF(:,:,newtime,iblock) = PGUESS(:,:,iblock)
 
      if (mix_pass == 2) then
         WORKX(:,:,iblock) = GRADPX(:,:,newtime,iblock)
         WORKY(:,:,iblock) = GRADPY(:,:,newtime,iblock)
      endif

   end do ! block loop

   !$OMP END PARALLEL DO

   call POP_HaloUpdate(RHS,POP_haloClinic,  &
            POP_gridHorzLocCenter,          &
            POP_fieldKindScalar, errorCode, &
            fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_BarotropicDriver: error updating RHS halo')
      return
   endif

!-----------------------------------------------------------------------
!
!  solve elliptic equation for surface pressure
!
!-----------------------------------------------------------------------

   call POP_SolversRun(PSURF(:,:,newtime,:), RHS, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_BarotropicDriver: error in solver')
      return
   endif

!-----------------------------------------------------------------------
!
!  calculate global sum of checkerboard nullspace
!
!-----------------------------------------------------------------------

   if (sfc_layer_type == sfc_layer_varthick) then

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic
         PCHECK(:,:,iblock) = PSURF(:,:,newtime,iblock)* &
                            CHECKER(:,:,iblock)
      end do
      !$OMP END PARALLEL DO

      xcheck = global_sum(PCHECK,distrb_clinic,field_loc_center)
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)  

!-----------------------------------------------------------------------
!
!     remove checkerboard nullspace from solution
!
!-----------------------------------------------------------------------

      if (sfc_layer_type == sfc_layer_varthick) then

         PSURF(:,:,newtime,iblock) = PSURF(:,:,newtime,iblock) +       &
                                   CONSTNT(:,:,iblock)*rcheck*xcheck - &
                                   CHECKER(:,:,iblock)*rconst*xcheck
      endif

!-----------------------------------------------------------------------
!
!     calculate gradient of PSURF(:,:,newtime)
!
!-----------------------------------------------------------------------

      call grad(1,GRADPX(:,:,newtime,iblock), &
                  GRADPY(:,:,newtime,iblock), &
                   PSURF(:,:,newtime,iblock),this_block)

!-----------------------------------------------------------------------
!
!     update new barotropic velocity, pressure, pressure gradient
!
!-----------------------------------------------------------------------

      if (leapfrogts) then               ! leapfrog
 
         UBTROP(:,:,newtime,iblock) = UH(:,:,iblock) -                &
                             beta*c2dtp*(GRADPX(:,:,newtime,iblock) - & 
                                         GRADPX(:,:,oldtime,iblock))
         VBTROP(:,:,newtime,iblock) = VH(:,:,iblock) -                &
                             beta*c2dtp*(GRADPY(:,:,newtime,iblock) - &
                                         GRADPY(:,:,oldtime,iblock))
 
      elseif (mix_pass == 1 .or. f_euler_ts) then  ! matsuno 1st pass
 
         UBTROP(:,:,newtime,iblock) = UH(:,:,iblock) -                &
                             beta*c2dtp*(GRADPX(:,:,newtime,iblock) - &
                                         GRADPX(:,:,curtime,iblock))
         VBTROP(:,:,newtime,iblock) = VH(:,:,iblock) -                &
                             beta*c2dtp*(GRADPY(:,:,newtime,iblock) - &
                                         GRADPY(:,:,curtime,iblock))
 
      else ! (mix_pass == 2)                ! matsuno 2nd pass
 
         UBTROP(:,:,newtime,iblock) = UH(:,:,iblock) -                &
                             beta*c2dtp*(GRADPX(:,:,newtime,iblock) - &
                                         WORKX(:,:,iblock))
         VBTROP(:,:,newtime,iblock) = VH(:,:,iblock) -                &
                             beta*c2dtp*(GRADPY(:,:,newtime,iblock) - &
                                         WORKY(:,:,iblock))
 
      endif

!-----------------------------------------------------------------------
!
!     accumulate tavg diagnostics for barotropic velocities
!
!-----------------------------------------------------------------------

      call accumulate_tavg_field(HU(:,:,iblock)*             &
                                 UBTROP(:,:,curtime,iblock), &
                                 tavg_SU, iblock, 1)

      call accumulate_tavg_field(HU(:,:,iblock)*             &
                                 VBTROP(:,:,curtime,iblock), &
                                 tavg_SV, iblock, 1)

   end do ! block loop

   !$OMP END PARALLEL DO

   call POP_HaloUpdate(PSURF(:,:,newtime,:),POP_haloClinic,   &
                             POP_gridHorzLocCenter,           &
                             POP_fieldKindScalar, errorCode,  &
                             fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_BarotropicDriver: error updating PSURF halo')
      return
   endif

   call POP_HaloUpdate(GRADPX(:,:,newtime,:),POP_haloClinic,  &
                             POP_gridHorzLocNECorner,         &
                             POP_fieldKindVector, errorCode,  &
                             fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_BarotropicDriver: error updating GRADPX halo')
      return
   endif

   call POP_HaloUpdate(GRADPY(:,:,newtime,:),POP_haloClinic,  &
                             POP_gridHorzLocNECorner,         &
                             POP_fieldKindVector, errorCode,  &
                             fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_BarotropicDriver: error updating GRADPY halo')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine barotropic_driver

!***********************************************************************

 end module barotropic

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
