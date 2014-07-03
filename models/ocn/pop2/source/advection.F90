!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module advection

!BOP
! !MODULE: advection
!
! !DESCRIPTION:
!  This module contains arrays and variables necessary for performing
!  advection of momentum and tracer quantities.  Currently, the
!  module supports leapfrog centered advection of momentum and
!  both leapfrog centered advection and third-order upwinding of
!  tracers.
!
! !REVISION HISTORY:
!  SVN:$Id: advection.F90 47361 2013-05-21 20:54:30Z mlevy@ucar.edu $

! !USES:
   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use kinds_mod, only: r4, r8, int_kind, char_len, log_kind, rtavg
   use constants, only: c0, c1, p5, p125, p25, blank_fmt, delim_fmt, ndelim_fmt, c2, &
       field_loc_center, field_type_scalar, &
       field_loc_Eface, field_type_vector
   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size
   use communicate, only: my_task, master_task
   use distribution, only: 
   use grid, only: dz, DXT, DYT, HUW, HUS, c2dz, KMT, HTE, UAREA_R, DZT,    &
       partial_bottom_cells, DYU, DZU, DXU, DZR, DZ2R, KMU, TAREA_R, HTN,   &
       sfc_layer_type, sfc_layer_varthick, FCORT, KMTE,                     &
       KMTW, KMTEE, KMTN, KMTS, KMTNN, ugrid_to_tgrid
   use domain, only: nblocks_clinic, blocks_clinic, distrb_clinic,          &
       POP_haloClinic
   use broadcast, only: broadcast_scalar, broadcast_array
   use diagnostics, only: cfl_advect
   use state_mod, only: state
   use operators, only: zcurl
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use io_types, only: nml_in, nml_filename, stdout
   use time_management, only: max_blocks_clinic, km, nt, mix_pass, c2dtt
   use timers, only: timer_start, timer_stop, get_timer
   use exit_mod, only: sigAbort, exit_pop, flushm
   use prognostic, only: UVEL, VVEL, curtime, tracer_d
   use passive_tracers, only: tadvect_ctype_passive_tracers
   use registry
   use overflows
   use overflow_type

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: init_advection, &
             comp_flux_vel_ghost, &
             advt,           &
             advu

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  choices for tracer advection
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      tadvect_centered = 1,  &! centered leap-frog tracer advection
      tadvect_upwind3  = 2,  &! 3rd-order upwind tracer advection
      tadvect_lw_lim   = 3    ! 1d Lax-Wendroff with 1d flux limiters

   integer (int_kind), dimension(nt) :: &
      tadvect_itype           ! users tracer advection choice

!-----------------------------------------------------------------------
!
!  coefficients for metric advection terms (KXU,KYU)
!
!-----------------------------------------------------------------------

   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) :: & 
      KXU,KYU

!-----------------------------------------------------------------------
!
!  geometric arrays necessary for non-centered advection
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable :: &
      AUX

!-----------------------------------------------------------------------
!
!  geometric arrays necessary for upwind advection
!
!-----------------------------------------------------------------------

   real (r8), dimension(:), allocatable :: &
      talfzp,tbetzp,tgamzp,                &
      talfzm,tbetzm,tdelzm

   real (r8), dimension(:,:,:), allocatable ::   &
      TALFXP,TBETXP,TGAMXP,TALFYP,TBETYP,TGAMYP, &
      TALFXM,TBETXM,TDELXM,TALFYM,TBETYM,TDELYM

!-----------------------------------------------------------------------
!
!     geometric arrays necessary for lw_lim advection
!
!-----------------------------------------------------------------------

   real (r8), dimension(:), allocatable :: &
      p5_dz_ph_r

   real (r8), dimension(:,:,:), allocatable :: &
      p5_DXT_ph_R,        &! 1/(DXT(i,j)+DXT(i+1,j))
      p5_DYT_ph_R,        &! 1/(DYT(i,j)+DYT(i,j+1))
      UTE_jbm2,           &! UTE for j==jb-2
      WTKB_jbm2,          &! WTKB for j==jb-2
      WTKB_jep2,          &! WTKB for j==je+2
      WTKB_ibm2,          &! WTKB for i==ib-2
      WTKB_iep2            ! WTKB for i==ie+2

   real (r8), dimension(:,:,:,:), allocatable :: &
      FLUX_VEL_prev,      &! flux velocities from prev k iteration
      UTE_to_UVEL_E,      &! converts UTE to UVEL_E
      VTN_to_VVEL_N        ! converts VTN to VVEL_E

!-----------------------------------------------------------------------
!
!  tavg ids for tavg diagnostics related to advection
!  north, east, zonal, merid refer here to logical space only
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_WVEL,         &! Vertical Velocity
      tavg_WVEL2,        &! Vertical Velocity Squared
      tavg_UEU,          &! flux of zonal momentum across east  face
      tavg_VNU,          &! flux of zonal momentum across north face
      tavg_WTU,          &! flux of zonal momentum across top   face
      tavg_UEV,          &! flux of merid momentum across east  face
      tavg_VNV,          &! flux of merid momentum across north face
      tavg_WTV,          &! flux of merid momentum across top   face
      tavg_PV,           &! potential vorticity
      tavg_Q,            &! z-derivative of pot density
      tavg_PD,           &! potential density 
      tavg_RHOU,         &! pot density times U velocity
      tavg_RHOV,         &! pot density times V velocity
      tavg_PVWM,         &! pot vorticity flux through bottom
      tavg_PVWP,         &! pot vorticity flux through top
      tavg_UPV,          &! pot vorticity flux through east  face
      tavg_VPV,          &! pot vorticity flux through north face
      tavg_URHO,         &! pot density   flux through east  face
      tavg_VRHO,         &! pot density   flux through north face
      tavg_WRHO,         &! pot density   flux through top   face
      tavg_UQ,           &! advection of Q across east  face
      tavg_VQ             ! advection of Q across north face

   integer (int_kind), dimension(nt) :: &
      tavg_ADV_TRACER,   &! vertical average of tracer advective tendency
      tavg_UE_TRACER,    &! flux of tracer across east  face
      tavg_VN_TRACER,    &! flux of tracer across north face
      tavg_WT_TRACER      ! flux of tracer across top   face

!-----------------------------------------------------------------------
!
!  advection timers
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      timer_advt,        &! timer for tracer   advection
      timer_advu          ! timer for momentum advection

!-----------------------------------------------------------------------
!
!  advection flags
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      luse_centered,       &! is centered used by any tracer
      luse_upwind3,        &! is upwind3 used by any tracer
      luse_lw_lim           ! is lw_lim used by any tracer

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_advection
! !INTERFACE:

 subroutine init_advection 

! !DESCRIPTION:
!  This subroutine initializes variables associated with advection
!  schemes.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  namelist input
!
!-----------------------------------------------------------------------

   character (*), parameter :: &
      docfmt = "( 3x, a16, 3x, a16 )"

   character (char_len) :: &
      tadvect_ctype        ! character string for tracer advect choice

   namelist /advect_nml/ tadvect_ctype

!-----------------------------------------------------------------------
!
!  local variables for setting up upwind coefficients
!
!-----------------------------------------------------------------------

   integer (int_kind) :: & 
      i,j,k,n,           &! dummy loop indices
      iblock,            &! local block number
      nattempts,         &! num of attempts to read namelist input
      nml_error           ! error flag for namelist read

   real (r8) ::            &
      dxc,dxcw,dxce,dxce2, &
      dyc,dycs,dycn,dycn2          

   real (r8), dimension(:), allocatable :: & 
      dzc

   real (r8), dimension(nx_block,ny_block) :: & 
     WORK1                ! local temp space

   type (block) ::        &
      this_block          ! block information for current block

!-----------------------------------------------------------------------
!
!  check that init_passive_tracers has been called
!
!-----------------------------------------------------------------------

   if (.not. registry_match('init_passive_tracers')) then
     call exit_POP(sigAbort,'ERROR: init_passive_tracers must be ' /&
                         &/ 'called before init_advection')
   end if

!-----------------------------------------------------------------------
!
!  read input namelist for choice of tracer advection method
!
!-----------------------------------------------------------------------

   tadvect_ctype = 'centered'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
        nml_error = -1
      else
        nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
        read(nml_in, nml=advect_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading advection namelist')
   endif
     
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Advection:'
      write(stdout,blank_fmt)
      write(stdout,*) ' advect_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,advect_nml)
      write(stdout,blank_fmt)

      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a17)') 'Advection options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      tadvect_itype = tadvect_ctype_to_tadvect_itype(tadvect_ctype)

      do n=3,nt
         if (tadvect_ctype_passive_tracers(n) == 'base_model') then
            tadvect_ctype_passive_tracers(n) = tadvect_ctype
         else
            tadvect_itype(n) = tadvect_ctype_to_tadvect_itype( &
               tadvect_ctype_passive_tracers(n))
         endif
      end do

      if (all(tadvect_itype == tadvect_itype(1))) then
         write(stdout,*) 'Using ' /&
            &/ trim(tadvect_ctype) /&
            &/ ' for all tracers.'
      else
         write(stdout,docfmt) 'Tracer Name', 'Advection Option'
         write(stdout,docfmt) 'TEMP', trim(tadvect_ctype)
         write(stdout,docfmt) 'SALT', trim(tadvect_ctype)
         do n=3,nt
            write(stdout,docfmt) trim(tracer_d(n)%short_name), &
                                 trim(tadvect_ctype_passive_tracers(n))
         end do
      end if
   endif

   call broadcast_array(tadvect_itype, master_task)

   if (minval(tadvect_itype) < 0) then
      call exit_POP(sigAbort,'ERROR: unknown tracer advection method')
   endif

!-----------------------------------------------------------------------
!
!  advection type must be implemented by overflows if on and interactive
!
!-----------------------------------------------------------------------

   if ( overflows_on .and. overflows_interactive ) then
      do n=1,nt
         if( tadvect_itype(n) == tadvect_centered ) then
            write(stdout,*) &
              'advection: ERROR: advection type not implemented ', &
                                'for interactive overflows'
            write(stdout,*) ' advection type = ',tadvect_itype(n), &
                            ' for tracer = ',n
            write(stdout,*) ' only upwind3 and lw_lim implemented'
            call exit_POP(sigAbort,'ERROR non-implemented advection with overflows')
         endif
      enddo
   endif

!-----------------------------------------------------------------------
!
!  set use flags
!
!-----------------------------------------------------------------------

   luse_centered = any(tadvect_itype == tadvect_centered)
   luse_upwind3  = any(tadvect_itype == tadvect_upwind3)
   luse_lw_lim   = any(tadvect_itype == tadvect_lw_lim)

!-----------------------------------------------------------------------
!
!  initialize metric advection coefficients
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)  

      KXU(:,:,iblock) = (eoshift(HUW(:,:,iblock),dim=1,shift=1) - & 
                                 HUW(:,:,iblock))*UAREA_R(:,:,iblock)
      KYU(:,:,iblock) = (eoshift(HUS(:,:,iblock),dim=2,shift=1) - &
                                 HUS(:,:,iblock))*UAREA_R(:,:,iblock)

   end do

   !*** KXU,KYU only needed in physical domain
   !*** no ghost cell update required assuming HUS,HUW,UAREA_R
   !*** were defined correctly in ghost cells

!-----------------------------------------------------------------------
!
!  allocate and initialize non-centered arrays if necessary
!
!-----------------------------------------------------------------------

   if (any(tadvect_itype /= tadvect_centered)) then

      allocate (AUX   (nx_block,ny_block,nt,nblocks_clinic))

   endif ! non-centered setup

!-----------------------------------------------------------------------
!
!  allocate and initialize upwinding grid arrays if necessary
!
!-----------------------------------------------------------------------

   if (luse_upwind3) then

      allocate (talfzp(km), &
                tbetzp(km), &
                tgamzp(km), &
                talfzm(km), &
                tbetzm(km), &
                tdelzm(km))

      allocate (TALFXP(nx_block,ny_block,nblocks_clinic), &
                TBETXP(nx_block,ny_block,nblocks_clinic), &
                TGAMXP(nx_block,ny_block,nblocks_clinic), &
                TALFYP(nx_block,ny_block,nblocks_clinic), &
                TBETYP(nx_block,ny_block,nblocks_clinic), &
                TGAMYP(nx_block,ny_block,nblocks_clinic), &
                TALFXM(nx_block,ny_block,nblocks_clinic), &
                TBETXM(nx_block,ny_block,nblocks_clinic), &
                TDELXM(nx_block,ny_block,nblocks_clinic), &
                TALFYM(nx_block,ny_block,nblocks_clinic), &
                TBETYM(nx_block,ny_block,nblocks_clinic), &
                TDELYM(nx_block,ny_block,nblocks_clinic))

      allocate (dzc(0:km+1))

      !***
      !*** vertical grid coefficients     
      !***

      dzc(0)    = dz(1)
      do k=1,km
         dzc(k)  = dz(k)
      enddo
      dzc(km+1) = dzc(km)

      do k=1,km-1
         talfzp(k) =  dz(k)*(c2*dz(k)+dzc(k-1))/ &
                      ((dz(k)+dz(k+1))*          &
                      (dzc(k-1)+c2*dz(k)+dz(k+1)))
         tbetzp(k) =  dz(k+1)*(c2*dz(k)+dzc(k-1))/ &
                      ((dz(k)+dz(k+1))*            &
                      (dz(k)+dzc(k-1)          ))
         tgamzp(k) =  -(dz(k)*dz(k+1))/ &
                      ((dz(k)+dzc(k-1))*           &
                      (dz(k+1)+dzc(k-1)+c2*dz(k)))
      enddo 
      tbetzp(1) = tbetzp(1) + tgamzp(1)
      tgamzp(1) = c0
      talfzp(km) = c0
      tbetzp(km) = c0
      tgamzp(km) = c0
 
      do k=1,km-1
         talfzm(k) =  dz(k)*(c2*dz(k+1)+dzc(k+2))/ &
                      ((dz(k)+dz(k+1))*            &
                      (dz(k+1)+dzc(k+2)            ))
         tbetzm(k) =  dz(k+1)*(c2*dz(k+1)+dzc(k+2))/ &
                      ((dz(k)+dz(k+1))*              &
                      (dz(k)+dzc(k+2)+c2*dz(k+1)))
         tdelzm(k) =  -(dz(k)*dz(k+1))/ &
                      ((dz(k+1)+dzc(k+2))*         &
                      (dz(k)+dzc(k+2)+c2*dz(k+1)))
      enddo    
      talfzm(km-1) = talfzm(km-1) + tdelzm(km-1)
      tdelzm(km-1) = c0
      talfzm(km) = c0
      tbetzm(km) = c0
      tdelzm(km) = c0

      deallocate (dzc)

      !***
      !*** horizontal grid coeffs
      !***

      do iblock = 1,nblocks_clinic

         this_block = get_block(blocks_clinic(iblock),iblock)  

         !***
         !*** zonal grid coefficients     
         !***

         do j=this_block%jb,this_block%je
         do i=this_block%ib-1,this_block%ie

            dxc   = DXT(i  ,j,iblock)
            dxcw  = DXT(i-1,j,iblock)
            dxce  = DXT(i+1,j,iblock)
            dxce2 = DXT(i+2,j,iblock)

            TALFXP(i,j,iblock) = dxc*(c2*dxc+dxcw)/ &
                                 ((dxc+dxce)*(dxcw+c2*dxc+dxce))
            TBETXP(i,j,iblock) = dxce*(c2*dxc+dxcw)/ &
                                 ((dxc+dxcw)*(        dxc+dxce))
            TGAMXP(i,j,iblock) =     -(   dxc*dxce)/ &
                                 ((dxc+dxcw)*(dxcw+c2*dxc+dxce))

            TALFXM(i,j,iblock) = dxc *(c2*dxce+dxce2)/ &
                                 ((dxc  +dxce)*(       dxce+dxce2))
            TBETXM(i,j,iblock) = dxce*(c2*dxce+dxce2)/ &
                                 ((dxc  +dxce)*(dxc+c2*dxce+dxce2))
            TDELXM(i,j,iblock) =     -(   dxc *dxce )/ &
                                 ((dxce2+dxce)*(dxc+c2*dxce+dxce2))

         end do
         end do

         !***
         !*** poloidal grid coefficients     
         !***

         do j=this_block%jb-1,this_block%je
         do i=this_block%ib,this_block%ie

            dyc   = DYT(i,j  ,iblock)
            dycs  = DYT(i,j-1,iblock)
            dycn  = DYT(i,j+1,iblock)
            dycn2 = DYT(i,j+2,iblock)

            TALFYP(i,j,iblock) = dyc *(c2*dyc+dycs)/ &
                                 ((dyc+dycn)*(dycs+c2*dyc+dycn))
            TBETYP(i,j,iblock) = dycn*(c2*dyc+dycs)/ &
                                 ((dyc+dycn)*(dycs+   dyc     ))
            TGAMYP(i,j,iblock) =     -(   dyc*dycn)/ &
                                 ((dyc+dycs)*(dycs+c2*dyc+dycn))
 
            TALFYM(i,j,iblock) = dyc *(c2*dycn+dycn2)/ &
                                 ((dyc+dycn)*(       dycn+dycn2))
            TBETYM(i,j,iblock) = dycn*(c2*dycn+dycn2)/ &
                                 ((dyc+dycn)*(dyc+c2*dycn+dycn2))
            TDELYM(i,j,iblock) =     -(   dyc *dycn )/ &
                                 ((dycn2+dycn)*(dyc+c2*dycn+dycn2))

         end do
         end do

      end do

      !*** assuming DXT,DYT were defined correctly in ghost cells
      !*** these are valid from (ib-1:ie,jb-1:je) and are only
      !*** accessed in that range so no halo update necessary

   endif ! 3rd order upwind setup

!-----------------------------------------------------------------------
!
!  allocate and initialize lw_lim grid arrays if necessary
!
!-----------------------------------------------------------------------

   if (luse_lw_lim) then

      allocate (p5_dz_ph_r(km), &
                p5_DXT_ph_R(nx_block,ny_block,nblocks_clinic), &
                p5_DYT_ph_R(nx_block,ny_block,nblocks_clinic), &
                UTE_jbm2(nx_block,km,nblocks_clinic), &
                WTKB_jbm2(nx_block,km,nblocks_clinic), &
                WTKB_jep2(nx_block,km,nblocks_clinic), &
                WTKB_ibm2(ny_block,km,nblocks_clinic), &
                WTKB_iep2(ny_block,km,nblocks_clinic), &
                FLUX_VEL_prev(nx_block,ny_block,5,nblocks_clinic) )

      if (partial_bottom_cells) then

         allocate (UTE_to_UVEL_E(nx_block,ny_block,km,nblocks_clinic), &
                   VTN_to_VVEL_N(nx_block,ny_block,km,nblocks_clinic))

      else

         allocate (UTE_to_UVEL_E(nx_block,ny_block,1,nblocks_clinic), &
                   VTN_to_VVEL_N(nx_block,ny_block,1,nblocks_clinic))

      endif

      p5_DXT_ph_R   = c0
      p5_DYT_ph_R   = c0
      UTE_jbm2      = c0
      WTKB_jbm2     = c0
      WTKB_jep2     = c0
      WTKB_ibm2     = c0
      WTKB_iep2     = c0
      UTE_to_UVEL_E = c0
      VTN_to_VVEL_N = c0

      p5_dz_ph_r(1:km-1) = c1 / (dz(1:km-1) + dz(2:km))
      p5_dz_ph_r(km)     = p5 / dz(km)

      !***
      !*** horizontal grid coeffs
      !***

      do iblock = 1,nblocks_clinic

         this_block = get_block(blocks_clinic(iblock),iblock)

         !***
         !*** zonal grid coefficients
         !***    compute in j halo because they are used
         !***

         do j=this_block%jb-2,this_block%je+2
         do i=this_block%ib-2,this_block%ie+1

            p5_DXT_ph_R(i,j,iblock) = c1 / (DXT(i,j,iblock) + &
                                            DXT(i+1,j,iblock))

         end do
         end do

         if (partial_bottom_cells) then

            do k=1,km
            do j=this_block%jb-2,this_block%je+2
            do i=this_block%ib-2,this_block%ie+1

               UTE_to_UVEL_E(i,j,k,iblock) = c1 / HTE(i,j,iblock) / &
                  min(DZT(i,j,k,iblock), DZT(i+1,j,k,iblock))

            end do
            end do
            end do

         else

            do j=this_block%jb-2,this_block%je+2
            do i=this_block%ib-2,this_block%ie+1

               UTE_to_UVEL_E(i,j,1,iblock) = c1 / HTE(i,j,iblock)

            end do
            end do

         endif

         !***
         !*** poloidal grid coefficients     
         !***

         do j=this_block%jb-2,this_block%je+1
         do i=this_block%ib,this_block%ie

            p5_DYT_ph_R(i,j,iblock) = c1 / (DYT(i,j,iblock) + &
                                            DYT(i,j+1,iblock))

         end do
         end do

         if (partial_bottom_cells) then

            do k=1,km
            do j=this_block%jb-2,this_block%je+1
            do i=this_block%ib,this_block%ie

               VTN_to_VVEL_N(i,j,k,iblock) = c1 / HTN(i,j,iblock) / &
                  min(DZT(i,j,k,iblock), DZT(i,j+1,k,iblock))

            end do
            end do
            end do

         else

            do j=this_block%jb-2,this_block%je+1
            do i=this_block%ib,this_block%ie

               VTN_to_VVEL_N(i,j,1,iblock) = c1 / HTN(i,j,iblock)

            end do
            end do

         endif

      end do

   endif ! lw_lim setup

!-----------------------------------------------------------------------
!
!  initialize timers
!
!-----------------------------------------------------------------------

   call get_timer(timer_advu,'ADVECTION_MOMENTUM', nblocks_clinic, &
                                                   distrb_clinic%nprocs)
 
   if (all(tadvect_itype == tadvect_centered)) then
      call get_timer(timer_advt,'ADVECTION_TRACER_CENTERED', &
                                 nblocks_clinic, distrb_clinic%nprocs)
   else if (all(tadvect_itype == tadvect_upwind3)) then
      call get_timer(timer_advt,'ADVECTION_TRACER_UPWIND3', &
                                 nblocks_clinic, distrb_clinic%nprocs)
   else if (all(tadvect_itype == tadvect_lw_lim)) then
      call get_timer(timer_advt,'ADVECTION_TRACER_LW_LIM', &
                                 nblocks_clinic, distrb_clinic%nprocs)
   else
      call get_timer(timer_advt,'ADVECTION_TRACER', &
                                 nblocks_clinic, distrb_clinic%nprocs)
   endif

!-----------------------------------------------------------------------
!
!  define tavg fields related to advection
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_WVEL,'WVEL',3,                          &
                          long_name='Vertical Velocity',               &
                          units='centimeter/s', grid_loc='3112',       &
                          coordinates='TLONG TLAT z_w time')

   call define_tavg_field(tavg_WVEL2,'WVEL2',3,                          &
                          long_name='Vertical Velocity**2',       &
                          units='centimeter^2/s^2', grid_loc='3112',       &
                          coordinates='TLONG TLAT z_w time')

   call define_tavg_field(tavg_UEU,'UEU',3,                            &
                          long_name='East Flux of Zonal Momentum',     &
                          units='cm/s^2', grid_loc='3321')

   call define_tavg_field(tavg_VNU,'VNU',3,                            &
                          long_name='North Flux of Zonal Momentum',    &
                          units='cm/s^2', grid_loc='3231')

   call define_tavg_field(tavg_UEV,'UEV',3,                            &
                          long_name='East Flux of Meridional Momentum',&
                          units='cm/s^2', grid_loc='3321')

   call define_tavg_field(tavg_VNV,'VNV',3,                            &
                        long_name='North Flux of Meridional Momentum', &
                          units='cm/s^2', grid_loc='3231')

   call define_tavg_field(tavg_WTU,'WTU',3,                            &
                          long_name='Top flux of Zonal Momentum',      &
                          units='cm/s^2', grid_loc='3222')

   call define_tavg_field(tavg_WTV,'WTV',3,                            &
                          long_name='Top flux of Meridional Momentum', &
                          units='cm/s^2', grid_loc='3222')

   call define_tavg_field(tavg_UE_TRACER(1),'UET',3,                   &
                          long_name='Flux of Heat in grid-x direction',&
                          units='degC/s', grid_loc='3211',             &
                          coordinates='ULONG TLAT z_t time' )

   call define_tavg_field(tavg_VN_TRACER(1),'VNT',3,                   &
                          long_name='Flux of Heat in grid-y direction',&
                          units='degC/s', grid_loc='3121',             &
                          coordinates='TLONG ULAT z_t time')

   call define_tavg_field(tavg_WT_TRACER(1),'WTT',3,                   &
                          long_name='Heat Flux Across Top Face',       &
                          units='degC/s', grid_loc='3112',             &
                          coordinates='TLONG TLAT z_w time' )

   call define_tavg_field(tavg_UE_TRACER(2),'UES',3,                   &
                          long_name='Salt Flux in grid-x direction',   &
                          scale_factor=1000.0_rtavg,                   &
                          units='gram/kilogram/s', grid_loc='3211',    &
                          coordinates='ULONG TLAT z_t time' )

   call define_tavg_field(tavg_VN_TRACER(2),'VNS',3,                   &
                          long_name='Salt Flux in grid-y direction',   &
                          scale_factor=1000.0_rtavg,                   &
                          units='gram/kilogram/s', grid_loc='3121',    &
                          coordinates='TLONG ULAT z_t time')

   call define_tavg_field(tavg_WT_TRACER(2),'WTS',3,                   &
                          long_name='Salt Flux Across Top Face',       &
                          scale_factor=1000.0_rtavg,                   &
                          units='gram/kilogram/s', grid_loc='3112',    &
                          coordinates='TLONG TLAT z_w time' )

   call define_tavg_field(tavg_ADV_TRACER(1),'ADVT',2,                     &
                    long_name='Vertically-Integrated T Advection Tendency',&
                          units='centimeter degC/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_ADV_TRACER(2),'ADVS',2,                        &
                    long_name='Vertically-Integrated S Advection Tendency',   &
                          scale_factor=1000.0_rtavg,                             &
                          units='centimeter gram/kilogram/s', grid_loc='2110',&
                          coordinates='TLONG TLAT time')

   do n=3,nt

      call define_tavg_field(tavg_UE_TRACER(n),                        &
                             'UE_' /&
                                    &/ trim(tracer_d(n)%short_name),3, &
                             long_name=trim(tracer_d(n)%short_name)   /&
                                    &/ ' Flux in grid-x direction',    &
                             units=trim(tracer_d(n)%tend_units),       &
                             scale_factor=tracer_d(n)%scale_factor,    &
                             grid_loc='3211',                          &
                             coordinates='ULONG TLAT z_t time' )

      call define_tavg_field(tavg_VN_TRACER(n),                        &
                             'VN_' /&
                                    &/ trim(tracer_d(n)%short_name),3, &
                             long_name=trim(tracer_d(n)%short_name)   /&
                                    &/ ' Flux in grid-y direction',    &
                             units=trim(tracer_d(n)%tend_units),       &
                             scale_factor=tracer_d(n)%scale_factor,    &
                             grid_loc='3121',                          &
                             coordinates='TLONG ULAT z_t time')

      call define_tavg_field(tavg_WT_TRACER(n),                        &
                             'WT_' /&
                                    &/ trim(tracer_d(n)%short_name),3, &
                             long_name=trim(tracer_d(n)%short_name)   /&
                                    &/ ' Flux Across Top Face',        &
                             units=trim(tracer_d(n)%tend_units),       &
                             scale_factor=tracer_d(n)%scale_factor,    &
                             grid_loc='3112',                          &
                             coordinates='TLONG TLAT z_w time' )

      call define_tavg_field(tavg_ADV_TRACER(n),                       &
                             'ADV_' /&
                                     &/ trim(tracer_d(n)%short_name),2,&
                             long_name='Vertically-Integrated '       /&
                                     &/trim(tracer_d(n)%short_name)   /&
                                     &/' Advection Tendency',          &
                             units=trim(tracer_d(n)%flux_units),       &
                             scale_factor=tracer_d(n)%scale_factor,    &
                             grid_loc='2110',                          &
                             coordinates='TLONG TLAT time' )

   end do

   call define_tavg_field(tavg_PV,'PV',3,                              &
                          long_name='Potential Vorticity',             &
                          units='1/s', grid_loc='3111',                &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_Q,'Q',3,                                &
                        long_name='Static Stability (d(rho(p_r))/dz)', &
                          units='gram/centimeter^4', grid_loc='3111',  &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_PD,'PD',3,                              &
                          long_name='Potential Density Ref to Surface',&
                          units='gram/centimeter^3', grid_loc='3111',  &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_RHOU,'RHOU',3,                          &
                          long_name='U times potential density',       &
                          units='g/cm^2/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_RHOV,'RHOV',3,                          &
                          long_name='V times potential density',       &
                          units='g/cm^2/s', grid_loc='3111',           &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_URHO,'URHO',3,                          &
                     long_name='flux of pot density across east face', &
                          units='g/cm^2/s', grid_loc='3321')

   call define_tavg_field(tavg_VRHO,'VRHO',3,                          &
                    long_name='flux of pot density across north face', &
                          units='g/cm^2/s', grid_loc='3231')

   call define_tavg_field(tavg_WRHO,'WRHO',3,                          &
                      long_name='flux of pot density across top face', &
                          units='g/cm^2/s', grid_loc='3112',           &
                          coordinates='TLONG TLAT z_w time' )

   call define_tavg_field(tavg_PVWM,'PVWM',3,                          &
                   long_name='flux of pot vorticity through top face', &
                          units='cm/s^2', grid_loc='3112',             &
                          coordinates='TLONG TLAT z_w time' )

   call define_tavg_field(tavg_PVWP,'PVWP',3,                          &
                long_name='flux of pot vorticity through bottom face', &
                          units='cm/s^2', grid_loc='3111',             &
                          coordinates='TLONG TLAT z_t time')

   call define_tavg_field(tavg_UPV,'UPV',3,                            &
                  long_name='flux of pot vorticity through east face', &
                          units='cm/s^2', grid_loc='3321')

   call define_tavg_field(tavg_VPV,'VPV',3,                            &
                 long_name='flux of pot vorticity through north face', &
                          units='cm/s^2', grid_loc='3231')

   call define_tavg_field(tavg_UQ,'UQ',3,                              &
                          long_name='flux of Q through east face',     &
                          units='g/cm^3/s', grid_loc='3321')

   call define_tavg_field(tavg_VQ,'VQ',3,                              &
                          long_name='flux of Q through north face',    &
                          units='g/cm^3/s', grid_loc='3231')

   call flushm (stdout)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_advection

!***********************************************************************
!BOP
! !IROUTINE: tadvect_ctype_to_tadvect_itype
! !INTERFACE:

 function tadvect_ctype_to_tadvect_itype(tadvect_ctype)

! !DESCRIPTION:
!  map an advection ctype to an itype
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(len=char_len), intent(in) :: &
      tadvect_ctype      ! character string for tracer advect choice

! !OUTPUT PARAMETERS:

   integer(kind=int_kind) :: tadvect_ctype_to_tadvect_itype

!-----------------------------------------------------------------------

   select case (tadvect_ctype(1:6))
   case ('center')
      tadvect_ctype_to_tadvect_itype = tadvect_centered
   case ('upwind')
      tadvect_ctype_to_tadvect_itype = tadvect_upwind3
   case ('lw_lim')
      tadvect_ctype_to_tadvect_itype = tadvect_lw_lim
   case default
      tadvect_ctype_to_tadvect_itype = -1000
   end select

!-----------------------------------------------------------------------
!EOC

 end function tadvect_ctype_to_tadvect_itype

!***********************************************************************
!BOP
! !IROUTINE: comp_flux_vel_ghost
! !INTERFACE:

 subroutine comp_flux_vel_ghost(DH, errorCode)

! !DESCRIPTION:
!  Compute and store tracer flux velocities in ghost cells.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(in) :: &
      DH                   ! change in surface height at T points

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!         
!  local variables:
!
!-----------------------------------------------------------------------

   type (block) ::        &
      this_block          ! block information for current block

   integer (int_kind) ::  &
      iblock,             &! local block number
      i,j,                &! dummy indices for horizontal directions
      k,n                  ! dummy indices for vertical level, tracer

   real (r8), dimension(nx_block,ny_block,nblocks_clinic) :: & 
      UTE,UTW,VTN,VTS,    &! tracer flux velocities across E,W,N,S faces
      WTK,                &! vertical velocity at bottom of T box
      WTKB                 ! vertical velocity at bottom of T box

!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. luse_lw_lim) return

   do k = 1,km

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)  
         if (k == 1) then
            WTK(:,:,iblock) = DH(:,:,iblock)
         else
            WTK(:,:,iblock) = WTKB(:,:,iblock)
         endif
         call comp_flux_vel(k,UVEL(:,:,:,curtime,iblock),&
                            VVEL(:,:,:,curtime,iblock),&
                            WTK(:,:,iblock),&
                            UTE(:,:,iblock),UTW(:,:,iblock),&
                            VTN(:,:,iblock),VTS(:,:,iblock),&
                            WTKB(:,:,iblock),this_block)
      end do

      !$OMP END PARALLEL DO

      call POP_HaloUpdate(UTE, POP_haloClinic,   &
                POP_gridHorzLocEFace,            &
                POP_fieldKindVector, errorCode,  &
                fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'comp_flux_vel_ghost: error updating halo for UTE')
         return
      endif

      call POP_HaloUpdate(WTKB, POP_haloClinic,   &
                POP_gridHorzLocCenter,            &
                POP_fieldKindScalar, errorCode,   &
                fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'comp_flux_vel_ghost: error updating halo for WTKB')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)  
         UTE_jbm2(:,k,iblock)  = UTE(:,this_block%jb-2,iblock)
         WTKB_jbm2(:,k,iblock) = WTKB(:,this_block%jb-2,iblock)
         WTKB_jep2(:,k,iblock) = WTKB(:,this_block%je+2,iblock)
         WTKB_ibm2(:,k,iblock) = WTKB(this_block%ib-2,:,iblock)
         WTKB_iep2(:,k,iblock) = WTKB(this_block%ie+2,:,iblock)
      end do

      !$OMP END PARALLEL DO

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_flux_vel_ghost

!***********************************************************************
!BOP
! !IROUTINE: advu
! !INTERFACE:

 subroutine advu(k,LUK,LVK,WUK,UUU,VVV,this_block)

! !DESCRIPTION:
!  This routine computes the $x,y$ components of the advection of 
!  momentum and metric terms given by:
!  \begin{equation}
!    L_U(u_x) + u_x u_y k_y - u_y^2 k_x
!  \end{equation}
!  \begin{equation}
!    L_U(u_y) + u_x u_y k_x - u_x^2 k_y
!  \end{equation}
!  where
!  \begin{equation}
!     L_U(\alpha) = {1\over{\Delta_y}}
!                \delta_x\left[
!                 \overline{\left(\overline{\Delta_y u_x}^y\right)}^{xy}
!                 \overline{\alpha}^x\right] + 
!                {1\over{\Delta_x}}
!                \delta_y\left[
!                 \overline{\left(\overline{\Delta_x u_y}^x\right)}^{xy}
!                 \overline{\alpha}^y\right] 
!           + \delta_z(w^U\overline{\alpha}^z),
!  \end{equation}
!  \begin{equation}
!  k_x = {1\over{\Delta_y}}\delta_x \Delta_y
!  \end{equation}
!  and
!  \begin{equation}
!  k_y = {1\over{\Delta_x}}\delta_y \Delta_x
!  \end{equation}
!
!  Comments:
!   \begin{itemize}
!     \item this routine must be called successively with k = 1,2,3,...
!
!     \item the vertical velocity $w^U$ in U columns is determined by 
!     integrating the continuity equation $L(1) = 0$ from the surface
!     down to level k.  In the rigid-lid formulation, the integration
!     starts with $w^U = 0$ at the surface.  In the free-surface 
!     formulation, the integration starts with $w^U =$ the area-weighted
!     average of $\partial\eta/\partial t$ at surrounding T points
!     ($\partial\eta/\partial t$ is the time change of the surface 
!     height, and satisfies the barotropic continuity equation 
!     $\partial\eta/\partial t + \nabla\cdot H {\rm\bf U}=q_w$
!     where ${\rm\bf U}$ is the barotropic velocity and $q_w$
!     is the surface fresh water flux.)
!   \end{itemize}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k ! depth level index

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: & 
      UUU,VVV             ! U,V velocity for this block 
                          ! at the current time

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout) :: & 
      WUK             ! on  input: flux velocity at top    of U box
                      ! on output: flux velocity at bottom of U box

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      LUK,               &! advection of U-momentum
      LVK                 ! advection of V-momentum

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,             &! loop indices
      ibeg,iend,         &! beginning and ending index of
      jbeg,jend,         &!  physical domain
      bid                 ! local block index

   real (r8) ::          &
      cc                  ! scalar central weight

   real (r8), dimension(nx_block,ny_block) :: &
      UUE,UUW,VUN,VUS, &! momen  flux velocities across E,W,N,S faces
      FVN,FUE,         &! stencil coeffs used by tavg
      WORK,            &! local temp space
      WUKB              ! vert velocity at bottom of level k U box

!-----------------------------------------------------------------------
!
!  advection fluxes for U box (4-point averages of adv fluxes
!  across T boxes).
!
!-----------------------------------------------------------------------

   bid  = this_block%local_id

   call timer_start(timer_advu, block_id=bid)

   ibeg = this_block%ib
   iend = this_block%ie
   jbeg = this_block%jb
   jend = this_block%je

   UUW = c0
   UUE = c0
   VUN = c0
   VUS = c0

   if (partial_bottom_cells) then

      do j=jbeg-1,jend+1
      do i=ibeg-1,iend+1

         UUW(i,j) = p25 *(UUU(i  ,j  ,k)*DYU(i  ,j    ,bid)*  &
                                         DZU(i  ,j  ,k,bid) + &
                          UUU(i-1,j  ,k)*DYU(i-1,j    ,bid)*  &
                                         DZU(i-1,j  ,k,bid))+ &
                    p125*(UUU(i  ,j-1,k)*DYU(i  ,j-1  ,bid)*  &
                                         DZU(i  ,j-1,k,bid) + &
                          UUU(i-1,j-1,k)*DYU(i-1,j-1  ,bid)*  &
                                         DZU(i-1,j-1,k,bid) + &
                          UUU(i  ,j+1,k)*DYU(i  ,j+1  ,bid)*  &
                                         DZU(i  ,j+1,k,bid) + &
                          UUU(i-1,j+1,k)*DYU(i-1,j+1  ,bid)*  &
                                         DZU(i-1,j+1,k,bid))

         UUE(i,j) = p25 *(UUU(i+1,j  ,k)*DYU(i+1,j    ,bid)*  &
                                         DZU(i+1,j  ,k,bid) + &
                          UUU(i  ,j  ,k)*DYU(i  ,j    ,bid)*  &
                                         DZU(i  ,j  ,k,bid))+ &
                    p125*(UUU(i+1,j-1,k)*DYU(i+1,j-1  ,bid)*  &
                                         DZU(i+1,j-1,k,bid) + &
                          UUU(i  ,j-1,k)*DYU(i  ,j-1  ,bid)*  &
                                         DZU(i  ,j-1,k,bid) + &
                          UUU(i+1,j+1,k)*DYU(i+1,j+1  ,bid)*  &
                                         DZU(i+1,j+1,k,bid) + &
                          UUU(i  ,j+1,k)*DYU(i  ,j+1  ,bid)*  &
                                         DZU(i  ,j+1,k,bid))
   
         VUS(i,j) = p25* (VVV(i  ,j  ,k)*DXU(i  ,j    ,bid)*  &
                                         DZU(i  ,j  ,k,bid) + &
                          VVV(i  ,j-1,k)*DXU(i  ,j-1  ,bid)*  &
                                         DZU(i  ,j-1,k,bid))+ &
                    p125*(VVV(i-1,j  ,k)*DXU(i-1,j    ,bid)*  &
                                         DZU(i-1,j  ,k,bid) + &
                          VVV(i-1,j-1,k)*DXU(i-1,j-1  ,bid)*  &
                                         DZU(i-1,j-1,k,bid) + &
                          VVV(i+1,j  ,k)*DXU(i+1,j    ,bid)*  &
                                         DZU(i+1,j  ,k,bid) + &
                          VVV(i+1,j-1,k)*DXU(i+1,j-1  ,bid)*  &
                                         DZU(i+1,j-1,k,bid))

         VUN(i,j) = p25* (VVV(i  ,j+1,k)*DXU(i  ,j+1  ,bid)*  &
                                         DZU(i  ,j+1,k,bid) + &
                          VVV(i  ,j  ,k)*DXU(i  ,j    ,bid)*  &
                                         DZU(i  ,j  ,k,bid))+ &
                    p125*(VVV(i-1,j+1,k)*DXU(i-1,j+1  ,bid)*  &
                                         DZU(i-1,j+1,k,bid) + &
                          VVV(i-1,j  ,k)*DXU(i-1,j    ,bid)*  &
                                         DZU(i-1,j  ,k,bid) + &
                          VVV(i+1,j+1,k)*DXU(i+1,j+1  ,bid)*  &
                                         DZU(i+1,j+1,k,bid) + &
                          VVV(i+1,j  ,k)*DXU(i+1,j    ,bid)*  &
                                         DZU(i+1,j  ,k,bid))

      end do
      end do

   else

      do j=jbeg-1,jend+1
      do i=ibeg-1,iend+1

         UUW(i,j) = p25 *(UUU(i  ,j  ,k)*DYU(i  ,j  ,bid) + &
                          UUU(i-1,j  ,k)*DYU(i-1,j  ,bid))+ &
                    p125*(UUU(i  ,j-1,k)*DYU(i  ,j-1,bid) + &
                          UUU(i-1,j-1,k)*DYU(i-1,j-1,bid) + &
                          UUU(i  ,j+1,k)*DYU(i  ,j+1,bid) + &
                          UUU(i-1,j+1,k)*DYU(i-1,j+1,bid))

         UUE(i,j) = p25 *(UUU(i+1,j  ,k)*DYU(i+1,j  ,bid) + &
                          UUU(i  ,j  ,k)*DYU(i  ,j  ,bid))+ &
                    p125*(UUU(i+1,j-1,k)*DYU(i+1,j-1,bid) + &
                          UUU(i  ,j-1,k)*DYU(i  ,j-1,bid) + &
                          UUU(i+1,j+1,k)*DYU(i+1,j+1,bid) + &
                          UUU(i  ,j+1,k)*DYU(i  ,j+1,bid))

         VUS(i,j) = p25* (VVV(i  ,j  ,k)*DXU(i  ,j  ,bid) + &
                          VVV(i  ,j-1,k)*DXU(i  ,j-1,bid))+ &
                    p125*(VVV(i-1,j  ,k)*DXU(i-1,j  ,bid) + &
                          VVV(i-1,j-1,k)*DXU(i-1,j-1,bid) + &
                          VVV(i+1,j  ,k)*DXU(i+1,j  ,bid) + &
                          VVV(i+1,j-1,k)*DXU(i+1,j-1,bid))

         VUN(i,j) = p25* (VVV(i  ,j+1,k)*DXU(i  ,j+1,bid) + &
                          VVV(i  ,j  ,k)*DXU(i  ,j  ,bid))+ &
                    p125*(VVV(i-1,j+1,k)*DXU(i-1,j+1,bid) + &
                          VVV(i-1,j  ,k)*DXU(i-1,j  ,bid) + &
                          VVV(i+1,j+1,k)*DXU(i+1,j+1,bid) + &
                          VVV(i+1,j  ,k)*DXU(i+1,j  ,bid))

      end do
      end do

   endif ! partial bottom cells


!-----------------------------------------------------------------------
!
!  calculate vertical velocity at bottom of kth level
!  (vertical velocity is nonzero at bottom of U columns
!  if topography is varying)
!
!-----------------------------------------------------------------------

   if (partial_bottom_cells) then
      WUKB = WUK + (VUN - VUS + UUE - UUW)*UAREA_R(:,:,bid)
   else
      WUKB = WUK + c2dz(k)*p5*(VUN - VUS + UUE - UUW)* &
                   UAREA_R(:,:,bid)
   endif

!#if drifter_particles
!#if drifter_particles != 2
!-----------------------------------------------------------------------
!
!     save the vertical velocity for advecting drifters
!
!-----------------------------------------------------------------------
!      WVEL(:,:,k) = WUKB
!#endif
!#endif

!-----------------------------------------------------------------------
!
!  advect momentum
!
!  horizontal advection
!
!-----------------------------------------------------------------------

   LUK = c0
   LVK = c0

   if (partial_bottom_cells) then

      do j=jbeg,jend
      do i=ibeg,iend

         cc = VUS(i,j+1) - VUS(i,j) + UUW(i+1,j) - UUW(i,j)

         LUK(i,j) = p5*(        cc*UUU(i  ,j  ,k) + &
                        VUS(i,j+1)*UUU(i  ,j+1,k) - &
                        VUS(i,j  )*UUU(i  ,j-1,k) + &
                        UUW(i+1,j)*UUU(i+1,j  ,k) - &
                        UUW(i  ,j)*UUU(i-1,j  ,k))* &
                        UAREA_R(i,j,bid)/    &
                        DZU(i,j,k,bid)

         LVK(i,j) = p5*(        cc*VVV(i  ,j  ,k) + &
                        VUS(i,j+1)*VVV(i  ,j+1,k) - &
                        VUS(i,j  )*VVV(i  ,j-1,k) + &
                        UUW(i+1,j)*VVV(i+1,j  ,k) - &
                        UUW(i  ,j)*VVV(i-1,j  ,k))* &
                        UAREA_R(i,j,bid)/    &
                        DZU(i,j,k,bid)

      end do
      end do

   else ! no partial bottom cells

      do j=jbeg,jend
      do i=ibeg,iend

         cc = VUS(i,j+1) - VUS(i,j) + UUW(i+1,j) - UUW(i,j)

         LUK(i,j) = p5*(        cc*UUU(i  ,j  ,k) + &
                        VUS(i,j+1)*UUU(i  ,j+1,k) - &
                        VUS(i,j  )*UUU(i  ,j-1,k) + &
                        UUW(i+1,j)*UUU(i+1,j  ,k) - &
                        UUW(i  ,j)*UUU(i-1,j  ,k))* &
                        UAREA_R(i,j,bid)

         LVK(i,j) = p5*(        cc*VVV(i  ,j  ,k) + &
                        VUS(i,j+1)*VVV(i  ,j+1,k) - &
                        VUS(i,j  )*VVV(i  ,j-1,k) + &
                        UUW(i+1,j)*VVV(i+1,j  ,k) - &
                        UUW(i  ,j)*VVV(i-1,j  ,k))* &
                        UAREA_R(i,j,bid)

      end do
      end do

   endif ! partial bottom cells

!-----------------------------------------------------------------------
!
!  vertical advection through top of U box
!
!-----------------------------------------------------------------------

   if (k == 1) then
      LUK = LUK + dzr(k)*WUK*UUU(:,:,k)
      LVK = LVK + dzr(k)*WUK*VVV(:,:,k)
   else
      if (partial_bottom_cells) then
         LUK = LUK + p5/DZU(:,:,k,bid)*WUK*(UUU(:,:,k-1) + &
                                                   UUU(:,:,k))
         LVK = LVK + p5/DZU(:,:,k,bid)*WUK*(VVV(:,:,k-1) + &
                                                   VVV(:,:,k))
      else
         LUK = LUK + dz2r(k)*WUK*(UUU(:,:,k-1) + UUU(:,:,k))
         LVK = LVK + dz2r(k)*WUK*(VVV(:,:,k-1) + VVV(:,:,k))
      endif
   endif

!-----------------------------------------------------------------------
!
!  vertical advection through bottom of U box
!  for k=km, UUU(k+1) is not defined, but WUKB=0
!
!-----------------------------------------------------------------------

   if (k < km) then  
      if (partial_bottom_cells) then
         LUK = LUK - p5/DZU(:,:,k,bid)*WUKB*(UUU(:,:,k) + &
                                                    UUU(:,:,k+1))
         LVK = LVK - p5/DZU(:,:,k,bid)*WUKB*(VVV(:,:,k) + &
                                                    VVV(:,:,k+1))
      else
         LUK = LUK - dz2r(k)*WUKB*(UUU(:,:,k) + UUU(:,:,k+1))
         LVK = LVK - dz2r(k)*WUKB*(VVV(:,:,k) + VVV(:,:,k+1))
      endif
   endif

!-----------------------------------------------------------------------
!
!  add metric terms, and zero fields at land points
!
!-----------------------------------------------------------------------

   do j=jbeg,jend
   do i=ibeg,iend
      if (k <= KMU(i,j,bid)) then
         LUK(i,j) = LUK(i,j) + UUU(i,j,k)*VVV(i,j,k)*KYU(i,j,bid) - &
                               VVV(i,j,k)**2*KXU(i,j,bid)
         LVK(i,j) = LVK(i,j) + UUU(i,j,k)*VVV(i,j,k)*KXU(i,j,bid) - &
                               UUU(i,j,k)**2*KYU(i,j,bid)
      else
         LUK(i,j) = c0
         LVK(i,j) = c0
      endif
   end do
   end do

!------------------------------------------------------------------------------------
!
!  accumulate time average if necessary. testing is internal to accumulate_tavg_field
!
!------------------------------------------------------------------------------------

   if (mix_pass /= 1) then

      if (partial_bottom_cells) then
         FUE =  UUE*p5*UAREA_R(:,:,bid)/DZU(:,:,k,bid)
      else
         FUE =  UUE*p5*UAREA_R(:,:,bid)
      endif

      do j=jbeg,jend
      do i=ibeg,iend
         WORK(i,j) = FUE(i,j)*(UUU(i,j,k) + UUU(i+1,j,k))
      end do
      end do
      call accumulate_tavg_field(WORK,tavg_UEU,bid,k)

      do j=jbeg,jend
      do i=ibeg,iend
         WORK(i,j) = FUE(i,j)*(VVV(i  ,j,k) + VVV(i+1,j,k))
      end do
      end do
      call accumulate_tavg_field(WORK,tavg_UEV,bid,k)

      if (partial_bottom_cells) then
         FVN =  VUN*p5*UAREA_R(:,:,bid)/DZU(:,:,k,bid)
      else
         FVN =  VUN*p5*UAREA_R(:,:,bid)
      endif

      do j=jbeg,jend
      do i=ibeg,iend
         WORK(i,j) = FVN(i,j)*(UUU(i,j  ,k) + UUU(i,j+1,k))
      end do
      end do
      call accumulate_tavg_field(WORK,tavg_VNU,bid,k)

      do j=jbeg,jend
      do i=ibeg,iend
         WORK(i,j) = FVN(i,j)*(VVV(i,j  ,k) + VVV(i,j+1,k))
      end do
      end do
      call accumulate_tavg_field(WORK,tavg_VNV,bid,k)

      if (k == 1) then
         WORK = dzr(k)*WUK*UUU(:,:,k)
      else
         WORK = dz2r(k)*WUK*(UUU(:,:,k) + UUU(:,:,k-1))
      endif
      call accumulate_tavg_field(WORK,tavg_WTU,bid,k)

      if (k == 1) then
         WORK = dzr(k)*WUK*VVV(:,:,k)
      else
         WORK = dz2r(k)*WUK*(VVV(:,:,k) + VVV(:,:,k-1))
      endif
      call accumulate_tavg_field(WORK,tavg_WTV,bid,k)

   endif ! mix_pass

!-----------------------------------------------------------------------
!
!  reset W for next k level
!
!-----------------------------------------------------------------------

   WUK = WUKB   ! bottom value becomes top value for next pass

   call timer_stop(timer_advu, block_id=bid)

!-----------------------------------------------------------------------
!EOC

 end subroutine advu

!***********************************************************************
!BOP
! !IROUTINE: advt
! !INTERFACE:

 subroutine advt(k,LTK,WTK,TMIX,TRCR,UUU,VVV,this_block)

! !DESCRIPTION:
!  Advection of tracers - this routine actually computes only
!  vertical velocities and then calls individual routines based
!  on the type of advection chosen
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k  ! depth level index

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX,              &! tracers for this block at mix     time
      TRCR                ! tracers for this block at current time

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU, VVV            ! U,V for this block at current time

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout) :: & 
      WTK            ! on  input flux velocity at top    of T box
                     ! on output flux velocity at bottom of T box

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(out) :: & 
      LTK             ! returned as L(T) for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local or common variables:
!
!-----------------------------------------------------------------------

   logical (log_kind), dimension(nt) :: &
     tr_mask              ! which tracers are using a particular advection scheme

   integer (int_kind) :: &
     i,j,n,              &! dummy loop indices
     ib, ie, jb, je,     &! domain limits
     bid                  ! local block address

   real (r8), dimension(nx_block,ny_block) :: & 
     UTE,UTW,VTN,VTS,  &! tracer flux velocities across E,W,N,S faces
     FVN,FUE,          &! north and east advective stencil weights
     WTKB,             &! vertical velocity at bottom of T box
     RHOK1,            &! pot density at k relative to k=1
     RHOK1M,           &! pot density at k-1
     WORK,             &! local temp space
     WORK1,WORK2,      &! local temp space
     WORK3,WORK4

   real (r8), dimension(nx_block,ny_block,nt) :: & 
     TRACER_E,         &! tracer value on east face
     TRACER_N,         &! tracer value on north face
     FLUX_T             ! tracer tendency due to flux through top   face, volume normalized

!-----------------------------------------------------------------------
!
!  advection fluxes for T box.
!
!-----------------------------------------------------------------------

   bid = this_block%local_id
   jb  = this_block%jb
   je  = this_block%je
   ib  = this_block%ib
   ie  = this_block%ie

   call timer_start(timer_advt, block_id=bid)

!-----------------------------------------------------------------------
!
!  calculate vertical velocity at bottom of kth level
!  (vertical velocity is zero at bottom of T columns)
!
!  use flux velocities from k+1 of previous iteration if available
!
!-----------------------------------------------------------------------

   if (k == 1 .or. .not. luse_lw_lim) then
      call comp_flux_vel(k,UUU,VVV,WTK, &
                         UTE,UTW,VTN,VTS,WTKB,this_block)
   else
      UTE  = FLUX_VEL_prev(:,:,1,bid)
      UTW  = FLUX_VEL_prev(:,:,2,bid)
      VTN  = FLUX_VEL_prev(:,:,3,bid)
      VTS  = FLUX_VEL_prev(:,:,4,bid)
      WTKB = FLUX_VEL_prev(:,:,5,bid)
   end if

!-----------------------------------------------------------------------
!
!  advect tracers
!
!-----------------------------------------------------------------------

   LTK = c0

   !*** lw_lim advection

   if (luse_lw_lim) then

      tr_mask = tadvect_itype == tadvect_lw_lim

!-----------------------------------------------------------------------
!     compute WTKB at bottom of level k+1 T box and store flux velocities
!     WTKBp1 is stored in FLUX_VEL_prev(:,:,5,bid)
!-----------------------------------------------------------------------

      call comp_flux_vel(k+1,UUU,VVV,WTKB, &
                         FLUX_VEL_prev(:,:,1,bid),&
                         FLUX_VEL_prev(:,:,2,bid),&
                         FLUX_VEL_prev(:,:,3,bid),&
                         FLUX_VEL_prev(:,:,4,bid),&
                         FLUX_VEL_prev(:,:,5,bid),this_block)

      call advt_lw_lim(k,LTK,TMIX,WTK,WTKB,FLUX_VEL_prev(:,:,5,bid),&
                       UTE,UTW,VTN,VTS,TRACER_E,TRACER_N,FLUX_T,&
                       tr_mask,this_block)

   end if !luse_lw_lim

   !*** upwind3 advection

   if (luse_upwind3) then

      tr_mask = tadvect_itype == tadvect_upwind3

      call advt_upwind3(k,LTK,TRCR,WTK,WTKB,UTE,UTW,VTN,VTS, &
                          TRACER_E,TRACER_N,FLUX_T,tr_mask,this_block)

   end if

   !*** centered advection

   if (luse_centered) then

      tr_mask = tadvect_itype == tadvect_centered

      call advt_centered(k,LTK,TRCR,WTK,WTKB,UTE,VTN,tr_mask,this_block)

   end if

   call timer_stop(timer_advt, block_id=bid)

!-----------------------------------------------------------------------
!
!  compute diagnostics if necessary. Note: testing is done internally
!     is accumulate_tavg_field
!
!-----------------------------------------------------------------------

   if (mix_pass /= 1 ) then

      if (partial_bottom_cells) then
         FVN =  p5*VTN*TAREA_R(:,:,bid)/DZT(:,:,k,bid)
         FUE =  p5*UTE*TAREA_R(:,:,bid)/DZT(:,:,k,bid)
      else
         FVN =  p5*VTN*TAREA_R(:,:,bid)
         FUE =  p5*UTE*TAREA_R(:,:,bid)
      endif

      call accumulate_tavg_field(WTK,tavg_WVEL,bid,k)
      call accumulate_tavg_field(WTK**2,tavg_WVEL2,bid,k)

      do n=1,nt
         if (tadvect_itype(n) == tadvect_centered) then

            WORK = FUE*(        TRCR(:,:,k,n) + eoshift(TRCR(:,:,k,n),dim=1,shift=1))
            call accumulate_tavg_field(WORK,tavg_UE_TRACER(n),bid,k)

            WORK = FVN*(        TRCR(:,:,k,n) +  eoshift(TRCR(:,:,k,n),dim=2,shift=1))
            call accumulate_tavg_field(WORK,tavg_VN_TRACER(n),bid,k)

            if (k == 1) then
               if (sfc_layer_type /= sfc_layer_varthick) then
                  WORK = dzr(k)*WTK*TRCR(:,:,k,n)
               else
                  WORK = c0
               endif
            else
               WORK = dz2r(k)*WTK*(TRCR(:,:,k  ,n) + TRCR(:,:,k-1,n))
            endif
            call accumulate_tavg_field(WORK,tavg_WT_TRACER(n),bid,k)

         else

            WORK = c2*FUE*TRACER_E(:,:,n)
            call accumulate_tavg_field(WORK,tavg_UE_TRACER(n),bid,k)

            WORK = c2*FVN*TRACER_N(:,:,n)
            call accumulate_tavg_field(WORK,tavg_VN_TRACER(n),bid,k)

            call accumulate_tavg_field(FLUX_T(:,:,n),tavg_WT_TRACER(n),bid,k)

         endif

         WORK = c0
         if (partial_bottom_cells) then
            do j=jb,je
            do i=ib,ie
               if (k <= KMT(i,j,bid)) then
                  WORK(i,j) = -DZT(i,j,k,bid)*LTK(i,j,n)
               endif
            end do
            end do
         else
            do j=jb,je
            do i=ib,ie
               if (k <= KMT(i,j,bid)) then
                  WORK(i,j) = -dz(k)*LTK(i,j,n)
               endif
            end do
            end do
         endif
         call accumulate_tavg_field(WORK,tavg_ADV_TRACER(n),bid,k)

      enddo

      !***
      !*** potential density referenced to k=1 needed for a variety of 
      !*** tavg fields
      !***

   if (accumulate_tavg_now(tavg_PD)   .or.  &
       accumulate_tavg_now(tavg_RHOU) .or.  &
       accumulate_tavg_now(tavg_RHOV) .or.  &
       accumulate_tavg_now(tavg_URHO) .or.  &
       accumulate_tavg_now(tavg_VRHO) .or.  &
       accumulate_tavg_now(tavg_WRHO)) then

       call state(k,1,TRCR(:,:,k,1), TRCR(:,:,k,2), this_block, RHOFULL=RHOK1)

       if (k == 1) then
          RHOK1M = RHOK1
       else
          call state(k-1,1,TRCR(:,:,k,1), TRCR(:,:,k,2), this_block, RHOFULL=RHOK1M)
       endif

       call accumulate_tavg_field(RHOK1,tavg_PD,bid,k)

       WORK = FUE*(RHOK1 + eoshift(RHOK1,dim=1,shift=1))
       call accumulate_tavg_field(WORK,tavg_URHO,bid,k)

       WORK = FVN*(RHOK1 + eoshift(RHOK1,dim=2,shift=1))
       call accumulate_tavg_field(WORK,tavg_VRHO,bid,k)
 
       WORK = dz2r(k)*WTK*(RHOK1 + RHOK1M)
       call accumulate_tavg_field(WORK,tavg_WRHO,bid,k)

       call ugrid_to_tgrid(WORK1,UUU(:,:,k),bid)
       WORK = RHOK1*WORK1
       call accumulate_tavg_field(WORK,tavg_RHOU,bid,k)

       call ugrid_to_tgrid(WORK1,VVV(:,:,k),bid)
       WORK = RHOK1*WORK1
       call accumulate_tavg_field(WORK,tavg_RHOV,bid,k)
   endif! accumulate_tavg_now

      !***
      !*** vertical density gradient and potential vorticity
      !***

   if (accumulate_tavg_now(tavg_Q)    .or.  &
       accumulate_tavg_now(tavg_PV)   .or.  &
       accumulate_tavg_now(tavg_PVWM) .or.  &
       accumulate_tavg_now(tavg_PVWP) .or.  &
       accumulate_tavg_now(tavg_UPV)  .or.  &
       accumulate_tavg_now(tavg_VPV)  .or.  &
       accumulate_tavg_now(tavg_UQ)   .or.  &
       accumulate_tavg_now(tavg_VQ) ) then

       call state(k,k,TRCR(:,:,k,1),TRCR(:,:,k,2), this_block,RHOOUT=WORK)

       if (k == 1 ) then
          WORK3 = WORK
       else
          call state(k-1,k,TRCR(:,:,k-1,1),TRCR(:,:,k-1,2), this_block,RHOOUT=WORK3)
          WORK3 = p5*(WORK3 + WORK)
       endif

       if (k == km) then
          WORK4 = WORK
       else
          call state(k+1,k,TRCR(:,:,k+1,1),TRCR(:,:,k+1,2),this_block,RHOOUT=WORK4)

          do j=jb,je
          do i=ib,ie
             if (k /= KMT(i,j,bid)) then
                WORK4(i,j) = p5*(WORK4(i,j) + WORK(i,j))
             else 
                WORK4(i,j) = WORK(i,j)
             endif
          end do
          end do
       endif

       call zero_ghost_cells(this_block,WORK1)

       do j=jb,je
       do i=ib,ie
          if (k <= KMT(i,j,bid)) then
             WORK1(i,j) = (WORK3(i,j) - WORK4(i,j))*dzr(k) ! drho/dz
          else
             WORK1(i,j) = c0
          endif
       end do
       end do

       call zcurl(k,WORK3,UUU(:,:,k),VVV(:,:,k),this_block)
       WORK2 = WORK1*(WORK3*TAREA_R(:,:,bid) + FCORT(:,:,bid)) ! PV = pot vorticity

       call accumulate_tavg_field(WORK1,tavg_Q,bid,k)
       call accumulate_tavg_field(WORK2,tavg_PV,bid,k)

      !***
      !*** advection of potential vorticity, Q
      !***

       call accumulate_tavg_field(WORK2*WTK,tavg_PVWM,bid,k)
       call accumulate_tavg_field(WORK2*WTKB,tavg_PVWP,bid,k)

       WORK = FUE*(WORK2 + eoshift(WORK2,dim=1,shift=1))
       call accumulate_tavg_field(WORK,tavg_UPV,bid,k)

       WORK = FVN*(WORK2 + eoshift(WORK2,dim=2,shift=1))
       call accumulate_tavg_field(WORK,tavg_VPV,bid,k)

       WORK = FUE*(WORK1 + eoshift(WORK1,dim=1,shift=1))
       call accumulate_tavg_field(WORK,tavg_UQ,bid,k)

       WORK = FVN*(WORK1 + eoshift(WORK1,dim=2,shift=1))
       call accumulate_tavg_field(WORK,tavg_VQ,bid,k)
   endif! accumulate_tavg_now

   endif ! mix_pass

   call cfl_advect(k,bid,UTE,UTW,VTN,VTS,WTK,WTKB,this_block)

!-----------------------------------------------------------------------
!
!  set bottom vertical velocity to top for next k-level
!
!-----------------------------------------------------------------------

   WTK = WTKB

!-----------------------------------------------------------------------
!EOC

 end subroutine advt

!***********************************************************************
!BOP
! !IROUTINE: comp_flux_vel
! !INTERFACE:

 subroutine comp_flux_vel(k,UUU,VVV,WTK, &
                          UTE,UTW,VTN,VTS,WTKB,this_block)

! !DESCRIPTION:
!  This routine computes the vertical velocity at bottom of a layer.
!  It also returns the tracer flux velocities across the lateral faces.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! depth level index

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UUU, VVV            ! U,V for this block at current time

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      WTK                 ! vertical velocity at top of T box

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      UTE,UTW,VTN,VTS,   &! tracer flux velocities across E,W,N,S faces
      WTKB                ! vertical velocity at bottom of T box

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! tracer loop index
      ibeg,iend,         &! loop limits
      jbeg,jend,         &! loop limits
      bid                 ! local block index

   real (r8), dimension(nx_block,ny_block) :: & 
      FC                  ! local temp space

!-----------------------------------------------------------------------

   call zero_ghost_cells(this_block,UTE)
   call zero_ghost_cells(this_block,UTW)
   call zero_ghost_cells(this_block,VTN)
   call zero_ghost_cells(this_block,VTS)

   if (k > km) then

      WTKB = c0
      return

   endif

   jbeg = this_block%jb-1
   jend = this_block%je+1
   ibeg = this_block%ib-1
   iend = this_block%ie+1

   if (luse_lw_lim) jend = this_block%je+2

   bid = this_block%local_id

   if (partial_bottom_cells) then

      do j=jbeg,jend
      do i=ibeg,iend

         UTE(i,j) = p5*(UUU(i  ,j  ,k)*DYU(i  ,j    ,bid)*  &
                                       DZU(i  ,j  ,k,bid) + &
                        UUU(i  ,j-1,k)*DYU(i  ,j-1  ,bid)*  &
                                       DZU(i  ,j-1,k,bid))
         UTW(i,j) = p5*(UUU(i-1,j  ,k)*DYU(i-1,j    ,bid)*  &
                                       DZU(i-1,j  ,k,bid) + &
                        UUU(i-1,j-1,k)*DYU(i-1,j-1  ,bid)*  &
                                       DZU(i-1,j-1,k,bid))

         VTN(i,j) = p5*(VVV(i  ,j  ,k)*DXU(i  ,j    ,bid)*  &
                                       DZU(i  ,j  ,k,bid) + &
                        VVV(i-1,j  ,k)*DXU(i-1,j    ,bid)*  &
                                       DZU(i-1,j  ,k,bid))
         VTS(i,j) = p5*(VVV(i  ,j-1,k)*DXU(i  ,j-1  ,bid)*  &
                                       DZU(i  ,j-1,k,bid) + &
                        VVV(i-1,j-1,k)*DXU(i-1,j-1  ,bid)*  &
                                       DZU(i-1,j-1,k,bid))

      end do
      end do

   else

      do j=jbeg,jend
      do i=ibeg,iend

         UTE(i,j) = p5*(UUU(i  ,j  ,k)*DYU(i  ,j  ,bid) + &
                        UUU(i  ,j-1,k)*DYU(i  ,j-1,bid))
         UTW(i,j) = p5*(UUU(i-1,j  ,k)*DYU(i-1,j  ,bid) + &
                        UUU(i-1,j-1,k)*DYU(i-1,j-1,bid))

         VTN(i,j) = p5*(VVV(i  ,j  ,k)*DXU(i  ,j  ,bid) + &
                        VVV(i-1,j  ,k)*DXU(i-1,j  ,bid))
         VTS(i,j) = p5*(VVV(i  ,j-1,k)*DXU(i  ,j-1,bid) + &
                        VVV(i-1,j-1,k)*DXU(i-1,j-1,bid))

      end do
      end do

   endif ! partial bottom cells

   if (luse_lw_lim) then
      UTE(:,this_block%jb-2)  = UTE_jbm2(:,k,bid)
      UTW(2:nx_block,this_block%jb-2) = &
         UTE(1:nx_block-1,this_block%jb-2)
      UTE(this_block%ib-2,:)  = UTW(this_block%ib-1,:)
      VTN(:,this_block%jb-2)  = VTS(:,this_block%jb-1)
   endif

!-----------------------------------------------------------------------
!
!  calculate vertical velocity at bottom of kth level
!  (vertical velocity is zero at bottom of T columns)
!
!-----------------------------------------------------------------------

   if (k < km) then

      FC = (VTN - VTS + UTE - UTW)*TAREA_R(:,:,bid)

      if ( overflows_on ) then
         WTKB = WTK+dz(k)*FC
         call ovf_wtkb_check(k,WTKB,this_block)
      endif

      if (partial_bottom_cells) then
         WTKB = merge(WTK+FC, c0, k < KMT(:,:,bid))
      else
         WTKB = merge(WTK+dz(k)*FC, c0, k < KMT(:,:,bid))
      endif

      if (luse_lw_lim) then
         WTKB(:,this_block%jb-2) = WTKB_jbm2(:,k,bid)
         WTKB(:,this_block%je+2) = WTKB_jep2(:,k,bid)
         WTKB(this_block%ib-2,:) = WTKB_ibm2(:,k,bid)
         WTKB(this_block%ie+2,:) = WTKB_iep2(:,k,bid)
      endif

   else

      WTKB = c0

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_flux_vel

!***********************************************************************
!BOP
! !IROUTINE: advt_centered
! !INTERFACE:

 subroutine advt_centered(k,LTK,TRCR,WTK,WTKB,UTE,VTN,tr_mask,this_block)

! !DESCRIPTION:
!  This routine computes the tracer advection tendencies using
!  a centered differencing for leapfrog advection given by:
!  \begin{equation}
!     L_T(\phi) = {1\over{\Delta_y}}
!                  \delta_x\left[
!                   \overline{\left(\overline{\Delta_y u_x}^y\right)}^{xy}
!                   \overline{\phi}^x\right] + 
!                  {1\over{\Delta_x}}
!                  \delta_y\left[
!                   \overline{\left(\overline{\Delta_x u_y}^x\right)}^{xy}
!                   \overline{\phi}^y\right] 
!             + \delta_z(w\overline{\phi}^z),
!  \end{equation}
!  where $\phi$ is the tracer concentration.
!
!  Comments:
!   \begin{itemize}
!     \item this routine must be called successively with k = 1,2,3,...
!
!     \item the vertical velocity $w$ in T columns is determined by 
!     integrating the continuity equation $L(1) = 0$ from the surface
!     down to level k.  In the rigid-lid formulation, the integration
!     starts with $w = 0$ at the surface.  In the free-surface 
!     formulation, the integration starts with 
!     $w = \partial\eta/\partial t$ 
!     ($\partial\eta/\partial t$ is the time change of the surface 
!     height, and satisfies the barotropic continuity equation 
!     $\partial\eta/\partial t + \nabla\cdot H {\rm\bf U}=q_w$
!     where ${\rm\bf U}$ is the barotropic velocity and $q_w$
!     is the surface fresh water flux.)
!   \end{itemize}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers for this block at current time

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      UTE,VTN,         &! tracer flux velocities across E,N faces
      WTK,             &! vert velocity at top    of level k T box
      WTKB              ! vert velocity at bottom of level k T box

   logical (log_kind), dimension(nt), intent(in) :: &
      tr_mask           ! true for tracers using centered

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      LTK                 ! returned as L(T) for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,             &! tracer loop index
      bid                 ! local block index

!-----------------------------------------------------------------------
!
!  advect tracers
!  horizontal advection
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

!-----------------------------------------------------------------------

   if (partial_bottom_cells) then 

      do n=1,nt
      if (.not. tr_mask(n)) cycle
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         LTK(i,j,n) = p5*((VTN(i,j)-VTN(i,j-1)+UTE(i,j)-UTE(i-1,j))  &
                                      *TRCR(i  ,j  ,k,n) +           &
                          VTN(i  ,j  )*TRCR(i  ,j+1,k,n) -           &
                          VTN(i  ,j-1)*TRCR(i  ,j-1,k,n) +           &
                          UTE(i  ,j  )*TRCR(i+1,j  ,k,n) -           &
                          UTE(i-1,j  )*TRCR(i-1,j  ,k,n))*           &
                          TAREA_R(i,j,bid)/                   &
                          DZT(i,j,k,bid)
      end do
      end do
      end do

   else ! no partial bottom cells

      do n=1,nt
      if (.not. tr_mask(n)) cycle
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         LTK(i,j,n) = p5*((VTN(i,j)-VTN(i,j-1)+UTE(i,j)-UTE(i-1,j))  &
                                      *TRCR(i  ,j  ,k,n) +           &
                          VTN(i  ,j  )*TRCR(i  ,j+1,k,n) -           &
                          VTN(i  ,j-1)*TRCR(i  ,j-1,k,n) +           &
                          UTE(i  ,j  )*TRCR(i+1,j  ,k,n) -           &
                          UTE(i-1,j  )*TRCR(i-1,j  ,k,n))*           &
                          TAREA_R(i,j,bid)
      end do
      end do
      end do

   endif ! partial bottom cells

!-----------------------------------------------------------------------
!
!  vertical advection through top and bottom of T box
!
!-----------------------------------------------------------------------

   do n = 1,nt

      if (.not. tr_mask(n)) cycle

      !*** no advection thru surface (k=1) in vs model

      if (k == 1) then
         if (sfc_layer_type /= sfc_layer_varthick) then
            LTK(:,:,n) = LTK(:,:,n) + &
                         dzr(k)*WTK*TRCR(:,:,k,n)
         endif
      else
         if (partial_bottom_cells) then
            LTK(:,:,n) = LTK(:,:,n) + p5/DZT(:,:,k,bid)*WTK*  &
                         (TRCR(:,:,k-1,n) + TRCR(:,:,k  ,n))
         else
            LTK(:,:,n) = LTK(:,:,n) + dz2r(k)*WTK*  &
                         (TRCR(:,:,k-1,n) + TRCR(:,:,k  ,n))
         endif
      endif

      !***
      !*** for k=km, TRACER(k+1) is not defined, but WTKB=0
      !***

      if (k < km) then
         if (partial_bottom_cells) then
            LTK(:,:,n) = LTK(:,:,n) - p5/DZT(:,:,k,bid)*WTKB* &
                         (TRCR(:,:,k,n) + TRCR(:,:,k+1,n))
         else
            LTK(:,:,n) = LTK(:,:,n) - dz2r(k)*WTKB* &
                         (TRCR(:,:,k,n) + TRCR(:,:,k+1,n))
         endif
      endif

   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine advt_centered

!***********************************************************************
!BOP
! !IROUTINE: advt_upwind3
! !INTERFACE:

 subroutine advt_upwind3(k,LTK,TRCR,WTK,WTKB,UTE,UTW,VTN,VTS, &
                           TRACER_E,TRACER_N,FLUX_T,tr_mask,this_block)

! !DESCRIPTION:
!  This routine computes the advection of tracers using a 3rd-order 
!  upwinding defined by the QUICKEST scheme in
!  Leonard, B.P. 1979, {\em Comp. Meth. Applied Math. and Engineering},
!        {\bf 19}, 59-98.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers for this block at current time

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      UTE,UTW,VTN,VTS, &! tracer flux velocities across E,W,N,S faces
      WTK,             &! vert velocity at top    of level k T box
      WTKB              ! vert velocity at bottom of level k T box

   logical (log_kind), dimension(nt), intent(in) :: &
      tr_mask           ! true for tracers using upwind3

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      LTK,              &! returned as L(T) for nth tracer at level k
      TRACER_E,         &! tracer value on east face
      TRACER_N,         &! tracer value on north face
      FLUX_T             ! tracer tendency due to flux through top   face, volume normalized

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: & 
      i,j,n,             &! loop index
      bid                 ! local block address for this block

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1,                                  &! local temp space
      FVN,FVS,FUE,FUW,                        &! flux velocities
      TPLUS,TMINUS,AZMINUS,DZMINUS

   real (r8), dimension(nx_block,ny_block,nt) :: &
      AUXB

!-----------------------------------------------------------------------
!
!  advect tracers
!  horizontal advection
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   if (partial_bottom_cells) then
      WORK1 = TAREA_R(:,:,bid) / DZT(:,:,k,bid)
      FVN =  VTN*WORK1
      FVS = -VTS*WORK1 ! note sign
      FUE =  UTE*WORK1
      FUW = -UTW*WORK1 ! note sign
   else
      FVN =  VTN*TAREA_R(:,:,bid)
      FVS = -VTS*TAREA_R(:,:,bid)  ! note sign change
      FUE =  UTE*TAREA_R(:,:,bid)
      FUW = -UTW*TAREA_R(:,:,bid)  ! note sign change
   endif

   call hupw3(k,LTK,TRCR,FVN,FVS,FUE,FUW,TRACER_E,TRACER_N, &
                    bid, tr_mask, this_block)

!-----------------------------------------------------------------------
!
!  vertical advection through top and bottom of T box
!
!-----------------------------------------------------------------------

   do j=1,ny_block
   do i=1,nx_block
      if ( k < KMT(i,j,bid)-1 ) then
         AZMINUS(i,j) = talfzm(k)
         DZMINUS(i,j) = tdelzm(k)
      else
         AZMINUS(i,j) = talfzm(k) + tdelzm(k)
         DZMINUS(i,j) = c0
      endif
   end do
   end do

   do n = 1,nt

      if (.not. tr_mask(n)) cycle

      if ( k < km-1 .and. k > 1) then
         TPLUS  = talfzp(k)*TRCR(:,:,k+1,n)  & 
                + tbetzp(k)*TRCR(:,:,k  ,n)  &
                + tgamzp(k)*TRCR(:,:,k-1,n)
         TMINUS = AZMINUS  *TRCR(:,:,k+1,n)  &
                + tbetzm(k)*TRCR(:,:,k  ,n)  &
                + DZMINUS  *TRCR(:,:,k+2,n)
         AUXB(:,:,n) = (WTKB-abs(WTKB))*TPLUS + &
                     (WTKB+abs(WTKB))*TMINUS
      else if ( k == 1) then
         TPLUS  = talfzp(k)*TRCR(:,:,k+1,n)  &
                + tbetzp(k)*TRCR(:,:,k  ,n)
         TMINUS = AZMINUS  *TRCR(:,:,k+1,n)  &
                + tbetzm(k)*TRCR(:,:,k  ,n)  &
                + DZMINUS  *TRCR(:,:,k+2,n)
         AUXB(:,:,n) = (WTKB-abs(WTKB))*TPLUS + &
                       (WTKB+abs(WTKB))*TMINUS
      else if ( k == km-1) then
         TPLUS  = talfzp(k)*TRCR(:,:,k+1,n)  &
                + tbetzp(k)*TRCR(:,:,k  ,n)  &
                + tgamzp(k)*TRCR(:,:,k-1,n)
         TMINUS = AZMINUS  *TRCR(:,:,k+1,n)  &
                + tbetzm(k)*TRCR(:,:,k  ,n)
         AUXB(:,:,n) = (WTKB-abs(WTKB))*TPLUS + &
                       (WTKB+abs(WTKB))*TMINUS
      else
         AUXB(:,:,n) = c0
      endif

      !*** no advection thru surface (k=1) in vs model

      if (k == 1) then
         if (sfc_layer_type /= sfc_layer_varthick) then
            FLUX_T(:,:,n) = dzr(k)*WTK*TRCR(:,:,k,n)
            LTK(:,:,n) = LTK(:,:,n)                            &
                       + FLUX_T(:,:,n)                         &
                       - dz2r(k)*AUXB(:,:,n)
         else
            FLUX_T(:,:,n) = c0
            LTK(:,:,n) = LTK(:,:,n)                            &
                       - dz2r(k)*AUXB(:,:,n)
         endif
      else
         if (partial_bottom_cells) then
            WORK1 = p5 / DZT(:,:,k,bid)
            FLUX_T(:,:,n) = WORK1*AUX(:,:,n,bid)
            LTK(:,:,n) = LTK(:,:,n) + WORK1* &
                         (AUX(:,:,n,bid) - AUXB(:,:,n))
         else
            FLUX_T(:,:,n) = dz2r(k)*AUX(:,:,n,bid)
            LTK(:,:,n) = LTK(:,:,n) + dz2r(k)* &
                         (AUX(:,:,n,bid) - AUXB(:,:,n))
         endif

      endif

      AUX(:,:,n,bid) = AUXB(:,:,n) ! set bottom to top for next k level

   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine advt_upwind3

!***********************************************************************
!BOP
! !IROUTINE: hupw3
! !INTERFACE:

 subroutine hupw3(k,XOUT,X,CN,CS,CE,CW,TRACER_E,TRACER_N, &
                         bid,tr_mask,this_block)

! !DESCRIPTION:
!  This routine is the horizontal stencil operator for 3rd-order upwind
!  advection.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      k,                 &! vertical level index
      bid                 ! local block address for this block

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: & 
      X                   ! input tracer array

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      CN,CS,CE,CW           ! stencil weights based on flux velocities

   logical (log_kind), dimension(nt), intent(in) :: &
      tr_mask           ! true for tracers using upwind3

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      XOUT,              &! output tracer advection
      TRACER_E,          &! tracer value on east face
      TRACER_N            ! tracer value on north face

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n               ! dummy loop indices

   real (r8) :: &
      am, bm, dm, ap, bp, gp, work

!-----------------------------------------------------------------------
!
!  3rd order upwinding stencil
!
!-----------------------------------------------------------------------

   do n=1,nt

      if (.not. tr_mask(n)) cycle

      call zero_ghost_cells(this_block,TRACER_E(:,:,n))
      call zero_ghost_cells(this_block,TRACER_N(:,:,n))

!-----------------------------------------------------------------------
!     Compute grid-x direction contribution
!-----------------------------------------------------------------------

      do j=this_block%jb,this_block%je
         do i=this_block%ib-1,this_block%ie

            if (k <= KMTE(i,j,bid)) then
               work  = TBETXP(i,j,bid)
               ap    = TALFXP(i,j,bid)
            else
               work  = TBETXP(i,j,bid) + TALFXP(i,j,bid)
               ap    = c0
            endif

            if (k <= KMTW(i,j,bid)) then
               bp   = work
               gp   = TGAMXP(i,j,bid)
            else
               bp   = work + TGAMXP(i,j,bid)
               gp   = c0
            endif

            if (k <= KMTEE(i,j,bid)) then
               am   = TALFXM(i,j,bid)
               dm   = TDELXM(i,j,bid)
            else
               am   = TALFXM(i,j,bid) + TDELXM(i,j,bid)
               dm   = c0
            endif
            bm   = TBETXM(i,j,bid)

            if (CE(i,j) > c0) then
               TRACER_E(i,j,n) = ap*X(i+1,j,k,n) + &
                  bp*X(i,j,k,n) + gp*X(i-1,j,k,n)
            else
               TRACER_E(i,j,n) = am*X(i+1,j,k,n) + &
                  bm*X(i,j,k,n) + dm*X(i+2,j,k,n)
            endif

         end do ! i loop for grid-x contribution

         !*** i loop broken up to avoid dependency on previous i iteration
         !*** and to allow different lower loop bounds

         end do ! j loop for grid-x TRACER_E evaluation
         !***  j loop broken up for overflow TRACER_E modification

         if ( overflows_on .and. overflows_interactive ) then
            call ovf_advt(k,TRACER_E(:,:,n),TRACER_N(:,:,n),n,this_block, &
                          CE,CW,CN,CS)
         endif

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie

            XOUT(i,j,n) = CE(i,j)*TRACER_E(i,j,n) + CW(i,j)*TRACER_E(i-1,j,n)

         end do ! i loop for grid-x contribution
      end do ! j loop for grid-x contribution

!-----------------------------------------------------------------------
!     Compute grid-y direction contribution
!-----------------------------------------------------------------------

      do j=this_block%jb-1,this_block%je
         do i=this_block%ib,this_block%ie

            if (k <= KMTN(i,j,bid)) then
               work  = TBETYP(i,j,bid)
               ap    = TALFYP(i,j,bid)
            else
               work  = TBETYP(i,j,bid) + TALFYP(i,j,bid)
               ap    = c0
            endif

            if (k <= KMTS(i,j,bid)) then
               bp   = work
               gp   = TGAMYP(i,j,bid)
            else
               bp   = work + TGAMYP(i,j,bid)
               gp   = c0
            endif

            if (k <= KMTNN(i,j,bid)) then
               am   = TALFYM(i,j,bid)
               dm   = TDELYM(i,j,bid)
            else
               am   = TALFYM(i,j,bid) + TDELYM(i,j,bid)
               dm   = c0
            endif
            bm   = TBETYM(i,j,bid)

            if (CN(i,j) > c0) then
              TRACER_N(i,j,n) = ap*X(i,j+1,k,n) + &
                bp*X(i,j,k,n) + gp*X(i,j-1,k,n)
            else
              TRACER_N(i,j,n) = am*X(i,j+1,k,n) + &
                bm*X(i,j,k,n) + dm*X(i,j+2,k,n)
            endif

         end do ! i loop for grid-y TRACER_N evaluation
         end do ! j loop for grid-y TRACER_N evaluation

         !***  i,j loops broken for overflow TRACER_N modification

         if ( overflows_on .and. overflows_interactive ) then
            call ovf_advt(k,TRACER_E(:,:,n),TRACER_N(:,:,n),n,this_block, &
                          CE,CW,CN,CS)
         endif

         do j=this_block%jb-1,this_block%je
         do i=this_block%ib,this_block%ie
            !*** The formula below is not correct for j==this_block%jb-1
            !*** because TRACER_N(:,this_block%jb-2,n) is not computed.
            !*** This is OK because XOUT(:,this_block%jb-1,n) is not used.

            XOUT(i,j,n) = XOUT(i,j,n) + &
              CN(i,j)*TRACER_N(i,j,n) + CS(i,j)*TRACER_N(i,j-1,n)

         end do ! i loop for grid-y contribution
      end do ! j loop for grid-y contribution
   end do ! tracer loop

!-----------------------------------------------------------------------
!EOC

 end subroutine hupw3 

!***********************************************************************
!BOP
! !IROUTINE: advt_lw_lim
! !INTERFACE:

 subroutine advt_lw_lim(k,LTK,TMIX,WTK,WTKB,WTKBp1,UTE,UTW,VTN,VTS, &
                          TRACER_E,TRACER_N,FLUX_T,tr_mask,this_block)

! !DESCRIPTION:
!  This routine computes the advection of tracers using a 2nd-order
!  forward in time scheme with one-dimensional flux limiters.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX              ! tracers for this block at mix     time

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      UTE,UTW,VTN,VTS, &! tracer flux velocities across E,W,N,S faces
      WTK,             &! vert velocity at top    of level k T box
      WTKB,            &! vert velocity at bottom of level k T box
      WTKBp1            ! vert velocity at bottom of level k+1 T box

   logical (log_kind), dimension(nt), intent(in) :: &
      tr_mask           ! true for tracers using lw_lim

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      LTK,              &! returned as L(T) for nth tracer at level k
      TRACER_E,         &! tracer value on east face
      TRACER_N,         &! tracer value on north face
      FLUX_T             ! tracer tendency due to flux through top   face, volume normalized

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: & 
      i,j,n,             &! loop index
      bid                 ! local block address for this block

   real (r8) :: &
      adv_dt              ! timestep for advection

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1,             &! local temp space
      WTK_EFF,           &! effective velocity at top of cell
      FVN,FVS,FUE,FUW,   &! flux velocities
      UVEL_E_dt,         &! dt * velocity on east face
      VVEL_N_dt,         &! dt * velocity on north face
      KMASKE,KMASKN       ! land mask on east and north faces

   real (r8), dimension(nx_block,ny_block,nt) :: &
      AUXB

!-----------------------------------------------------------------------
!
!  set up scalars & arrays for call to lw_lim
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   adv_dt = c2dtt(k)

   if (partial_bottom_cells) then
      UVEL_E_dt = adv_dt * UTE * UTE_to_UVEL_E(:,:,k,bid)
      VVEL_N_dt = adv_dt * VTN * VTN_to_VVEL_N(:,:,k,bid)

      WORK1 = TAREA_R(:,:,bid) / DZT(:,:,k,bid)
      FVN =  VTN*WORK1
      FVS = -VTS*WORK1 ! note sign
      FUE =  UTE*WORK1
      FUW = -UTW*WORK1 ! note sign
   else
      UVEL_E_dt = adv_dt * UTE * UTE_to_UVEL_E(:,:,1,bid)
      VVEL_N_dt = adv_dt * VTN * VTN_to_VVEL_N(:,:,1,bid)

      FVN =  VTN*TAREA_R(:,:,bid)
      FVS = -VTS*TAREA_R(:,:,bid)  ! note sign change
      FUE =  UTE*TAREA_R(:,:,bid)
      FUW = -UTW*TAREA_R(:,:,bid)  ! note sign change
   endif

   KMASKE = merge(c1, c0, k <= KMT(:,:,bid) .and. k <= KMTE(:,:,bid))
   KMASKN = merge(c1, c0, k <= KMT(:,:,bid) .and. k <= KMTN(:,:,bid))

!-----------------------------------------------------------------------
!
!  set flux at surface (k=1)
!  store top of cell volume averaged flux
!
!-----------------------------------------------------------------------

   if (k == 1 .and. sfc_layer_type == sfc_layer_varthick) then
      WTK_EFF = c0
   else
      WTK_EFF = WTK
   endif

   if (partial_bottom_cells) then
      do n = 1,nt
         if (.not. tr_mask(n)) cycle
         if (k == 1) AUX(:,:,n,bid) = WTK_EFF * TMIX(:,:,k,n)
         FLUX_T(:,:,n) = AUX(:,:,n,bid) / DZT(:,:,k,bid)
      enddo
   else
      do n = 1,nt
         if (.not. tr_mask(n)) cycle
         if (k == 1) AUX(:,:,n,bid) = WTK_EFF * TMIX(:,:,k,n)
         FLUX_T(:,:,n) = AUX(:,:,n,bid)*dzr(k)
      enddo
   endif

   call lw_lim(k,adv_dt,UVEL_E_dt,VVEL_N_dt,LTK,TMIX, &
                 WTK_EFF,WTKB,WTKBp1,AUXB, &
                 FVN,FVS,FUE,FUW,KMASKE,KMASKN, &
                 TRACER_E,TRACER_N, &
                 bid, tr_mask, this_block)

!-----------------------------------------------------------------------
!
!  copy bottom flux to top for next k level
!
!-----------------------------------------------------------------------

   do n = 1,nt
      if (tr_mask(n)) AUX(:,:,n,bid) = AUXB(:,:,n)
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine advt_lw_lim

!***********************************************************************
!BOP
! !IROUTINE: lw_lim
! !INTERFACE:

 subroutine lw_lim(k,adv_dt,UVEL_E_dt,VVEL_N_dt,XOUT,X, &
                      WTK,WTKB,WTKBp1,AUXB, &
                      CN,CS,CE,CW,KMASKE,KMASKN, &
                      TRACER_E,TRACER_N,bid,tr_mask,this_block)

! !DESCRIPTION:
!  This routine computes the advection of tracers using a 2nd-order
!  forward in time scheme with one-dimensional flux limiters.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      k,                 &! vertical level index
      bid                 ! local block address for this block

   real (r8), intent(in) :: &
      adv_dt              ! advection timestep

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      X                   ! input tracer array

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      UVEL_E_dt,         &! dt * velocity on east face
      VVEL_N_dt,         &! dt * velocity on north face
      WTK,WTKB,WTKBp1,   &! vert velocity at top of T box for levels k,k+1,k+2
      CN,CS,CE,CW,       &! stencil weights based on flux velocities
      KMASKE,KMASKN       ! land mask on east and north faces

   logical (log_kind), dimension(nt), intent(in) :: &
      tr_mask             ! true for tracers using lw_lim

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout)  :: &
      AUXB,              &! tracer flux through bottom of cell
      XOUT,              &! output tracer advection
      TRACER_E,          &! tracer value on east face
      TRACER_N            ! tracer value on north face

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nu                ,&! i/o unit number
      reclength         ,&! record length
      i,j,n               ! dummy loop indices

   real (r8) :: &
      adv_dt_r,          &
      dzk_div_dt,        &
      dzkp1_div_dt,      &
      dTRm1, dTR, dTRp1, &! tracer differences at cell faces
      psi,               &! fraction of tracer difference appearing in TRACER_E
      psi_dTR,           &! psi * dTR
      work1, work2, work3 ! temporary work space

   real (r8), dimension(nx_block,ny_block)  :: &
      DIV,               &! divergence of velocity field
      XSTAR,             &! tracer w/ advective tendency terms added
      MU_z,LW_z,         &! coefficients for z contribution
      MU_x,LW_x,         &! coefficients for x contribution
      MU_y,LW_y           ! coefficients for y contribution

!-----------------------------------------------------------------------
!  compute tracer independent coefficients
!-----------------------------------------------------------------------

   adv_dt_r = c1 / adv_dt

   if (partial_bottom_cells) then

      do j=this_block%jb-2,this_block%je+2
         do i=this_block%ib-2,this_block%ie+2
            if (WTKB(i,j) > c0) then
               LW_z(i,j) = (DZT(i,j,k+1,bid) - adv_dt * WTKB(i,j)) / &
                           (DZT(i,j,k,bid) + DZT(i,j,k+1,bid))
               if (WTKBp1(i,j) > c0) then
                  MU_z(i,j) = (DZT(i,j,k+1,bid) * adv_dt_r - WTKBp1(i,j)) / WTKB(i,j)
               else if (WTKBp1(i,j) < c0) then
                  MU_z(i,j) = -WTKBp1(i,j) / WTKB(i,j) * &
                     (DZT(i,j,k+1,bid) + adv_dt * WTKBp1(i,j)) / &
                     (DZT(i,j,k+1,bid) + DZT(i,j,k+2,bid))
               else
                  MU_z(i,j) = c0
               end if
            else if (WTKB(i,j) < c0) then
               LW_z(i,j) = (DZT(i,j,k,bid) + adv_dt * WTKB(i,j)) / &
                           (DZT(i,j,k,bid) + DZT(i,j,k+1,bid))
               if (WTK(i,j) < c0) then
                  MU_z(i,j) = -(DZT(i,j,k,bid) * adv_dt_r + WTK(i,j)) / WTKB(i,j)
               else if (WTK(i,j) > c0) then
                  MU_z(i,j) =  -WTK(i,j) / WTKB(i,j) * &
                     (DZT(i,j,k,bid) - adv_dt * WTK(i,j)) / &
                     (DZT(i,j,k-1,bid) + DZT(i,j,k,bid))
               else
                  MU_z(i,j) = c0
               end if
               if (WTK(i,j) > c0) then
               end if
            end if
         end do
      end do

   else

      dzk_div_dt = dz(k) * adv_dt_r
      if (k < km) dzkp1_div_dt = dz(k+1) * adv_dt_r

      work1 = dz(k) * p5_dz_ph_r(k)
      if (k < km) work2 = dz(k+1) * p5_dz_ph_r(k)
      work3 = adv_dt * p5_dz_ph_r(k)

      do j=this_block%jb-2,this_block%je+2
         do i=this_block%ib-2,this_block%ie+2
            if (WTKB(i,j) > c0) then
               LW_z(i,j) = work2 - work3 * WTKB(i,j)
               if (WTKBp1(i,j) > c0) then
                  MU_z(i,j) = (dzkp1_div_dt - WTKBp1(i,j)) / WTKB(i,j)
               else if (WTKBp1(i,j) < c0) then
                  MU_z(i,j) = -WTKBp1(i,j) / WTKB(i,j) * &
                     (dz(k+1) + adv_dt * WTKBp1(i,j)) * p5_dz_ph_r(k+1)
               else
                  MU_z(i,j) = c0
               end if
            else if (WTKB(i,j) < c0) then
               LW_z(i,j) = work1 + work3 * WTKB(i,j)
               if (WTK(i,j) < c0) then
                  MU_z(i,j) = -(dzk_div_dt + WTK(i,j)) / WTKB(i,j)
               else if (WTK(i,j) > c0) then
                  MU_z(i,j) =  -WTK(i,j) / WTKB(i,j) * &
                     (dz(k) - adv_dt * WTK(i,j)) * p5_dz_ph_r(k-1)
               else
                  MU_z(i,j) = c0
               end if
            end if
         end do
      end do

   endif

   do j=this_block%jb-2,this_block%je+2
      i=this_block%ib-2
      if (UVEL_E_dt(i,j) < c0) then
         LW_x(i,j) = (DXT(i+1,j,bid) + UVEL_E_dt(i,j)) * p5_DXT_ph_R(i,j,bid)
      end if
      do i=this_block%ib-1,this_block%ie
         if (UVEL_E_dt(i,j) > c0) then
            LW_x(i,j) = (DXT(i,j,bid) - UVEL_E_dt(i,j)) * p5_DXT_ph_R(i,j,bid)
            if (UVEL_E_dt(i-1,j) > c0) then
               MU_x(i,j) = (DXT(i,j,bid) - UVEL_E_dt(i-1,j)) / UVEL_E_dt(i,j)
            else
               MU_x(i,j) = c0
            end if
         else if (UVEL_E_dt(i,j) < c0) then
            LW_x(i,j) = (DXT(i+1,j,bid) + UVEL_E_dt(i,j)) * p5_DXT_ph_R(i,j,bid)
            if (UVEL_E_dt(i+1,j) < c0) then
               MU_x(i,j) = -(DXT(i+1,j,bid) + UVEL_E_dt(i+1,j)) / UVEL_E_dt(i,j)
            else
               MU_x(i,j) = c0
            end if
         else
            LW_x(i,j) = DXT(i,j,bid) * p5_DXT_ph_R(i,j,bid)
         end if
      end do
      i=this_block%ie+1
      if (UVEL_E_dt(i,j) > c0) then
         LW_x(i,j) = (DXT(i,j,bid) - UVEL_E_dt(i,j)) * p5_DXT_ph_R(i,j,bid)
      end if
      do i=this_block%ib-1,this_block%ie
         if ((UVEL_E_dt(i,j) > c0) .and. (UVEL_E_dt(i-1,j) < c0)) then
            MU_x(i,j) = -UVEL_E_dt(i-1,j)/UVEL_E_dt(i,j)*LW_x(i-1,j)
         end if
         if ((UVEL_E_dt(i,j) < c0) .and. (UVEL_E_dt(i+1,j) > c0)) then
            MU_x(i,j) = -UVEL_E_dt(i+1,j)/UVEL_E_dt(i,j)*LW_x(i+1,j)
         end if
      end do
   end do

   j=this_block%jb-2
   do i=this_block%ib,this_block%ie
      if (VVEL_N_dt(i,j) < c0) then
         LW_y(i,j) = (DYT(i,j+1,bid) + VVEL_N_dt(i,j)) * p5_DYT_ph_R(i,j,bid)
      end if
   end do
   do j=this_block%jb-1,this_block%je
      do i=this_block%ib,this_block%ie
         if (VVEL_N_dt(i,j) > c0) then
            LW_y(i,j) = (DYT(i,j,bid) - VVEL_N_dt(i,j)) * p5_DYT_ph_R(i,j,bid)
            if (VVEL_N_dt(i,j-1) > c0) then
               MU_y(i,j) = (DYT(i,j,bid) - VVEL_N_dt(i,j-1)) / VVEL_N_dt(i,j)
            else
               MU_y(i,j) = c0
            end if
         else if (VVEL_N_dt(i,j) < c0) then
            LW_y(i,j) = (DYT(i,j+1,bid) + VVEL_N_dt(i,j)) * p5_DYT_ph_R(i,j,bid)
            if (VVEL_N_dt(i,j+1) < c0) then
               MU_y(i,j) = -(DYT(i,j+1,bid) + VVEL_N_dt(i,j+1)) / VVEL_N_dt(i,j)
            else
               MU_y(i,j) = c0
            end if
         else
            LW_y(i,j) = DYT(i,j,bid) * p5_DYT_ph_R(i,j,bid)
         end if
      end do
   end do
   j=this_block%je+1
   do i=this_block%ib,this_block%ie
      if (VVEL_N_dt(i,j) > c0) then
         LW_y(i,j) = (DYT(i,j,bid) - VVEL_N_dt(i,j)) * p5_DYT_ph_R(i,j,bid)
      end if
   end do
   do j=this_block%jb-1,this_block%je
      do i=this_block%ib,this_block%ie
         if ((VVEL_N_dt(i,j) > c0) .and. (VVEL_N_dt(i,j-1) < c0)) then
            MU_y(i,j) = -VVEL_N_dt(i,j-1)/VVEL_N_dt(i,j)*LW_y(i,j-1)
         end if
         if ((VVEL_N_dt(i,j) < c0) .and. (VVEL_N_dt(i,j+1) > c0)) then
            MU_y(i,j) = -VVEL_N_dt(i,j+1)/VVEL_N_dt(i,j)*LW_y(i,j+1)
         end if
      end do
   end do

   if (partial_bottom_cells) then
      DIV = (WTK - WTKB) / DZT(:,:,k,bid) + CE + CW + CN + CS
   else
      DIV = (WTK - WTKB) * dzr(k) + CE + CW + CN + CS
   end if

!-----------------------------------------------------------------------
!  loop over tracers
!-----------------------------------------------------------------------

   do n=1,nt
      if (.not. tr_mask(n)) cycle

      call zero_ghost_cells(this_block,TRACER_E(:,:,n))
      call zero_ghost_cells(this_block,TRACER_N(:,:,n))
      call zero_ghost_cells(this_block,AUXB(:,:,n))

!-----------------------------------------------------------------------
!     Compute vertical contribution
!     Computations are done in the grix-x & grid-y halo regions
!
!     Compute grid-x direction contribution
!     Computations are done in the grid-y halo regions because the
!        grid-y computations, which use the grid-y halo region,
!        rely on the grid-x results.
!
!     j loops for z and x contributions fused to improve performance
!
!     j loop for y contribution cannot easily be fused
!        because of data dependencies
!-----------------------------------------------------------------------

      do j=this_block%jb-2,this_block%je+2
         do i=this_block%ib-2,this_block%ie+2
            if (k+1 <= KMT(i,j,bid)) then
               dTR = (X(i,j,k+1,n) - X(i,j,k,n))
               if (WTKB(i,j) > c0) then
                  AUXB(i,j,n) = WTKB(i,j) * X(i,j,k+1,n)
                  if (k+2 <= KMT(i,j,bid)) then
                     dTRp1 = (X(i,j,k+2,n) - X(i,j,k+1,n))
!!!                  if (dTR * dTRp1 > c0) then
!!!                     psi = min(LW_z(i,j), MU_z(i,j) * dTRp1 / dTR)
!!!                     AUXB(i,j,n) = WTKB(i,j) * (X(i,j,k+1,n) - psi * dTR)
!!!                  end if
                     if (dTR > c0 .and. dTRp1 > c0) then
                        psi_dTR = min(LW_z(i,j) * dTR, MU_z(i,j) * dTRp1)
                        AUXB(i,j,n) = WTKB(i,j) * (X(i,j,k+1,n) - psi_dTR)
                     else if (dTR < c0 .and. dTRp1 < c0) then
                        psi_dTR = max(LW_z(i,j) * dTR, MU_z(i,j) * dTRp1)
                        AUXB(i,j,n) = WTKB(i,j) * (X(i,j,k+1,n) - psi_dTR)
                     end if
                  end if
               else if (WTKB(i,j) < c0) then
                  AUXB(i,j,n) = WTKB(i,j) * X(i,j,k,n)
                  if (k > 1) then
                     dTRm1 = (X(i,j,k,n) - X(i,j,k-1,n))
!!!                  if (dTR * dTRm1 > c0) then
!!!                     psi = min(LW_z(i,j), MU_z(i,j) * dTRm1 / dTR)
!!!                     AUXB(i,j,n) = WTKB(i,j) * (X(i,j,k,n) + psi * dTR)
!!!                  end if
                     if (dTR > c0 .and. dTRm1 > c0) then
                        psi_dTR = min(LW_z(i,j) * dTR, MU_z(i,j) * dTRm1)
                        AUXB(i,j,n) = WTKB(i,j) * (X(i,j,k,n) + psi_dTR)
                     else if (dTR < c0 .and. dTRm1 < c0) then
                        psi_dTR = max(LW_z(i,j) * dTR, MU_z(i,j) * dTRm1)
                        AUXB(i,j,n) = WTKB(i,j) * (X(i,j,k,n) + psi_dTR)
                     end if
                  end if
               else
                  AUXB(i,j,n) = c0
               end if
            else ! (k+1>KMT)
               AUXB(i,j,n) = c0
            end if

            if (partial_bottom_cells) then
                XOUT(i,j,n) = ( AUX(i,j,n,bid) - AUXB(i,j,n) - &
                   (WTK(i,j) - WTKB(i,j)) * X(i,j,k,n) ) / DZT(i,j,k,bid)
            else
                XOUT(i,j,n) = ( AUX(i,j,n,bid) - AUXB(i,j,n) - &
                   (WTK(i,j) - WTKB(i,j)) * X(i,j,k,n) ) * dzr(k)
            end if
            XSTAR(i,j) = X(i,j,k,n) - adv_dt * XOUT(i,j,n)
         end do ! i loop for z contribution

         do i=this_block%ib-1,this_block%ie

             dTR = (XSTAR(i+1,j) - XSTAR(i,j)) * KMASKE(i,j)
             if (CE(i,j) > c0) then
                dTRm1 = (XSTAR(i,j) - XSTAR(i-1,j)) * KMASKE(i-1,j)
!!!             if (dTR * dTRm1 > c0) then
!!!                psi = min(LW_x(i,j), MU_x(i,j) * dTRm1 / dTR)
!!!                TRACER_E(i,j,n) = XSTAR(i,j) + psi * dTR
                if (dTR > c0 .and. dTRm1 > c0) then
                   psi_dTR = min(LW_x(i,j) * dTR, MU_x(i,j) * dTRm1)
                   TRACER_E(i,j,n) = XSTAR(i,j) + psi_dTR
                else if (dTR < c0 .and. dTRm1 < c0) then
                   psi_dTR = max(LW_x(i,j) * dTR, MU_x(i,j) * dTRm1)
                   TRACER_E(i,j,n) = XSTAR(i,j) + psi_dTR
                else
                   TRACER_E(i,j,n) = XSTAR(i,j)
                end if
             else if (CE(i,j) < c0) then
                dTRp1 = (XSTAR(i+2,j) - XSTAR(i+1,j)) * KMASKE(i+1,j)
!!!             if (dTR * dTRp1 > c0) then
!!!                psi = min(LW_x(i,j), MU_x(i,j) * dTRp1 / dTR)
!!!                TRACER_E(i,j,n) = XSTAR(i+1,j) - psi * dTR
                if (dTR > c0 .and. dTRp1 > c0) then
                   psi_dTR = min(LW_x(i,j) * dTR, MU_x(i,j) * dTRp1)
                   TRACER_E(i,j,n) = XSTAR(i+1,j) - psi_dTR
                else if (dTR < c0 .and. dTRp1 < c0) then
                   psi_dTR = max(LW_x(i,j) * dTR, MU_x(i,j) * dTRp1)
                   TRACER_E(i,j,n) = XSTAR(i+1,j) - psi_dTR
                else
                   TRACER_E(i,j,n) = XSTAR(i+1,j)
                end if
             else
                TRACER_E(i,j,n) = XSTAR(i,j) + LW_x(i,j) * dTR
             end if
          end do ! i loop for grid-x contribution

          !*** i loop broken up to avoid dependency on previous i
          !*** iteration and to avoid using an XSTAR value that
          !*** contains grid-x direction contribution

          end do ! j loop for z and grid-x contributions

          !***  j loop broken up for overflow TRACER_E modification

          if ( overflows_on .and. overflows_interactive ) then
             call ovf_advt(k,TRACER_E(:,:,n),TRACER_N(:,:,n),n,this_block, &
                           CE,CW,CN,CS)
          endif

          do j=this_block%jb-2,this_block%je+2
          do i=this_block%ib,this_block%ie
             work1 = &
                CE(i,j) * TRACER_E(i,j,n) + CW(i,j) * TRACER_E(i-1,j,n) - &
                (CE(i,j) + CW(i,j)) * X(i,j,k,n)
             XOUT(i,j,n) = XOUT(i,j,n) + work1
             XSTAR(i,j) = XSTAR(i,j) - adv_dt * work1
         end do ! i loop for grid-x contribution
      end do ! j loop for z and grid-x contributions

!-----------------------------------------------------------------------
!     Compute grid-y direction contribution
!     Add divergence term
!-----------------------------------------------------------------------

      do j=this_block%jb-1,this_block%je
         do i=this_block%ib,this_block%ie

             dTR = (XSTAR(i,j+1) - XSTAR(i,j)) * KMASKN(i,j)
             if (CN(i,j) > c0) then
                dTRm1 = (XSTAR(i,j) - XSTAR(i,j-1)) * KMASKN(i,j-1)
!!!             if (dTR * dTRm1 > c0) then
!!!                psi = min(LW_y(i,j), MU_y(i,j) * dTRm1 / dTR)
!!!                TRACER_N(i,j,n) = XSTAR(i,j) + psi * dTR
                if (dTR > c0 .and. dTRm1 > c0) then
                   psi_dTR = min(LW_y(i,j) * dTR, MU_y(i,j) * dTRm1)
                   TRACER_N(i,j,n) = XSTAR(i,j) + psi_dTR
                else if (dTR < c0 .and. dTRm1 < c0) then
                   psi_dTR = max(LW_y(i,j) * dTR, MU_y(i,j) * dTRm1)
                   TRACER_N(i,j,n) = XSTAR(i,j) + psi_dTR
                else
                   TRACER_N(i,j,n) = XSTAR(i,j)
                end if
             else if (CN(i,j) < c0) then
                dTRp1 = (XSTAR(i,j+2) - XSTAR(i,j+1)) * KMASKN(i,j+1)
!!!             if (dTR * dTRp1 > c0) then
!!!                psi = min(LW_y(i,j), MU_y(i,j) * dTRp1 / dTR)
!!!                TRACER_N(i,j,n) = XSTAR(i,j+1) - psi * dTR
                if (dTR > c0 .and. dTRp1 > c0) then
                   psi_dTR = min(LW_y(i,j) * dTR, MU_y(i,j) * dTRp1)
                   TRACER_N(i,j,n) = XSTAR(i,j+1) - psi_dTR
                else if (dTR < c0 .and. dTRp1 < c0) then
                   psi_dTR = max(LW_y(i,j) * dTR, MU_y(i,j) * dTRp1)
                   TRACER_N(i,j,n) = XSTAR(i,j+1) - psi_dTR
                else
                   TRACER_N(i,j,n) = XSTAR(i,j+1)
                end if
             else
                TRACER_N(i,j,n) = XSTAR(i,j) + LW_y(i,j) * dTR
             end if
         end do ! i loop for grid-y contribution
         end do ! j loop for grid-y contribution

        !***  i,j loops broken for overflow TRACER_N modification

         if ( overflows_on .and. overflows_interactive ) then
            call ovf_advt(k,TRACER_E(:,:,n),TRACER_N(:,:,n),n,this_block, &
                          CE,CW,CN,CS)
         endif

         do j=this_block%jb-1,this_block%je
         do i=this_block%ib,this_block%ie

             !*** The formula below is not correct for j==this_block%jb-1
             !*** because TRACER_N(:,this_block%jb-2,n) is not computed.
             !*** This is OK because XOUT(:,this_block%jb-1,n) is not used.

             XOUT(i,j,n) = XOUT(i,j,n) + &
               CN(i,j) * TRACER_N(i,j,n) + CS(i,j) * TRACER_N(i,j-1,n) - &
               (CN(i,j) + CS(i,j) - DIV(i,j)) * X(i,j,k,n)

             !*** do not update XSTAR in order to avoid using an XSTAR
             !*** value that contains grid-y direction contribution

         end do ! i loop for grid-y contribution
      end do ! j loop for grid-y contribution

   end do ! tracer loop

!-----------------------------------------------------------------------

 end subroutine lw_lim

!-----------------------------------------------------------------------
!EOC

!***********************************************************************
!BOP
! !IROUTINE: zero_ghost_cells
! !INTERFACE:

 subroutine zero_ghost_cells(this_block,FIELD)

! !DESCRIPTION:
!  This routine sets ghost cell values to zero.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block          ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(inout)  :: &
      FIELD               ! FIELD to be masked

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j                 ! dummy loop indices

!-----------------------------------------------------------------------

   do j=1,this_block%jb-1
      do i=1,nx_block
         FIELD(i,j) = c0
      end do
   end do

   do j=this_block%jb,this_block%je
      do i=1,this_block%ib-1
         FIELD(i,j) = c0
      end do
      do i=this_block%ie+1,nx_block
         FIELD(i,j) = c0
      end do
   end do

   do j=this_block%je+1,ny_block
      do i=1,nx_block
         FIELD(i,j) = c0
      end do
   end do

 end subroutine zero_ghost_cells

!-----------------------------------------------------------------------
!EOC

!***********************************************************************

 end module advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
