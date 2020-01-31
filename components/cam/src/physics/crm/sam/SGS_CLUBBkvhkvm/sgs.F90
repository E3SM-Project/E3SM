module sgs

! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
! Marat Khairoutdinov, 2012

use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s 
use params, only: dosgs
use vars, only: tke2, tk2
#ifdef CLUBB_CRM
use clubbvars, only: khzt, khzm
use params, only: doclubb
#endif
implicit none

!----------------------------------------------------------------------
! Required definitions:

!!! prognostic scalar (need to be advected arround the grid):

integer, parameter :: nsgs_fields = 1   ! total number of prognostic sgs vars

real sgs_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields)

!!! sgs diagnostic variables that need to exchange boundary information (via MPI):

#ifndef CLUBB_CRM
integer, parameter :: nsgs_fields_diag = 2   ! total number of diagnostic sgs vars
#else
integer, parameter :: nsgs_fields_diag = 4   ! total number of diagnostic sgs vars
#endif

! diagnostic fields' boundaries:
integer, parameter :: dimx1_d=0, dimx2_d=nxp1, dimy1_d=1-YES3D, dimy2_d=nyp1

real sgs_field_diag(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag)

logical:: advect_sgs = .false. ! advect prognostics or not, default - not (Smagorinsky)
logical, parameter:: do_sgsdiag_bound = .true.  ! exchange boundaries for diagnostics fields

! SGS fields that output by default (if =1).
integer, parameter :: flag_sgs3Dout(nsgs_fields) = (/0/)
#ifndef CLUBB_CRM
integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0/)
#else
integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0,0,0/)
#endif

real fluxbsgs (nx, ny, 1:nsgs_fields) ! surface fluxes 
real fluxtsgs (nx, ny, 1:nsgs_fields) ! top boundary fluxes 

!!! these arrays may be needed for output statistics:

real sgswle(nz,1:nsgs_fields)  ! resolved vertical flux
real sgswsb(nz,1:nsgs_fields)  ! SGS vertical flux
real sgsadv(nz,1:nsgs_fields)  ! tendency due to vertical advection
real sgslsadv(nz,1:nsgs_fields)  ! tendency due to large-scale vertical advection
real sgsdiff(nz,1:nsgs_fields)  ! tendency due to vertical diffusion

!------------------------------------------------------------------
! internal (optional) definitions:

! make aliases for prognostic variables:

real tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
equivalence (tke(dimx1_s,dimy1_s,1),sgs_field(dimx1_s,dimy1_s,1,1))

! make aliases for diagnostic variables:

real tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
real tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
equivalence (tk(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,1))
equivalence (tkh(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,2))
#ifdef CLUBB_CRM
real tk_clubb  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
real tkh_clubb (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
equivalence (tk_clubb(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,3))
equivalence (tkh_clubb(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,4))
#endif

real grdf_x(nzm)! grid factor for eddy diffusion in x
real grdf_y(nzm)! grid factor for eddy diffusion in y
real grdf_z(nzm)! grid factor for eddy diffusion in z

logical:: dosmagor   ! if true, then use Smagorinsky closure

! Local diagnostics:

real tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz)

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm (namelist) file

subroutine sgs_setparm()

  use grid, only: case
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  !======================================================================
  ! UW ADDITION
  NAMELIST /SGS_TKE/ &
       dosmagor ! Diagnostic Smagorinsky closure

  NAMELIST /BNCUIODSBJCB/ place_holder

  dosmagor = .true. ! default 

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  !open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

  !read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  !rewind(55) !note that one must rewind before searching for new namelists

  !read (55,SGS_TKE,IOSTAT=ios)

  advect_sgs = .not.dosmagor

  !if (ios.ne.0) then
  !   !namelist error checking
  !   if(ios.ne.ios_missing_namelist) then
  !      write(*,*) '****** ERROR: bad specification in SGS_TKE namelist'
  !      call task_abort()
  !   end if
  !end if
  !close(55)

  ! END UW ADDITION
  !======================================================================

end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


subroutine sgs_init()

  use grid, only: nrestart, dx, dy, dz, adz, masterproc
  use params, only: LES
#ifdef CLUBB_CRM
  use params, only: doclubb
#endif
  integer k

  if(nrestart.eq.0) then

     sgs_field = 0.
     sgs_field_diag = 0.

     fluxbsgs = 0.
     fluxtsgs = 0.

  end if

!  if(masterproc) then
!     if(dosmagor) then
!        write(*,*) 'Smagorinsky SGS Closure'
!     else
!        write(*,*) 'Prognostic TKE 1.5-order SGS Closure'
!     end if
!#ifdef CLUBB_CRM
!      if ( doclubb ) then
!        write(*,*) 'CLUBB Parameterization'
!      end if
!#endif
!  end if

  if(LES) then
    do k=1,nzm
       grdf_x(k) = dx**2/(adz(k)*dz)**2
       grdf_y(k) = dy**2/(adz(k)*dz)**2
       grdf_z(k) = 1.
    end do
  else
    do k=1,nzm
       grdf_x(k) = min(16.,dx**2/(adz(k)*dz)**2)
       grdf_y(k) = min(16.,dy**2/(adz(k)*dz)**2)
       grdf_z(k) = 1.
    end do
  end if

  sgswle = 0.
  sgswsb = 0.
  sgsadv = 0.
  sgsdiff = 0.
  sgslsadv = 0.


end subroutine sgs_init

!----------------------------------------------------------------------
!!! make some initial noise in sgs:
!
subroutine setperturb_sgs(ptype)

use vars, only: q0, z
integer, intent(in) :: ptype
integer i,j,k

select case (ptype)

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(k.le.4.and..not.dosmagor) then
            tke(i,j,k)=0.04*(5-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(2)

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.0.5e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.150..and..not.dosmagor) then
            tke(i,j,k)=0.15*(1.-z(k)/150.)
         endif
       end do
      end do
     end do


  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.3000..and..not.dosmagor) then
            tke(i,j,k)=1.-z(k)/3000.
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case default

end select

end subroutine setperturb_sgs

!----------------------------------------------------------------------
!!! Estimate Courant number limit for SGS
!

subroutine kurant_sgs(cfl)

use grid, only: dt, dx, dy, dz, adz, adzw
implicit none

real, intent(out) :: cfl

integer k
real tkhmax(nz)

do k = 1,nzm
 tkhmax(k) = maxval(tkh(1:nx,1:ny,k))
end do

cfl = 0.
do k=1,nzm
  cfl = max(cfl,        &
     0.5*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2, &
     0.5*tkhmax(k)*grdf_x(k)*dt/dx**2, &
     YES3D*0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
end do

end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
subroutine sgs_mom()
#ifdef CLUBB_CRM
  use params, only: doclubb
  use clubb_sgs, only: apply_clubb_sgs_tndcy_mom
  use vars, only: dudt, dvdt
#endif

#ifdef CLUBB_CRM
     if ( doclubb ) then
!          call apply_clubb_sgs_tndcy_mom &
!               ( dudt, dvdt ) ! in/out
     endif
#endif /*CLUBB_CRM*/

   call diffuse_mom()

end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
subroutine sgs_scalars()

  use vars
  use microphysics
  use crmtracers
  use params, only: dotracers, doclubb, doclubb_sfc_fluxes, doclubbnoninter, docam_sfc_fluxes
#ifdef CLUBB_CRM
  use clubbvars, only: edsclr_dim, sclr_dim
  use clubb_sgs, only: total_energy
  use clubb_sgs, only: apply_clubb_sgs_tndcy_scalars
  use grid, only: dtn
  use clubb_precision, only: time_precision
#endif  /*CLUBB_CRM*/
  implicit none

    real dummy(nz)
    real f2lediff_xy(nz), f2lediss_xy(nz), fwlediff_xy(nz)
    real f2lediff_z(nz), f2lediss_z(nz), fwlediff_z(nz)
    real sdiff_xy(nz), sdiff_z(nz)
    real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
    integer k


#ifdef CLUBB_CRM
      total_energy_evap = total_energy_evap - total_energy(t)
#endif

!  Update for t, qv, qcl from clubb_sgs
#ifdef CLUBB_CRM
     if ( doclubb ) then

      ! Recalculate q, qv, qcl based on new micro_fields (updated by horizontal
      ! diffusion)
       call micro_update()

      ! Then Re-compute q/qv/qcl based on values computed in CLUBB
       call apply_clubb_sgs_tndcy_scalars &
            ( real( dtn, kind=time_precision), & ! in
              t, qv, qcl) ! in/out

       call micro_adjust( qv, qcl ) ! in
     end if
#endif /*CLUBB_CRM*/

      f2lediff_xy = 0.0
      f2lediss_xy = 0.0
      fwlediff_xy = 0.0

!      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
!                           t2lediff,t2lediss,twlediff,.true.)
      call diffuse_scalar_xy(t,fluxbt,fluxtt,tdiff_xy,twsb, &
                           f2lediff_xy,f2lediss_xy,fwlediff_xy,.true.)
      f2lediff_z =0.0
      f2lediss_z =0.0
      fwlediff_z =0.0
#ifdef CLUBB_CRM     
      ! Diffuse moist static energy in the vertical only if CLUBB is not being
      ! called
      if ( .not. doclubb ) then
         call diffuse_scalar_z(t,fluxbt,fluxtt,tdiff_z,twsb, &
                             f2lediff_z,f2lediss_z,fwlediff_z,.true.)
      else  ! doclubb
         if(doclubb_sfc_fluxes .or. docam_sfc_fluxes) then
           ! The flux will be applied in advance_clubb_core, so the 2nd argument
           ! is zero.
            call fluxes_scalar_z(t,fzero,fluxtt,tdiff_z,twsb, &
                             f2lediff_z,f2lediss_z,fwlediff_z,.true.)
         else
            call fluxes_scalar_z(t,fluxbt,fluxtt,tdiff_z,twsb, &
                             f2lediff_z,f2lediss_z,fwlediff_z,.true.)
         end if
      end if
#else
      call diffuse_scalar_z(t,fluxbt,fluxtt,tdiff_z,twsb, &
                           f2lediff_z,f2lediss_z,fwlediff_z,.true.)
#endif
 
      tdiff = tdiff_xy + tdiff_z

      t2lediff = f2lediff_xy + f2lediff_z
      t2lediss = f2lediss_xy + f2lediss_z 
      twlediff = fwlediff_xy + fwlediff_z

#ifdef CLUBB_CRM
      total_energy_evap = total_energy_evap + total_energy(t)
#endif
    
      if(advect_sgs) then
!         call diffuse_scalar(tke,fzero,fzero,dummy,sgswsb, &
!                                    dummy,dummy,dummy,.false.)
         call diffuse_scalar_xy(tke,fzero,fzero,dummy,sgswsb, &
                                    dummy,dummy,dummy,.false.)
         call diffuse_scalar_z(tke,fzero,fzero,dummy,sgswsb, &
                                    dummy,dummy,dummy,.false.)
      end if


!
!    diffusion of microphysics prognostics:
!
      call micro_flux()

      total_water_evap = total_water_evap - total_water()

      do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
#ifdef CLUBB_CRM
        .or. ( docloud.or.doclubb.or.doclubbnoninter ).and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#else
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#endif

         .or. doprecip.and.flag_precip(k).eq.1 ) then

           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
           sdiff_xy = 0.0
           sdiff_z = 0.0

!           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
!                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
           call diffuse_scalar_xy(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                sdiff_xy,mkwsb(:,k), dummy,dummy,dummy,.false.)
           if(k.ne.index_water_vapor) then
              call diffuse_scalar_z(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                 sdiff_z,mkwsb(:,k), dummy,dummy,dummy,.false.)
           else  ! k==index_water_vapor 
             if(.not. doclubb) then
                call diffuse_scalar_z(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                   sdiff_z,mkwsb(:,k), dummy,dummy,dummy,.false.)
             else   ! doclubb
                call fluxes_scalar_z(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                   sdiff_z,mkwsb(:,k), dummy,dummy,dummy,.false.)
             end if
           end if
           mkdiff(:, k) = sdiff_xy + sdiff_z
       end if
      end do

      total_water_evap = total_water_evap + total_water()

 ! diffusion of tracers:

      if(dotracers) then

        call tracers_flux()

        do k = 1,ntracers

#ifdef CLUBB_CRM
          ! If CLUBB is using the high-order or eddy diffusivity scalars, then
          ! we should apply the flux within advance_clubb_core when
          ! doclubb_sfc_fluxes is set to true. -dschanen UWM 2 Mar 2010
          if ( ( edsclr_dim > 0 .or. sclr_dim > 0 ) .and. (doclubb_sfc_fluxes .or. docam_sfc_fluxes)) then
            fluxbtmp = 0. ! Apply surface flux in CLUBB
          else
            fluxbtmp = fluxbtr(:,:,k)
          end if
#else
          fluxbtmp = fluxbtr(:,:,k)
#endif /*CLUBB_CRM*/
          fluxttmp = fluxttr(:,:,k)
!          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
!               trdiff(:,k),trwsb(:,k), &
!               dummy,dummy,dummy,.false.)
          call diffuse_scalar_xy(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)

#ifdef CLUBB_CRM
          ! Only diffuse the tracers if CLUBB is either disabled or using the
          ! eddy scalars code to diffuse them.
          if ( .not. doclubb .or. ( doclubb .and. edsclr_dim < 1 .and. sclr_dim < 1 ) ) then
            call diffuse_scalar_z(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
                 trdiff(:,k),trwsb(:,k), &
                 dummy,dummy,dummy,.false.)
          end if
#else
          call diffuse_scalar_z(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)
#endif
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)

        end do

      end if



end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc()

   use grid, only: nstep,dt,icycle
   use params, only: dosmoke
#ifdef CLUBB_CRM
   use clubbvars, only: khzt, khzm
   use microphysics 
   use params, only: doclubb, doclubbnoninter, nclubb
   use grid, only: dtn, time, dt
   use vars, only: u, v, w, rho, rhow, wsub, qpl, qci, qpi, t, qv, qcl
   use clubb_precision, only: time_precision
   use clubb_sgs, only: advance_clubb_sgs
#endif

! SGS CLUBB
#ifdef CLUBB_CRM
     if ( doclubb .or. doclubbnoninter ) then
        ! In case of ice fall, we recompute qci here for the 
        ! single-moment scheme.  Also, subsidence, diffusion and advection have
        ! been applied to micro_field but not qv/qcl so they must be updated.
        call micro_update()

        ! We call CLUBB here because adjustments to the wind
        ! must occur prior to adams() -dschanen 26 Aug 2008
        ! Here we call clubb only if nstep divides the current timestep,
        ! or we're on the very first timestep

! in the case with m2005, clubb is only called in the first subscycle (icycle=1))
        if ( ((nstep == 1 .or. mod( nstep, nclubb ) == 0) .and.  &
          (icycle == 1)).and.(nclubb .ne. 1) ) then  ! call every CRM step, so dt is used 
          call advance_clubb_sgs &
               ( real( dt*real( nclubb ), kind=time_precision), & ! in
                 real( 0., kind=time_precision ),         & ! in
                 real( time, kind=time_precision ),       & ! in
                 rho, rhow, wsub, u, v, w, qpl, qci, qpi, & ! in
                 t, qv, qcl ) ! in
        else if(nclubb.eq.1) then   ! call every icycle, so dtn is used 
          call advance_clubb_sgs &
               ( real( dtn*real( nclubb ), kind=time_precision), & ! in
                 real( 0., kind=time_precision ),         & ! in
                 real( time, kind=time_precision ),       & ! in
                 rho, rhow, wsub, u, v, w, qpl, qci, qpi, & ! in
                 t, qv, qcl ) ! in
        end if ! nstep == 1 .or. mod( nstep, nclubb) == 0

      end if ! doclubb .or. doclubbnoninter
#endif 

!    SGS TKE equation:

     if(dosgs) call tke_full()

     tke2 = tke
     tk2 = tk

#ifdef CLUBB_CRM
     if(doclubb) then
!       tk = khzt
!       tkh = khzt

!       tk_clubb = khzt
!       tkh_clubb = khzt
       tk_clubb = khzm
       tkh_clubb = khzm
     end if
#endif


end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
! None 

end subroutine sgs_diagnose

!----------------------------------------------------------------------
! called when stepout() called

subroutine sgs_print()

 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)

end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics 
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,sgscount

character*8 name
character*80 longname
character*10 units

#ifdef CLUBB
if (doclubb) then
name = 'TKCLUBB'
longname = 'Eddy diffusivity from CLUBB'
units = 'm2/s'
call add_to_namelist(count,sgscount,name,longname,units,0)

name = 'TKHCLUBB'
longname = 'Eddy diffusivity from CLUBB'
units = 'm2/s'
call add_to_namelist(count,sgscount,name,longname,units,0)
end if
#endif

end subroutine sgs_hbuf_init


end module sgs



