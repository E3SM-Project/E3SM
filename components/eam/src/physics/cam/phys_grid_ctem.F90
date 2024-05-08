!----------------------------------------------------------------------------------
! circulation diagnostics -- terms of the Transformed Eulerian Mean (TEM) equation
!
! This module makes use of the zonal_mean_mod methods to produce terms needed 
! for the TEM budget from the physics grid state. Without this online 
! calculation of terms, high-frequency 3D history output would be required to 
! calculate the terms offline. The resulting zonal mean quantities and eddy
! covariance terms take up much less disk space and can be output monthly
! 
! History output variables:
!   Uzm      Zonal-Mean zonal wind
!   Vzm      Zonal-Mean meridional wind
!   Wzm      Zonal-Mean vertical wind
!   THzm     Zonal-Mean potential temp
!   VTHzm    Meridional Heat Flux
!   WTHzm    Vertical Heat Flux
!   UVzm     Meridional Flux of Zonal Momentum
!   UWzm     Vertical Flux of Zonal Momentum
!   THphys   Potential temp
!----------------------------------------------------------------------------------
module phys_grid_ctem
use shr_kind_mod,  only: r8 => shr_kind_r8
use ppgrid,        only: begchunk, endchunk, pcols, pver
use physics_types, only: physics_state
use cam_history,   only: addfld, outfld, horiz_only
use zonal_mean_mod,only: ZonalAverage_t, ZonalMean_t
use physconst,     only: pi
use cam_logfile,   only: iulog
use cam_abortutils,only: endrun, handle_allocate_error
use namelist_utils,only: find_group_name
use spmd_utils,    only: masterproc, mpi_integer, masterprocid, mpicom
use time_manager,  only: get_nstep

use shr_const_mod, only: rgas => shr_const_rgas ! J/K/kmole
use shr_const_mod, only: grav => shr_const_g ! m/s2
use string_utils,  only: int2str

implicit none

private
public :: phys_grid_ctem_readnl
public :: phys_grid_ctem_reg
public :: phys_grid_ctem_init
public :: phys_grid_ctem_diags
public :: phys_grid_ctem_final

type(ZonalMean_t) :: ZMobj
type(ZonalAverage_t) :: ZAobj

integer :: nzalat = -huge(1)
integer :: nzmbas = -huge(1)

integer :: ntimesteps = -huge(1) ! number of time steps bewteen TEM calculations

logical :: do_tem_diags = .false.

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine phys_grid_ctem_readnl(nlfile,dtime)
   character(len=*), intent(in) :: nlfile
   integer,          intent(in) :: dtime   ! Step time in seconds
   integer :: ierr, unitn

   character(len=*), parameter :: prefix = 'phys_grid_ctem_readnl: '

   integer :: phys_grid_ctem_zm_nbas
   integer :: phys_grid_ctem_za_nlat
   integer :: phys_grid_ctem_nfreq

   namelist /phys_grid_ctem_opts/ phys_grid_ctem_zm_nbas, phys_grid_ctem_za_nlat, phys_grid_ctem_nfreq

   phys_grid_ctem_zm_nbas = 0
   phys_grid_ctem_za_nlat = 0
   phys_grid_ctem_nfreq = 0

   ! Read in namelist values
   !------------------------
   if(masterproc) then
      open(newunit=unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'phys_grid_ctem_opts', status=ierr)
      if(ierr == 0) then
         read(unitn,phys_grid_ctem_opts,iostat=ierr)
         if(ierr /= 0) then
            call endrun(prefix//'ERROR reading namelist')
         end if
      end if
      close(unitn)
   end if

   call MPI_bcast(phys_grid_ctem_zm_nbas, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(prefix//'FATAL: mpi_bcast: phys_grid_ctem_zm_nbas')
   call MPI_bcast(phys_grid_ctem_za_nlat, 1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(prefix//'FATAL: mpi_bcast: phys_grid_ctem_za_nlat')
   call MPI_bcast(phys_grid_ctem_nfreq,   1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(prefix//'FATAL: mpi_bcast: phys_grid_ctem_nfreq')

   do_tem_diags = .false.
   if (phys_grid_ctem_nfreq/=0) then
      if (.not.(phys_grid_ctem_zm_nbas>0 .and. phys_grid_ctem_za_nlat>0)) then
         call endrun(prefix//'inconsistent phys_grid_ctem namelist settings -- phys_grid_ctem_zm_nbas=' &
                     //int2str(phys_grid_ctem_zm_nbas)//', phys_grid_ctem_za_nlat='//int2str(phys_grid_ctem_za_nlat))
      end if
      if (phys_grid_ctem_nfreq>0) then
         ntimesteps = phys_grid_ctem_nfreq
      else
         ntimesteps = nint( -phys_grid_ctem_nfreq*3600._r8/dtime )
      end if
      if (ntimesteps<1) then
         call endrun(prefix//'invalid ntimesteps -- phys_grid_ctem_nfreq needs to be a larger negative value ' &
                           //'or the model time step needs to be shorter')
      end if
      do_tem_diags = .true.
   end if

   if (masterproc) then
      if (do_tem_diags) then
         write(iulog,*) 'TEM diagnostics will be calculated every ',ntimesteps,' time steps'
         write(iulog,*) ' phys_grid_ctem_zm_nbas = ', phys_grid_ctem_zm_nbas
         write(iulog,*) ' phys_grid_ctem_za_nlat = ', phys_grid_ctem_za_nlat
         write(iulog,*) ' phys_grid_ctem_nfreq = ', phys_grid_ctem_nfreq
      else
         write(iulog,*) 'TEM diagnostics will not be performed'
      end if
   endif

   if (do_tem_diags) then
      nzalat = phys_grid_ctem_za_nlat
      nzmbas = phys_grid_ctem_zm_nbas
   end if

end subroutine phys_grid_ctem_readnl

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine phys_grid_ctem_reg

   use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register

   type(horiz_coord_t), pointer :: zalon_coord
   type(horiz_coord_t), pointer :: zalat_coord
   integer(iMap),       pointer :: grid_map(:,:)

   real(r8) :: zalats(nzalat)
   real(r8) :: area(nzalat)
   real(r8) :: zalons(1)
   real(r8) :: dlatrad, dlatdeg, lat1, lat2
   real(r8) :: total_area
   real(r8) :: total_wght
   integer :: j, astat

   real(r8), parameter :: latdeg0 = -90._r8
   real(r8), parameter :: latrad0 = -pi*0.5_r8
   real(r8), parameter :: fourpi = pi*4._r8

   integer, parameter :: ctem_zavg_phys_decomp = 333 ! Must be unique within CAM

   if (.not.do_tem_diags) return

   nullify(zalat_coord)
   nullify(zalon_coord)
   nullify(grid_map)

   zalons(1) = 0._r8

   dlatrad = pi/real(nzalat,kind=r8)
   dlatdeg = 180._r8/real(nzalat,kind=r8)
   total_area = 0._r8
   total_wght = 0._r8

   ! calculate latitudes and areas of zonal average grid boxes
   do j = 1,nzalat
      zalats(j) = latdeg0 + (real(j,kind=r8)-0.5_r8)*dlatdeg
      lat1 = latrad0 + real(j-1,kind=r8)*dlatrad
      lat2 = latrad0 + real(j  ,kind=r8)*dlatrad
      area(j) = 2._r8*pi*(sin(lat2)-sin(lat1))
      total_area = total_area + area(j)
      total_wght = total_wght + 0.5_r8*(sin(lat2)-sin(lat1))
   end do

   ! sanity check
   if ( abs(1._r8-total_wght)>1.e-12_r8 .or. abs(fourpi-total_area)>1.e-12_r8 ) then
      call endrun('phys_grid_ctem_reg: problem with area/wght calc')
   end if

   ! initialize zonal-average and zonal-mean utility objects
   call ZAobj%init(zalats,area,nzalat,GEN_GAUSSLATS=.false.)
   call ZMobj%init(nzmbas)

   ! Zonal average grid for history fields

   zalat_coord => horiz_coord_create('zalat', '', nzalat, 'latitude', 'degrees_north', 1, nzalat, zalats)
   zalon_coord => horiz_coord_create('zalon', '', 1, 'longitude', 'degrees_east', 1, 1, zalons)

   ! grid decomposition map
   allocate(grid_map(4,nzalat), stat=astat)
   call handle_allocate_error(astat, 'phys_grid_ctem_reg', 'grid_map')

   do j = 1,nzalat
      grid_map(1,j) = 1
      grid_map(2,j) = j
      if (masterproc) then
         grid_map(3,j) = 1
         grid_map(4,j) = j
      else
         grid_map(3,j) = 0
         grid_map(4,j) = 0
      end if
   end do

   ! register the zonal average grid
   call cam_grid_register('ctem_zavg_phys', ctem_zavg_phys_decomp, zalat_coord, zalon_coord, grid_map, &
                          unstruct=.false., zonal_grid=.true.)

end subroutine phys_grid_ctem_reg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine phys_grid_ctem_init

   if (.not.do_tem_diags) return

   call addfld ('PSzm',horiz_only, 'A','m s-1',  'Zonal-Mean surface pressure', gridname='ctem_zavg_phys' )
   call addfld ('Uzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean zonal wind', gridname='ctem_zavg_phys' )
   call addfld ('Vzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean meridional wind', gridname='ctem_zavg_phys' )
   call addfld ('Wzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean vertical wind', gridname='ctem_zavg_phys' )
   call addfld ('THzm', (/'lev'/), 'A','K',      'Zonal-Mean potential temp', gridname='ctem_zavg_phys' )
   call addfld ('VTHzm',(/'lev'/), 'A','K m s-1','Meridional Heat Flux:', gridname='ctem_zavg_phys')
   call addfld ('WTHzm',(/'lev'/), 'A','K m s-1','Vertical Heat Flux:', gridname='ctem_zavg_phys')
   call addfld ('UVzm', (/'lev'/), 'A','m2 s-2', 'Meridional Flux of Zonal Momentum', gridname='ctem_zavg_phys')
   call addfld ('UWzm', (/'lev'/), 'A','m2 s-2', 'Vertical Flux of Zonal Momentum', gridname='ctem_zavg_phys')
   call addfld ('THphys',(/'lev'/), 'A', 'K',    'Potential temp', gridname='physgrid' )

end subroutine phys_grid_ctem_init

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine phys_grid_ctem_diags(phys_state)
   use physconst,     only: rair, cpair

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

   character(len=*), parameter :: prefix = 'phys_grid_ctem_diags: '

   real(r8) :: ps(pcols,begchunk:endchunk)

   real(r8) :: u(pcols,pver,begchunk:endchunk)
   real(r8) :: v(pcols,pver,begchunk:endchunk)
   real(r8) :: w(pcols,pver,begchunk:endchunk)

   real(r8) :: pszm(pcols,begchunk:endchunk)

   real(r8) :: uzm(pcols,pver,begchunk:endchunk)
   real(r8) :: vzm(pcols,pver,begchunk:endchunk)
   real(r8) :: wzm(pcols,pver,begchunk:endchunk)

   real(r8) :: ud(pcols,pver,begchunk:endchunk)
   real(r8) :: vd(pcols,pver,begchunk:endchunk)
   real(r8) :: wd(pcols,pver,begchunk:endchunk)
   real(r8) :: thd(pcols,pver,begchunk:endchunk)

   real(r8) :: uvp(pcols,pver,begchunk:endchunk)
   real(r8) :: uwp(pcols,pver,begchunk:endchunk)
   real(r8) :: vthp(pcols,pver,begchunk:endchunk)
   real(r8) :: wthp(pcols,pver,begchunk:endchunk)

   integer  :: lchnk, ncol, j, k

   ! potential temperature
   real(r8) :: theta(pcols,pver,begchunk:endchunk)
   real(r8) :: thzm(pcols,pver,begchunk:endchunk)

   real(r8) :: uvza(nzalat,pver)
   real(r8) :: uwza(nzalat,pver)
   real(r8) :: vthza(nzalat,pver)
   real(r8) :: wthza(nzalat,pver)

   real(r8) :: psza(nzalat)
   real(r8) :: uza(nzalat,pver)
   real(r8) :: vza(nzalat,pver)
   real(r8) :: wza(nzalat,pver)
   real(r8) :: thza(nzalat,pver)

   ! In CESM/WACCM the variable mbarv is provided by the "air_composition" 
   ! module, which is not in E3SM, so we just use a rough approximation
   real(r8) :: mbarv = 28.97 ! molecular weight of dry air (g/mol)

   if (.not.do_calc()) return

   do lchnk = begchunk,endchunk
      ncol = phys_state(lchnk)%ncol
      ! surface pressure
      ps(:ncol,lchnk) = phys_state(lchnk)%ps(:ncol)
      ! potential temperature
      theta(:ncol,:,lchnk) = phys_state(lchnk)%t(:ncol,:) * ( 1000e2 / phys_state(lchnk)%pmid(:ncol,:) )**(rair/cpair)
      ! vertical pressure velocity
      w(:ncol,:,lchnk) = phys_state(lchnk)%omega(:ncol,:)
      ! horizontal velocity
      u(:ncol,:,lchnk) =  phys_state(lchnk)%u(:ncol,:)
      v(:ncol,:,lchnk) =  phys_state(lchnk)%v(:ncol,:)
   end do

   ! zonal means evaluated on the physics grid (3D) to be used in the deviations calculation below
   pszm(:,:)   = zmean_fld_2D(ps(:,:))
   uzm(:,:,:)  = zmean_fld_3D(u(:,:,:))
   vzm(:,:,:)  = zmean_fld_3D(v(:,:,:))
   wzm(:,:,:)  = zmean_fld_3D(w(:,:,:))
   thzm(:,:,:) = zmean_fld_3D(theta(:,:,:))

   ! diagnostic output
   do lchnk = begchunk, endchunk
      call outfld( 'THphys', theta(:,:,lchnk), pcols, lchnk)
   end do

   do lchnk = begchunk,endchunk
      ncol = phys_state(lchnk)%ncol
      do k = 1,pver
         ! zonal deviations
         thd(:ncol,k,lchnk) = theta(:ncol,k,lchnk) - thzm(:ncol,k,lchnk)
         ud(:ncol,k,lchnk)  = u(:ncol,k,lchnk)     - uzm(:ncol,k,lchnk)
         vd(:ncol,k,lchnk)  = v(:ncol,k,lchnk)     - vzm(:ncol,k,lchnk)
         wd(:ncol,k,lchnk)  = w(:ncol,k,lchnk)     - wzm(:ncol,k,lchnk)
         ! fluxes
         uvp(:ncol,k,lchnk)  = ud(:ncol,k,lchnk) * vd(:ncol,k,lchnk)
         uwp(:ncol,k,lchnk)  = ud(:ncol,k,lchnk) * wd(:ncol,k,lchnk)
         vthp(:ncol,k,lchnk) = vd(:ncol,k,lchnk) * thd(:ncol,k,lchnk)
         wthp(:ncol,k,lchnk) = wd(:ncol,k,lchnk) * thd(:ncol,k,lchnk)
      end do
   end do

   ! evaluate and output fluxes on the zonal-average grid
   call ZAobj%binAvg(uvp,  uvza)
   call ZAobj%binAvg(uwp,  uwza)
   call ZAobj%binAvg(vthp, vthza)
   call ZAobj%binAvg(wthp, wthza)


   if (any(abs(uvza)>1.e20_r8))  call endrun(prefix//'bad values in uvza')
   if (any(abs(uwza)>1.e20_r8))  call endrun(prefix//'bad values in uwza')
   if (any(abs(vthza)>1.e20_r8)) call endrun(prefix//'bad values in vthza')
   if (any(abs(wthza)>1.e20_r8)) call endrun(prefix//'bad values in wthza')

   call ZAobj%binAvg(pszm, psza)
   call ZAobj%binAvg(uzm,  uza)
   call ZAobj%binAvg(vzm,  vza)
   call ZAobj%binAvg(wzm,  wza)
   call ZAobj%binAvg(thzm, thza)

   if (any(abs(psza)>1.e20_r8)) call endrun(prefix//'bad values in psza')
   if (any(abs(uza)>1.e20_r8))  call endrun(prefix//'bad values in uza')
   if (any(abs(vza)>1.e20_r8))  call endrun(prefix//'bad values in vza')
   if (any(abs(wza)>1.e20_r8))  call endrun(prefix//'bad values in wza')
   if (any(abs(thza)>1.e20_r8)) call endrun(prefix//'bad values in thza')

   ! diagnostic output
   do j = 1,nzalat
      call outfld('PSzm',psza(j),1,j)
      call outfld('Uzm',uza(j,:),1,j)
      call outfld('Vzm',vza(j,:),1,j)
      call outfld('Wzm',wza(j,:),1,j)
      call outfld('THzm',thza(j,:),1,j)
      call outfld('UVzm',uvza(j,:),1,j)
      call outfld('UWzm',uwza(j,:),1,j)
      call outfld('VTHzm',vthza(j,:),1,j)
      call outfld('WTHzm',wthza(j,:),1,j)
   end do

   contains

   !------------------------------------------------------------------------------
   ! utility function for evaluating 3D zonal mean fields
   !------------------------------------------------------------------------------
   function zmean_fld_3D( fld ) result(fldzm)
      real(r8), intent(in) :: fld(pcols,pver,begchunk:endchunk)
      real(r8) :: fldzm(pcols,pver,begchunk:endchunk)
      real(r8) :: Zonal_Bamp3d(nzmbas,pver)
      call ZMobj%calc_amps(fld,Zonal_Bamp3d)
      call ZMobj%eval_grid(Zonal_Bamp3d,fldzm)
   end function zmean_fld_3D

   !------------------------------------------------------------------------------
   ! utility function for evaluating 2D zonal mean fields
   !------------------------------------------------------------------------------
   function zmean_fld_2D( fld ) result(fldzm)
      real(r8), intent(in) :: fld(pcols,begchunk:endchunk)
      real(r8) :: fldzm(pcols,begchunk:endchunk)
      real(r8) :: Zonal_Bamp2d(nzmbas)
      call ZMobj%calc_amps(fld,Zonal_Bamp2d)
      call ZMobj%eval_grid(Zonal_Bamp2d,fldzm)
   end function zmean_fld_2D

   !------------------------------------------------------------------------------
   ! utility function returns TRUE when time to update TEM diags
   !------------------------------------------------------------------------------
   logical function do_calc()
      integer :: nstep
      nstep = get_nstep()
      do_calc = do_tem_diags .and. mod(nstep,ntimesteps) == 0
   end function do_calc

end subroutine phys_grid_ctem_diags

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine phys_grid_ctem_final
   call ZAobj%final()
   call ZMobj%final()
end subroutine phys_grid_ctem_final

end module phys_grid_ctem
