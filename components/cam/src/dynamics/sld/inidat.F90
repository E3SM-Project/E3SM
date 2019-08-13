module inidat
!----------------------------------------------------------------------- 
! 
! Purpose: Read initial dataset and process fields as appropriate
!
! Method: Read and process one field at a time
! 
! Author: J. Olson  May 2004
!
! Modified: A. Gettelman and C. Craig Nov 2010 - put micro/macro physics into separate routines
! 
!-----------------------------------------------------------------------

   use pmgrid,          only: plon, plat, plev, beglat, endlat, plnlv, plevp
   use rgrid,           only: nlon
   use prognostics,     only : ps, u3, v3, t3, tm, tl, ed1, div, ps, phis, &
        phism, phisl, q3, qm, dpsm, dpsl, ql
   use ncdio_atm,       only: infld
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use cam_abortutils  ,    only: endrun
   use phys_grid,       only: get_ncols_p
   use spmd_utils,      only: masterproc, npes
   use cam_control_mod, only: ideal_phys, aqua_planet, moist_physics, adiabatic
#if ( defined SPMD )
   use mpishorthand,    only: mpir8, mpicom
#endif
   use cam_logfile,     only: iulog
   use pio,             only : file_desc_t, pio_inq_attlen, pio_noerr, pio_inq_varid, &
                               pio_get_att, pio_seterrorhandling,pio_bcast_error, pio_internal_error

   implicit none

PRIVATE

   real(r8), allocatable :: ps_tmp  (:,:  )
   real(r8), allocatable :: phis_tmp(:,:  )
   real(r8), allocatable :: q3_tmp  (:,:,:)
   real(r8), allocatable :: t3_tmp  (:,:,:)
   real(r8), allocatable :: arr3d_a (:,:,:)
   real(r8), allocatable :: arr3d_b (:,:,:)

   logical readvar            ! inquiry flag:  true => variable exists on netCDF file

   public :: read_inidat


contains




  subroutine global_int
!
!-----------------------------------------------------------------------
!
! Purpose:
! Compute global integrals of mass, moisture and geopotential height
! and fix mass of atmosphere
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
    use commap
    use physconst,    only: gravit
    use hycoef, only : hyai, ps0
#if ( defined SPMD )
    use mpishorthand
    use spmd_dyn, only: compute_gsfactors
#endif
    use sld_control_mod, only: pdela, tmassf, fixmas, qmass1, qmass2, qmassf, zgsint, tmass0
!
!---------------------------Local workspace-----------------------------
!
    integer i,k,lat,ihem,irow  ! grid indices
    real(r8) pdelb(plon,plev)  ! pressure diff between interfaces
                               ! using "B" part of hybrid grid only
    real(r8) pssum             ! surface pressure sum
    real(r8) dotproda          ! dot product
    real(r8) dotprodb          ! dot product
    real(r8) zgssum            ! partial sums of phis
    real(r8) hyad (plev)       ! del (A)
    real(r8) tmassf_tmp        ! Global mass integral
    real(r8) qmass1_tmp        ! Partial Global moisture mass integral
    real(r8) qmass2_tmp        ! Partial Global moisture mass integral
    real(r8) qmassf_tmp        ! Global moisture mass integral
    real(r8) zgsint_tmp        ! Geopotential integral

#if ( defined SPMD )
    integer :: numperlat         ! number of values per latitude band
    integer :: numsend(0:npes-1) ! number of items to be sent
    integer :: numrecv           ! number of items to be received
    integer :: displs(0:npes-1)  ! displacement array
#endif
!
!-----------------------------------------------------------------------
!
    if(masterproc) then
!        
! Initialize mass and moisture integrals for summation
! in a third calculation loop (assures bit-for-bit compare
! with non-random history tape).
!
       tmassf_tmp = 0._r8
       qmass1_tmp = 0._r8
       qmass2_tmp = 0._r8
       zgsint_tmp = 0._r8
!
! Compute pdel from "A" portion of hybrid vertical grid for later use in global integrals
!
       do k = 1,plev
          hyad(k) = hyai(k+1) - hyai(k)
       end do
       do k = 1,plev
          do i = 1,plon
             pdela(i,k) = hyad(k)*ps0
          end do
       end do
!
! Compute integrals of mass, moisture, and geopotential height
!
       do irow = 1,plat/2
          do ihem = 1,2
             if (ihem.eq.1) then
                lat = irow
             else
                lat = plat - irow + 1
             end if
!              
! Accumulate average mass of atmosphere
!
             call pdelb0 (ps_tmp(1,lat),pdelb   ,nlon(lat))
             pssum  = 0._r8
             zgssum = 0._r8
             do i = 1,nlon(lat)
                pssum  = pssum  + ps_tmp  (i,lat)
                zgssum = zgssum + phis_tmp(i,lat)
             end do
             tmassf_tmp = tmassf_tmp + w(irow)*pssum/nlon(lat)
             zgsint_tmp = zgsint_tmp + w(irow)*zgssum/nlon(lat)
!
! Calculate global integrals needed for water vapor adjustment
!
             do k = 1,plev
                dotproda = 0._r8
                dotprodb = 0._r8
                do i = 1,nlon(lat)
                   dotproda = dotproda + q3_tmp(i,k,lat)*pdela(i,k)
                   dotprodb = dotprodb + q3_tmp(i,k,lat)*pdelb(i,k)
                end do
                qmass1_tmp = qmass1_tmp + w(irow)*dotproda/nlon(lat)
                qmass2_tmp = qmass2_tmp + w(irow)*dotprodb/nlon(lat)
             end do
          end do
       end do                  ! end of latitude loop
!
! Normalize average mass, height
!
       tmassf_tmp = tmassf_tmp*.5_r8/gravit
       qmass1_tmp = qmass1_tmp*.5_r8/gravit
       qmass2_tmp = qmass2_tmp*.5_r8/gravit
       zgsint_tmp = zgsint_tmp*.5_r8/gravit
       qmassf_tmp = qmass1_tmp + qmass2_tmp
!
! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
!
       tmass0 = 98222._r8/gravit
       if (adiabatic)   tmass0 =  tmassf_tmp
       if (ideal_phys ) tmass0 =  100000._r8/gravit
       if (aqua_planet) tmass0 = (101325._r8-245._r8)/gravit
       write(iulog,800) tmassf_tmp,tmass0,qmassf_tmp
       write(iulog,810) zgsint_tmp
800    format(/72('*')//'INIDAT: Mass of initial data before correction = ' &
              ,1p,e20.10,/,' Dry mass will be held at = ',e20.10,/, &
              ' Mass of moisture after removal of negatives = ',e20.10) 
810    format('INIDAT: Globally averaged geopotential ', &
              'height = ',f16.10,' meters'//72('*')/)
!
! Compute and apply an initial mass fix factor which preserves horizontal
! gradients of ln(ps).
!
       if (.not. moist_physics) then
          fixmas = tmass0/tmassf_tmp
       else
          fixmas = (tmass0 + qmass1_tmp)/(tmassf_tmp - qmass2_tmp)
       end if
       do lat = 1,plat
          do i = 1,nlon(lat)
             ps_tmp(i,lat) = ps_tmp(i,lat)*fixmas
          end do
       end do
!
! Global integerals
!
       tmassf = tmassf_tmp
       qmass1 = qmass1_tmp
       qmass2 = qmass2_tmp
       qmassf = qmassf_tmp
       zgsint = zgsint_tmp

    end if   ! end of if-masterproc

#if ( defined SPMD )
    call mpibcast (tmass0,1,mpir8,0,mpicom)
    call mpibcast (tmassf,1,mpir8,0,mpicom)
    call mpibcast (qmass1,1,mpir8,0,mpicom)
    call mpibcast (qmass2,1,mpir8,0,mpicom)
    call mpibcast (qmassf,1,mpir8,0,mpicom)
    call mpibcast (zgsint,1,mpir8,0,mpicom)
#endif

#if ( defined SPMD )
    numperlat = plon
    call compute_gsfactors (numperlat, numrecv, numsend, displs)
    ! Don't do this for empty arrays. (Satisfy NAG checks.)
    if (beglat <= endlat) then
       call mpiscatterv (ps_tmp    ,numsend, displs, mpir8,ps    (1,beglat,1) ,numrecv, mpir8,0,mpicom)
    end if
#else
!$omp parallel do private(lat)
    do lat = 1,plat
       ps(:nlon(lat),lat,1) = ps_tmp(:nlon(lat),lat)
    end do
#endif
    return

  end subroutine global_int

  subroutine copytimelevels()
   use pmgrid,       only: plon, plev, plevp, beglat, endlat
   use prognostics,  only: ps, u3, v3, t3, q3, qm1, div, ptimelevels

!---------------------------Local variables-----------------------------
!
   integer n                    ! index
!
! Make all time levels of prognostics contain identical data.
! Fields to be convectively adjusted only *require* n3 time
! level since copy gets done in linems.
!
   do n=2,ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(:,:,:,:,n) = q3(:,:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
   end do
   qm1(:,:,:,:) = q3(:,:,:,:,1)

  end subroutine copytimelevels

  subroutine read_inidat( ncid_ini, ncid_topo, dyn_in, datatype)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Read initial dataset and spectrally truncate as appropriate.
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
    use dyn_comp, only : dyn_import_t
    use ppgrid,       only: pverp
    use phys_grid,        only: scatter_field_to_chunk, gather_chunk_to_field
    
    use constituents,     only: pcnst, cnst_name, cnst_read_iv, cnst_get_ind
!
! Arguments
!
   type(file_desc_t), intent(inout) :: ncid_ini, ncid_topo
   type(dyn_import_t)            :: dyn_in   ! not used in eul dycore
   integer, optional, intent(in) :: datatype !  in:  model or analysis dataset

!
!---------------------------Local workspace-----------------------------
!
    integer i,c,m,n                         ! indices
    integer ncol

    real(r8), pointer, dimension(:,:,:)   :: convptr_2d
    real(r8), pointer, dimension(:,:,:,:) :: convptr_3d
    real(r8), pointer, dimension(:,:,:,:) :: cldptr
    real(r8), pointer, dimension(:,:    ) :: arr2d_tmp
    real(r8), pointer, dimension(:,:    ) :: arr2d
    character*16 fieldname                  ! field name

    character*16 :: subname='READ_INIDAT'   ! subroutine name
!
!-----------------------------------------------------------------------
!     May 2004 revision described below (Olson)
!-----------------------------------------------------------------------
!
! This routine reads and processes fields one at a time to minimize 
! memory usage.
!
!   State fields (including PHIS) are read into a global array on 
!     masterproc, processed, and scattered to all processors on the
!     appropriate grid 
!
!   Physics fields are read in and scattered immediately to all
!     processors on the physics grid.
!
!-----------------------------------------------------------------------
!

!
!-------------------------------------
! Allocate memory for temporary arrays
!-------------------------------------
!
! Note if not masterproc still might need to allocate array for spmd case
! since each processor calls MPI_scatter 
!
    allocate ( ps_tmp  (plon,plat     ) )
    allocate ( phis_tmp(plon,plat     ) )
    allocate ( q3_tmp  (plon,plev,plat) )
    allocate ( t3_tmp  (plon,plev,plat) )
!
!---------------------
! Read required fields
!---------------------

!
!-----------
! 3-D fields
!-----------
!
    allocate ( arr3d_a (plon,plev,plat) )
    allocate ( arr3d_b (plon,plev,plat) )

    fieldname = 'U'
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                     arr3d_a, readvar, grid_map='global')
    if(.not. readvar) call endrun()

    fieldname = 'V'
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                     arr3d_b, readvar, grid_map='global')
    if(.not. readvar) call endrun()

    call process_inidat_pio('UV')
    

    fieldname = 'T'
    call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                      t3_tmp, readvar, grid_map='global')
    if(.not. readvar) call endrun()
    call process_inidat_pio('T')

    ! Constituents (read and process one at a time)

    do m = 1,pcnst

       readvar   = .false.
       fieldname = cnst_name(m)
       if(cnst_read_iv(m)) &
            call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, plon, 1, plev, 1, plat, &
                                                           arr3d_a, readvar, grid_map='global')
       call process_inidat_pio('CONSTS', m_cnst=m, ncid=ncid_ini)

    end do

    deallocate ( arr3d_a  )
    deallocate ( arr3d_b  )
!
!-----------
! 2-D fields
!-----------
!
    
    fieldname = 'PHIS'
    readvar   = .false.
    if (ideal_phys .or. aqua_planet) then
       phis_tmp(:,:) = 0._r8
    else
       call infld(fieldname, ncid_topo, 'lon', 'lat', 1, plon, 1, plat, &
            phis_tmp, readvar, grid_map='global')
       if(.not. readvar) call endrun()
    end if
    
    call process_inidat_pio('PHIS', ncid=ncid_topo)

    fieldname = 'PS'
    call infld(fieldname, ncid_ini, 'lon', 'lat', 1, plon, 1, plat, &
         ps_tmp  , readvar, grid_map='global')
    if(.not. readvar) call endrun()

    call process_inidat_pio('PS')

!
! Integrals of mass, moisture and geopotential height
! (fix mass of moisture as well)
!
    call global_int

    deallocate ( ps_tmp   )
    deallocate ( phis_tmp )
!
! Initialization of other misc. required fields (this could be moved elsewhere)
!
    ed1(:,:,:) = 0._r8
!

    deallocate ( q3_tmp  )
    deallocate ( t3_tmp  )

    call copytimelevels()

    return

  end subroutine read_inidat

!*********************************************************************

  subroutine process_inidat_pio(fieldname, m_cnst, ncid)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Post-process input fields
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
!
    use pspect,       only: psp
    use commap
    use comspe
    use spetru,       only: spetru_phis, spetru_ps, spetru_3d_scalar, spetru_uv
    use physconst   , only: rair
    use constituents, only: cnst_name, qmin
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use carma_intr,   only: carma_implements_cnst, carma_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers , only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use stratiform,   only: stratiform_implements_cnst, stratiform_init_cnst
    use microp_driver, only: microp_driver_implements_cnst, microp_driver_init_cnst
    use phys_control,  only: phys_getopts
    use co2_cycle   , only: co2_implements_cnst, co2_init_cnst
    use unicon_cam,   only: unicon_implements_cnst, unicon_init_cnst
    use cam_control_mod, only: pertlim
#if ( defined SPMD )
    use spmd_dyn, only: compute_gsfactors
#endif

!
! Input arguments
!
    character(len=*), intent(in)           :: fieldname ! fields to be processed
    integer         , intent(in), optional :: m_cnst    ! constituent index
    type(file_desc_t)         , intent(inout), optional :: ncid
!
!---------------------------Local workspace-----------------------------
!
    integer i,ii,j,k,n,lat,irow            ! grid and constituent indices
    real(r8) pertval                       ! perturbation value
    real(r8) tmp1                          ! tmp space
    real(r8) phi(2,psp/2)                  ! used in spectral truncation of phis
!                                          ! using "B" part of hybrid grid only
    integer  varid                         ! netCDF variable id
    integer  ret, attlen                   ! netcdf return values
    logical  phis_hires                    ! true => PHIS came from hi res topo
    character*256 text
    character*256 trunits                  ! tracer untis

    real(r8), pointer, dimension(:,:,:) :: q_tmp
    real(r8), pointer, dimension(:,:,:) :: tmp3d_a, tmp3d_b, tmp3d_c, tmp3d_extend
    real(r8), pointer, dimension(:,:  ) :: tmp2d_a, tmp2d_b

#if ( defined SPMD )
    integer :: numperlat                   ! number of values per latitude band
    integer :: numsend(0:npes-1)           ! number of items to be sent
    integer :: numrecv                     ! number of items to be received
    integer :: displs(0:npes-1)            ! displacement array
#endif
    integer, allocatable :: gcid(:)
    character*16 :: subname='PROCESS_INIDAT' ! subroutine name

    logical :: nonempty_arrays             ! Test to see if arrays are empty.
                                           ! (used to satisfy NAG checks)

    nonempty_arrays = (beglat <= endlat)

    select case (fieldname)

!------------
! Process U/V
!------------

    case ('UV')

       allocate ( tmp3d_a(plon,plev,plat) )
!
! Spectral truncation
!
       call spetru_uv(arr3d_a ,arr3d_b ,div=tmp3d_a)

#if ( defined SPMD )
       numperlat = plnlv
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       if (nonempty_arrays) then
          call mpiscatterv (arr3d_a ,numsend, displs, mpir8,u3 (1,1,beglat,1)  ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (arr3d_b ,numsend, displs, mpir8,v3 (1,1,beglat,1)  ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (tmp3d_a ,numsend, displs, mpir8,div(1,1,beglat,1)  ,numrecv, mpir8,0,mpicom)
       end if
#else
       u3    (1:plon,:,1:plat,1) = arr3d_a(:plon,:plev,:plat)
       v3    (1:plon,:,1:plat,1) = arr3d_b(:plon,:plev,:plat)
       div   (:plon,:,:,1) = tmp3d_a(:plon,:,:)
#endif

       deallocate ( tmp3d_a )

!----------
! Process T
!----------

    case ('T')

       allocate ( tmp3d_a(plon,plev,plat) )
       allocate ( tmp3d_b(plon,plev,plat) )
!
! Add random perturbation to temperature if required
!
       if(masterproc) then	
          if (pertlim.ne.0.0_r8) then
             write(iulog,*)trim(subname), ':  Adding random perturbation bounded by +/-', &
                  pertlim,' to initial temperature field'
             do lat = 1,plat
                do k = 1,plev
                   do i = 1,nlon(lat)
                      call random_number (pertval)
                      pertval = 2._r8*pertlim*(0.5_r8 - pertval)
                      t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1._r8 + pertval)
                   end do
                end do
             end do
          end if
       end if
!
! Spectral truncation
!
       call spetru_3d_scalar(t3_tmp ,dl=tmp3d_a ,dm=tmp3d_b)

#if ( defined SPMD )
       numperlat = plnlv
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       if (nonempty_arrays) then
          call mpiscatterv (t3_tmp  ,numsend, displs, mpir8,t3(1,1,beglat,1) ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (tmp3d_a ,numsend, displs, mpir8,tl( 1,1,beglat) ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (tmp3d_b ,numsend, displs, mpir8,tm( 1,1,beglat) ,numrecv, mpir8,0,mpicom)
       end if
#else
       t3    (:plon,:,:,1) = t3_tmp(:plon,:,:)
       tl    (:,:,:) = tmp3d_a(:plon,:,:)
       tm    (:,:,:) = tmp3d_b (:plon,:,:)
#endif

       deallocate ( tmp3d_a )
       deallocate ( tmp3d_b )

!---------------------
! Process Constituents
!---------------------

    case ('CONSTS')

       if (.not. present(m_cnst)) then
          call endrun('  '//trim(subname)//' Error:  m_cnst needs to be present in the'// &
                      ' argument list')
       end if

       allocate ( tmp3d_extend(plon,plev,beglat:endlat) )
!
! If "Q", then allocate extra space for spectral truncation
!
       if(m_cnst == 1) then
          allocate ( tmp3d_a(plon,plev,plat) )
          allocate ( tmp3d_b(plon,plev,plat) )
          allocate ( tmp3d_c(plon,plev,plat) )
       end if
!
! Check that all tracer units are in mass mixing ratios
!
       if(readvar) then
          ret = pio_inq_varid    (NCID,cnst_name(m_cnst), varid)
          ret = pio_get_att (NCID,varid,'units',trunits)
          if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
             call endrun('  '//trim(subname)//' Error:  Units for tracer ' &
                  //trim(cnst_name(m_cnst))//' must be in KG/KG')
          end if
!
! Constituents not read from initial file are initialized by the package that implements them.
!
       else
          if(m_cnst == 1 .and. moist_physics) then
             call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
          end if

        arr3d_a(:,:,:) = 0._r8
        allocate(gcid(plon))
        do j=1,plat
           gcid(:) = j
           if(masterproc) write(iulog,*) 'Warning:  Not reading ',cnst_name(m_cnst), ' from IC file.'
           if (microp_driver_implements_cnst(cnst_name(m_cnst))) then
              call microp_driver_init_cnst(cnst_name(m_cnst),arr3d_a(:,:,j) , gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "microp_driver_init_cnst"'
           else if (stratiform_implements_cnst(cnst_name(m_cnst))) then
              call stratiform_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "stratiform_init_cnst"'
           else if (chem_implements_cnst(cnst_name(m_cnst))) then
              call chem_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "chem_init_cnst"'
           else if (tracers_implements_cnst(cnst_name(m_cnst))) then
              call tracers_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "tracers_init_cnst"'
           else if (aoa_tracers_implements_cnst(cnst_name(m_cnst))) then
              call aoa_tracers_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "aoa_tracers_init_cnst"'
           else if (carma_implements_cnst(cnst_name(m_cnst))) then
              call carma_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "carma_init_cnst"'
           else if (co2_implements_cnst(cnst_name(m_cnst))) then
              call co2_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "co2_init_cnst"'
           else if (unicon_implements_cnst(cnst_name(m_cnst))) then
              call unicon_init_cnst(cnst_name(m_cnst), arr3d_a(:,:,j), gcid)
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' initialized by "unicon_init_cnst"'
           else
              if(masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' set to 0.'
           end if
        end do
        
        deallocate ( gcid )
       end if

!$omp parallel do private(lat)
       do lat = 1,plat
          call qneg3(trim(subname), lat   ,nlon(lat),plon   ,plev    , &
               m_cnst, m_cnst, qmin(m_cnst) ,arr3d_a(1,1,lat))
       end do
!
! if "Q", "CLDLIQ", or "CLDICE", save off for later use
!
       if(m_cnst == 1       ) q3_tmp (:plon,:,:) = arr3d_a(:plon,:,:)
!
! Spectral truncation of "Q" (only to get spectral derivatives)
!
       if(m_cnst == 1) then
          tmp3d_a(:plon,:,:) = arr3d_a(:plon,:,:)
          call spetru_3d_scalar(tmp3d_a ,dl=tmp3d_b, dm=tmp3d_c)
       end if


#if ( defined SPMD )
       numperlat = plnlv
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       if (nonempty_arrays) then
          call mpiscatterv (arr3d_a  , numsend, displs, mpir8, tmp3d_extend(1,1,beglat) ,numrecv, mpir8,0,mpicom)
          q3(:plon,:,m_cnst,beglat:endlat,1) = tmp3d_extend(:plon,:,beglat:endlat)
          if(m_cnst == 1) then
             call mpiscatterv (tmp3d_b, numsend, displs, mpir8, ql( 1,1  ,beglat  ) ,numrecv, mpir8,0,mpicom)
             call mpiscatterv (tmp3d_c, numsend, displs, mpir8, qm( 1,1  ,beglat  ) ,numrecv, mpir8,0,mpicom)
          end if
       end if
#else
       q3(:plon,:plev,m_cnst,:plat,1) = arr3d_a(:plon,:plev,:plat)
       if(m_cnst == 1) then
          ql(:,:,:) = tmp3d_b(:plon,:,:)
          qm(:,:,:) = tmp3d_c(:plon,:,:)
       end if
#endif
       deallocate ( tmp3d_extend )
       if(m_cnst == 1) then
          deallocate ( tmp3d_a )
          deallocate ( tmp3d_b )
          deallocate ( tmp3d_c )
       end if

!-----------
! Process PS
!-----------

    case ('PS')

       allocate ( tmp2d_a(plon,plat) )
       allocate ( tmp2d_b(plon,plat) )
!
! Spectral truncation
!
       call spetru_ps  (ps_tmp, tmp2d_a, tmp2d_b)

#if ( defined SPMD )
       numperlat = plon
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       if (nonempty_arrays) then
          call mpiscatterv (tmp2d_a ,numsend, displs, mpir8,dpsl(1,beglat) ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (tmp2d_b ,numsend, displs, mpir8,dpsm(1,beglat) ,numrecv, mpir8,0,mpicom)
       end if
#else
       dpsl(:,:) = tmp2d_a(:plon,:)
       dpsm(:,:) = tmp2d_b(:plon,:)
#endif

       deallocate ( tmp2d_a )
       deallocate ( tmp2d_b )

!-------------
! Process PHIS
!-------------

    case ('PHIS')

       allocate ( tmp2d_a(plon,plat) )
       allocate ( tmp2d_b(plon,plat) )

!
! Check for presence of 'from_hires' attribute to decide whether to filter
!
       if(readvar) then
          ret = pio_inq_varid (ncid, 'PHIS', varid)
          call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)

          ret = pio_inq_attlen (ncid, varid, 'from_hires', attlen)
          if (ret.eq.PIO_NOERR .and. attlen.gt.256) then
             call endrun('  '//trim(subname)//' Error:  from_hires attribute length is too long')
          end if
          ret = pio_get_att (ncid, varid, 'from_hires', text)
          if (ret.eq.PIO_NOERR .and. text(1:4).eq.'true') then
             phis_hires = .true.
             if(masterproc) write(iulog,*) trim(subname), ': Will filter input PHIS: attribute from_hires is true'
          else
             phis_hires = .false.
             if(masterproc) write(iulog,*)trim(subname), ': Will not filter input PHIS: attribute ', &
                  'from_hires is either false or not present'
          end if
          call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

       else
          phis_hires = .false.
       end if
!
! Spectral truncation
!
       call spetru_phis  (phis_tmp, phis_hires, phisl=tmp2d_a, phism=tmp2d_b, phi_out=phi)
!
! Compute ln(Ps*) (Ritchie & Tanguay, 1995) in spectral space
!
       tmp1 = 1._r8/(rair*t0(plev))
       do ii = 1,psp/2
          i = 2*ii - 1
          lnpstar(i  ) = -phi(1,ii)*tmp1
          lnpstar(i+1) = -phi(2,ii)*tmp1
       end do
          
#if ( defined SPMD )
       call mpibcast   (lnpstar   ,psp,mpir8,0 , mpicom)
       numperlat = plon
       call compute_gsfactors (numperlat, numrecv, numsend, displs)
       if (nonempty_arrays) then
          call mpiscatterv (phis_tmp  ,numsend, displs, mpir8,phis  (1,beglat) ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (tmp2d_a   ,numsend, displs, mpir8,phisl (1,beglat) ,numrecv, mpir8,0,mpicom)
          call mpiscatterv (tmp2d_b   ,numsend, displs, mpir8,phism (1,beglat) ,numrecv, mpir8,0,mpicom)
       end if
#else
!$omp parallel do private(lat)
       do lat = 1,plat
          phis (:nlon(lat),lat) = phis_tmp(:nlon(lat),lat)
          phisl(:nlon(lat),lat) = tmp2d_a (:nlon(lat),lat)
          phism(:nlon(lat),lat) = tmp2d_b (:nlon(lat),lat)
       end do
#endif

       deallocate ( tmp2d_a )
       deallocate ( tmp2d_b )

    end select

    return
    
  end subroutine process_inidat_pio

end module inidat
