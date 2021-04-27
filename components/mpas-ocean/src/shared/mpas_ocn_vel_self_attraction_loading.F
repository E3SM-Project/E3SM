! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_self_attraction_loading
!
!> \brief MPAS ocean module
!> \author Kristin Barton 
!> \date   August 2020
!> \details
!>  This module contains routines for the tidal potential forcing.
!>  Design document located in :
!>    MPAS-Model/docs/ocean/design_docs
!
!-----------------------------------------------------------------------

module ocn_vel_self_attraction_loading

    use mpas_kind_types
    use mpas_constants
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_timekeeping
    use mpas_timer
    use mpas_threading, only : mpas_threading_get_thread_num
    use ocn_constants
    use ocn_diagnostics
    use ocn_diagnostics_variables
    use ocn_config
    use ocn_mesh
    use netcdf
#ifdef _MPI
    use mpi
#endif
    use iso_c_binding
    use iso_fortran_env, only: real64
   
    implicit none
    private
    save
    
    include 'shtns.f03'
    
    !--------------------------------------------------------------------
    !
    ! Public parameters
    !
    !--------------------------------------------------------------------
    
    !--------------------------------------------------------------------
    !
    ! Public member functions
    !
    !--------------------------------------------------------------------
    
    public :: ocn_compute_self_attraction_loading, &
              ocn_vel_self_attraction_loading_init
    
    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------
    
    logical :: selfAttractionLoadingOn
    
    ! Interpolation weights variables
    integer, dimension(:), allocatable :: toRowValues, toColValues
    integer, dimension(:), allocatable :: fromRowValues, fromColValues
    real, dimension(:), allocatable :: toSValues, fromSValues
    integer:: nMpas, nGrid
    
    ! SHTns routine variables
    integer, parameter :: dp=real64
    type(c_ptr) :: shtns_c
    type(shtns_info), pointer :: shtns
    real, dimension(:), pointer :: scaling
    complex(dp), allocatable :: Slm(:)
    
    ! MPI variables
    integer :: nCellsGlobal
    integer, dimension(:), allocatable :: nCellsDisplacement, indexToCellIDGathered
    integer, dimension(:), allocatable :: nCellsPerProc

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_compute_self_attraction_loading
!
!> \brief   Computes self-attraction and loading 
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine computes the sea surface height perturbation due to
!>     self-attraction and loading. 
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_self_attraction_loading(domain, forcingPool, dminfo, ssh, err)!{{{

        !-----------------------------------------------------------------
        !
        ! input variables
        !
        !-----------------------------------------------------------------
        type (dm_info), intent(in) :: dminfo
        real (kind=RKIND), dimension(:), intent(in) :: ssh
        
        !-----------------------------------------------------------------
        !
        ! input/output variables
        !
        !-----------------------------------------------------------------
        type (domain_type), intent(inout) :: domain !< Input/Output: domain information
        type (mpas_pool_type), intent(inout) :: forcingPool 
        
        !-----------------------------------------------------------------
        !
        ! output variables
        !
        !-----------------------------------------------------------------
        
        integer, intent(out) :: err !< Output: Error flag
        
        !-----------------------------------------------------------------
        !
        ! local variables
        !
        !-----------------------------------------------------------------
        
        ! SAL variables
        real (kind=RKIND), dimension(:), pointer :: ssh_gg
        ! MPI Variables
        integer :: iCell, ilm, curProc
        real (kind=RKIND), dimension(:), allocatable :: globalArray, gatheredArray

        call mpas_timer_start('SAL Calculation')

        err = 0
        if (.not. selfAttractionLoadingOn) then
            call mpas_timer_stop('SAL Calculation')
            return
        endif

        ! Allocate globalArray and gatheredArray only on process 0
        call MPI_COMM_RANK( dminfo % comm, curProc, err)
        if (curProc.eq.0) then
            allocate(globalArray(nCellsGlobal), gatheredArray(nCellsGlobal))
        endif

        ! Gather only the nCellsOwned from ssh (does not include Halos)
        call MPI_GATHERV(ssh, nCellsOwned, MPI_DOUBLE, gatheredArray, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, dminfo % comm, err)

        ! Perform SAL calculation only on process 0
        if (curProc.eq.0) then
            allocate(ssh_gg(nGrid))
            ssh_gg(:) = 0.0_dp
            Slm(:)=(0.0_dp, 0.0_dp)

            ! Rearrange ssh into CellID order
            do iCell = 1,nCellsGlobal
               globalArray(indexToCellIDGathered(iCell)) = gatheredArray(iCell)
            enddo
           
            ! Interpolate ssh onto Gaussian Grid
            call interpolate( toColValues, toRowValues, toSValues, globalArray, ssh_gg)
            ! Perform spherical harmonic transform
            call Spat_to_SH(shtns_c, ssh_gg, Slm)
            ! Multiply each harmonic coefficient by the scaling factor
            do ilm = 1,shtns%nlm
                Slm(ilm) = Slm(ilm) * scaling(ilm)
            enddo
            ! Perform inverse spherical harmonic transform
            call SH_to_spat(shtns_c, Slm, ssh_gg)
            ! Interpolate back to MPAS mesh
            call interpolate( fromColValues, fromRowValues, fromSValues, ssh_gg, globalArray)

            ! Rearrange back to index order
            do iCell = 1,nCellsGlobal
                gatheredArray(iCell) = globalArray(indexToCellIDGathered(iCell))
            enddo
        endif

        ! Scatter back to ssh_sal
!!mrp check nCellsAll on next line
        call MPI_SCATTERV(gatheredArray, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                    ssh_sal, nCellsAll, MPI_DOUBLE, 0, dminfo % comm, err) 

        if (curProc.eq.0) then
            deallocate(globalArray, gatheredArray, ssh_gg)
        endif

        ! Perform Halo exchange update
        call mpas_dmpar_field_halo_exch(domain,'ssh_sal')

        call mpas_timer_stop('SAL Calculation')

    end subroutine ocn_compute_self_attraction_loading!}}}


!***********************************************************************
!
!  routine ocn_vel_self_attraction_loading_init
!
!> \brief   Initializes ocean tidal protential forcing module.
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine initializes the ocean self-attraction and loading module 
!>  and variables.
!
!-----------------------------------------------------------------------

    subroutine ocn_vel_self_attraction_loading_init(domain,err)!{{{

        type (domain_type), intent(inout) :: domain
        integer, intent(out) :: err !< Output: error flag
       
        ! Config variables
        type (block_type), pointer :: block_ptr
        type (mpas_pool_type), pointer :: forcingPool, meshPool
       
        ! NetCDF and weights file variables
        integer :: toNcId, toNsDimId, toRowId, toColId, toSId
        integer :: fromNcId, fromNsDimId, fromRowId, fromColId, fromSId
        integer:: nMpasDimId, nGridDimId, toNsLen, fromNsLen 
        character (len = NF90_MAX_NAME) :: toNsName, fromNsName, nMpasName, nGridName
        integer, pointer :: n_s
        character(len=StrKIND) :: mpasToGridFile, gridToMpasFile

        ! SHTns routine variables
        character(len=StrKIND) :: loadLoveFile
        integer :: lmax, mmax, mres, layout, norm
        integer, pointer :: nphi, nlat
        real(dp),parameter :: pi=acos(-1.0_dp)
        real(dp), pointer :: cosTheta(:), sinTheta(:)
        real(dp) :: eps_polar
        real(dp), allocatable :: Sh(:,:), Sh_init(:,:)
        integer, allocatable :: lvec(:)
       
        ! MPI variables
        integer :: curProc
        integer, pointer ::  nCells
        integer, dimension(:), pointer :: indexToCellID
        integer :: iProc, l, ilm, nProcs

        ! Alarm variables
        type (MPAS_Time_Type) :: alarmTime
        type (MPAS_TimeInterval_type) :: alarmTimeStep
        
        call mpas_timer_start('SAL Init')

        err = 0
        
        selfAttractionLoadingOn= .false.
        if (config_use_self_attraction_loading) then
            selfAttractionLoadingOn = .true.
        else
            call mpas_timer_stop('SAL Init')
            return
        endif
   
        call mpas_pool_get_config(ocnConfigs, 'config_load_Love_file', config_load_Love_file)
        call mpas_pool_get_config(ocnConfigs, 'config_mpas_to_grid_weights_file', config_mpas_to_grid_weights_file)
        call mpas_pool_get_config(ocnConfigs, 'config_grid_to_mpas_weights_file', config_grid_to_mpas_weights_file)
        call mpas_pool_get_config(ocnConfigs, 'config_self_attraction_loading_compute_interval', config_self_attraction_loading_compute_interval)
        call mpas_pool_get_config(ocnConfigs, 'config_nLatitude', nlat)
        call mpas_pool_get_config(ocnConfigs, 'config_nLongitude', nphi)

        block_ptr => domain % blocklist
        do while ( associated( block_ptr ) )
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_dimension(meshPool, 'nCells', nCells) 
            block_ptr => block_ptr % next
        end do
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

        ! Setup Alarm for SAL
        alarmTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr=err)
        call mpas_set_timeInterval(alarmTimeStep, timeString = &
              config_self_attraction_loading_compute_interval, ierr = err)
        call mpas_add_clock_alarm(domain % clock, 'salComputeAlarm', alarmTime, &
              alarmTimeInterval=alarmTimeStep, ierr=err)

        ! Begin MPI portion
        call MPI_COMM_RANK( domain % dminfo % comm, curProc, err)
        call MPI_COMM_SIZE( domain % dminfo % comm, nProcs, err)
        
        if (curProc.eq.0) then
            allocate(nCellsPerProc(nProcs))
            allocate(nCellsDisplacement(nProcs))
        endif
       
        ! Gather nCellsOwned
        call MPI_GATHER( nCellsOwned, 1, MPI_INTEGER, nCellsPerProc, 1, MPI_INTEGER, &
                         0, domain % dminfo % comm, err)
       
        ! Set Displacement variable for GATHERV command
        if (curProc.eq.0) then
            nCellsGlobal = sum(nCellsPerProc)
            allocate(indexToCellIDGathered(nCellsGlobal))
            nCellsDisplacement(1) = 0
            if (nProcs > 1) then
                do iProc=2,nProcs
                    nCellsDisplacement(iProc) = nCellsDisplacement(iProc-1) + nCellsPerProc(iProc-1)
                enddo
            endif
        endif

        ! Gather indexToCellID
        call MPI_GATHERV( indexToCellID, nCellsOwned, MPI_INTEGER, indexToCellIDGathered, &
                nCellsPerProc, nCellsDisplacement, MPI_INTEGER, 0, domain % dminfo % comm, err)
     
        ! SAL calculations only need to be initialized on process 0
        if (curProc.eq.0) then
            ! Set variables for file locations
            loadLoveFile = trim(config_load_Love_file)
            mpasToGridFile = trim(config_mpas_to_grid_weights_file)
            gridToMpasFile = trim(config_grid_to_mpas_weights_file)

            ! Open netcdf weights files
            call check( nf90_open(path = mpasToGridFile, mode = nf90_nowrite, ncid = toNcId) ,err)
            call check( nf90_open(path = gridToMpasFile, mode = nf90_nowrite, ncid = fromNcId) ,err)
            
            ! Get dimension ID
            call check( nf90_inq_dimid(toNcId, "n_s", toNsDimId) ,err)
            call check( nf90_inq_dimid(fromNcId, "n_s", fromNsDimId) ,err)
            call check( nf90_inq_dimid(toNcId, "n_a", nMpasDimId) ,err)
            call check( nf90_inq_dimid(fromNcId, "n_a", nGridDimId) ,err)
            
            ! Get Variable IDs
            call check( nf90_inq_varid(toNcId, "row", toRowId) ,err)
            call check( nf90_inq_varid(toNcId, "col", toColId) ,err)
            call check( nf90_inq_varid(toNcId, "S", toSId) ,err)
            call check( nf90_inq_varid(fromNcId, "row", fromRowId) ,err)
            call check( nf90_inq_varid(fromNcId, "col", fromColId) ,err)
            call check( nf90_inq_varid(fromNcId, "S", fromSId) ,err)
            
            ! Get Dimension Length
            call check( nf90_inquire_dimension(toNcId, toNsDimId, toNsName, toNsLen) ,err)
            call check( nf90_inquire_dimension(fromNcId, fromNsDimId, fromNsName, fromNsLen) ,err)
            call check( nf90_inquire_dimension(toNcId, nMpasDimId, nMpasName, nMpas) ,err)
            call check( nf90_inquire_dimension(fromNcId, nGridDimId, nGridName, nGrid) ,err)
            
            ! Allocate matrices to read data into
            allocate ( toRowValues (toNsLen) )
            allocate ( toColValues (toNsLen) )
            allocate ( toSValues (toNsLen) )
            allocate ( fromRowValues (fromNsLen) )
            allocate ( fromColValues (fromNsLen) )
            allocate ( fromSValues (fromNsLen) )
            
            ! Retrieve data
            call check( nf90_get_var(toNcId, toColId, toColValues(:) ) ,err)
            call check( nf90_get_var(toNcId, toRowId, toRowValues(:) ) ,err)
            call check( nf90_get_var(toNcId, toSId, toSValues(:) ) ,err)
            call check( nf90_get_var(fromNcId, fromColId, fromColValues(:) ) ,err)
            call check( nf90_get_var(fromNcId, fromRowId, fromRowValues(:) ) ,err)
            call check( nf90_get_var(fromNcId, fromSId, fromSValues(:) ) ,err)
            
            ! Initialize SHTns routine
            lmax = nlat-1
            mmax = lmax
            mres = 1
            eps_polar = 1.0e-10_dp
            norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
            layout = SHT_GAUSS + SHT_PHI_CONTIGUOUS
            call shtns_verbose(0)
            shtns_c = shtns_create(lmax,mmax,mres,norm)
            call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)
            call c_f_pointer(cptr=shtns_c, fptr=shtns)
            call c_f_pointer(cptr=shtns%ct, fptr=cosTheta, shape=[shtns%nlat])
            call c_f_pointer(cptr=shtns%st, fptr=sinTheta, shape=[shtns%nlat])
            allocate( Slm(1:shtns%nlm), scaling(1:shtns%nlm), lvec(1:shtns%nlm) )
            
            ! List l values in lvec (for use in reading load Love numbers)
            lvec = 0
            do ilm = 1,shtns%nlm
                l = shtns_lm2l(shtns_c, ilm)
                lvec(ilm) = l
            enddo
            
            call getloadLoveNums(shtns%nlm, lvec, loadLoveFile)
        endif

        call mpas_timer_stop('SAL Init')

   end subroutine ocn_vel_self_attraction_loading_init!}}}

!***********************************************************************
!
!  routine check
!
!> \brief   Check status of netcdf operations
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine checks to status of the netcdf file
!
!-----------------------------------------------------------------------

    subroutine check(status,err) !{{{
        integer, intent ( in) :: status
        integer, intent(inout) :: err
        
        if(status /= nf90_noerr) then
            err = 1
        end if
    end subroutine!}}}

!***********************************************************************
!
!  routine interpolate
!
!> \brief   Perform interpolation
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine contains the sparse matrix multiplication
!>    algorithm to interpolate between MPAS and Gaussian Grid
!
!-----------------------------------------------------------------------

    subroutine interpolate(colValues, rowValues, sValues, dataIn, dataOut) !{{{

        integer, dimension(:) :: rowValues, colValues
        real (kind=RKIND), dimension(:) :: sValues, dataOut, dataIn
        real (kind=RKIND) :: rhs = 0
        integer :: n_S, n, nRow, nCol
        
        n_S = size(sValues)
        n = 1
        
        do while (n .LE. n_S)
            nRow = rowValues(n)
            do while ( rowValues(n) .EQ. nRow )
                nCol = colValues(n)
                rhs = rhs + dataIn(nCol) * sValues(n)
                n = n + 1
            end do
            dataOut(nRow) = rhs
            rhs = 0
        end do

    end subroutine!}}}

!***********************************************************************
!
!  routine getloadLoveNums
!
!> \brief   Get the load Love numbers
!> \author  Kristin Barton
!> \date    February 2021
!> \details
!>  Modified by Kristin Barton from routine Written by K. Quinn (March 
!>     2010) and modifed by M. Schindelegger (May 2017) 
!> 	This program contains load Love numbers Wang et al. (2012):
!>  https://www.sciencedirect.com/science/article/pii/S0098300412002245
!>      Note:  these load love number are in the CM reference frame,
!>      i.e. center of mass of total Earth system (solid+hydro+atmos).
!>      To convert to other reference frames, see Blewitt, JGR, 2003,
!>      https://doi.org/10.1029/2002JB002082 
!>      For CE (center of mass solid Earth), ke0 = 0, he1 = hm1 + 1,
!>      le1 = lm1 + 1, ke1 = km1 + 1.  For CF (center of figure ~ CE, within 2%),
!>      kf0 = 0, hf1 = 2/3(hm1 - lm1), lf1 = -1/3(hm1 - lm1),
!>      kf1 = -1/3*hm1 -2/3*lm1 - 1.
!>      All love numbers for degree >= 2 are the same for all ref frames.
!
!-----------------------------------------------------------------------

    subroutine getloadLoveNums(nlm, lvec, filename) !{{{

        integer, dimension(:) :: lvec
        real :: n(1441), H(1441), L(1441), K(1441)
        real :: H1, L1, K1
        integer :: i, j, nlm
        real(dp), parameter :: rhoE=5517d0 ! Average density of Earth (kg/m^3)
        real(dp), parameter :: rhoW=1035.0d0 ! Density of water (kg/m^3)
        character(len=StrKIND) :: filename 
        
        ! Format: LovDat = (n, H, L, K)
        open(action = 'read', file=filename, unit=11)
        
        do i=1,1441
            read (11, *) n(i), H(i), L(i), K(i) 
        enddo
        close(11)

        ! Convert from center of mass (CM) to center of frame (CF)
        H1 = H(2)
        L1 = L(2)
        K1 = K(2)
        H(2) = 2.0d0/3.0d0*(H1 - L1)
        L(2) = -1.0d0/3.0d0*(H1 - L1)
        K(2) = -1.0d0/3.0d0*H1 - 2.0d0/3.0d0*L1 - 1.0d0

        ! Get scaling factors from load Love numbers 
        do i=1,nlm
            j = lvec(i)
            scaling(i) = (1.0d0 + K(j+1) - H(j+1)) / ((2*j)+1)
        enddo
        scaling = scaling * 3.0d0 * rhoW / rhoE

    end subroutine!}}}

end module ocn_vel_self_attraction_loading!}}}


!***********************************************************************

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
