!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_read_nc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.1 $
!     created:   $Date: 2009/05/22 22:20:11 $

!=============================================================================== 
! rrtmg_read_nc.f90
!
! Description: This program reads all of the RRTM shortwave data from a NetCDF
!              file in band by band subroutines, as a replacement for the 
!              rrtmg_sw_k_g.f90 data statements.
!
! Written By: Patrick Hofmann
! Last Update: 4/3/2009
!===============================================================================

!*******************************************************************************
subroutine sw_kgb16
	use rrsw_kg16, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no16
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"

	integer(kind=im), parameter :: bandNumber = 1, numGPoints = no16
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
	 
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
													 
	status(14) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 16 variables from file" 
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeignlower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
	rayl = ncrayl(1)
	
end subroutine sw_kgb16
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb17
        use rrsw_kg17, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no17
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"

   	integer(kind=im), parameter :: bandNumber = 2
	integer(kind=im), parameter :: numGPoints = no17
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionUpperAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keyupper,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo(1:3,:), start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
													 
	status(14) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsUpperAtmos",varID)
	status(15) = nf90_get_var(ncid, varID, forrefo(4,:), start = (/2,1,bandNumber,gPointSetNumber/), &
                     count = (/1,numGPoints,1,1/))	
	
	status(14) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 17 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, keyupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, keyupper*tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, 4*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
	rayl = ncrayl(1)
	
end subroutine sw_kgb17
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb18
	use rrsw_kg18, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no18
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"

    integer(kind=im), parameter :: bandNumber = 3
	integer(kind=im), parameter :: numGPoints = no18
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
	
	status(14) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 18 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeignlower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
	rayl = ncrayl(1)
	
end subroutine sw_kgb18
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb19
        use rrsw_kg19, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no19
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
		
    integer(kind=im), parameter :: bandNumber = 4
	integer(kind=im), parameter :: numGPoints = no19
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
	
	status(14) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 19 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeignlower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)

	rayl = ncrayl(1)
	
end subroutine sw_kgb19
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb20
        use rrsw_kg20, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, absch4o, no20
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
    
    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 5
	integer(kind=im), parameter :: numGPoints = no20
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo(1:3,:), start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
													 
	status(14) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsUpperAtmos",varID)
	status(15) = nf90_get_var(ncid, varID, forrefo(4,:), start = (/2,1,bandNumber,gPointSetNumber/), &
                     count = (/1,numGPoints,1,1/))												 
													 
	!Get absorber index for CH4
	call getAbsorberIndex('CH4',ab)
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, absch4o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
																				 
	status(16) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 20 variables from file"																		 
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, (tforeignlower+1)*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( absch4o, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
	rayl = ncrayl(1)
	
end subroutine sw_kgb20
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb21
        use rrsw_kg21, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no21
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
	
    integer(kind=im), parameter :: bandNumber = 6
	integer(kind=im), parameter :: numGPoints = no21
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo(1:3,:), start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
													 
	status(14) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsUpperAtmos",varID)
	status(15) = nf90_get_var(ncid, varID, forrefo(4,:), start = (/2,1,bandNumber,gPointSetNumber/), &
                     count = (/1,numGPoints,1,1/))			
													 
	status(14) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 21 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, keyupper*tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, (tforeignlower+1)*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)

	rayl = ncrayl(1)
	
end subroutine sw_kgb21
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb22	
        use rrsw_kg22, only: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no22
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    	
   	integer(kind=im), parameter :: bandNumber = 7
	integer(kind=im), parameter :: numGPoints = no22
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
	
	status(14) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 22 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeignlower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
	rayl = ncrayl(1)
	
end subroutine sw_kgb22
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb23		
        use rrsw_kg23, only: sfluxrefo, kao, selfrefo, forrefo, raylo, no23
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"

    integer(kind=im), parameter :: bandNumber = 8
	integer(kind=im), parameter :: numGPoints = no23
        integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
        status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, raylo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	 
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(9)  = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(11) = nf90_get_var(ncid, varID, forrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
	
	status(12) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 23 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( raylo, numgpoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeignlower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine sw_kgb23
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb24	
        use rrsw_kg24, only: sfluxrefo, kao, kbo, selfrefo, forrefo, &
						 raylao, raylbo, abso3ao, abso3bo, no24
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    	
    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 9
	integer(kind=im), parameter :: numGPoints = no24
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, raylao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	 
	status(6)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsUpperAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, raylbo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
		
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(11) = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(13) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(14) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(15) = nf90_get_var(ncid, varID, forrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
	
	!Get absorber index for O3
	call getAbsorberIndex('O3',ab)
	status(16) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(17) = nf90_get_var(ncid, varID, abso3ao, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/1,1,numGPoints,1,1,1/))
	status(18) = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(19) = nf90_get_var(ncid, varID, abso3bo, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/1,1,numGPoints,1,1,1/))
											
	status(20) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 24 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( raylao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( raylbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeignlower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( abso3ao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( abso3bo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine sw_kgb24
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb25		
        use rrsw_kg25, only: sfluxrefo, kao, raylo, abso3ao, abso3bo, no25
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
    	
    integer(kind=im) :: ab	
   	integer(kind=im), parameter :: bandNumber = 10
	integer(kind=im), parameter :: numGPoints = no25
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, raylo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	 
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	!Get absorber index for O3
	call getAbsorberIndex('O3',ab)
	status(8)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, abso3ao, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/1,1,numGPoints,1,1,1/))
	status(10) = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(11) = nf90_get_var(ncid, varID, abso3bo, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/1,1,numGPoints,1,1,1/))
											
	status(12) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 25 variables from file"
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( raylo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( abso3ao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( abso3bo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine sw_kgb25
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb26
        use rrsw_kg26, only: sfluxrefo, raylo, no26
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
	
	integer(kind=im), parameter :: bandNumber = 11
	integer(kind=im), parameter :: numGPoints = no26
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, raylo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	 
	status(6)  = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 26 variables from file"										
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( raylo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine sw_kgb26
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb27
        use rrsw_kg27, only: sfluxrefo, kao, kbo, raylo, no27
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
		
	integer(kind=im), parameter :: bandNumber = 12
	integer(kind=im), parameter :: numGPoints = no27
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, raylo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 27 variables from file"	
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( raylo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine sw_kgb27
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb28		
        use rrsw_kg28, only: sfluxrefo, kao, kbo, rayl, no28
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"

	integer(kind=im), parameter :: bandNumber = 13
	integer(kind=im), parameter :: numGPoints = no28  
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionUpperAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keyupper,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))
				
	status(10) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 28 variables from file"										
  endif	

    CALL MPI_BCAST( sfluxrefo, keyupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, keyupper*tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
	rayl = ncrayl(1)
	
end subroutine sw_kgb28
!*******************************************************************************

!*******************************************************************************
subroutine sw_kgb29	
        use rrsw_kg29, only: sfluxrefo, kao, kbo, selfrefo, forrefo, &
						 absh2oo, absco2o, rayl, no29
	use rrsw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
	
	integer(kind=im) :: ab
	integer(kind=im), parameter :: bandNumber = 14
	integer(kind=im), parameter :: numGPoints = no29
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	real(kind=rb) :: ncrayl(1)
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_sw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"SolarSourceFunctionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, sfluxrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"RayleighExtinctionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, ncrayl, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,1,1,1/))
		
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H2OSelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsLowerAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo(1:3,:), start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeignlower,numGPoints,1,1/))
													 
	status(14) = nf90_inq_varid(ncid,"H2OForeignAbsorptionCoefficientsUpperAtmos",varID)
	status(15) = nf90_get_var(ncid, varID, forrefo(4,:), start = (/2,1,bandNumber,gPointSetNumber/), &
                     count = (/1,numGPoints,1,1/))	
	
	!Get absorber index for H2O
	call getAbsorberIndex('H2O',ab)
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, absh2oo, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
							
	!Get absorber index for CO2
	call getAbsorberIndex('CO2',ab)
	status(18)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(19)  = nf90_get_var(ncid, varID, absco2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
	
	status(18) = nf90_close(ncid)
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading band 29 variables from file"	
  endif	

    CALL MPI_BCAST( sfluxrefo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ncrayl, 1, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, (tforeignlower+1)*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( absh2oo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( absco2o, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)

	rayl = ncrayl(1)
	
end subroutine sw_kgb29
!*******************************************************************************
