!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_read_nc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.1 $
!     created:   $Date: 2009/05/22 21:02:13 $
!

!=============================================================================== 
! rrtmg_lw_read_nc.f90
!
! Description: This program reads all of the RRTM longwave data from a NetCDF
!              file in band by band subroutines, as a replacement for the 
!              rrtmg_lw_k_g.f90 data statements.
!
! Written By: Patrick Hofmann
! Last Update: 1/23/2009
!===============================================================================

!*******************************************************************************
subroutine lw_kgb01
	use rrlw_kg01, only : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
						selfrefo, forrefo, no1
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
	
	integer(kind=im) :: ab
	integer(kind=im), parameter :: bandNumber = 1, numGPoints = no1
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr

	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for N2
	call getAbsorberIndex('N2',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mn2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, kbo_mn2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mn2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mn2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)

end subroutine lw_kgb01
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb02
        use rrlw_kg02, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no2
        use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
		
   	integer(kind=im), parameter :: bandNumber = 2
	integer(kind=im), parameter :: numGPoints = no2
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb02
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb03
        use rrlw_kg03, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                          kbo_mn2o, selfrefo, forrefo, no3
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
	
	integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 3
	integer(kind=im), parameter :: numGPoints = no3
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keyupper,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for N2
	call getAbsorberIndex('N2O',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))
	
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, kbo_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,T,numGPoints,1,1,1/))

	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, keyupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, keyupper*tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mn2o, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mn2o, keyupper*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb03
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb04
        use rrlw_kg04, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no4
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
		
    integer(kind=im), parameter :: bandNumber = 4
	integer(kind=im), parameter :: numGPoints = no4
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo(:,1:5), &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keyupper,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, keyupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, keyupper*tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb04
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb05
        use rrlw_kg05, only : fracrefao, fracrefbo, kao, kbo, kao_mo3, &
                          selfrefo, forrefo, ccl4o, no5
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    
    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 5
	integer(kind=im), parameter :: numGPoints = no5
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keyupper,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, & 
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keyupper,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for O3
	call getAbsorberIndex('O3',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mo3, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))
	
	!Get absorber index for CCL4
	call getAbsorberIndex('CCL4',ab)
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, ccl4o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, keyupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, keyupper*tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mo3, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ccl4o, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb05
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb06
        use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
                          cfc11adjo, cfc12o, no6
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
	
	integer(kind=im) :: ab
        integer(kind=im), parameter :: bandNumber = 6
	integer(kind=im), parameter :: numGPoints = no6
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(9)   = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

	!Get absorber index for CO2
	call getAbsorberIndex('CO2',ab)
	status(12)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(13)  = nf90_get_var(ncid, varID, kao_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	!Get absorber index for CFC11
	call getAbsorberIndex('CFC11',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, cfc11adjo, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
	
	!Get absorber index for CFC12
	call getAbsorberIndex('CFC12',ab)
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, cfc12o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mco2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( cfc11adjo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( cfc12o, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb06
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb07	
        use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
                          kbo_mco2, selfrefo, forrefo, no7
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
    	
    integer(kind=im) :: ab
   	integer(kind=im), parameter :: bandNumber = 7
	integer(kind=im), parameter :: numGPoints = no7
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, & 
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for CO2
	call getAbsorberIndex('CO2',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))
	
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, kbo_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mco2, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mco2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb07
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb08		
        use rrlw_kg08, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                          kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                          cfc12o, cfc22adjo, no8
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    	
    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 8
	integer(kind=im), parameter :: numGPoints = no8
    integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
        status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for O3
	call getAbsorberIndex('O3',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mo3, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	!Get absorber index for CO2
	call getAbsorberIndex('CO2',ab)
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, kao_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	status(18)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(19)  = nf90_get_var(ncid, varID, kbo_mco2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	!Get absorber index for N2O
	call getAbsorberIndex('N2O',ab)
	status(20)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(21)  = nf90_get_var(ncid, varID, kao_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	status(22)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(23)  = nf90_get_var(ncid, varID, kbo_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	!Get absorber index for CFC12
	call getAbsorberIndex('CFC12',ab)
	status(24)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(25)  = nf90_get_var(ncid, varID, cfc12o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
	
	!Get absorber index for CFC22
	call getAbsorberIndex('CFC22',ab)
	status(26)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(27)  = nf90_get_var(ncid, varID, cfc22adjo, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,1,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mo3, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mco2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mco2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mn2o, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mn2o, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( cfc12o, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( cfc22adjo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb08
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb09	
        use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
                            kbo_mn2o, selfrefo, forrefo, no9
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    	
    integer(kind=im) :: ab
    integer(kind=im), parameter :: bandNumber = 9
	integer(kind=im), parameter :: numGPoints = no9
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,keylower,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))

	!Get absorber index for N2O
	call getAbsorberIndex('N2O',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/keylower,T,numGPoints,1,1,1/))
	
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, kbo_mn2o, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))

	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mn2o, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mn2o, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb09
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb10		
        use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no10
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    	
   	integer(kind=im), parameter :: bandNumber = 10
	integer(kind=im), parameter :: numGPoints = no10
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb10
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb11
        use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
                          kbo_mo2, selfrefo, forrefo, no11
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
	
	integer(kind=im) :: ab
	integer(kind=im), parameter :: bandNumber = 11
	integer(kind=im), parameter :: numGPoints = no11
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)   = nf90_NoErr
  if (masterproc) then	
	status(1)   = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)   = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)   = nf90_get_var(ncid, varID, fracrefao, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(4)   = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)   = nf90_get_var(ncid, varID, fracrefbo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/numGPoints,1,1,1/))
	
	status(6)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)   = nf90_get_var(ncid, varID, kao, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,plower,numGPoints,1,1/))
	
	status(8)   = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)   = nf90_get_var(ncid, varID, kbo, &
                      start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                      count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11)  = nf90_get_var(ncid, varID, selfrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tself,numGPoints,1,1/))
	
	status(12)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13)  = nf90_get_var(ncid, varID, forrefo, &
                      start = (/1,1,bandNumber,gPointSetNumber/), &
                      count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for O2
	call getAbsorberIndex('O2',ab)
	status(14)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15)  = nf90_get_var(ncid, varID, kao_mo2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	status(16)  = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(17)  = nf90_get_var(ncid, varID, kbo_mo2, &
                      start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                      count = (/1,T,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mo2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mo2, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)

end subroutine lw_kgb11
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb12
        use rrlw_kg12, only : fracrefao, kao, selfrefo, forrefo, no12
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
		
	integer(kind=im), parameter :: bandNumber = 12
	integer(kind=im), parameter :: numGPoints = no12
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:) = nf90_NoErr
  if (masterproc) then	
	status(1) = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2) = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3) = nf90_get_var(ncid, varID, fracrefao, &
                    start = (/1,1,bandNumber,gPointSetNumber/), &
                    count = (/numGPoints,keylower,1,1/))
	
	status(4) = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(5) = nf90_get_var(ncid, varID, kao, &
                    start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                    count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(6) = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(7) = nf90_get_var(ncid, varID, selfrefo, &
                    start = (/1,1,bandNumber,gPointSetNumber/), &
                    count = (/Tself,numGPoints,1,1/))
	
	status(8) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(9) = nf90_get_var(ncid, varID, forrefo, &
                    start = (/1,1,bandNumber,gPointSetNumber/), &
                    count = (/Tforeign,numGPoints,1,1/))

	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb12
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb13		
        use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                          kbo_mo3, selfrefo, forrefo, no13
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"
    
    integer(kind=im) :: ab
	integer(kind=im), parameter :: bandNumber = 13
	integer(kind=im), parameter :: numGPoints = no13  
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, fracrefbo, &
                     start = (/1,1,bandNumber,gPointSetNumber/),  &
                     count = (/numGPoints,1,1,1/))
	
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(9)  = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for O3
	call getAbsorberIndex('O3',ab)
	status(12) = nf90_inq_varid(ncid,"AbsorptionCoefficientsUpperAtmos",varID)
	status(13) = nf90_get_var(ncid, varID, kbo_mo3, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/1,T,numGPoints,1,1,1/))
	
	!Get absorber index for CO2
	call getAbsorberIndex('CO2',ab)
	status(14) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(15) = nf90_get_var(ncid, varID, kao_mco2, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/keylower,T,numGPoints,1,1,1/))
	
	!Get absorber index for CO
	call getAbsorberIndex('CO',ab)
	status(16) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
        status(17) = nf90_get_var(ncid, varID, kao_mco, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/keylower,T,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo_mo3, t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mco2, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mco, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb13
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb14	
        use rrlw_kg14, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no14
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 
	
    include "mpif.h"

	integer(kind=im), parameter :: bandNumber = 14
	integer(kind=im), parameter :: numGPoints = no14
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, fracrefbo, &
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
	
	status(10) = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))

	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb14
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb15	
        use rrlw_kg15, only : fracrefao, kao, kao_mn2, selfrefo, forrefo, no15
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
    
    integer(kind=im) :: ab
 	integer(kind=im), parameter :: bandNumber = 15
	integer(kind=im), parameter :: numGPoints = no15
	integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
	status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(6)  = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(7)  = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(9)  = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))
	
	!Get absorber index for N2
	call getAbsorberIndex('N2',ab)
	status(10) = nf90_inq_varid(ncid,"AbsorptionCoefficientsLowerAtmos",varID)
	status(11) = nf90_get_var(ncid, varID, kao_mn2, &
                     start = (/1,1,1,ab,bandNumber,gPointSetNumber/), &
                     count = (/keylower,T,numGPoints,1,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao_mn2, keylower*t*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
	
end subroutine lw_kgb15
!*******************************************************************************

!*******************************************************************************
subroutine lw_kgb16		
        use rrlw_kg16, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no16
	use rrlw_ncpar
	use rrtm_grid, only : masterproc
	use netcdf
	
	implicit none
	save 

    include "mpif.h"
	  
    integer(kind=im), parameter :: bandNumber = 16
 	integer(kind=im), parameter :: numGPoints = no16
    integer(kind=im), parameter :: gPointSetNumber = 1
	integer(kind=im) :: ncid, varID, ierr
	
        status(:)  = nf90_NoErr
  if (masterproc) then	
	status(1)  = nf90_open('RUNDATA/rrtmg_lw.nc',nf90_nowrite,ncid)
	
	status(2)  = nf90_inq_varid(ncid,"PlanckFractionLowerAtmos",varID)
	status(3)  = nf90_get_var(ncid, varID, fracrefao, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,keylower,1,1/))
	
	status(4)  = nf90_inq_varid(ncid,"PlanckFractionUpperAtmos",varID)
	status(5)  = nf90_get_var(ncid, varID, fracrefbo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/numGPoints,1,1,1/))
	
	status(6)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsLowerAtmos",varID)
	status(7)  = nf90_get_var(ncid, varID, kao, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/keylower,Tdiff,plower,numGPoints,1,1/))
	
	status(8)  = nf90_inq_varid(ncid,"KeySpeciesAbsorptionCoefficientsUpperAtmos",varID)
	status(9)  = nf90_get_var(ncid, varID, kbo, &
                     start = (/1,1,1,1,bandNumber,gPointSetNumber/), &
                     count = (/1,Tdiff,pupper,numGPoints,1,1/))
	
	status(10) = nf90_inq_varid(ncid,"H20SelfAbsorptionCoefficients",varID)
	status(11) = nf90_get_var(ncid, varID, selfrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tself,numGPoints,1,1/))
	
	status(12) = nf90_inq_varid(ncid,"H20ForeignAbsorptionCoefficients",varID)
	status(13) = nf90_get_var(ncid, varID, forrefo, &
                     start = (/1,1,bandNumber,gPointSetNumber/), &
                     count = (/Tforeign,numGPoints,1,1/))
	
	if(any(status(:) /= nf90_NoErr)) stop  "Error reading variables from file" 
	
	status(1) = nf90_close(ncid)
  endif	

    CALL MPI_BCAST( fracrefao, keylower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( fracrefbo, numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kao, keylower*tdiff*plower*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( kbo, tdiff*pupper*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( selfrefo, tself*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( forrefo, tforeign*numGPoints, MPI_DOUBLE_PRECISiON, 0, MPI_COMM_WORLD, ierr)

end subroutine lw_kgb16
!*******************************************************************************
