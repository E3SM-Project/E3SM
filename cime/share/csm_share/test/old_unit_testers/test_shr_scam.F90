program test_shr_scam

  use shr_kind_mod, only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_scam_mod
  use shr_mpi_mod
  use shr_sys_mod
  use shr_ncread_mod
  use test_mod
  use netcdf
  use pio
  implicit none
#include <mpif.h>

  real(r8) :: targetLat, targetLon    ! target latitude/longitude
  real(r8) :: closeLat,  closeLon     ! close latitude/longitude
  real(r8) :: expect(2)               ! lat lon of expected
  integer :: closeLatIdx, closeLonIdx ! indices of returned points
  integer :: rc                       ! return code
  integer :: ncid                     ! NetCDF id
  integer :: npes, mype               ! number of processors and my processor rank
  character(len=CL) :: filename       ! Filename to read
  character(len=CL) :: badfilename    ! bad Filename to read
  character(len=CL) :: csmdata        ! directory to inputdata
  type(file_desc_t) :: pioid          ! pio file ID
  type (iosystem_desc_t), pointer :: piosystems
  logical :: found                    ! if found or NOT

  call test_init( 22 )

  ! Test simple valid tests
  csmdata = "/fs/cgd/csm/inputdata"
  filename = trim(csmdata)//"/lnd/clm2/surfdata/surfdata_1.9x2.5_simyr2000_c100505.nc"
  write(6,*) "Test file: "//trim(filename)
  targetLat = 45.0
  targetLon = 180.0
  expect = (/ 44.5263157894736d00, targetLon /)
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found )
  write(6,*) "closest values to target of : ", targetLat, targetLon, " is: ", &
              closeLat, closeLon
  call test_is( found, "Test that a a simple call with filename works" )
  call test_close( expect, (/ closeLat, closeLon /), 1.e-13_r8, "Test lat/lon found correct" )
  expect = (/ closeLat, closeLon /)
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx )
  call test_is( expect, (/ closeLat, closeLon /), "Test OK without found" )
  rc = nf90_open( filename, NF90_NOWRITE, ncid )
  if ( rc /= NF90_NOERR ) call shr_sys_abort( "NetCDF error opening file: "//trim(filename) )
  call shr_scam_getCloseLatLon( ncid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found )
  
  call test_is( found, "Test that a a simple call to NetCDF id works" )
  call test_is( expect, (/ closeLat, closeLon /), "Test lat/lon found correct" )
  call shr_scam_getCloseLatLon( ncid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx )
  call test_is( expect, (/ closeLat, closeLon /), "Test OK without found" )
  
  if ( nf90_close( ncid ) /= NF90_NOERR ) call shr_sys_abort( "NetCDF error closing file" )
  write(6,*) "init mpi"
  call shr_mpi_init( )
  call shr_mpi_commsize( MPI_COMM_WORLD, npes )
  call shr_mpi_commrank( MPI_COMM_WORLD, mype )
  write(6,*) "init PIO"
  allocate( piosystems )
  call PIO_init(mype, MPI_COMM_WORLD, npes, 1, 1, pio_rearr_box, piosystems, base=0)

  rc  = pio_openfile(piosystems, pioid, iotype_netcdf, filename, pio_nowrite)
  if(rc/= PIO_NOERR) call shr_sys_abort( "PIO error opening file: "//trim(filename) )
  write(6,*) "PIO open on file"
  call shr_scam_getCloseLatLon( pioid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found )

  call test_is( found, "Test that a a simple call to the PIO interface works" )
  call test_is( expect, (/ closeLat, closeLon /), "Test lat/lon found correct" )
  call shr_scam_getCloseLatLon( pioid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx )
  call test_is( expect, (/ closeLat, closeLon /), "Test OK without found" )
  call pio_closefile(pioid)

  ! Test that can find periodic longitudes 
  targetLat = 1.0
  targetLon = 842.0
  expect = (/ 0.947368421052549d00, 122.5d00 /)
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( found, "Test that periodic longitude targets returns" )
  write(6,*) "closest values to target of : ", targetLat, targetLon, " is: ", &
              closeLat, closeLon
  call test_close( expect, (/ closeLat, closeLon /), 1.e-13_r8, "Test lat/lon found correct" )
  expect = (/ closeLat, closeLon /)
  filename = trim(csmdata)// &
  "/lnd/clm2/initdata/clmi.BCN.2000-01-01_1.9x2.5_gx1v6_simyr2000_c100309.nc"
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( found, "Test that can find targets for clmi file" )
  call test_close( expect, (/ closeLat, closeLon /), 1.d-13, &
                  "Test that clmi targets same as other file" )
  ! Test abort tests
  ! non-existant filename
  call shr_ncread_setAbort( .false. )
  badfilename = "ZZTop.nc"
  call shr_scam_getCloseLatLon( badfilename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that non existant file returns NOT found" )
  call shr_scam_getCloseLatLon( ncid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that non existant NetCDF ID returns NOT found" )
  call shr_scam_getCloseLatLon( pioid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that non existant PIO ID returns NOT found" )
  ! Test that targets outside of global lat/lons return not found
  targetLat = -91.0
  targetLon = 0.0
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that bad negative lat returns NOT found" )
  if ( found ) then
     write(6,*) "closest values to target of : ", targetLat, targetLon, " is: ", &
                 closeLat, closeLon
  end if
  targetLat = +91.0
  targetLon = 0.0
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that bad positive lat returns NOT found" )
  if ( found ) then
     write(6,*) "closest values to target of : ", targetLat, targetLon, " is: ", &
                 closeLat, closeLon
  end if
  targetLat = 45.
  targetLon = 180.
  filename = trim(csmdata)// &
  "/lnd/clm2/snicardata/snicar_optics_5bnd_c090915.nc"
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that can NOT find targets for snicar optics file" )
  filename = trim(csmdata)// &
  "/lnd/clm2/pftdata/pft-physiology.c110425.nc"
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that can NOT find targets for pft-phys file" )
  filename = trim(csmdata)// &
  "/lnd/clm2/mappingdata/maps/10x15/map_0.5x0.5_landuse_to_10x15_aave_da_110307.nc"
  call shr_scam_getCloseLatLon( filename, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call test_is( .not. found, "Test that can NOT find targets for mapping file" )
  rc  = pio_openfile(piosystems, pioid, iotype_netcdf, filename, pio_nowrite)
  if(rc/= PIO_NOERR) call shr_sys_abort( "PIO error opening file: "//trim(filename) )
  call shr_scam_getCloseLatLon( pioid, targetLat,  targetLon, closeLat, closeLon, &
                                closeLatIdx, closeLonIdx, found=found, rc=rc )
  call pio_closefile(pioid)
  call test_is( .not. found, "Test that can NOT find targets for PIO clmi file" )

  call test_final()

end
