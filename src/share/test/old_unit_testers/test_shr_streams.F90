module streams_exp
   use shr_kind_mod,    only : SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_CX
   use shr_sys_mod,     only : shr_sys_abort
   use shr_file_mod,    only : shr_file_getUnit, shr_file_freeUnit
   use shr_stream_mod

   implicit none

   private

   public streams_exp_init
   public streams_exp_set
   public streams_exp_write_strm_txt
   public is_streams_expected

   public streams_exp_data

   integer, public,  parameter :: maxFiles = 2000

   type streams_exp_data
      character(SHR_KIND_CL) :: dataSource
      character(SHR_KIND_CL) :: filePath
      character(SHR_KIND_CX) :: fldListFile
      character(SHR_KIND_CX) :: fldListModel
      character(SHR_KIND_CL) :: domFilePath
      character(SHR_KIND_CL) :: domFileName
      character(SHR_KIND_CL) :: domTvarName
      character(SHR_KIND_CL) :: domXvarName
      character(SHR_KIND_CL) :: domYvarName
      character(SHR_KIND_CL) :: domAreaName
      character(SHR_KIND_CL) :: domMaskName
      integer                :: nfiles
      character(SHR_KIND_CL) :: filenames(maxFiles)
   end type streams_exp_data

contains

subroutine streams_exp_init( streams_exp )
 implicit none
 type(streams_exp_data), intent(OUT) :: streams_exp

 integer :: i

 streams_exp%dataSource   = "dataSource"
 streams_exp%filePath     = "filePath/"
 streams_exp%fldListFile  = "T:U"
 streams_exp%fldListModel = "Temp:Wind_u"
 streams_exp%domFilePath  = "domFilePath/"
 streams_exp%domFileName  = "domFileName"
 streams_exp%domTvarName  = "time"
 streams_exp%domXvarName  = "xc"
 streams_exp%domYvarName  = "yc"
 streams_exp%domAreaName  = "area"
 streams_exp%domMaskName  = "mask"
 streams_exp%nfiles       = 1
 do i = 1, streams_exp%nfiles
   write(streams_exp%filenames(i), '(a,i2.2)') "filename", i
 end do
end subroutine streams_exp_init

subroutine streams_exp_set( streams_exp, datasource, filePath, fldListfile,      &
                            fldListModel, domFilePath, domFileName, domTvarName, &
                            domXvarName, domYvarName, domAreaName, domMaskName,  &
                            nfiles, filenames )
 implicit none
 type(streams_exp_data), intent(INOUT) :: streams_exp
 character(*), intent(IN), optional :: dataSource
 character(*), intent(IN), optional :: filePath
 character(*), intent(IN), optional :: fldListFile
 character(*), intent(IN), optional :: fldListModel
 character(*), intent(IN), optional :: domFilePath
 character(*), intent(IN), optional :: domFileName
 character(*), intent(IN), optional :: domTvarName
 character(*), intent(IN), optional :: domXvarName
 character(*), intent(IN), optional :: domYvarName
 character(*), intent(IN), optional :: domAreaName
 character(*), intent(IN), optional :: domMaskName
 integer     , intent(IN), optional :: nfiles
 character(*), intent(IN), optional :: filenames(:)

 integer :: i

 if ( present(dataSource)   ) streams_exp%dataSource   = datasource
 if ( present(filePath)     ) streams_exp%filePath     = filePath
 if ( present(fldListFile)  ) streams_exp%fldListFile  = fldListFile
 if ( present(fldListModel) ) streams_exp%fldListModel = fldListModel
 if ( present(domFilePath)  ) streams_exp%domFilePath  = domFilePath
 if ( present(domFileName)  ) streams_exp%domFileName  = domFileName
 if ( present(domTvarName)  ) streams_exp%domTvarName  = domTvarName
 if ( present(domXvarName)  ) streams_exp%domXvarName  = domXvarName
 if ( present(domYvarName)  ) streams_exp%domYvarName  = domYvarName
 if ( present(domAreaName)  ) streams_exp%domAreaName  = domAreaName
 if ( present(domMaskName)  ) streams_exp%domMaskName  = domMaskName
 if ( present(nfiles) .and. present(filenames) )then
    streams_exp%nfiles       = nfiles
    do i = 1, streams_exp%nfiles
       streams_exp%filenames(i) = filenames(i)
    end do
 end if

end subroutine streams_exp_set


subroutine streams_exp_write_strm_txt( stream_filename, streams_exp )
 use shr_string_mod, only : shr_string_listGetNum, shr_string_listGetName
 use shr_sys_mod,    only : shr_sys_system
 implicit none
 character(SHR_KIND_CL), intent(IN) :: stream_filename
 type(streams_exp_data), intent(IN) :: streams_exp

 integer :: unit, n, rcode, nfModel, nfFile
 character(SHR_KIND_CS) :: varModel, varFile
 character(*), parameter :: sub = "write_streams_txt"

   unit = shr_file_getUnit( )
   write(*,*) "Write streams text file out to: ", trim(stream_filename)
   open( unit, file=stream_filename, status="unknown")

   write(unit,*) "<dataSource>"
   write(unit,*) "   ", trim(streams_exp%dataSource)
   write(unit,*) "</dataSource>"
   write(unit,*) "<domainInfo>"
   write(unit,*) "  <variableNames>"
   write(unit,*) "     ", trim(streams_exp%domTvarName), "  time"
   write(unit,*) "     ", trim(streams_exp%domXvarName), "  lon"
   write(unit,*) "     ", trim(streams_exp%domYvarName), "  lat"
   write(unit,*) "     ", trim(streams_exp%domAreaName), "  area"
   write(unit,*) "     ", trim(streams_exp%domMaskName), "  mask"
   write(unit,*) "  </variableNames>"
   write(unit,*) "  <filePath>"
   write(unit,*) "     ", trim(streams_exp%domFilePath)
   write(unit,*) "  </filePath>"
   write(unit,*) "  <fileNames>"
   write(unit,*) "     ", trim(streams_exp%domFileName)
   write(unit,*) "  </fileNames>"
   write(unit,*) "</domainInfo>"
   write(unit,*) "<fieldInfo>"
   write(unit,*) "  <variableNames>"
   nfModel = shr_string_listGetNum( streams_exp%fldListModel )
   nfFile  = shr_string_listGetNum( streams_exp%fldListFile  )
   do n = 1, max( nfModel, nfFile )
       if ( n > nfFile ) then
          varFile = " "
       else
          call shr_string_listGetName(streams_exp%fldListFile,  n, varFile  )
       end if
       if ( n > nfModel ) then
          varModel = " "
       else
          call shr_string_listGetName(streams_exp%fldListModel, n, varModel )
       end if
       write(unit,*) &
                 "     ", trim(varFile), "  ", &
                 "     ", trim(varModel)
   end do
   write(unit,*) "  </variableNames>"
   write(unit,*) "  <filePath>"
   write(unit,'(A,A)') "     ", trim(streams_exp%FilePath)
   write(unit,*) "  </filePath>"
   write(unit,*) "  <fileNames>"
   do n = 1, streams_exp%nfiles
       write(unit,*) &
                 "     ", trim(streams_exp%filenames(n))
   end do
   write(unit,*) "  </fileNames>"
   write(unit,*) "</fieldInfo>"
   close(unit)
   call shr_file_freeUnit(unit)
   call shr_sys_system( "cat "//trim(stream_filename), rcode )

end subroutine streams_exp_write_strm_txt

logical function is_streams_expected( stream, streams_exp )
   implicit none
   type(shr_stream_streamType)  ,intent(in)  :: stream  ! stream in question
   type(streams_exp_data), intent(IN) :: streams_exp

      character(SHR_KIND_CL) :: dataSource
      character(SHR_KIND_CL) :: filePath
      character(SHR_KIND_CX) :: fldListFile
      character(SHR_KIND_CX) :: fldListModel
      character(SHR_KIND_CL) :: domFilePath
      character(SHR_KIND_CL) :: domFileName
      character(SHR_KIND_CL) :: domTvarName
      character(SHR_KIND_CL) :: domXvarName
      character(SHR_KIND_CL) :: domYvarName
      character(SHR_KIND_CL) :: domAreaName
      character(SHR_KIND_CL) :: domMaskName
      character(SHR_KIND_CL) :: filen, file_next, file_first
      integer :: n

   is_streams_expected = .true.

   call shr_stream_getFileFieldList(  stream, fldlistFile  )
   call shr_stream_getModelFieldList( stream, fldlistModel )
   call shr_stream_getFilePath(       stream, filePath     )
   call shr_stream_getDataSource(     stream, dataSource   )
   call shr_stream_getDomainInfo(     stream, domFilePath, domfileName,      &
                                      domTvarName, domXvarName, domYvarName, &
                                      dommaskName, domareaName)
   if ( trim(fldListFile) /= trim(streams_exp%fldListFile) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "fldListFile different"
   if ( .not. is_streams_expected )then
      write(*,*) trim(fldListFile)
      write(*,*) trim(streams_exp%fldListFile)
   end if
   if ( trim(fldListModel) /= trim(streams_exp%fldListModel) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "fldListModel different"
   if ( trim(filePath) /= trim(streams_exp%filePath) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "filePath different"
   if ( trim(dataSource) /= trim(streams_exp%dataSource) ) &
          is_streams_expected = .false.
   if ( trim(domFilePath) /= trim(streams_exp%domFilePath) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domfilePath different"
   if ( trim(domFileName) /= trim(streams_exp%domFileName) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domfileName different"
   if ( trim(domTvarName) /= trim(streams_exp%domTvarName) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domTvarName different"
   if ( trim(domXvarName) /= trim(streams_exp%domXvarName) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domXvarName different"
   if ( trim(domYvarName) /= trim(streams_exp%domYvarName) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domYvarName different"
   if ( trim(domAreaName) /= trim(streams_exp%domAreaName) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domAreaName different"
   if ( trim(domMaskName) /= trim(streams_exp%domMaskName) ) &
          is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "domMaskName different"
   n = 1
   call shr_stream_getFirstFileName( stream, filen )
   file_first = filen
   if ( trim(filen) /= trim(streams_exp%filenames(1)) ) is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "first file different"
   do while( n < streams_exp%nfiles )
       n = n + 1
       call shr_stream_getNextFileName( stream, filen, file_next )
       if ( trim(file_next) /= trim(streams_exp%filenames(n)) ) &
                                                   is_streams_expected = .false.
       if ( .not. is_streams_expected ) write(*,*) "next file different"
       if ( trim(file_next) == trim(file_first)  ) is_streams_expected = .false.
       if ( .not. is_streams_expected ) write(*,*) "Too few files"
       filen = file_next
   end do
   call shr_stream_getNextFileName( stream, filen, file_next )
   if ( trim(file_next) /= trim(file_first) ) is_streams_expected = .false.
   if ( .not. is_streams_expected ) write(*,*) "too many files"

end function is_streams_expected

end module streams_exp

program test_shr_streams

  use shr_kind_mod
  use shr_string_mod
  use shr_sys_mod
  use shr_stream_mod
  use streams_exp
  use test_mod

  implicit none

  type(shr_stream_streamType), pointer  :: streams(:)  ! stream in question
  type(shr_stream_streamType), pointer  :: streams2(:)  ! stream in question
  integer :: yearFirst, yearLast, yearAlign
  character(SHR_KIND_CL) :: stream_filename = "sfile.txt"
  character(SHR_KIND_CL) ::   rest_filename = "sfile_rest.nc"
  character(SHR_KIND_CL) :: test_descrip, filenames1(maxFiles)
  type(streams_exp_data) :: stream_exp  ! stream in question
  integer :: series, n, i
  integer, pointer :: expected(:), value(:)
  character(SHR_KIND_CS) :: clmncep(12) = (/ &
                                            "clmforc.Qian.c2006.T62.Solr.2003-01.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-02.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-03.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-04.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-05.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-06.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-07.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-08.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-09.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-10.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-11.nc", &
                                            "clmforc.Qian.c2006.T62.Solr.2003-12.nc"  &
                                          /)
  character(SHR_KIND_CS) :: clmncepTPQW(12) = (/ &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-01.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-02.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-03.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-04.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-05.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-06.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-07.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-08.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-09.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-10.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-11.nc", &
                                            "clmforc.Qian.c2006.T62.TPQW.2003-12.nc"  &
                                          /)
   character(SHR_KIND_CS) :: filenames2(12)
   integer :: mDateIn, SecIn, year, month, rcode, exp_int, nfiles
   integer :: mDateLB, dDateLB, secLB, n_lb
   integer :: mDateUB, dDateUB, secUB, n_ub
   character(SHR_KIND_CL) :: fileLB, fileUB
   integer :: num_series, num_fail
   integer, parameter :: bogus_TEST         = 1, &
                         CLMNCEP_TEST       = 2, &
                         CLMNCEP_ALOGO_TEST = 3, &
                         GISS_TEST          = 4, &
                         CAMHIST_TEST       = 5

#ifdef LINUX
  num_series = CLMNCEP_ALOGO_TEST
#else
  num_series = CAMHIST_TEST
#endif
  num_fail = 3 + 12
  call test_init( 2 + (num_series-1)*3 + num_fail )
  do series = 2, num_series
     yearAlign = 1
     yearFirst = 1
     yearLast  = 1
     allocate( streams(1)  )
     allocate( streams2(1) )
     write(*,*) "Initialize expected streams"
     call streams_exp_init( stream_exp )
     if (      series  == bogus_TEST )then
        test_descrip = "bogus"
     else if ( series  == CLMNCEP_TEST )then
        test_descrip = "CLMNCEP"
        call streams_exp_set( stream_exp, datasource="CLMNCEP", &
                              fldListfile ="FSDS",      &
                              fldListModel="fsds",      &
                              filepath= &
  "/fs/cgd/csm/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727/Solar6Hrly",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
                              domfilename="domain.T62.050609.nc", &
                              nfiles=12, filenames=clmncep(1:12) )
        yearAlign = 2003
        yearFirst = 2003
        yearLast  = 2003
     else if ( series  == CLMNCEP_ALOGO_TEST )then
        test_descrip = "CLMNCEP-ALOGO"
        call streams_exp_set( stream_exp, datasource="CLMNCEP", &
                              fldListfile ="TBOT:QBOT:WIND:PSRF",      &
                              fldListModel="tbot:qbot:wind:psrf",      &
                              filepath=&
  "/fs/cgd/csm/inputdata/atm/datm7/atm_forcing.datm7.Qian.T62.c080727/TmpPrsHumWnd3Hrly",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
                              domfilename="domain.T62.050609.nc", &
                              nfiles=12, filenames=clmncepTPQW(1:12) )
        yearAlign = 1
        yearFirst = 2003
        yearLast  = 2003
#ifndef LINUX
     else if ( series  == GISS_TEST )then
        test_descrip = "GISS"
        call streams_exp_set( stream_exp, datasource="GISS", &
                              fldListfile = "lwdn:swdn:swup", &
                              fldListModel= "lwdn:swdn:swup", &
                              filepath="/fs/cgd/csm/inputdata/atm/datm7/TN460/",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/TN460/", &
                              domXvarName="lon", &
                              domYvarName="lat", &
                              domfilename="tn460nyf.giss.T62.051007.nc", &
                              nfiles=1, filenames=(/ "tn460nyf.giss.T62.051007.nc" /) )
     else if ( series  == CAMHIST_TEST )then
        test_descrip = "CAMHIST"
        yearAlign = 5
        yearFirst = 5
        yearLast  = 6
        call streams_exp_set( stream_exp, datasource="CAMHIST", &
                              fldListfile = &
           "FSNS:PRECC:PRECL:PRECSC:PRECSL:PS:PSL:QBOT:SOLL:SOLLD:SOLS:SOLSD:SRFRAD:FSNS:TBOT:UBOT:VBOT:ZBOT", &
                              fldListModel= &
           "swnet:precc:precl:snowc:snowl:ps:pslv:shum:swndr:swndf:swvdr:swvdf:srfrad:swnet:tbot:u:v:z", &
                              filepath="/fs/cgd/csm/inputdata/atm/datm7/CAMHIST/",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/CAMHIST/",        &
                              domfilename="domain.T42.050516.nc", &
                              nfiles=2, filenames=(/  &
                   "eul64x128_datm6.01.cam2.h1.0005-01-01-00000.nc", &
                   "eul64x128_datm6.01.cam2.h1.0006-01-01-00000.nc"  &
                    /) )
#endif
     end if
     write(*,*) "Write streams out to file"
     call streams_exp_write_strm_txt( stream_filename, stream_exp )
     write(*,*) "Initialize shr_streams"
     call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, yearAlign )
     if ( series  > 1 )then
        write(*,*) "Get time bounds..."
        secIn = 0
        write(*,*) "mDateIn, SecIn, mDateLB,mDateUB, dDateLB,dDateUB, secLB, secUB"
        allocate( expected((yearLast-yearFirst+3)*12) )
        allocate( value((yearLast-yearFirst+3)*12) )
        n = 0
        do year = yearAlign-1, yearAlign+1+(yearLast-yearFirst)
           do month = 1, 12
              n = n + 1
              mDateIn = year * 10000 + month*100 + 1
              call shr_stream_findBounds(streams(1),mDateIn,        secIn,       &
                                        &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
                                        &   mDateUB,dDateUB,secUB,n_ub,fileUB )
              if (      year < yearFirst )then
                 expected(n) = yearLast * 10000 + month*100 + 1
              else if ( year > yearLast  )then
                 expected(n) = yearFirst * 10000 + month*100 + 1
              else
                 expected(n) = year * 10000 + month*100 + 1
              end if
              if ( series == CAMHIST_TEST ) expected(n) = expected(n) + 1
              value(n) = dDateUB
              write(6,'(8i9)') mDateIn, SecIn, mDateLB,mDateUB, dDateLB,dDateUB, &
                               secLB, secUB
           end do
        end do
        call test_is( value, expected, " test if expected values")
        deallocate( expected )
        deallocate( value )
     end if
     call shr_stream_dataDump( streams(1) )
     write(*,*) "Check if it is as expected..."
     call test_is( is_streams_expected( streams(1), stream_exp ),  &
            "test if initialization is what expected "//trim(test_descrip) )
     write(*,*) "Write restart file out"
     call shr_stream_restWrite( streams, rest_filename, caseName="clmrun", &
                                caseDesc="clmrun description" )
     write(*,*) "Read that file into a different stream"
     call shr_stream_init( streams2(1), stream_filename, yearFirst, yearLast, yearAlign )
     call shr_stream_restRead(  streams2, rest_filename )
     write(*,*) "Check if read restart is as expected..."
     call test_is( is_streams_expected( streams2(1), stream_exp ),  &
                   "test after read restart "//trim(test_descrip) )
     deallocate( streams  )
     deallocate( streams2 )
     call shr_sys_system( "/bin/rm -f "//trim(stream_filename), rcode )
     call shr_sys_system( "/bin/rm -f "//trim(rest_filename),   rcode )
  end do

  ! Fail tests
  call shr_stream_setAbort( .false. )
  call shr_string_setAbort( .false. )
  allocate( streams(1)  )
  allocate( streams2(1)  )

  write(*,*) "Try to write uninitialized stream out"
  call shr_stream_restWrite( streams, rest_filename, caseName="clmrun", &
                                caseDesc="clmrun description", rc=rcode )
  call test_is( rcode, 1, "test that writing uninitialized stream fails" )

  write(*,*) "Try to read uninitialized stream in"
  call shr_stream_restRead(  streams2, rest_filename, rc=rCode )
  call test_is( rcode, 1, "test that reading uninitialized stream fails" )

  mDateIn = 20000101
  write(*,*) "Try to find bounds on uninitialized stream"
  call shr_stream_findBounds(streams(1),mDateIn,        secIn,      &
                           &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
                           &   mDateUB,dDateUB,secUB,n_ub,fileUB, rc=rCode )
  call test_is( rcode, 1, "test that find bounds of uninitialized stream fails" )


  do series = 1, 99
     yearAlign = 1
     yearFirst = 1
     yearLast  = 1
     call streams_exp_init( stream_exp )
     if (      series == 1 )then
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, yearAlign )
        test_descrip = "Try to read restart file that does not exist"
        call shr_sys_system( "/bin/rm -f "//trim(rest_filename),   rcode )
        call shr_stream_restRead(  streams, rest_filename, rc=rCode )
        exp_int = 2
     else if ( series == 2 )then
        test_descrip = "Try to initialize streams with too many files"
        nfiles = 1001
        do i = 1, nfiles
          write(filenames1(i),'("filename",i4.4,".nc")' ) i
        end do
        call streams_exp_set( stream_exp, datasource="CAMHIST", &
                              fldListfile = &
           "FSNS:PRECC:PRECL:PRECSC:PRECSL:PS:PSL:QBOT:SOLL:SOLLD:SOLS:SOLSD:SRFRAD:FSNS:TBOT:UBOT:VBOT:ZBOT", &
                              fldListModel= &
           "swnet:precc:precl:snowc:snowl:ps:pslv:shum:swndr:swndf:swvdr:swvdf:srfrad:swnet:tbot:u:v:z", &
                              filepath="/fs/cgd/csm/inputdata/atm/datm7/CAMHIST/",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/CAMHIST/",        &
                              domfilename="domain.T42.050516.nc", &
                              nfiles=nfiles, filenames=filenames1(1:nfiles) )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 3 )then
        test_descrip = "variable name lists do not have same number of values"
        call streams_exp_set( stream_exp, datasource="CLMNCEP", &
                              fldListfile ="TBOT:QBOT:WIND:PRECTmms:FSDS:PSRF",      &
                              fldListModel="tbot:qbot:wind:prectMMS",                &
                              filepath="/fs/cgd/csm/inputdata/atm/datm7/CLMNCEP/",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
                              domfilename="domain.T62.050609.nc", &
                              nfiles=12, filenames=clmncep(1:12) )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 4 )then
        test_descrip = "Mask name set to blank"
        call streams_exp_set( stream_exp, domMaskName=" " )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 5 )then
        test_descrip = "Area name set to blank"
        call streams_exp_set( stream_exp, domAreaName=" " )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 6 )then
        test_descrip = "Yvar name set to blank"
        call streams_exp_set( stream_exp, domYVarName=" " )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 7 )then
        test_descrip = "Xvar name set to blank"
        call streams_exp_set( stream_exp, domXVarName=" " )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 8 )then
        test_descrip = "tvar name set to blank"
        call streams_exp_set( stream_exp, domTVarName=" " )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 9 )then
        test_descrip = "no filenames"
        call streams_exp_set( stream_exp, nfiles=0, filenames=(/" "/) )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 10 )then
        test_descrip = "no fieldnames"
        call streams_exp_set( stream_exp, fldListfile ="", fldListModel="" )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        exp_int = 1
     else if ( series == 11 )then
        test_descrip = "Dates are out of range"
        call streams_exp_set( stream_exp, datasource="CLMNCEP", &
                              fldListfile ="TBOT:QBOT:WIND:PRECTmms:FSDS:PSRF",      &
                              fldListModel="tbot:qbot:wind:prectMMS:fsds:psrf",      &
                              filepath="/fs/cgd/csm/inputdata/atm/datm7/CLMNCEP/",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
                              domfilename="domain.T62.050609.nc", &
                              nfiles=12, filenames=clmncep(1:12) )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        yearAlign = 1948
        yearFirst = 1952
        yearLast  = 1952
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign )
        secIn = 0
        mDateIn = yearAlign * 10000 + 12*100 + 1
        call shr_stream_findBounds(streams(1),mDateIn,        secIn,       &
                               &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
                               &   mDateUB,dDateUB,secUB,n_ub,fileUB, rc=rCode )
        exp_int = 1
     else if ( series == 12 )then
        test_descrip = "One file is out of sequence"
        filenames2     = clmncep
        filenames2(2)  = clmncep(4)
        filenames2(4)  = clmncep(2)
        call streams_exp_set( stream_exp, datasource="CLMNCEP", &
                              fldListfile ="TBOT:QBOT:WIND:PRECTmms:FSDS:PSRF",      &
                              fldListModel="tbot:qbot:wind:prectMMS:fsds:psrf",      &
                              filepath="/fs/cgd/csm/inputdata/atm/datm7/CLMNCEP/",   &
                              domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
                              domfilename="domain.T62.050609.nc", &
                              nfiles=12, filenames=filenames2 )
        call streams_exp_write_strm_txt( stream_filename, stream_exp )
        yearAlign = 1948
        yearFirst = 1948
        yearLast  = 1948
        call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, &
                              yearAlign, rc=rCode )
        secIn = 0
        mDateIn = yearAlign * 10000 + 12*100 + 1
        call shr_stream_findBounds(streams(1),mDateIn,        secIn,       &
                               &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
                               &   mDateUB,dDateUB,secUB,n_ub,fileUB, rc=rCode )
        exp_int = 1
!    else if ( series == 12 )then
!       test_descrip = "year range is out of bounds"
!       call streams_exp_set( stream_exp, datasource="CLMNCEP", &
!                             fldListfile ="TBOT:QBOT:WIND:PRECTmms:FSDS:PSRF",      &
!                             fldListModel="tbot:qbot:wind:prectMMS:fsds:psrf",      &
!                             filepath="/fs/cgd/csm/inputdata/atm/datm7/CLMNCEP/",   &
!                             domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
!                             domfilename="domain.T62.050609.nc", &
!                             nfiles=12, filenames=clmncep(1:12) )
!       yearAlign = 1948
!       yearFirst = 1948
!       yearLast  = 1972
!       call streams_exp_write_strm_txt( stream_filename, stream_exp )
!       call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, yearAlign, rCode )
!       secIn = 0
!       mDateIn = yearAlign * 10000 + 12*100 + 1
!       call shr_stream_findBounds(streams(1),mDateIn,        secIn,       &
!                              &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
!                              &   mDateUB,dDateUB,secUB,n_ub,fileUB )
!       exp_int = 1
!    else if ( series == 13 )then
!       test_descrip = "Dates are backwards"
!       call streams_exp_set( stream_exp, datasource="CLMNCEP", &
!                             fldListfile ="TBOT:QBOT:WIND:PRECTmms:FSDS:PSRF",      &
!                             fldListModel="tbot:qbot:wind:prectMMS:fsds:psrf",      &
!                             filepath="/fs/cgd/csm/inputdata/atm/datm7/CLMNCEP/",   &
!                             domfilepath="/fs/cgd/csm/inputdata/atm/datm7/",        &
!                             domfilename="domain.T62.050609.nc", &
!                             nfiles=12, filenames=clmncep(12:1:-1) )
!       call streams_exp_write_strm_txt( stream_filename, stream_exp )
!       yearAlign = 1948
!       yearFirst = 1948
!       yearLast  = 1948
!       call shr_stream_init( streams(1), stream_filename, yearFirst, yearLast, yearAlign, rCode )
!       secIn = 0
!       mDateIn = yearAlign * 10000 + 12*100 + 1
!       call shr_stream_findBounds(streams(1),mDateIn,        secIn,       &
!                              &   mDateLB,dDateLB,secLB,n_lb,fileLB,   &
!                              &   mDateUB,dDateUB,secUB,n_ub,fileUB, rc=rCode )
!       exp_int = 1
     else
        exit
     end if
     write(*,*) trim(test_descrip)
     call test_is( rcode, exp_int, "test that "//trim(test_descrip)//" fails" )
  end do

  call shr_sys_system( "/bin/rm -f "//trim(stream_filename), rcode )
  deallocate( streams  )
  deallocate( streams2 )

  call test_final()

end program test_shr_streams

