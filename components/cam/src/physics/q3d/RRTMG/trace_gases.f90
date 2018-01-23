MODULE trace_gases
    
      USE shr_kind_mod,   only: dbl_kind => shr_kind_r8, real_kind => shr_kind_r4

      IMPLICIT NONE
      PRIVATE

      include "mpif.h"

      PUBLIC :: trace_gas_input
      
      REAL (KIND=real_kind), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: &
          o3, co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4

      INTEGER, PUBLIC :: mpierr

      CONTAINS 
      
!==============================================================================
      SUBROUTINE TRACE_GAS_INPUT(IM, JM, KM, PL, P0)
!==============================================================================      
      
      USE rrtm_params, only: ggr
      USE rrtm_grid,   only: masterproc

      INTEGER, INTENT(IN) :: &
          IM, JM, KM    ! Number of grid points in x-, y-, and z-directions
      
      REAL (KIND = dbl_kind), INTENT(IN), DIMENSION(km) :: &
          PL     ! model layer pressure (mb)

      REAL (KIND = dbl_kind), INTENT(IN), DIMENSION(1:km+1) :: &
          P0     ! model interface pressure (mb)
      
      INTEGER :: I, J, K, L, M
          
      INTEGER, PARAMETER :: &
          nzls = 59, &         ! Number of levels in fort.91 datafile
          nTraceGases = 9, &   ! Number of trace gases
          nPress = nzls        ! Number of pressure levels in raw trace gas profiles

! Trace gas profiles from fort.91

      REAL (KIND = real_kind), DIMENSION(nzls) :: &
          pMLS, o3MLS, co2MLS, n2oMLS, ch4MLS, o2MLS, &
          cfc11MLS, cfc12MLS, cfc22MLS, ccl4MLS

      REAL (KIND = real_kind), DIMENSION(km+1) :: tmppres, tmpTrace
      REAL (KIND = real_kind), DIMENSION(km+2, nTraceGases) :: trpath
      REAL (KIND = real_kind), DIMENSION(nTraceGases,nPress) :: trace
      REAL (KIND = real_kind), DIMENSION(km+2) :: tmppresi
      REAL (KIND = dbl_kind) :: plow, pupp, pmid, wgtlow, wgtupp, godp

      CHARACTER(LEN = 5), DIMENSION(nTraceGases), PARAMETER :: &
          TraceGasNameOrder = (/        &
     				'O3   ',  &
     				'CO2  ',  &
     				'CH4  ',  &
     				'N2O  ',  & 
     				'O2   ',  &
     				'CFC11',  &
     				'CFC12',  &
     				'CFC22',  &
     				'CCL4 '  /)

!------------------------------------------------------------------------------

      IF (.NOT. (ALLOCATED(o3)) )    ALLOCATE(o3(im,jm,km))
      IF (.NOT. (ALLOCATED(co2)) )   ALLOCATE(co2(im,jm,km))
      IF (.NOT. (ALLOCATED(n2o)) )   ALLOCATE(n2o(im,jm,km))
      IF (.NOT. (ALLOCATED(ch4)) )   ALLOCATE(ch4(im,jm,km))
      IF (.NOT. (ALLOCATED(o2)) )    ALLOCATE(o2(im,jm,km))
      IF (.NOT. (ALLOCATED(cfc11)) ) ALLOCATE(cfc11(im,jm,km))
      IF (.NOT. (ALLOCATED(cfc12)) ) ALLOCATE(cfc12(im,jm,km))
      IF (.NOT. (ALLOCATED(cfc22)) ) ALLOCATE(cfc22(im,jm,km))
      IF (.NOT. (ALLOCATED(ccl4)) )  ALLOCATE(ccl4(im,jm,km))

!------------------------------------------------------------------------------
! Read in trace gas profiles, provided by the rrtmg_lw.nc datafile supplied in the 
! RRTMG code distribution.

  if(masterproc) then
 	  read(91,*)
 	  read(91,*)

	  do 10 k=nzls,1,-1
        read(91,*) pMLS(k), o3MLS(k), co2MLS(k), ch4MLS(k), n2oMLS(k), o2MLS(k), &
                     cfc11MLS(k), cfc12MLS(k), cfc22MLS(k), ccl4MLS(k)
  10  continue
      close (91)
  999 format(f8.2,9a12)
  
      print*,' '
      print*,'RRTMG rrtmg_lw.nc trace gas profiles: '
      write(*,1001) 'k','p (mb)', TraceGasNameOrder
      do 1010 k = nPress,1,-1
        write(*,1002) k, pMLS(k),o3MLS(k),co2MLS(k),ch4MLS(k),n2oMLS(k),o2MLS(k),&
                     cfc11MLS(k),cfc12MLS(k),cfc22MLS(k),ccl4MLS(k)
 1010 continue
   endif
   
      CALL MPI_BCAST( pmls,     nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( o3MLS,    nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( co2MLS,   nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ch4MLS,   nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( n2oMLS,   nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( o2MLS,    nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( cfc11MLS, nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( cfc12MLS, nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( cfc22MLS, nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ccl4MLS,  nzls, MPI_REAL, 0, MPI_COMM_WORLD, mpierr)
  
 1001 format(A4,A8, 9A12)
 1002 format(i4,f8.2,9e12.4)
 
      trace(1,:) = o3MLS(:)
      trace(2,:) = co2MLS(:)
      trace(3,:) = ch4MLS(:)
      trace(4,:) = n2oMLS(:)
      trace(5,:) = o2MLS(:)
      trace(6,:) = cfc11MLS(:)
      trace(7,:) = cfc12MLS(:)
      trace(8,:) = cfc22MLS(:)
      trace(9,:) = ccl4MLS(:)

!===================================================================================================
!bloss: modify routine to compute trace gas paths from surface to
! top of supplied sounding.  Then, interpolate these paths onto the
!  interface pressure levels of the model grid, with an extra level
!  at the top for the overlying atmosphere.  Differencing these
!   paths and dividing by dp/g will give the mean mass concentration
!    in that level.

!  This procedure has the advantage that the total trace gas path
!   will be invariant to changes in the vertical grid.

! pressure sounding
      tmppres(1:km)    = pl(1:km)  ! pressure at model levels (mb)
      tmppresi(1:km+1) = p0(1:km+1)   ! pressure at model interfaces (mb)
    
! add a level for the overlying atmosphere.
      tmppres(km+1)  = 0.5_dbl_kind*p0(km+1) ! half of pressure at top of model
      tmppresi(km+2) = MIN(1.e-4_real_kind,0.25_dbl_kind*tmppres(km+1)) ! near-zero pressure at top of extra laye

! trace gas paths at surface are zero.
      trpath(1,:) = 0.0_dbl_kind

! start with trace path at interface below

!------------------------------------------------------------------------------
      do 100 k = 2,km+2
        trpath(k,:) = trpath(k-1,:)

! if pressure greater than sounding, assume concentration at bottom.

! dp/g * tr
        if (tmppresi(k-1).gt.pMLS(1)) then
          trpath(k,:) = trpath(k,:) &
               + (tmppresi(k-1) - MAX(tmppresi(k),pMLS(1)))/ggr &
               *trace(:,1)
        end if

! limit pMLS(m:m-1) so that they are within the model level
!  tmppresi(k-1:k).

        do 110 m = 2,nPress
          plow = MIN(tmppresi(k-1),MAX(tmppresi(k),pMLS(m-1)))
          pupp = MIN(tmppresi(k-1),MAX(tmppresi(k),pMLS(m)))

          if(plow.gt.pupp) then
            pmid = 0.5_dbl_kind*(plow+pupp)

            wgtlow = (pmid-pMLS(m))/(pMLS(m-1)-pMLS(m))
            wgtupp = (pMLS(m-1)-pmid)/(pMLS(m-1)-pMLS(m))
!!$          write(*,*) pMLS(m-1),pmid,pMLS(m),wgtlow,wgtupp

! include this level of the sounding in the trace gas path
            trpath(k,:) = trpath(k,:) &
                 + (plow - pupp)/ggr*(wgtlow*trace(:,m-1)+wgtupp*trace(:,m)) ! dp/g*tr
          end if
  110   continue

! if pressure is off top of trace gas sounding, assume
!  concentration at top
        if (tmppresi(k).lt.pMLS(nPress)) then
          trpath(k,:) = trpath(k,:) &
               + (MIN(tmppresi(k-1),pMLS(nPress)) - tmppresi(k))/ggr & ! dp/g
               *trace(:,nPress)                               ! *tr
        end if
  100 continue

!------------------------------------------------------------------------------

      do 200 m = 1,nTraceGases
        do 210 k = 1,km+1
          godp = ggr/(tmppresi(k) - tmppresi(k+1))
          tmpTrace(k) = (trpath(k+1,m) - trpath(k,m))*godp
  210   continue
        if(TRIM(TraceGasNameOrder(m))=='O3') then
          o3(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='CO2') then
          co2(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='CH4') then
          ch4(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='N2O') then
          n2o(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='O2') then
          o2(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='CFC11') then
          cfc11(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='CFC12') then
          cfc12(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='CFC22') then
          cfc22(1,1,1:km) = tmpTrace(1:km)
        elseif(TRIM(TraceGasNameOrder(m))=='CCL4') then
          ccl4(1,1,1:km) = tmpTrace(1:km)
        end if
200 continue

  if(masterproc) then
      print*,' '
      print*,'Interpolated trace gas vertical profiles (mass mixing ratio -- g/g):'
      write(*,101) 'p (hPa)', (TraceGasNameOrder(m),m=1,nTraceGases)
      do 300 k=1,km
        write(*,102) tmppres(k),o3(1,1,k),co2(1,1,k),ch4(1,1,k),n2o(1,1,k),o2(1,1,k), &
             cfc11(1,1,k),cfc12(1,1,k), cfc22(1,1,k),ccl4(1,1,k)
  300 continue
   endif
   
  101 FORMAT(A8, 9A18)
  102 FORMAT(F10.2, 9E18.8)

      do 400 i=1,im
      do 400 j=1,jm
      do 400 k=1,km
        o3(i,j,k)    = o3(1,1,k)
        n2o(i,j,k)   = n2o(1,1,k)
        ch4(i,j,k)   = ch4(1,1,k)
        cfc11(i,j,k) = cfc11(1,1,k)
        cfc12(i,j,k) = cfc12(1,1,k)
        co2(i,j,k)   = co2(1,1,k)
        o2(i,j,k)    = o2(1,1,k)
        cfc22(i,j,k) = cfc22(1,1,k)
        ccl4(i,j,k)  = ccl4(1,1,k)
  400 continue
  
      END SUBROUTINE trace_gas_input
      
      END MODULE trace_gases
