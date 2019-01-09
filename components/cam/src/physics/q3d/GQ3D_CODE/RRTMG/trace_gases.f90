MODULE trace_gases
    
      USE shr_kind_mod, only: r8 => shr_kind_r8
      USE rrtm_grid,    only: nx, ny, nzm
      
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: trace_gas_input

!JUNG      REAL, DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: &
!JUNG          o3, co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4      
     
      REAL, DIMENSION(nx,ny,nzm), PUBLIC :: &
          o3, co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4

      CONTAINS 
      
!==============================================================================
      SUBROUTINE TRACE_GAS_INPUT(IM, JM, KM, PL, P0)
!==============================================================================      
      USE rrtm_params, only: ggr

      INTEGER, INTENT(IN) :: IM, JM, KM    ! Number of grid points in x-, y-, and z-directions
      
      REAL (KIND = r8), INTENT(IN), DIMENSION(km) :: PL     ! model layer pressure (mb)
      REAL (KIND = r8), INTENT(IN), DIMENSION(1:km+1) :: P0 ! model interface pressure (mb)
      
      INTEGER :: I, J, K, L, M
          
      INTEGER, PARAMETER :: &
          nzls = 59,        &  ! Number of levels of data
          nTraceGases = 9,  &  ! Number of trace gases
          nPress = nzls        ! Number of pressure levels in raw trace gas profiles

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
     				
      REAL, DIMENSION(nzls) :: pMLS, o3MLS, co2MLS, n2oMLS, ch4MLS, o2MLS, &
                               cfc11MLS, cfc12MLS, cfc22MLS, ccl4MLS

      REAL, DIMENSION(km+1) :: tmppres, tmpTrace
      REAL, DIMENSION(km+2, nTraceGases) :: trpath
      REAL, DIMENSION(nTraceGases,nPress) :: trace
      REAL, DIMENSION(km+2) :: tmppresi
      
      REAL (KIND = r8) :: plow, pupp, pmid, wgtlow, wgtupp, godp

!JUNG      IF (.NOT. (ALLOCATED(o3)))    ALLOCATE(o3(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(co2)))   ALLOCATE(co2(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(n2o)))   ALLOCATE(n2o(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(ch4)))   ALLOCATE(ch4(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(o2)))    ALLOCATE(o2(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(cfc11))) ALLOCATE(cfc11(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(cfc12))) ALLOCATE(cfc12(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(cfc22))) ALLOCATE(cfc22(im,jm,km))
!JUNG      IF (.NOT. (ALLOCATED(ccl4)))  ALLOCATE(ccl4(im,jm,km))

!JUNG: Specify the data, instead of reading.
!------------------------------------------------------------------------------
! Read in trace gas profiles, provided by the rrtmg_lw.nc datafile supplied in the 
! RRTMG code distribution.
!
! 	  read(91,*)
! 	  read(91,*)
!
!	  do 10 k=nzls,1,-1
!        read(91,*) pMLS(k), o3MLS(k), co2MLS(k), ch4MLS(k), n2oMLS(k), o2MLS(k), &
!                     cfc11MLS(k), cfc12MLS(k), cfc22MLS(k), ccl4MLS(k)
!  10  continue
!      close (91)
!  
!      print*,' '
!      print*,'RRTMG rrtmg_lw.nc trace gas profiles: '
!      write(*,1001) 'k','p (mb)', TraceGasNameOrder
!      do 1010 k = nPress,1,-1
!        write(*,1002) k, pMLS(k),o3MLS(k),co2MLS(k),ch4MLS(k),n2oMLS(k),o2MLS(k),&
!                     cfc11MLS(k),cfc12MLS(k),cfc22MLS(k),ccl4MLS(k)
! 1010 continue
!  
! 1001 format(A4,A8, 9A12)
! 1002 format(i4,f8.2,9e12.4)
!------------------------------------------------------------------------------
      pMLS(:) = (/ &
      1053.63, 862.64, 706.27, 578.25, 473.43, 387.61, 317.35, 259.82, 212.72, 174.16, &
       142.59, 116.75,  95.58,  78.26,  64.07,  52.46,  42.95,  35.16,  28.79,  23.57, &
        19.30,  15.80,  12.94,  10.59,   8.67,   7.10,   5.81,   4.76,   3.90,   3.19, &
         2.61,   2.14,   1.75,   1.43,   1.17,   0.96,   0.79,   0.64,   0.53,   0.43, &
         0.35,   0.29,   0.24,   0.19,   0.16,   0.13,   0.11,   0.09,   0.07,   0.06, &
         0.05,   0.04,   0.03,   0.03,   0.02,   0.02,   0.01,   0.01,   0.01/)
      
      o3MLS(:) = (/ &
      0.5000E-07, 0.5754E-07, 0.7039E-07, 0.8743E-07, 0.1109E-06, 0.1444E-06, &
      0.1888E-06, 0.2598E-06, 0.3611E-06, 0.5376E-06, 0.7721E-06, 0.9413E-06, &
      0.1153E-05, 0.1854E-05, 0.2920E-05, 0.3856E-05, 0.4901E-05, 0.6064E-05, &
      0.7614E-05, 0.8814E-05, 0.9879E-05, 0.1079E-04, 0.1171E-04, 0.1275E-04, &
      0.1368E-04, 0.1443E-04, 0.1464E-04, 0.1444E-04, 0.1341E-04, 0.1215E-04, &
      0.1046E-04, 0.8894E-05, 0.7429E-05, 0.6362E-05, 0.5440E-05, 0.4679E-05, &
      0.4127E-05, 0.3587E-05, 0.3047E-05, 0.2754E-05, 0.2494E-05, 0.2235E-05, &
      0.1984E-05, 0.1737E-05, 0.1490E-05, 0.1265E-05, 0.1083E-05, 0.9018E-06, &
      0.7201E-06, 0.6035E-06, 0.5169E-06, 0.4303E-06, 0.3437E-06, 0.3173E-06, &
      0.3209E-06, 0.3245E-06, 0.3281E-06, 0.1265E-05, 0.1083E-05/)
      
      co2MLS(:) = (/ &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, &
      0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5393E-03, 0.5388E-03, &
      0.5381E-03, 0.5375E-03, 0.5368E-03, 0.5393E-03, 0.5393E-03/)
      
      ch4MLS(:) = (/ &
      0.9390E-06, 0.9390E-06, 0.9390E-06, 0.9337E-06, 0.9209E-06, 0.9032E-06, &
      0.8892E-06, 0.8611E-06, 0.8352E-06, 0.8143E-06, 0.7946E-06, 0.7734E-06, &
      0.7497E-06, 0.7253E-06, 0.6911E-06, 0.6445E-06, 0.5829E-06, 0.5153E-06, &
      0.4486E-06, 0.4156E-06, 0.3907E-06, 0.3686E-06, 0.3465E-06, 0.3237E-06, &
      0.3025E-06, 0.2844E-06, 0.2663E-06, 0.2483E-06, 0.2304E-06, 0.2124E-06, &
      0.1944E-06, 0.1764E-06, 0.1582E-06, 0.1408E-06, 0.1241E-06, 0.1082E-06, &
      0.1011E-06, 0.9439E-07, 0.8773E-07, 0.8582E-07, 0.8453E-07, 0.8325E-07, &
      0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, &
      0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, &
      0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07, 0.8286E-07/)
      
      n2oMLS(:) = (/ &
      0.4861E-06, 0.4861E-06, 0.4861E-06, 0.4861E-06, 0.4861E-06, 0.4856E-06, &
      0.4790E-06, 0.4615E-06, 0.4469E-06, 0.4328E-06, 0.4203E-06, 0.4021E-06, &
      0.3689E-06, 0.3183E-06, 0.2612E-06, 0.2089E-06, 0.1721E-06, 0.1524E-06, &
      0.1387E-06, 0.1298E-06, 0.1221E-06, 0.1115E-06, 0.1002E-06, 0.8512E-07, &
      0.7153E-07, 0.6073E-07, 0.5010E-07, 0.3959E-07, 0.3200E-07, 0.2520E-07, &
      0.1977E-07, 0.1533E-07, 0.1158E-07, 0.9290E-08, 0.7090E-08, 0.4991E-08, &
      0.4327E-08, 0.3740E-08, 0.3153E-08, 0.2818E-08, 0.2517E-08, 0.2215E-08, &
      0.2004E-08, 0.1834E-08, 0.1663E-08, 0.1516E-08, 0.1415E-08, 0.1314E-08, &
      0.1212E-08, 0.1141E-08, 0.1084E-08, 0.1027E-08, 0.9695E-09, 0.9263E-09, &
      0.8901E-09, 0.8540E-09, 0.8178E-09, 0.1516E-08, 0.1415E-08/)
      
      o2MLS(:) = (/ &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, &
      0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00, 0.2309E+00/)
      
      cfc11MLS(:) = (/ &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/)
      
      cfc12MLS(:) = (/ &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/)      
      
      cfc22MLS(:) = (/ &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/)
      
      ccl4MLS(:) = (/ &   
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/)       
!------------------------------------------------------------------------------
 
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
      tmppres(km+1)  = 0.5*p0(km+1) ! half of pressure at top of model
      tmppresi(km+2) = MIN(1.e-4,0.25*tmppres(km+1)) ! near-zero pressure at top of extra laye

! trace gas paths at surface are zero.
      trpath(1,:) = 0.0

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
            pmid = 0.5_r8*(plow+pupp)

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

!      print*,' '
!      print*,'Interpolated trace gas vertical profiles (mass mixing ratio -- g/g):'
!      write(*,101) 'p (hPa)', (TraceGasNameOrder(m),m=1,nTraceGases)
!      do 300 k=1,km
!        write(*,102) tmppres(k),o3(1,1,k),co2(1,1,k),ch4(1,1,k),n2o(1,1,k),o2(1,1,k), &
!             cfc11(1,1,k),cfc12(1,1,k), cfc22(1,1,k),ccl4(1,1,k)
!  300 continue
!  101 FORMAT(A8, 9A18)
!  102 FORMAT(F10.2, 9E18.8)

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
