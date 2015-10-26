!
MODULE WRM_start_op_year 
! Description: module to provide initialization information for the WRM mmodel
! 
! Developed by Nathalie Voisin 
! REVISION HISTORY: 2012
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use WRM_type_mod, only : ctlSubwWRM, WRMUnit, StorWater
  
  implicit none
  private

  public WRM_init_StOp_FC

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
  
  subroutine WRM_init_StOp_FC
     ! !DESCRIPTION: define start of operation year - define irrigation releases pattern
     implicit none

     integer nn,nd, iunit, j , mth, match, nsub, mth1, mth2, mth3              ! local loop indices
     integer nio,ierror                 ! unit number of a file, flag number of IO status
     real(r8) :: peak                   ! peak value to define the start of operationalyr
     integer :: sgn,curr_sgn, nsc, ct, ct_mx, mth_op       ! number of sign change
  
     ! initialize start of the operationnal year based on long term simulation
     ! multiple hydrograph - 1 peak, 2 peaks, multiple small peaks
     print*,"find sign"
     do iunit=1,ctlSubwWRM%NDam
        WRMUnit%MthStOp(iunit) = 6
        print*, iunit
        peak = WRMUnit%MeanMthFlow(iunit,1)
        match = 1
        do j=2,12
           if ( WRMUnit%MeanMthFlow(iunit,j) > peak ) then
              match = j
              peak = WRMUnit%MeanMthFlow(iunit,j)
              !initialize else issue with very low flows
              WRMUnit%MthStOp(iunit) = j
           end if
        end do
  
        nsc=12 ! way to keep track of problematic reservoir
        ct=1
        ct_mx=1
        mth_op=match
        sgn = 1
        mth = 1

        if ( WRMUnit%MeanMthFlow(iunit,13) > 0.05_r8 ) then
           nsc = 0

           if ( abs(peak - WRMUnit%MeanMthFlow(iunit,13)) > 0.01_r8 ) then
              sgn =  (peak - WRMUnit%MeanMthFlow(iunit,13) ) / abs(peak - WRMUnit%MeanMthFlow(iunit,13))
           endif
           curr_sgn = sgn
           do j=1,12
              mth=match+j
              if (mth > 12) then
                 mth = mth-12
              end if
              if ( abs(WRMUnit%MeanMthFlow(iunit,mth) - WRMUnit%MeanMthFlow(iunit,13)) > 0.01_r8 ) then
                 curr_sgn = (WRMUnit%MeanMthFlow(iunit,mth) - WRMUnit%MeanMthFlow(iunit,13) ) / abs(WRMUnit%MeanMthFlow(iunit,mth) - WRMUnit%MeanMthFlow(iunit,13))
              else
                 curr_sgn = sgn
              endif
              if ( curr_sgn .ne. sgn ) then
                 nsc = nsc + 1
                 if ( curr_sgn > 0 .and. nsc > 0 .and. ct > ct_mx) then
                    ct_mx = ct
                    WRMUnit%MthStOp(iunit) = mth_op
                 end if
                 mth_op = mth
                 ct=1
              else
                 ct = ct+1
              end if
              sgn = curr_sgn   
           end do  
        endif ! condition on minimum flow of 0.05 
        print*,iunit,WRMUnit%DamName(iunit)," final start of op year is ",WRMUnit%MthStOp(iunit),ct_mx,nsc,WRMUnit%MeanMthFlow(iunit,13)

        ! FC part
        !only for flow larger than 1ms
        if (WRMUnit%MeanMthFlow(iunit,13) > 1._r8 ) then 
           j=0
           match = 0
           do while (j< 8) 
              j = j + 1
              mth =  WRMUnit%MthStOp(iunit) - j
              if ( mth < 1 ) then
                 mth = mth + 12
              endif
              mth1 = WRMUnit%MthStOp(iunit) - j + 1
              if ( mth1 < 1 ) then
                 mth1 = mth1 + 12 
              endif
              mth2 = WRMUnit%MthStOp(iunit) - j - 1
              if ( mth2 < 1 ) then
                 mth2 = mth2 + 12
              endif
              !print*,  WRMUnit%MthStOp(iunit), mth, mth1, mth2   
              !print*, WRMUnit%MeanMthFlow(iunit,13), WRMUnit%MeanMthFlow(iunit,mth), WRMUnit%MeanMthFlow(iunit,mth1), WRMUnit%MeanMthFlow(iunit,mth2)
              if ( (WRMUnit%MeanMthFlow(iunit,mth) >= WRMUnit%MeanMthFlow(iunit,13)) .and. (WRMUnit%MeanMthFlow(iunit,mth2) <= WRMUnit%MeanMthFlow(iunit,13)).and. (match == 0 )) then
                 WRMUnit%MthNdFC(iunit) = mth
                 match = 1
              endif
              if ( (WRMUnit%MeanMthFlow(iunit,mth) <= WRMUnit%MeanMthFlow(iunit,mth1)) .and. (WRMUnit%MeanMthFlow(iunit,mth) <= WRMUnit%MeanMthFlow(iunit,mth2)) .and. (WRMUnit%MeanMthFlow(iunit,mth) <= WRMUnit%MeanMthFlow(iunit,13))) then
                 WRMUnit%MthStFC(iunit) = mth
                 j=12 !get out of the loop
                 ! twist for hydropower - need to maintain flow
                 if ( WRMUnit%use_Elec(iunit) > 0 ) then
                    mth3 =  mth2 - 1
                    if ( mth3 < 1 ) then
                       mth3 = mth3 + 12
                    endif
                  
                    !if (WRMUnit%MeanMthFlow(iunit,mth2) <= WRMUnit%MeanMthFlow(iunit,13) ) then
                    !  WRMUnit%MthStFC(iunit) = mth2
                       !if (WRMUnit%MeanMthFlow(iunit,mth3) <= WRMUnit%MeanMthFlow(iunit,13) ) then
                       !  WRMUnit%MthStFC(iunit) = mth3
                    !endif
                    !endif
                 endif
                    
              endif
           enddo  
           print*, iunit, "final start of FC op is ", WRMUnit%MthStFC(iunit) , WRMUnit%MthNdFC(iunit), WRMUnit%MthStOp(iunit)
           ! ENForce the FC targets
           print*,"calibration? ",WRMUnit%StorageCalibFlag(iunit), WRMUnit%MinStorTarget(iunit), WRMUnit%MaxStorTarget(iunit)
           if ( WRMUnit%use_FCon(iunit) > 0 .and. WRMUnit%MthStFC(iunit).eq.0) then
              print*,iunit, "DOUBLE CHECK start of FC - run on the river or mostly irrigation and a little FC"
              mth = WRMUnit%MthNdFC(iunit)-2
              if (mth<0) then
                 mth = mth + 12
              endif
              WRMUnit%MthStFC(iunit) =  mth
              print*,iunit, "DOUBLE CHECK start of FC - run on the river or mostly irrigation and a little FC", WRMUnit%MthStFC(iunit)
              !print* WRMUnit%MeanMthFlow(iunit, 13)
              !stop
           end if
                
        end if
        !print*, peak, (WRMUnit%MeanMthFlow(iunit,mth), mth=1,13)
     end do
     ! one peak - nsc=2, two peaks - nsc = 4. Choose th eone with the longest period under mean annual flow

     ! Bieman and Haddeland's way - issue when multiple peaks in mountainous regions like the NW in tributaries, or in tropics

 end subroutine WRM_init_StOp_FC

end MODULE WRM_start_op_year

