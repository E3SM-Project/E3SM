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
  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit, StorWater
  use RtmVar        , only : iulog
  use shr_sys_mod   , only : shr_sys_flush
  
  implicit none
  private

  public WRM_init_StOp_FC

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
  
  subroutine WRM_init_StOp_FC
     ! !DESCRIPTION: define start of operation year - define irrigation releases pattern
     implicit none

     integer nn,nd, idam, j , mth, match, nsub, mth1, mth2, mth3              ! local loop indices
     integer nio,ierror                 ! unit number of a file, flag number of IO status
     real(r8) :: peak                   ! peak value to define the start of operationalyr
     integer :: sgn,curr_sgn, nsc, ct, ct_mx, mth_op       ! number of sign change
     character(len=*),parameter :: subname='WRM_init_StOp_FC'
  
     ! initialize start of the operationnal year based on long term simulation
     ! multiple hydrograph - 1 peak, 2 peaks, multiple small peaks
     !write(iulog,*) subname,"find sign"

     do idam=1,ctlSubwWRM%localNumDam
        WRMUnit%MthStOp(idam) = 6
        !write(iulog,*) subname, idam
        peak = WRMUnit%MeanMthFlow(idam,1)
        match = 1
        do j=2,12
           if ( WRMUnit%MeanMthFlow(idam,j) > peak ) then
              match = j
              peak = WRMUnit%MeanMthFlow(idam,j)
              !initialize else issue with very low flows
              WRMUnit%MthStOp(idam) = j
           end if
        end do
  
        nsc=12 ! way to keep track of problematic reservoir
        ct=1
        ct_mx=1
        mth_op=match
        sgn = 1
        mth = 1

        if ( WRMUnit%MeanMthFlow(idam,13) > 0.05_r8 ) then
           nsc = 0

           if ( abs(peak - WRMUnit%MeanMthFlow(idam,13)) > 0.01_r8 ) then
              sgn =  (peak - WRMUnit%MeanMthFlow(idam,13) ) / abs(peak - WRMUnit%MeanMthFlow(idam,13))
           endif
           curr_sgn = sgn
           do j=1,12
              mth=match+j
              if (mth > 12) then
                 mth = mth-12
              end if
              if ( abs(WRMUnit%MeanMthFlow(idam,mth) - WRMUnit%MeanMthFlow(idam,13)) > 0.01_r8 ) then
                 curr_sgn = (WRMUnit%MeanMthFlow(idam,mth) - WRMUnit%MeanMthFlow(idam,13) ) / abs(WRMUnit%MeanMthFlow(idam,mth) - WRMUnit%MeanMthFlow(idam,13))
              else
                 curr_sgn = sgn
              endif
              if ( curr_sgn .ne. sgn ) then
                 nsc = nsc + 1
                 if ( curr_sgn > 0 .and. nsc > 0 .and. ct > ct_mx) then
                    ct_mx = ct
                    WRMUnit%MthStOp(idam) = mth_op
                 end if
                 mth_op = mth
                 ct=1
              else
                 ct = ct+1
              end if
              sgn = curr_sgn   
           end do  
        endif ! condition on minimum flow of 0.05 
        !write(iulog,*) subname,idam,trim(WRMUnit%DamName(idam))," final start of op year is ",WRMUnit%MthStOp(idam),ct_mx,nsc,WRMUnit%MeanMthFlow(idam,13)

        ! FC part
        !only for flow larger than 1ms
        if (WRMUnit%MeanMthFlow(idam,13) > 1._r8 ) then 
           j=0
           match = 0
           do while (j< 8) 
              j = j + 1
              mth =  WRMUnit%MthStOp(idam) - j
              if ( mth < 1 ) then
                 mth = mth + 12
              endif
              mth1 = WRMUnit%MthStOp(idam) - j + 1
              if ( mth1 < 1 ) then
                 mth1 = mth1 + 12 
              endif
              mth2 = WRMUnit%MthStOp(idam) - j - 1
              if ( mth2 < 1 ) then
                 mth2 = mth2 + 12
              endif
              !write(iulog,*) subname,  WRMUnit%MthStOp(idam), mth, mth1, mth2   
              !write(iulog,*) subname, WRMUnit%MeanMthFlow(idam,13), WRMUnit%MeanMthFlow(idam,mth), WRMUnit%MeanMthFlow(idam,mth1), WRMUnit%MeanMthFlow(idam,mth2)
              call shr_sys_flush(iulog)
              if ( (WRMUnit%MeanMthFlow(idam,mth) >= WRMUnit%MeanMthFlow(idam,13)) .and. (WRMUnit%MeanMthFlow(idam,mth2) <= WRMUnit%MeanMthFlow(idam,13)).and. (match == 0 )) then
                 WRMUnit%MthNdFC(idam) = mth
                 match = 1
              endif
              if ( (WRMUnit%MeanMthFlow(idam,mth) <= WRMUnit%MeanMthFlow(idam,mth1)) .and. (WRMUnit%MeanMthFlow(idam,mth) <= WRMUnit%MeanMthFlow(idam,mth2)) .and. (WRMUnit%MeanMthFlow(idam,mth) <= WRMUnit%MeanMthFlow(idam,13))) then
                 WRMUnit%MthStFC(idam) = mth
                 j=12 !get out of the loop
                 ! twist for hydropower - need to maintain flow
                 if ( WRMUnit%use_Elec(idam) > 0 ) then
                    mth3 =  mth2 - 1
                    if ( mth3 < 1 ) then
                       mth3 = mth3 + 12
                    endif
                  
                    !if (WRMUnit%MeanMthFlow(idam,mth2) <= WRMUnit%MeanMthFlow(idam,13) ) then
                    !  WRMUnit%MthStFC(idam) = mth2
                       !if (WRMUnit%MeanMthFlow(idam,mth3) <= WRMUnit%MeanMthFlow(idam,13) ) then
                       !  WRMUnit%MthStFC(idam) = mth3
                    !endif
                    !endif
                 endif
                    
              endif
           enddo  
           write(iulog,*) subname, idam, "final start of FC op is ", WRMUnit%MthStFC(idam) , WRMUnit%MthNdFC(idam), WRMUnit%MthStOp(idam)
           write(iulog,*) subname, WRMUnit%StorCap(idam)
           ! ENForce the FC targets
           !write(iulog,*) subname,"calibration? ",WRMUnit%StorageCalibFlag(idam), WRMUnit%MinStorTarget(idam), WRMUnit%MaxStorTarget(idam)
           if ( WRMUnit%use_FCon(idam) > 0 .and. WRMUnit%MthStFC(idam).eq.0) then
              !write(iulog,*) subname,idam, "DOUBLE CHECK start of FC - run on the river or mostly irrigation and a little FC"
              mth = WRMUnit%MthNdFC(idam)-2
              if (mth<0) then
                 mth = mth + 12
              endif
              WRMUnit%MthStFC(idam) =  mth
              !write(iulog,*) subname,idam, "DOUBLE CHECK start of FC - run on the river or mostly irrigation and a little FC", WRMUnit%MthStFC(idam)
              !write(iulog,*) subname, WRMUnit%MeanMthFlow(idam, 13)
              !call shr_sys_abort(subname//' ERROR start of FC')
           end if
                
        end if
        !write(iulog,*) subname, peak, (WRMUnit%MeanMthFlow(idam,mth), mth=1,13)
     end do
     !call shr_sys_flush(iulog)
     write(iulog,*) subname,' Done'
     call shr_sys_flush(iulog)
     ! one peak - nsc=2, two peaks - nsc = 4. Choose th eone with the longest period under mean annual flow

     ! Bieman and Haddeland's way - issue when multiple peaks in mountainous regions like the NW in tributaries, or in tropics

 end subroutine WRM_init_StOp_FC

end MODULE WRM_start_op_year

