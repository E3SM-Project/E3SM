!
MODULE WRM_returnflow
! Description: core code of the WRM. 
! 
! Developed by Nathalie Voisin Feb 2012
! REVISION HISTORY:
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use RunoffMod, only : Trunoff, Tctl, Tunit
  use rof_cpl_indices, only : nt_nliq
  use WRM_type_mod, only : TWRMctl => ctlSubwWRM, WRMUnit, StorWater

  implicit none
  private

  public insert_returnflow_channel
  
! !PUBLIC MEMBER FUNCTIONS:

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
!________________________________________________________________________________________________
!
!_______________________________________________________________________________________________
  subroutine estimate_returnflow_deficit()
  ! !DESCRIPTION: subnetwork channel routing irrigation extraction
     implicit none    
     integer :: iunit      ! local index
     real(r8) :: frac, firr, fnonirr , tempfloatI, tempfloatNI ! flow in cubic meter rather than cms
               
     StorWater%ReturnIrrig = 0._r8
     StorWater%SuppIrrig = 0._r8
     StorWater%ReturnNonIrrig = 0._r8
     StorWater%SuppNonIrrig = 0._r8
     ! previous baseline was consomptive use, no explicit return flow. If No return flow, expect consmtpive us eonly
     if ( TWRMctl%ReturnFlowFlag .eq. 0) then
        if  (TWRMctl%TotalDemandFlag > 0 ) then
           do iunit=1,Tctl%NUnit
              if ( StorWater%ConDemNonIrrig(iunit) > 0._r8 .and. StorWater%supply(iunit) > 0._r8 ) then
                 frac = StorWater%supply(iunit)/StorWater%ConDemNonIrrig(iunit) 
                 if ( frac >= 1._r8) then
                    fnonirr = 1._r8
                 else
                    fnonirr = frac
                 endif
                 StorWater%SuppNonIrrig(iunit) = fnonirr * StorWater%ConDemNonIrrig(iunit)
                 StorWater%SuppIrrig(iunit) = StorWater%supply(iunit) - StorWater%SuppNonIrrig(iunit)
                 StorWater%ConDemNonIrrig(iunit) =  StorWater%ConDemNonIrrig(iunit) - StorWater%SuppNonIrrig(iunit)
                 StorWater%ConDemIrrig(iunit) = StorWater%ConDemIrrig(iunit) - StorWater%SuppIrrig(iunit)
              endif 
           end do
           if ( (TWRMctl%GroundwaterFlag > 0) ) then
              do iunit=1,Tctl%NUnit
                 StorWater%SuppNonIrrig(iunit) = StorWater%SuppNonIrrig(iunit) + StorWater%GWShareNonIrrig(iunit)*StorWater%ConDemNonIrrig(iunit)
                 StorWater%SuppIrrig(iunit) = StorWater%SuppIrrig(iunit) +  StorWater%GWShareIrrig(iunit)* StorWater%ConDemIrrig(iunit)
              enddo
           endif
        else !! no total demand
           if ( (TWRMctl%GroundwaterFlag > 0) ) then
              do iunit=1,Tctl%NUnit
                 StorWater%supply(iunit) = StorWater%supply(iunit) +  StorWater%GWShareIrrig(iunit) * StorWater%demand(iunit)
              enddo
           endif
        endif 
     else ! have return flow option
        if ( (TWRMctl%TotalDemandFlag > 0) ) then

           ! initialize return flow and supply with gw component
           if (TWRMctl%GroundwaterFlag > 0 ) then
              do iunit=1,Tctl%NUnit
                 StorWater%ReturnIrrig(iunit)= StorWater%GWShareIrrig(iunit) * (StorWater%WithDemIrrig(iunit)   - StorWater%ConDemIrrig(iunit) )
                 StorWater%ReturnNonIrrig(iunit)= StorWater%GWShareNonIrrig(iunit) * (StorWater%WithDemNonIrrig(iunit) - StorWater%ConDemNonIrrig(iunit))
                 StorWater%SuppIrrig(iunit) = StorWater%GWShareIrrig(iunit)*StorWater%WithDemIrrig(iunit)
                 StorWater%SuppNonIrrig(iunit) = StorWater%GWShareNonIrrig(iunit)*StorWater%WithDemNonIrrig(iunit)
              enddo
           endif

           do iunit=1,Tctl%NUnit
              tempfloatI = 0._r8
              tempfloatNI = 0._r8
              if ( StorWater%WithDemNonIrrig(iunit) > 0._r8 .and. StorWater%supply(iunit) > 0._r8 ) then
                 frac = StorWater%supply(iunit)/StorWater%WithDemNonIrrig(iunit)
                 if ( frac >= 1._r8) then
                    fnonirr = 1._r8
                 else
                    fnonirr = frac
                 endif 
                 StorWater%ReturnNonIrrig(iunit)=  StorWater%ReturnNonIrrig(iunit) + fnonirr * ( StorWater%WithDemNonIrrig(iunit) - StorWater%ConDemNonIrrig(iunit))
                 !StorWater%SuppNonIrrig(iunit) =  fnonirr * StorWater%WithDemNonIrrig(iunit)
                 tempfloatNI = fnonirr * StorWater%WithDemNonIrrig(iunit)
                 StorWater%SuppNonIrrig(iunit) = StorWater%SuppNonIrrig(iunit) + tempfloatNI
                 StorWater%WithDemNonIrrig(iunit) =  StorWater%WithDemNonIrrig(iunit) - tempfloatNI
              endif
              !StorWater%SuppIrrig(iunit) = StorWater%supply(iunit) - StorWater%SuppNonIrrig(iunit)
              StorWater%SuppIrrig(iunit) = StorWater%SuppIrrig(iunit) + StorWater%supply(iunit) - tempfloatNI
              tempfloatI = StorWater%supply(iunit) - tempfloatNI

              !if ( StorWater%WithDemIrrig(iunit) > 0._r8 .and. StorWater%SuppIrrig(iunit) > 0._r8 ) then
              if ( StorWater%WithDemIrrig(iunit) > 0._r8 .and. tempfloatI  > 0._r8 ) then
                 !frac = StorWater%SuppIrrig(iunit)/StorWater%WithDemIrrig(iunit)
                 frac = tempfloatI /StorWater%WithDemIrrig(iunit)
                 if ( frac > 1._r8) then
                    firr = 1._r8
                    print*,"Error, supply larger than demand", iunit, frac, StorWater%supply(iunit), StorWater%demand(iunit), StorWater%WithDemIrrig(iunit), StorWater%WithDemNonIrrig(iunit), StorWater%SuppNonIrrig(iunit), StorWater%SuppIrrig(iunit)
                    !stop
                 else
                    firr = frac
                 endif
                 StorWater%ReturnIrrig(iunit)= StorWater%ReturnIrrig(iunit) + firr * (StorWater%WithDemIrrig(iunit) - StorWater%ConDemIrrig(iunit)) 
                 !StorWater%WithDemIrrig(iunit) = StorWater%WithDemIrrig(iunit) - StorWater%SuppIrrig(iunit)
                 StorWater%WithDemIrrig(iunit) = StorWater%WithDemIrrig(iunit) -  tempfloatI
              end if
           end do
        else  ! total demand flag
           do iunit=1,Tctl%NUnit
              if ( StorWater%supply(iunit) > 0._r8  ) then
                 frac = StorWater%supply(iunit) / StorWater%demand(iunit)
                 if ( frac > 1._r8) then
                    firr = 1._r8
                    print*,"Error, supply larger than demand", iunit, frac, StorWater%supply(iunit), StorWater%demand(iunit), StorWater%WithDemIrrig(iunit)
                    !stop
                 else
                    firr = frac
                 endif
                 StorWater%ReturnIrrig(iunit)= StorWater%ReturnIrrig(iunit) + firr * (StorWater%WithDemIrrig(iunit) - StorWater%ConDemIrrig(iunit))
                 StorWater%SuppIrrig(iunit) = StorWater%SuppIrrig(iunit) + StorWater%supply(iunit) 
                 StorWater%WithDemIrrig(iunit) = StorWater%WithDemIrrig(iunit) - StorWater%Supply(iunit)
              end if
           end do
        end if 
     end if 

  end subroutine estimate_returnflow_deficit

!__________________________________________________________________________________________________
!
!_________________________________________________________________________________________________
   subroutine insert_returnflow_soilcolumn()
      implicit none
      integer :: iunit      ! local index
      real(r8) :: frac, firr, fnonirr, temp  ! flow in cubic meter rather than cms

      !irrigation demand by default goes into baseflow 
      ! version 1.0: no lag between application and return
      do iunit=1,Tctl%NUnit
         temp = StorWater%ReturnIrrig(iunit) / ( TUnit%area(iunit) * TUnit%frac(iunit) * Tctl%DATAH ) !m/seconds! 
         Trunoff%qsub(iunit,nt_nliq) = Trunoff%qsub(iunit,nt_nliq) + temp
         StorWater%ReturnIrrig(iunit) = 0._r8
      end do
      ! oput non irrigation back into main channel instead
      !if ( (ctlSubwWRM%TotalDemandFlag > 0) ) then
      !  do iunit=1,Tctl%NUnit
      !    Trunoff%qsur(iunit,nt_nliq) = Trunoff%qsur(iunit,nt_nliq) + StorWater%ReturnNonIrrig(iunit) / (TUnit%area(iunit) * TUnit%frac(iunit))
      !    StorWater%ReturnNonIrrig(iunit) = 0._r8
      !  end do
      !endif
  end subroutine insert_returnflow_soilcolumn
!__________________________________________________________________________________________________
!
!_________________________________________________________________________________________________
  subroutine insert_returnflow_channel(iunit, TheDeltaT)
     implicit none
     integer, intent(in) :: iunit
     real(r8), intent(in) :: theDeltaT
     real(r8) :: flow_vol         ! flow in cubic meter rather than cms
     ! version 1.0: no lag between application and return
     !flow_vol = (-Trunoff%erout(iunit,nt_nliq) + (StorWater%ReturnNonIrrig(iunit) / Tctl%DATAH) ) * theDeltaT ! m3 into m3/s
     !FIX WR Trunoff%erout(iunit,nt_nliq) = Trunoff%erout(iunit,nt_nliq) - StorWater%ReturnNonIrrig(iunit)/Tctl%DATAH
     Trunoff%wr(iunit,nt_nliq) = Trunoff%wr(iunit,nt_nliq) + StorWater%ReturnNonIrrig(iunit)/Tctl%DATAH
     !Trunoff%erout(iunit,nt_nliq) = -flow_vol / (theDeltaT)
  end subroutine insert_returnflow_channel


end MODULE WRM_returnflow
