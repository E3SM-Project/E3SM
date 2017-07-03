module Mphase_pckr_module

  use Material_module
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private
#include "petsc/finclude/petscsys.h"
  PetscReal, private, parameter :: pckr_sat_water_cut = 1.D0 - 5.D-7
  
  public ::  pckrNH_noderiv, pckrHY_noderiv     
  contains 

! ************************************************************************** !

subroutine pflow_pckr(ipckrtype,pckr_swir,pckr_lambda,pckr_alpha,&
              pckr_m ,pckr_pcmax,sg,pc,pc_s,kr,kr_s,pckr_beta,pckr_pwr) 
       
      implicit none
     
      PetscInt :: ipckrtype
      PetscReal :: sg
      PetscReal :: pckr_swir,pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax,pckr_pwr
      PetscReal :: pc(1:2),pc_s(1:2),kr(1:2),kr_s(1:2)
      PetscReal :: pckr_beta
     
      PetscReal :: sw,se,swir,sw0,lam,ala,um,un,upc,upc_s
      PetscReal :: temp,pcmax,ser
      PetscReal :: uum,pckr_betac,betac,st

    ! if (present(pckr_beta))
      pckr_betac=pckr_beta
      sw=1.D0-sg
      swir=pckr_swir
      sw0=1.d0
      pcmax=pckr_pcmax
       
!     print *,'pflow_pckr: ',ipckrtype,sg,pckr_pcmax

      select case(ipckrtype)

      case(1) ! van Genuchten
      
        ala=pckr_alpha
        um=pckr_m
        un=1.D0/(1.D0-um)

        if (sw>pckr_sat_water_cut) then
          upc=0.D0; upc_s=0.D0
          kr(1)=1.d0; kr(2)=0.d0;
          kr_s(1)=0.d0; kr_s(2)=0.d0;
        elseif (sw > (1.05D0*swir)) then
          se=(sw-swir)/(1.D0-swir)
          temp=se**(-1.D0/um)
          upc=(temp-1.D0)**(1.d0/un)/ala
          upc_s=-1.D0/um/un*upc*(se**(-1.D0-1.D0/um))/(se**(-1.D0/um)-1.d0)
          temp=1.D0/temp
          kr(1)=sqrt(se)*(1.D0-(1.D0-temp)**um)**2.D0
          kr_s(1)=0.5d0*kr(1)/se+2.D0*sqrt(se)*(1.d0-(1.d0-temp)**um)* &
                  (1.d0-temp)**(um-1.d0)*temp/se
          kr(2)=1.D0-kr(1)
          kr_s(2)= -kr_s(1)
!         print *,'in pckr  ',um, sw , ala, se,upc,kr 
        else  ! use linear extropolation
          se=(0.05d0*swir)/(1.D0-swir)
          temp=se**(-1.D0/um)
          upc=(temp-1.D0)**(1.d0/un)/ala
          upc_s=-1.D0/um/un*upc*(se**(-1.D0-1.D0/um))/(se**(-1.D0/um)-1.d0)
          if (sw>swir) then
            temp=1.D0/temp
            kr(1)=sqrt(se)*(1.D0-(1.D0-temp)**um)**2.D0
            kr_s(1)=0.5d0*kr(1)/se+2.d0*sqrt(se)*(1.D0-(1.D0-temp)**um)* &
                     (1.D0-temp)**(um-1.d0)*temp/se
            ser=(sw-swir)/(1.D0-swir)
            upc=upc+(ser-se)*upc_s
!           print *,se,ser,kr(1),kr_s(1)
            kr(1)=kr(1)+(ser-se)*kr_s(1)
            kr(2)=1.D0-kr(1)
            kr_s(2)= -kr_s(1)
          else
            upc=upc-upc_s*se
            upc_s=0.D0
            kr(1)=0.d0
            kr_s(1)=0.d0
            kr(2)=1.d0
            kr_s(2)=0.d0
          end if
        end if

        pc(1)=upc; pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) / temp

      case(2) !Brooks-Corey
      
        lam=pckr_lambda
        ala=pckr_alpha

        if (sw > (1.05D0*swir)) then
          se=(sw-swir)/(sw0-swir)
          upc=se**(-1.D0/lam)/ala
          upc_s=-upc/se/lam
          kr(1)=se**(2.d0/lam+3.d0)
          kr_s(1)=(2.d0/lam+3.d0)*kr(1)/se
          kr(2)=1.D0-kr(1)
          kr_s(1)= -kr_s(1)
        else   ! use linear extropolation
          se=(0.05d0*swir)/(1.D0-swir)
          upc=se**(-1.D0/lam)/ala
          upc_s=-upc/se/lam
          if (sw > swir) then
            kr(1)=se**(2.d0/lam+3.d0)
            kr_s(1)=(2.d0/lam+3.d0)*kr(1)/se
            ser=(sw-swir)/(1.D0-swir)
            upc=upc+(ser-se)*upc_s
            kr(1)=kr(1)+(ser-se)*kr_s(1)
            kr(2)=1.D0-kr(1)
            kr_s(2)=-kr_s(1)
          else
            upc=upc-upc_s*se
            kr(1)=0.D0
            kr_s(1)=0.D0
            kr(2)=1.D0
            kr_s(2)=0.D0
          end if
        end if
        pc(1)=upc; pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) / temp

      case(3) !linear interpolation, need pcmax, assign krmax=1.
         
        if (sw>swir)then
          se=(sw-swir)/(sw0-swir)
          upc=pcmax*(1.D0-se)
          upc_s= - pcmax
          kr(1)=se
          kr(2)=1.D0 - kr(1)
          kr_s(1)=1.D0
          kr_s(2)=-1.D0
        else
          upc=pcmax
          upc_s=0.D0
          kr(1)=0.d0
          kr(2)=1.d0
          kr_s(1)=0.d0
          kr_s(2)=0.d0
        end if
       
        pc(1)=upc; pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) / temp
       
       
      case(4)  ! po model with gas phase residual
        
        if (sw>1.D0) sw=1.D0
        if (sw<0.D-0) sw= 0.D-0
      
        if (sw> pckr_sat_water_cut)then
          upc=0.D0;  kr(1)=1.0;  kr(2)=0.D0;
          upc_s=0.D0; kr_s(1)=0.D0; kr(2)=0.D0;
        else
          ala=pckr_alpha
          betac=pckr_betac 
          um = pckr_m
          uum = 1.D0/um 
          un=1.D0/(1.D0-um)
     !    print *,'pckr:',  ala, betac, pckr_m, pckr_pwr, pcmax,swir
          se=sw  
          temp=se**uum
          upc=(1.D0/temp - 1.D0)**(1.D0 - um) / ala / betac
          upc_s= -((1.D0/temp -1.D0)**(-um))*(se**(-uum-1.D0))*uum*(1.D0-um)&
                    /ala/betac
          
          se=(sw-swir)/(1.D0-swir)
          st= 1.D0
          if (sw>=swir)then
            kr(1)= sqrt(se)*(1.D0-(1.D0-se**uum)**um)**2.D0
!           kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um )**7.D0
            kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um )**pckr_pwr
       
            kr_s(1) = 0.5D0 / sqrt(se)*(1.D0 - (1.D0 - se**uum)**um)**2.D0 &
                + 2.D0 * sqrt(se) * (1.D0- (1.D0-se**uum)**um) &
                * ((1.D0 - se**uum)**(um-1.D0)) * (Se**(uum-1.d0))
!               kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0 -  se**uum)**(7.D0*um) &
!                     - 7.D0 * sqrt(1.D0- se) * (1.D0 - se**uum)**(7.D0*um-1.D0) &
!                     * se**(uum -1.D0)    
            kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0-se**uum)**(pckr_pwr*um) &
                - pckr_pwr * sqrt(1.D0-se) * (1.D0-se**uum)**(pckr_pwr*um-1.D0) &
                * se**(uum-1.D0)    
          else  
            kr(1)=0.D0
            kr_s(1)=0.D0

!           kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um)**7.D0
            kr(2)= sqrt(st-se)* ((1.D0-se**uum)**um)**pckr_pwr
!               kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0 - se**uum)**(7.D0*um) &
!               - 7.D0 * sqrt(1.D0- se) * (1.D0 - se**uum)**(7.D0*um-1.D0) &
!               * se**(uum -1.D0)    
            kr_s(2) = -0.5D0 / sqrt(1.D0-se) * (1.D0 - se**uum)**(pckr_pwr*um) &
                - pckr_pwr * sqrt(1.D0- se) * (1.D0 - se**uum)**(pckr_pwr*um-1.D0) &
                * se**(uum -1.D0)    
          endif
        endif
     
        if (kr(1)<0.D0) kr(1)=0.D0; kr_s(1)=0.D0!; kr(2)=1.D0; kr_s(2)=0.D0
        if (kr(1)>1.D0) kr(1)=1.D0; kr_s(1)=0.D0!;kr(2)=0.D0; kr_s(2)=0.D0
        if (kr(2)<0.D0) kr(2)=0.D0; kr_s(2)=0.D0!;kr(1)=1.D0; kr_s(1)=0.D0
        if (kr(2)>1.D0) kr(2)=1.D0; kr_s(2)=0.D0!;kr(1)=0.D0; kr_s(1)=0.D0

     
        pc(1)=upc; pc(2)=0.d0
        pc_s(1)=upc_s 
        pc_s(2)=0.d0

! since our primary variable for saturation is sg, and all of the
! derivative here is to se based on sw, so need transfer
        temp=sw0-swir
        kr_s(:)=-kr_s(:) / temp
        pc_s(1)=-pc_s(1) 

    !   if (pc(1)>pcmax) print *, 'pckr4: ',sg,pc,kr,pc_s,kr_s
      ! if (sw<pckr_sat_water_cut) print *, sg,pc,kr,pc_s,kr_s
       
      end select

    return

end subroutine pflow_pckr

! ************************************************************************** !
!subroutine pflow_pckr_noderiv_exec(ipckrtype,pckr_sir,pckr_lambda, &
!     pckr_alpha,pckr_m,pckr_pcmax,sg,pc,kr,pckr_beta,pckr_pwr) 

subroutine pflow_pckr_noderiv_exec(ipckrtype,ikrtype,pckr_sir,kr0, &
      pckr_lambda, pckr_alpha,pckr_m,pckr_pcmax,sg,pc,kr,pckr_beta,pckr_pwr) 
  ! 
  ! pckrNH_noderiv: Non-hysteric S-Pc-kr relation excuting routine
  ! Copied from pflotran_orig
  ! 
  ! Author: Chuan Lu
  ! Date: 05/12/08
  ! 
  use Saturation_Function_module
  implicit none 

  
  PetscInt, intent(in) :: ipckrtype, ikrtype
  PetscReal, intent(in) :: pckr_sir(:)
  PetscReal, intent(in) :: kr0(:)
  PetscReal, intent(in) :: pckr_lambda,pckr_alpha,pckr_m,pckr_pcmax
  PetscReal, intent(in) :: pckr_beta,pckr_pwr
  PetscReal, intent(in) :: sg
  PetscReal, intent(out) :: pc(1:2)
  PetscReal, intent(out) :: kr(1:2)
       
  PetscReal :: se,swir,sgir,sw0,lam,ala,um,un,upc,upc_s,kr_s,krg_s
  PetscReal :: temp,ser,pcmax,sw
  PetscReal :: uum,pckr_betac,betac,st
  PetscReal :: se0,upc0,upc_s0
     
      pckr_betac = pckr_beta
      sw = 1.D0 - sg
      if (sw > 1.D0) sw = 1.D0
      if (sw < 0.D0) sw = 0.D0
     
      sw0 = 1.d0
      pcmax = pckr_pcmax
      swir = pckr_sir(1)
      sgir = pckr_sir(2)
      upc = 0.d0

      select case(ipckrtype)

      case(0) ! kr = 1
      
        if (sw >= 1.01D0*swir) then
          kr(1) = 1.D0
        elseif (sw <= swir) then
          kr(1) = 0.D0
        else
          kr(1) = (sw - swir)/swir*1.D2
        endif
       
        if (sg >= 1.01D0*sgir) then
          kr(2) = 1.D0
        elseif (sg <= sgir) then
          kr(2) = 0.D0
        else
          kr(2) = (sg - sgir)/sgir*1.D2
        endif

        kr = 1.D0
       
        upc = 0.D0

      case(VAN_GENUCHTEN) ! van Genuchten

        ala = pckr_alpha
        um = pckr_m
        un = 1.D0/(1.D0 - um)
        if (sw > pckr_sat_water_cut) then
          upc = 0.D0; kr(1) = 1.d0; kr(2) = 0.d0;
        elseif (sw > (1.05D0*swir)) then
          if (sw <= 0.99D0) then
            se = (sw - swir)/(1.D0 - swir)
            temp = se**(-1.D0/um)
            upc = (temp - 1.D0)**(1.d0/un)/ala
            kr(1) = sqrt(se)*(1.D0 - (1.D0 - 1.D0/temp)**um)**2.d0
          ! kr(2) = 1.D0-kr(1)
            kr(2) = sqrt(1.D0 - se)*((1.D0 - se**(1.D0/um))**um)**2.D0
       !    print *,'in pckr nond ',sw,se,upc,kr
          else
            se = (sw - swir)/(1.D0 - swir)
            temp = se**(-1.D0/um)
            kr(1) = sqrt(se)*(1.D0 - (1.D0 - 1.D0/temp)**um)**2.d0 
            kr(2) = sqrt(1.D0 - se)*((1.D0 - se**(1.D0/um))**um)**2.D0
            se = (0.99D0 - swir)/(1.D0 - swir)
            temp = se**(-1.D0/um)
            upc = (temp - 1.D0)**(1.d0/un)/ala
             ! kr(1)=sqrt(se)*(1.D0-(1.D0-1.D0/temp)**um)**2.d0
             ! kr(2)=sqrt(1.D0 - se)*((1.D0-se**(1.D0/um))**um)**2.D0
            upc = upc*(1.D0 - sw)/1.D-2
             ! kr(1) = kr(1) * (1.D0-sw)/5D-3
             ! kr(2) = kr(2) * (1.D0-sw)/5D-3
          endif  
          
        else  ! use linear extropolation
          se0 = (0.05D0*swir)/(1.D0 - swir)
          temp = se0**(-1.D0/um)
          upc0 = (temp - 1.D0)**(1.d0/un)/ala
          upc_s0 = -1.D0/um/un*upc0*(se0**(-1.D0 - 1.D0/um))/(se0**(-1.D0/um) - 1.d0)
          upc_s0 = upc_s0 /(1.D0 - swir)
          if (sw > swir) then
            se = (sw - swir)/(1.D0 - swir)
            temp = se**(-1.D0/um)
            kr(1) = sqrt(se)*(1.D0 - (1.D0 - 1.D0/temp)**um)**2.D0
            kr(2) = sqrt(1.D0 - se)*((1.D0 - se**(1.D0/um))**um)**2.D0
             ! temp=1.D0/temp
             ! kr_s=0.5d0*kr(1)/se+2.d0*sqrt(se)*(1.d0-(1.d0-temp)**um)* &
             !        (1.d0-temp)**(um-1.d0)*temp/se
             ! krg_s=-0.5D0/(1.D0-se)*kr(2) -2.D0*sqrt(1.D0-se)*((1.D0-se**(1.D0/um))**um)&
             !  *((1.D0-se**(1.D0/um))**(um-1.D0)) * (se**(1.D0/um-1.D0))
             ! ser=(sw-swir)/(1.D0-swir)
            upc = upc0 + (sw - 1.05D0 * swir) * upc_s0   
             ! kr(1)=kr(1)+(ser-se)*kr_s
             ! kr(2)=kr(2)+ (ser-se)*krg_s
              
        !  kr(2)=1.D0-kr(1)
          else
            upc = upc0 + (sw - 1.05D0 * swir) * upc_s0
            kr(1) = 0.D0
            kr(2) = 1.D0
          end if
        end if

      case(BROOKS_COREY) !Brooks-Corey
       
        lam = pckr_lambda
        ala = pckr_alpha
       !  swir=pckr_swir
       
        if (sw > (1.05D0*swir)) then
          se = (sw - swir)/(sw0 - swir)
          upc = se**(-1.D0/lam)/ala
          kr(1) = se**(2.d0/lam+3.d0)
            !kr(2) = 1.D0- kr(1)
          kr(2) = (1.D0 - se)**2.D0 * (1.D0 - se**(2.D0/lam + 1.D0)) 
        else   ! use linear extropolation
          se0 = (0.05d0*swir)/(1.D0-swir)
          upc0 = se0**(-1.D0/lam)/ala
          upc_s0 = -upc0/se0/lam
          upc_s0 = upc_s0 /(1.D0 - swir) 
          if (sw > swir) then
            se = (sw - swir)/(1.D0 - swir)
            kr(1) = se**(2.D0/lam+3.d0)
            kr(2) = (1.D0 - se)**2.D0 * (1.D0 - se**(2.D0/lam + 1.D0)) 
            upc = upc0 + (sw - 1.05D0 * swir) * upc_s0

          ! kr_s=(2.d0/lam+3.d0)*kr(1)/se
          ! krg_s = -2.D0*kr(2)/(1.D0-se) -(2.D0+lam)/lam*(1.D0-se)**2.D0*(se**(2.D0/lam))
          ! ser=(sw-swir)/(1.D0-swir)
            
          ! kr(1)=kr(1)+(ser-se)*kr_s
        !   kr(2)=kr(2)+(ser-se)*krg_s
          ! kr(2)=1.D0-kr(1)
          else
            upc = upc0 + (sw - 1.05D0 * swir) * upc_s0
            kr(1) = 0.D0
            kr(2) = 1.D0
          end if
        end if

           
      case(THOMEER_COREY) !linear interpolation, need pcmax, assign krmax=1.
        
        if (sw > swir) then
          se = (sw - swir)/(sw0 - swir)
          upc = pcmax * (1.D0 - se)
          kr(1) = se*se
            !  kr(2)=1.D0 - kr(1)
        else
          upc = pcmax
          kr(1) = 0.d0
             ! kr(2)=1.d0
        end if

        if (sg > sgir) then
          se = (sg - sgir)/(1.D0 - sgir)
        ! upc=pcmax*(1.D0-se)
        ! kr(1)=se
          kr(2) = se*se
        else
        ! upc=pcmax
        ! kr(1)=0.d0
          kr(2) = 0.d0
        end if

      case(NMT_EXP)  ! po model with gas phase residual (NMT)
        
        if (sw > 1.D0) sw = 1.D0
        if (sw < 0.D-0) sw = 1.D-5
        if (sw > pckr_sat_water_cut) then
          upc = 0.D0; kr(1) = 1.0d0; kr(2) = 0.D0;
        else
          ala = pckr_alpha
          betac = pckr_betac 
          um = pckr_m
          uum = 1.D0/um 
          un = 1.D0/(1.D0-um)
            !print *,'pckr:',  ala, betac, pckr_m, pckr_pwr, pcmax,swir
          se = sw
          temp = se**uum
          if (sw < 0.95D0) then 
            upc = (1.D0/temp - 1.D0)**(1.D0 - um) / ala / betac
          else
            temp = 0.95D0**uum
            upc = (1.D0/temp - 1.D0)**(1.D0 - um) / ala / betac
            upc = upc/0.05D0 * (1.D0 - se) 
          endif     
       
          if (upc > pcmax) upc = pcmax
          se = (sw - swir)/(1.D0 - swir)
          st = 1.D0
          if (sw >= swir) then
            kr(1) = sqrt(se)*(1.D0 - (1.D0 - se**uum)**um)**2.D0
!           kr(2)= sqrt(st-se)*((1.D0-se**uum)**um)**7.D0
            kr(2) = sqrt(st-se)*((1.D0-se**uum)**um)**pckr_pwr
          else
!         if (se <= 0.D0) se = 1.D-7
            kr(1) = 0.D0
!         kr(2) = sqrt(st-se)*((1.D0-se**uum)**um)**7.D0
            kr(2) = sqrt(st-se)!*((1.D0-se**uum)**um)**pckr_pwr
          endif
        endif
        
        if (kr(1) < 0.D0) kr(1)=0.D0!; kr(2)=1.D0
        if (kr(1) > 1.D0) kr(1)=1.D0!; kr(2)=0.D0
        if (kr(2) < 0.D0) kr(2)=0.D0!; kr(1)=1.D0
        if (kr(2) > 1.D0) kr(2)=1.D0!; kr(1)=0.D0

      case(PRUESS_1) !linear interpolation, need pcmax, assign krmax=1 (Pruess_1).
        
        if (sw > swir) then
          se = (sw - swir)/(sw0 - swir)
          upc = pcmax * (1.D0 - se)
          kr(1) = se
        else
          upc = pcmax
          kr(1) = 0.d0
        end if

        if (sg > sgir) then
          se = (sg - sgir)/(1.D0 - sgir)
          kr(2) = se
          if (kr(2) > 1.D0) kr(2) = 1.D0
        else
          kr(2) = 0.d0
        end if

      case(VAN_GENUCHTEN_PARKER) ! van Genuchten-Parker
       
        lam = pckr_lambda
        ala = pckr_alpha

!       Water phase using van Genuchten
        um = pckr_m
        un = 1.D0/(1.D0 - um)

        se = (sw - swir)/(sw0 - swir)

        if (sw > swir) then
          temp = se**(-1.D0/um)
!         if (temp < 1.D0+1e-6) temp = 1.D0+1e-6
          upc = (temp - 1.D0)**(1.d0/un)/ala
          if (upc > pcmax) upc = pcmax

!         Mualem rel. perm.
          kr(1) = sqrt(se)*(1.D0 - (1.D0-1.D0/temp)**um)**2.d0
        else
          upc = pcmax
          se = 0.D0
          kr(1) = 0.D0
        endif

! Gas phase using BC

        se = (sw - swir)/(1.D0 - swir - sgir)

        if (sw < swir) then
          se = 0.D0
          kr(2) = 1.0
        elseif (se >= 1.D0-1d-6) then
          kr(2) = 0.D0
        else
          kr(2) = (1.D0 - se)**0.5D0 * (1.D0 - se**(1.D0/um))**(2.D0*um)
!         kr(2) = (1.D0 - se)**0.33333333D0 * (1.D0 - se**(1.D0/um))**(2.D0*um)
        endif

      case(VAN_GENUCHTEN_DOUGHTY) ! Doughty (2007) - Van Genutchen Mualem model adjusted for sgir/=0
        ! Modeling geologic storage of carbon dioxide: Comparison of non-hysteretic and
        ! hysteretic characteristic curves, Doughty (2007)
        ala = pckr_alpha
        um = pckr_m
        un = 1.D0/(1.D0 - um)
        if (sw > pckr_sat_water_cut) then
          upc = 0.D0; kr(1) = 1.d0; kr(2) = 0.d0;
        elseif (sw > (1.05D0*swir)) then
          select case(ikrtype)
            case(2) 
              ! Krl Van Genuchten-Mualem
              se = (sw - swir)/(1.D0 - swir)
              temp = se**(-1.D0/um)
              kr(1) = sqrt(se)*(1.D0 - (1.D0 - 1.D0/temp)**um)**2.d0
              ! Krg Corey - TOUGH2 
              se = (sw - swir)/(1.D0 - sgir - swir)
              kr(2) = (1.D0 - se)**2.0D0 * (1.D0 - se**(2.D0))  
          end select    
          !! capillary pressure 
          if (sw <= ( 0.99*(1 - sgir) )) then   
            se = (sw - swir)/(1.D0 - sgir - swir)
            temp = se**(-1.D0/um)
            upc = (temp - 1.D0)**(1.d0/un)/ala
          else ! for sw > (1 - sgir) linear interpolation to pc = 0
            se = ( 0.99*(1.0D0 - sgir) - swir)/(1.D0 - sgir - swir)
            temp = se**(-1.D0/um)
            upc = (temp - 1.D0)**(1.d0/un)/ala
            upc = upc*(1.D0 - sw)/(1.0D0 - 0.99*(1.0D0 - sgir))
          endif  
  
        else  ! use linear extropolation
          se0 = (0.05D0*swir)/(1.D0 - sgir - swir)
          temp = se0**(-1.D0/um)
          upc0 = (temp - 1.D0)**(1.d0/un)/ala
          upc_s0 = -1.D0/um/un*upc0*(se0**(-1.D0 - 1.D0/um))/(se0**(-1.D0/um) - 1.d0)
          upc_s0 = upc_s0 /(1.D0 - sgir - swir)
          if (sw > swir) then
            select case(ikrtype)
              case(2)
                ! Krl Van Genutchen-Mualem
                se = (sw - swir)/(1.D0 - swir)
                temp = se**(-1.D0/um)
                kr(1) = sqrt(se)*(1.D0 - (1.D0 - 1.D0/temp)**um)**2.d0
                ! Krg Corey - TOUGH2
                se = (sw - swir)/(1.D0 - sgir - swir) 
                kr(2) = (1.D0 - se)**2.0D0 * (1.D0 - se**(2.D0))  
            end select    
            !! pressure 
            upc = upc0 + (sw - 1.05D0 * swir) * upc_s0   
          else
            upc = upc0 + (sw - 1.05D0 * swir) * upc_s0
            kr(1) = 0.D0
            kr(2) = 1.D0
          end if
        end if
        !if (upc > pcmax) upc = pcmax
    end select

    ! scaling kr with end points
    kr(1) = kr(1) * Kr0(1)
    kr(2) = kr(2) * Kr0(2)

    pc(1) = upc; pc(2) = 0.d0;

  return

end subroutine pflow_pckr_noderiv_exec

! ************************************************************************** !

subroutine pckrNH_noderiv(sat, pc, kr, saturation_function, option)
  !
  ! pckrHY_noderiv: Hysteric S-Pc-kr relation driver
  ! 
  ! Author: Chuan Lu
  ! Date: 05/12/08
  ! 

  use Saturation_Function_module

  implicit none
  
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option
  PetscReal :: sat(option%nphase),pc(option%nphase),kr(option%nphase)

  PetscReal :: pckr_sir(option%nphase)
  PetscReal :: Kr0(option%nphase)
  PetscReal :: pckr_lambda, &
       pckr_alpha,pckr_m,pckr_pcmax,sg ,pckr_beta,pckr_pwr
  
  
  pckr_sir(:) = saturation_function%Sr(:)
  Kr0(:) = saturation_function%Kr0(:)
  pckr_m = saturation_function%m
  pckr_lambda = saturation_function%lambda
  pckr_alpha = saturation_function%alpha
  pckr_pcmax = saturation_function%pcwmax
  pckr_beta = saturation_function%betac
  pckr_pwr  = saturation_function%power
  
  sg = sat(2)
  
 ! call pflow_pckr_noderiv_exec(saturation_function%saturation_function_itype,&  
 !      pckr_sir,pckr_lambda, pckr_alpha,pckr_m,pckr_pcmax,sg,pc,kr,pckr_beta,pckr_pwr) 

  call pflow_pckr_noderiv_exec(saturation_function%saturation_function_itype,&
       saturation_function%permeability_function_itype, pckr_sir, Kr0, & 
       pckr_lambda, pckr_alpha,pckr_m,pckr_pcmax, sg,pc,kr,pckr_beta,pckr_pwr)
      
end subroutine pckrNH_noderiv

! ************************************************************************** !

subroutine pckrHY_noderiv(sat, hysdat, pc, kr, saturation_function, option)
  ! 
  ! Hysteric S-Pc-kr relation driver
  ! 
  ! Author: Chuan Lu
  ! Date: 05/12/08
  ! 

  use Saturation_Function_module
  
  implicit none
  
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option
  PetscReal :: sat(option%nphase),pc(option%nphase),kr(option%nphase)
  PetscReal :: hysdat(:)
  
  pc=0.D0
  kr=1.D0  

end subroutine pckrHY_noderiv



end module Mphase_pckr_module
