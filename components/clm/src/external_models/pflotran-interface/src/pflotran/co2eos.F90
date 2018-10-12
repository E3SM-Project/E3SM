module co2eos_module

  private

#include "petsc/finclude/petscsys.h"

  public HENRY_co2_noderiv,VISCO2,duanco2,denmix,Henry_duan_sun, &
         Henry_duan_sun_0NaCl,CO2

contains

! ************************************************************************** !

subroutine CO2(TX,PCX,DC,FC,PHI,HC)
  ! 
  ! VERSION/REVISION HISTORY
  ! $Id: co2eos.F90,v 1.1.1.1 2004/07/30 21:49:42 lichtner Exp $
  ! $Log: co2eos.F90,v $
  ! Revision 1.1.1.1  2004/07/30 21:49:42  lichtner
  ! initial import
  ! Revision 1.2  2004/01/10 18:32:06  lichtner
  ! Began work on 2 phase capability.
  ! Revision 1.1.1.1  2003/11/23 20:12:46  lichtner
  ! initial entry
  ! Revision 1.2  2003/05/09 15:22:41  lichtner
  ! commented out icall statements
  ! Revision 1.1.1.1  2003/03/03 01:33:27  lichtner
  ! PFLOTRAN initial implementation
  ! Revision 1.6  2002/09/28 17:25:49  lichtner
  ! Improved fit of dissolved CO2.
  ! Revision 1.5  2002/05/19 18:53:01  lichtner
  ! Added documentation of CO2 EOS
  ! Revision 1.4  2002/05/19 00:21:46  lichtner
  ! Modified Crovetto (1991) fit to Henry's law to be consistent with
  ! Duan et al. CO2 EOS.
  ! Revision 1.3  2002/05/07 03:14:39  lichtner
  ! Modified Henry's law subroutine to output the Poynting term.
  ! Revision 1.2  2002/05/04 18:28:34  lichtner
  ! Added Duan and Weare CO2 eos and new Henry's law.
  ! Revision 1.1  2002/04/12 19:03:10  lichtner
  ! Initial entry
  ! 
      use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT    

      implicit none
      
      PetscReal :: AT,T,TOL,TX,PCX,DC,DV,FC,PHI,H,HC,XMWC,R,B,V, &
                Y,Y2,Y3,Z
      
      PetscInt :: KOUNT
!
!     This subroutine calculates the specific density, fugacity, and 
!     specific enthalpy of gaseous and supercritical CO2 as a function of
!     the pressure of CO2 (PCX) and temperature (TX) using a Modified 
!     Redlich-Kwong (MRK) equation of state (EOS) and standard thermo-
!     dynamic equations.  This formulation of the MRK EOS is based on the 
!     work of Kerrick and Jacobs (1981) and Weir et al. (1996). Weir et al.
!     extended the MRK EOS of Kerrick and Jacobs to low temperatures.
!     Accuracy is suspect outside the temperature and pressure ranges of
!     50 < T < 350 deg C and 0.1 < PCX < 45 MPa, respectively.
!
!       Input:
!               TX   = Temperature in degrees C
!               PCX  = Pressure of CO2 in Pa
!
!       Output:
!               DC   = Specific density of CO2 in kg/m3
!               HC   = Specific enthalpy of CO2 in J/kg
!               FC   = Fugacity of CO2 in Pa
!              PHI   = Fugacity coefficient of CO2
!
!       Constants:
!               XMWC = Molecular weight of CO2 in Kg/mol
!               R    = Universal gas constant in m3Pa/molK
!               B    = Covolume in m3/mol (value from K&J)
!
!       Variables:
!               T    = Temperature in K
!               V    = Molar volume of CO2 in m3/mol    
!               Y    = Dimensionless variable (B/4V)
!               Z    = Compressibility factor [-]
!               PHI  = Fugacity coefficient [-] 
!               H    = Molar enthalpy in J/mol
!

      PARAMETER(XMWC = 4.40098D-02)
      PARAMETER(R    = IDEAL_GAS_CONSTANT)
      PARAMETER(B    = 5.8D-05)
!      
!
      if (PCX.LE..1D6)THEN
        DC=0.D0
        FC=0.D0
        HC=0.D0
        RETURN
      END IF
!    Convert temperature from degrees C to K:
      T = TX + 2.7315D+02
!
!    First calculate V as a function of T and PCX using Newton Iteration 
!    with tolerance TOL:
      TOL = 1.0D-06
!
!    Initial guess of V, DV, and Y from ideal gas law:
      V = (R*T) / PCX
      DV = V
      Y = B / (4.D0*V)
!
!    Initialize attractive term (AT) of MRK EOS:
      AT  = 0.D0
!
!    Newton Iteration for V as a function of T and PCX:
      KOUNT = 0
!--> use the following as a crude method to set up file incon (revise 
!--> appropriate boundary cells to have correct pressures); uncomment,
!--> run model till it stops, then use the SAVE file:
      DO WHILE(DABS(DV/V).GT.TOL)
        CALL MRK(Y,T,PCX,V,DV,AT)
        V = V - DV
        Y = B / (4.D0*V)
        KOUNT = KOUNT + 1
        if (KOUNT.GT.25000) GO TO 5
      END DO

!    Calculate density (DC) in kg/m3 from V in m3/mol:
      DC = XMWC / V

!    Calculate Y to the 2nd and 3rd powers for later use:
      Y2 = Y * Y
      Y3 = Y2 * Y
 
!    Calculate compressibility factor (Z) by substituting MRK EOS
!    into Z=PV/RT:
      Z = ((1.D0+Y+Y2-Y3)/((1.D0-Y)**3.D0)) - (AT/(R*T*DSQRT(T)*(V+B)))      

!    Initialize fugacity coefficient (PHI):
      PHI = 0.D0

!    Calculate fugacity (FC):
      CALL FUGACITY(Y,T,V,Z,PHI)
      FC = PHI * PCX

!    Initialize molar enthalpy (H):
      H = 0.D0

!    Calculate specific enthalpy (HC):
      CALL ENTHALPY(T,V,Z,H)
      HC = (H/XMWC)+8.1858447D+05
      RETURN

!    Come here when no convergence:
5     CONTINUE
      PRINT 6
6     FORMAT('NO CONVERGENCE IN SUBROUTINE CO2')
      print*, PCX,T,V,Y,DV
      
      RETURN
end subroutine CO2

! ************************************************************************** !

subroutine MRK(Y,T,PCX,V,DV,AT)
      
      use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT    

      implicit none

!     MODIFIED REDLICH-KWONG EQUATION OF STATE FOR CO2
!     This subroutine is called from subroutine CO2 during the Newton
!     Iteration for the molar volume (V) of CO2 as function of temperature
!     (T) and pressure of CO2 (PCX).  This subroutine calculates  
!     the V for which the MRK EOS is 0 at the given T and PCX, and the 
!     value of the derivative of the MRK EOS wrt V for the calculated V.
!
!       Input:
!               Y   = Dimensionless variable (B/4V)
!               T   = Temperature in K
!               PCX = Pressure of CO2 in Pa
!               V   = Prev. estimate of molar volume of CO2 in m3/mol   
!
!       Output:
!               DV  = Change in molar volume of CO2 in m3/mol 
!
!       Constants:
!               R   = Universal gas constant in m3Pa/molK
!               B   = Covolume in m3/mol (value from K&J)
!       Ci thru Fi= Coefficients of the MRK EOS (i=1,2,3)   
!                   Values from Weir et al. (1996)
!
!       Variables:
!       CT thru FT= Temperature-dependent functions for evaluating
!                     attractive term of MRK EOS
!               AT  = Attractive term of MRK EOS        
!               FV  = V at which MRK EOS is 0 for T and PCX
!               DV  = -FV / Value of derivative wrt V of MRK EOS

      PetscReal :: Y,T,PCX,V,DV,AT
      PetscReal :: R,B,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3
      PetscReal :: T2,V2,V3,V4,Y2,Y3
      PetscReal :: B2,B3,CT,DT,ET,FT,FV

      PARAMETER(R   =  IDEAL_GAS_CONSTANT)
      PARAMETER(B   =  5.8D-05)
      PARAMETER(C1  =  2.39534D+01)
      PARAMETER(C2  = -4.55309D-02)
      PARAMETER(C3  =  3.65168D-05)
      PARAMETER(D1  = -4.09844D-03)
      PARAMETER(D2  =  1.23158D-05)
      PARAMETER(D3  = -8.99791D-09)
      PARAMETER(E1  =  2.89224D-07)
      PARAMETER(E2  = -8.02594D-10)
      PARAMETER(E3  =  7.30975D-13)
      PARAMETER(F1  = -6.43556D-12)
      PARAMETER(F2  =  2.01284D-14)
      PARAMETER(F3  = -2.17304D-17)

!    Calculate T squared for later use:
      T2 = T * T
!
!    Calculate V to the 2nd, 3rd, and 4th powers for later use:
      V2 = V * V
      V3 = V2 * V
      V4 = V3 * V
!
!    Calculate Y to the 2nd and 3rd powers for later use:
      Y2 = Y * Y
      Y3 = Y2 * Y
!
!    Calculate B to the 2nd and 3rd powers for later use:
      B2 = B * B
      B3 = B2 * B
!
!    Calculate temperature-dependent functions for evaluating attractive
!    term in MRK EOS:
      CT = C1 + (C2*T) + (C3*T2)
      DT = D1 + (D2*T) + (D3*T2)
      ET = E1 + (E2*T) + (E3*T2)
      FT = F1 + (F2*T) + (F3*T2)
!
!    Calculate attractive term in MRK EOS: 
      AT = CT + (DT/V) + (ET/V2) + (FT/V3)
!
!    Calculate V at which MRK EOS equals 0:
      FV = PCX - (((R*T*(1.D0+Y+Y2-Y3))/(V*((1.D0-Y)**3.D0)))- &
           (AT/(DSQRT(T)*V*(V+B))))     
!
!    Calculate -FV / value of derivative wrt V of MRK EOS 
      DV = -FV / (((-3.D0*B*R*T*(1.D0+Y+Y2-Y3))/(4.D0*V3* &
           ((1.D0-Y)**4.D0)))-((R*T*(1.D0+Y+Y2-Y3))/ &
           (V2*((1.D0-Y)**3.D0)))+((R*T*(((3.D0*B3)/(64.D0*V4)) &
           -(B2/(8.D0*V3))-(B/(4.D0*V2))))/(V*((1.D0-Y)**3.D0))) &
           -(AT/(DSQRT(T)*V*(V+B))))
      RETURN
end subroutine MRK

! ************************************************************************** !

subroutine FUGACITY(Y,T,V,Z,PHI)
      
      use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT    

      implicit none

!
!     This subroutine is called from subroutine CO2 during the       
!     calculation of fugacity of CO2 as function of temperature (T),
!     pressure of CO2 (PCX), and molar volume of CO2 (V).  This
!     subroutine calculates the fugacity coefficient of CO2 (PHI) by   
!     substituting the MRK EOS into RTln(PHI)=Integral from V to infinity
!     of (PCX-RT/V)dV - RTln(Z) + RT(Z-1).  This expression comes from
!     Prausnitz (1969).
!
!       Input:
!               Y   = Dimensionless variable (B/4V)
!               T   = Temperature in K
!               V   = Molar volume of CO2 in m3/mol
!               Z   = Compressibility factor of CO2 [-] 
!
!       Output:
!               PHI = Fugacity coefficient of CO2 [-]
!
!       Constants:
!               R   = Universal gas constant in m3Pa/molK
!               B   = Covolume in m3/mol (value from K&J)
!       Ci thru Fi= Coefficients of the MRK EOS (i=1,2,3)   
!                     Values from Weir et al. (1996)
!
!       Variables:
!       CT thru FT= Temperature-dependent functions for evaluating
!                     attractive term of MRK EOS
!
      PetscReal :: Y,T,V,Z,PHI
      PetscReal :: R,B,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3
      PetscReal :: T2,V2,V3,B2,B3,B4,CT,DT,ET,FT
      
      PARAMETER(R   =  IDEAL_GAS_CONSTANT)
      PARAMETER(B   =  5.8D-05)
      PARAMETER(C1  =  2.39534D+01)
      PARAMETER(C2  = -4.55309D-02)
      PARAMETER(C3  =  3.65168D-05)
      PARAMETER(D1  = -4.09844D-03)
      PARAMETER(D2  =  1.23158D-05)
      PARAMETER(D3  = -8.99791D-09)
      PARAMETER(E1  =  2.89224D-07)
      PARAMETER(E2  = -8.02594D-10)
      PARAMETER(E3  =  7.30975D-13)
      PARAMETER(F1  = -6.43556D-12)
      PARAMETER(F2  =  2.01284D-14)
      PARAMETER(F3  = -2.17304D-17)

!
!    Calculate T to the 2nd power for later use:
      T2 = T * T
!
!    Calculate V to the 2nd, 3rd, and 4th powers for later use:
      V2 = V * V
      V3 = V2 * V
!
!    Calculate B to the 2nd, 3rd, and 4th powers for later use:
      B2 = B * B
      B3 = B2 * B
      B4 = B3 * B
!
!    Calculate temperature dependent functions for evaluating attractive
!    term in MRK EOS:
      CT = C1 + (C2*T) + (C3*T2)
      DT = D1 + (D2*T) + (D3*T2)
      ET = E1 + (E2*T) + (E3*T2)
      FT = F1 + (F2*T) + (F3*T2)
!      
!    Calculate fugacity coefficient:
      PHI = Y * (8.D0 + Y * (-9.D0 + 3.D0 * Y))/(1.D0-Y)**3.D0 &
           - DLOG(Z) &
           - CT / (R * T * DSQRT(T) * (V + B)) &
           - DT / (R * T * DSQRT(T) * V * (V + B)) &
           - ET / (R * T * DSQRT(T) * V2 * (V + B)) &
           - FT / (R * T * DSQRT(T) * V3 * (V + B)) &
           + CT * DLOG(V / (V + B))/ (R * T * DSQRT(T) * B) &
           - DT / (R * T * DSQRT(T) * B * V) &
           + DT * DLOG((V + B) / V) / (R * T * DSQRT(T) * B2) &
           - ET / (R * T * DSQRT(T) * 2.D0 * B * V2) &
           + ET / (R * T * DSQRT(T) * B2 * V) &
           - ET * DLOG((V + B) / V) / (R * T * DSQRT(T) * B3) &
           - FT / (R * T * DSQRT(T) * 3.D0 * B * V3) &
           + FT / (R * T * DSQRT(T) * 2.D0 * B2 * V2) &
           - FT / (R * T * DSQRT(T) * B3 * V) &
           - FT * DLOG(V / (V + B)) / (R * T * DSQRT(T) * B4)
      PHI = DEXP(PHI)

      RETURN
end subroutine FUGACITY

! ************************************************************************** !

subroutine ENTHALPY(T,V,Z,H)

      use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT    

      implicit none
!
!     This subroutine is called from subroutine CO2 during the       
!     calculation of the specific enthalpy of CO2 as function of
!     temperature (T), pressure of CO2 (PCX), and molar volume 
!     of CO2 (V).  This subroutine calculates the molar enthalpy of CO2
!     using residual properties.  A residual property is defined as the
!     difference between the real fluid property and the perfect gas   
!     state property.  Following Patel and Eubank (1988) for molar enthalpy:
!     H-H'ref = H(T,rho) - H'(Tref,Pref/RTref), where ' indicates the 
!     perfect gas state.  Integration is done along the path
!     H(T,rho)-->H'(T,0)-->H'(Tref,0)-->H'(Tref,Pref/RTref).
!
!     Determine residual internal energy first: (U-U'ref)/RT = 1/T integral
!     from 0 to rho of dZ/d(1/T) drho/rho + 1/T integral from Tref to T of
!     Cv/R dT, where Cv is molar heat capacity in J/(mol K).  Then determine
!     residual enthalpy: (H-H'ref)/RT = (U-U'ref)/RT + Z - Tref/T.  Using
!     Tref=273.16 K and Pref=1000 Pa, H'ref=0 (from Patel and Eubank).
!
!       Input:
!               T    = Temperature in K
!               V    = Molar volume of CO2 in m3/mol
!               Z    = Compressibility factor of CO2 [-] 
!
!       Output:
!               H    = Molar enthalpy of CO2 in J/mol
!
!       Constants:
!               R    = Universal gas constant in m3Pa/molK
!               B    = Covolume in m3/mol (value from K&J)
!       Ci thru Fi = Coefficients of the MRK EOS (i=1,2,3)   
!                      Values from Weir et al. (1996)
!           Gi   = Coefficients of molar heat capacity
!                  Values from Angus et al. (1976)
!               Tref = Reference temperature in K (value from P&E)
!
!       Variables:
!               RHO  = Molar density of CO2 in mol/m3
!               XI1  = First Integral (see above)
!               XI2  = Second Integral (see above)
!               URES = Residual internal energy
!
      PetscReal :: T,V,Z,H
      PetscReal :: BETA,R,TREF,B,C1,C2,C3,D1,D2,D3,E1,E2,E3,F1,F2,F3, &
                   G0,G1,G2,G3,G4,G5,G6,G7, &
                   RHO,RHO2,BETA2,BETA3,BETA4,BETA5,BETA6,BETA7, &
                   TREF2,TREF3,TREF4,TREF5,TREF6,T2,T3,T4,T5,T6, &
                   B2,B3,B4,XI1,XI2,URES
      
      PARAMETER(BETA=  304.21D0)
      PARAMETER(R   =  IDEAL_GAS_CONSTANT)
      PARAMETER(TREF=  2.7316D+02)
      PARAMETER(B   =  5.8D-05)
      PARAMETER(C1  =  2.39534D+01)
      PARAMETER(C2  = -4.55309D-02)
      PARAMETER(C3  =  3.65168D-05)
      PARAMETER(D1  = -4.09844D-03)
      PARAMETER(D2  =  1.23158D-05)
      PARAMETER(D3  = -8.99791D-09)
      PARAMETER(E1  =  2.89224D-07)
      PARAMETER(E2  = -8.02594D-10)
      PARAMETER(E3  =  7.30975D-13)
      PARAMETER(F1  = -6.43556D-12)
      PARAMETER(F2  =  2.01284D-14)
      PARAMETER(F3  = -2.17304D-17)
      PARAMETER(G0  =  0.769441246D+01)
      PARAMETER(G1  = -0.249610766D+00)
      PARAMETER(G2  = -0.254000397D+02)
      PARAMETER(G3  =  0.651102201D+02)
      PARAMETER(G4  = -0.820863624D+02)
      PARAMETER(G5  =  0.574148450D+02)
      PARAMETER(G6  = -0.212184243D+02)
      PARAMETER(G7  =  0.323362153D+01)

!     SAVE ICALL
!     DATA ICALL/0/
!     ICALL=ICALL+1
!     if (ICALL.EQ.1) WRITE(11,899)
! 899 FORMAT(6X,'ENTHALPY 2.0      29 JUNE      1999',6X,
!    X'CALCULATE MOLAR ENTHALPY OF SEPARATE PHASE CO2')
!
!    Calculate molar density (RHO):
      RHO = 1.D0 / V
!
!    Calculate rho to the 2nd and 3rd powers for later use:
      RHO2 = RHO * RHO
!
!    Calculate beta to the 2nd thru 7th powers for later use:
      BETA2 = BETA * BETA
      BETA3 = BETA2 * BETA
      BETA4 = BETA3 * BETA
      BETA5 = BETA4 * BETA
      BETA6 = BETA5 * BETA
      BETA7 = BETA6 * BETA
!
!    Calculate Tref to the 2nd thru 6th powers for later use:
      TREF2 = TREF * TREF
      TREF3 = TREF2 * TREF
      TREF4 = TREF3 * TREF
      TREF5 = TREF4 * TREF
      TREF6 = TREF5 * TREF
!
!    Calculate T to the 2nd thru 6th powers for later use:
      T2 = T * T
      T3 = T2 * T
      T4 = T3 * T
      T5 = T4 * T
      T6 = T5 * T
!
!    Calculate B to the 2nd, 3rd, and 4th powers for later use:
      B2 = B * B
      B3 = B2 * B
      B4 = B3 * B
!
!    Calculate 1/T times the integral from 0 to rho of dZ/d(1/T) drho/rho:
      XI1 = (B*RHO*(-6.D0*(3.D0*F1+T*(F2-F3*T))-B2*(18.D0*D1 &
            +6.D0*T*(D2-D3*T)+3.D0*(3.D0*E1+T*(E2-E3*T))*RHO &
            +2.D0*(3.D0*F1+T*(F2-F3*T))*RHO2)+3.D0*B*(6.D0*E1 &
            +3.D0*F1*RHO+T*(2.D0*E2-2.D0*E3*T+F2*RHO-F3*T*RHO))) &
            +6.D0*(3.D0*F1+T*(F2-F3*T)+B3*(-3.D0*C1+T*(-C2+C3*T)) &
            +B2*(3.D0*D1+T*(D2-D3*T))+B*(-3.D0*E1+T*(-E2+E3*T))) &
            *DLOG(1+B*RHO))/(12*B4*R*T**1.5)
! 
!    Calculate 1/T times the integral from Tref to T of Cv/R dT, where Cv
!    is molar heat capacity in J/(mol K).  The expression for Cv is
!    derived an expression from Angus et al. (1976) for molar
!    heat capacity at constant pressure: 
      XI2 = G0-1.D0+((TREF/T)*(1.D0-G0))+(((BETA*G1)/T)*DLOG &
            (T/TREF))+(((BETA2*G2)/T)*((1.D0/TREF)-(1.D0/T))) &
            +(((BETA3*G3)/(2.D0*T))*((1.D0/TREF2)-(1.D0/T2))) &
            +(((BETA4*G4)/(3.D0*T))*((1.D0/TREF3)-(1.D0/T3))) &
            +(((BETA5*G5)/(4.D0*T))*((1.D0/TREF4)-(1.D0/T4))) &
            +(((BETA6*G6)/(5.D0*T))*((1.D0/TREF5)-(1.D0/T5))) &
            +(((BETA7*G7)/(6.D0*T))*((1.D0/TREF6)-(1.D0/T6)))
!
!    Calculate residual internal energy (URES):
      URES = XI1+XI2
!
!    Calculate molar enthalpy (H):
      H = (URES+Z-(TREF/T)) * R * T + 8.1858447D5*.044D0
      RETURN
end subroutine ENTHALPY

! ************************************************************************** !

subroutine duanco2 (tt,p,dc,fc,phi)
  ! 
  ! Subroutine: duanco2.f
  ! Input: tt   [C]        temperature
  ! p    [Pa]       CO2 partial pressure
  ! Output: fc   [bars]    CO2 fugacity
  ! phi  [-]       CO2 fugacity coefficient
  ! dc   [g/cm^3]  CO2 density
  ! Duan, Z., Moller, N., and Weare, J.H. (1992) An equation of state for
  ! the CH4-CO2-H2O system: I. Pure systems from 0 to 1000 oC and 0 to
  ! 8000 bar. Geochimica Cosmochimica Acta, 56, 2605-2617.
  ! 

      use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT    

      implicit none
      
      PetscReal :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,alpha,beta,gamma
      PetscReal :: t,tt,tc,pc,tr,tr1,tr2,tr3,pr,p,a,b,c,d,e,rgas,vc
      PetscReal :: epsilon,videal,v,v1,xmwc,vr,vr1,vr2,gamvr,expgam,z,zeq
      PetscReal :: f1,ov,df1,phi,fc,dc
      
      PetscInt :: itmax,iter
      
      data xmwc /4.40098d-02/

      data a1    / 8.99288497d-2/
      data a2    /-4.94783127d-1/
      data a3    / 4.77922245d-2/
      data a4    / 1.03808883d-2/
      data a5    /-2.82516861d-2/
      data a6    / 9.49887563d-2/
      data a7    / 5.20600880d-4/
      data a8    /-2.93540971d-4/
      data a9    /-1.77265112d-3/
      data a10   /-2.51101973d-5/
      data a11   / 8.93353441d-5/
      data a12   / 7.88998563d-5/
      data alpha /-1.66727022d-2/
      data beta  / 1.398d0/
      data gamma / 2.96d-2/

      t = tt + 273.15d0
      

      tc = 31.05d0 + 273.15d0
      pc = 73.825d0

      tr = t/tc
      tr1 = 1.d0/tr
      tr2 = tr1*tr1
      tr3 = tr1*tr2
      pr = p /pc
      a = alpha*tr3
      b = a1 + (a2 + a3*tr1)*tr2
      c = a4 + (a5 + a6*tr1)*tr2
      d = a7 + (a8 + a9*tr1)*tr2
      e = a10 + (a11 + a12*tr1)*tr2

      rgas = IDEAL_GAS_CONSTANT * 1.d-2 !bar dm**3 / k
      vc = rgas*tc/pc

!-----solve
!     f(v) = z[p, t, v] - zeq[t, v] = 0

!-----newton-raphson
!     f(x) = 0
!     f'(x1) (x-x1) = -f(x1)
!     x = x1-f/f'

!-----tolerance parameters
      epsilon = 1.d-8
      itmax = 150

!-----set initial guess from ideal gas law
      videal = rgas*t/p
      v = videal ! ideal gas law

      iter = 0
   10 continue
      iter = iter + 1
      if (iter .gt. itmax) goto 20

      vr = v/vc
      vr1 = 1.d0/vr
      vr2 = vr1*vr1
      gamvr = gamma*vr2
      expgam = exp(-gamvr)

      z = pr*vr/tr

      zeq = 1.d0 + (b + (c + (d + e*vr1)*vr2)*vr1)*vr1 + &
            a*vr2*(beta + gamvr)*expgam

      f1 = z-zeq

!-----check convergence
      if (dabs(f1).lt.epsilon .and. iter.gt.2) goto 30

      ov = 1.d0/v
      df1 = pr/vc/tr + (b + (2.d0*c + (4.d0*d + 5.d0*e*vr1)*vr2 &
           + 2.d0*a*((beta + gamvr)*(1.d0-gamvr) + gamvr) &
           * expgam)*vr1)*ov*vr1

      v1 = v
      v = v1-f1/df1
      
!     Protect against negative volume solutions
      if (v < 0) v = v1 /2.d0

!     write(*,*) 'v,v1,f,df: ',iter,v,v1,f1,df1

      goto 10

   20 continue

      write(*,*) 'duanco2: iter > itmax:  stop', tt,p
      stop

   30 continue

!-----fugacity coefficient
      phi = exp(z-1.d0-log(z) + (b + (0.5d0*c + (0.25d0*d + &
      0.2d0*e*vr1)*vr2)*vr1)*vr1 + a/(2.d0*gamma)* &
      (beta + 1.d0 - (beta + 1.d0 + gamvr)*expgam))

!-----fugacity
      fc = phi * p

!-----density
      dc = xmwc / v ! g/cm^3

      return

end subroutine duanco2

! ************************************************************************** !

subroutine Henry_duan_sun_0NaCl (p,tc,henry)

  implicit none

  PetscReal :: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
  PetscReal :: p,t,tc,lnt,muco2, henry

  c1  = 28.9447706d0
  c2  = -0.0354581768d0
  c3  = -4770.67077d0
  c4  =  1.02782768d-5
  c5  = 33.8126098d0
  c6  =  9.04037140d-3
  c7  = -1.14934031d-3
  c8  = -0.307405726d0
  c9  = -0.0907301486d0
  c10 =  9.32713393d-4

  t = tc + 273.15d0
  lnt = log(t)
  muco2 = c1 + c2*t + c3/t + c4*t*t + c5/(630.d0-t) + c6*p + c7*p*lnt &
    + c8*p/t + c9*p/(630.d0-t) + c10*p*p/(630.d0-t)**2

  henry = exp(-muco2)

  return
end subroutine Henry_duan_sun_0NaCl

! ************************************************************************** !

subroutine Henry_duan_sun(tc,p,keqco2,lngamco2,mc,ma,co2_aq_actcoef)
  
! t[c], p[bar], mco2[mol/Kg-H2O], mc[cation: mol/kg-H2O], 
! ma[anion: mol/kg-H2O], psat[bars]

  implicit none
  PetscReal, save :: coef(3,11)
  PetscReal :: tc,p,keqco2,mc,ma,t
  PetscReal, optional :: co2_aq_actcoef
  PetscReal :: temparray(11)

  PetscReal :: lngamco2, tmp, mu0, lamc, lamca

  data coef / 28.9447706, -0.411370585, 3.36389723e-4, &
  -0.0354581768,  6.07632013e-4,  -1.98298980e-5, &
  -4770.67077,    97.5347708,     0., &
  1.02782768e-5,  0.,             0., &
  33.8126098,     0.,             0., &
  9.04037140e-3,  0.,             0., &
  -1.14934031e-3, 0.,             0., &
  -0.307405726,  -0.0237622469,   2.12220830e-3, & 
  -0.0907301486,  0.0170656236,  -5.24873303e-3, & 
  9.32713393e-4,  0.,             0., &
  0.,             1.41335834e-5,  0./


  t=tc+273.15

  ! adding temparray to improve efficiency (and remove Intel warning) - geh
  temparray = coef(1,:)
  call duan_sun_param(t,p,temparray,mu0) ! mu0/RT
  temparray = coef(2,:)
  call duan_sun_param(t,p,temparray,lamc) ! lambda_CO2-Na Pitzer 2nd order int. param.
  temparray = coef(3,:)
  call duan_sun_param(t,p,temparray,lamca) ! zeta_CO2-Na-Cl Pitzer 3rd order int. param.
  
  !activity coef. aqueous co2
  lngamco2 = 2.d0*lamc*mc + lamca*mc*ma ! = log(gam(jco2))
  if (present(co2_aq_actcoef)) then
    co2_aq_actcoef = exp(lngamco2)
  endif 
  tmp = mu0 + lngamco2 !- ln(phico2) [ln y P/m = mu0 + ln gamma - ln phico2]
  
  ! yco2 = (p-psat)/p ! mole fraction of CO2 in supercritical CO2 phase: based on
                    ! assumption that mole fraction of H2O in SC CO2 = Psat(T)/P.
  
  ! mco2 = molality co2
  ! mco2 = phico2 * yco2 * p * exp(-mu0) / gamma
  
  keqco2 = exp(-tmp) ! = K_co2 / gamco2
  !print *, 'keqco2: ', mu0,lngamco2,tc,p,ma,mc
  return
end subroutine Henry_duan_sun

! ************************************************************************** !

subroutine duan_sun_param(t,p,c,par)
  
  implicit none
  
  PetscReal :: t,p,par,fac,c(11)
  
  fac = 1.d0/(630.d0-t)
  
  par = c(1) + c(2) * t + c(3) / t + c(4) * t * t &
  + c(5) * fac + c(6) * p + c(7) * p * log(t) &
  + c(8) * p / t + c(9) * p * fac &
  + C(10) * p * p * fac * fac + c(11) * t * log(p)

  return
end subroutine duan_sun_param

! ************************************************************************** !

subroutine HENRY_co2_noderiv(xmole,x1m,tx,pcx,xphi,rkh,poyn)
  ! 
  ! Subroutine: henry.f
  ! Input: tx   [C]       temperature
  ! pco2 [Pa]      CO2 partial pressure
  ! psys [Pa]      total system pressure ***not implemented***
  ! phi  [-]       CO2 fugacity coefficient
  ! Output: xmole [-]     mole fraction CO2 in aqueous phase
  ! x1m   [-]     mass fraction CO2
  ! rkh   [Pa]  Henry constant
  ! poyn  [-]     Poynting factor
  ! Crovetto, R. (1991) Evaluation of solubility data of the system CO2-H2O
  ! from 273�Z K to the critical point of water. Journal of Physical and
  ! Chemical Reference Data, 20(3), 575-589.
  ! 

!     input:
!     tx   [C]  temperature
!     pcx  [Pa] CO2 partial pressure
!     xphi [-]  fugacity coef.

!     output:
!     rkh   [Pa] Henry's constant
!     poyn  [-]  Poynting correction
!     xmole [-]  mole fraction CO2 in liquid water at t and p
!     x1m   [-]  mass fraction CO2
      
      implicit none
      
      PetscReal, intent(in) :: tx,xphi,pcx
      PetscReal, intent(out) :: xmole,x1m,rkh,poyn
      PetscReal :: ams,ama,tk,otk

      ams = 18.01534D-3 !h2o
      ama = 44.0098D-3  !co2

!     pcx = pco2! *1.d5 !suspious, pco2 should be in bar as well

!     COMMON/GASLAW/R,AMS,AMA,CVNCG
!
!     SAVE ICALL
!     DATA ICALL/0/
!     ICALL=ICALL+1
!     if (ICALL.EQ.1) WRITE(11,899)
! 899 FORMAT(6X,'HENRY    1.0      25 FEBRUARY  1993',6X,
! 899 FORMAT(6X,'HENRY    1.1      30 October   1997',6X,
!    X'Henry',1H','s law for dissolved CO2 mass',
!    X' fraction as function of partial pressure'/
!    x47X,'adapted from Battistelli et al. formulation in EWASG')
!
!     T2=TX*TX
!     T3=T2*TX
!     T4=T2*T2
!
!     CO2: CRAMER, 1982.
!        RKH=(783.666E0+19.6025E0*TX+0.820574E0*T2-7.40674E-03*T3+
!    X       2.18380E-05*T4-2.20999E-08*T2*T3)*1.E5
!
      tk = tx + 273.15d0
      otk = 1.d0/tk

      if (tx.gt.372.d0) then
        rkH = (-272998.463142d0*otk + 568.443483d0)*otk - &
        0.00625205d0*tk + 1.099025d0*Log(tk)
        write(*,*) "Henry's law out of range: t= ",tx,rkH
      else

!-------crovetto form
!       rkH = 2.45621d0 + (948.4254679d0-458199.88483d0*otk)*otk + (1008.599d0*
!    .  (1.d0 - 0.001545294d0*tk)**0.333333333d0)*otk

!-------new fit 1 - fixed Vco2 = 35 cm^3/mol
!       rkH = 9.38406d0 - (6192.46 + 132377.d0*otk)*otk + &
!       (5885.39d0*(1.d0 - 0.001545294d0*tk)**0.333333333d0)*otk
        
!       rkH = 9.38406d0 - (6192.46d0 + 132377.d0*otk - &
!       5885.39d0*(1.d0 - 0.001545294d0*tk)**0.333333333d0)*otk


        rkH = (-272998.463142d0*otk + 568.443483d0)*otk - 0.00625205d0*tk + 1.099025d0*Log(tk)
      endif

      rkH = 10.d0**rkH*1.d5/xphi ! units of rkH in [Pa]

      xmole = PCX/RKH ! note: pcx in [Pa]

!-----Poynting term
!     vid = 3.5d-5
!     rgas = IDEAL_GAS_CONSTANT
!     call sat(tx,ps)
 
!     poyn = exp(-(VID*(PCX-PS))/(Rgas*Tk))
      poyn  = 1.d0
      
!     xmole = xmole * poyn

      x1m = xmole*ama/(xmole*ama+(1.d0-xmole)*ams)
      
!     print *,'henry_co2_noderiv: ',pcx,tx,xmole,x1m,xphi,rkh

      RETURN
end subroutine HENRY_CO2_NODERIV

! ************************************************************************** !

subroutine HENRY_sullivan (TX,PCX,PS,FC,X1M,XCO2,HP)

      use PFLOTRAN_Constants_module, only : IDEAL_GAS_CONSTANT    

      implicit none

!     This subroutine calculates the mass fraction of CO2 in the liquid 
!     phase using an extended Henry's Law relationship from Reid et al.
!     (1987). The relationship is ln(FC/XCO2) = ln(HP) + VID(PCX-PS)/RT. 
!     See below for variable definitions.
!
!     The expression for Henry's Constant is from O'Sullivan et al. (1985).
!     The expression was created using a piece-wise quadratic fit to data
!     published by Ellis and Goulding (1963), Malinin (1959), Takenouchi 
!     and Kennedy (1964), and Gibb and Van Ness (1971). 
!
!     The value for the the partial molar volume of CO2 at infinite 
!     dilution is assumed to be constant at 30E-6 from the work of 
!     Takenouchi and Kennedy (1964) (and others). This assumption is 
!     reasonable at temperatures below 150 C.
!
!       Input:
!               TX   = Temperature in degrees C
!               PCX  = Pressure of CO2 in Pa
!               PS   = Saturation pressure of water in Pa
!               FC   = Fugacity of CO2 in Pa
!
!       Output:
!               X1M  = Mass fraction of CO2 in liquid phase [-]
!                HP  = Henry's contant [Pa] <pcl>
!
!       Constants:
!               XMWC = Molecular weight of CO2 in Kg/mol
!               XMWW = Molecular weight of H2O in Kg/mol
!               R    = Universal gas constant in m3Pa/molK
!               VID  = Partial molar volume of CO2 at infinite
!                      dilution in m3/mol (value from T&K)
!
!       Variables:
!               T    = Temperature in K
!               TAU  = Temperature variable used in calculation
!                      of Henry's Coefficient in degrees C
!               HP   = Henry's Coefficient in bars, then Pa
!               XCO2 = Mole fraction CO2 in liquid phase [-]

      PetscReal :: TX,PCX,PS,FC,X1M,XCO2,HP
      PetscReal :: XMWC,XMWW,R,VID,TAU,TAU2,TAU4,T
      
      PARAMETER(XMWC = 4.40098D-02)
      PARAMETER(XMWW = 1.801534D-02)
      PARAMETER(R    = IDEAL_GAS_CONSTANT)
      PARAMETER(VID  = 3.0D-05)

!    Calculate TAU:
      TAU = (TX-1.7D+02) / 1.0D+02

!    Calculate TAU to the 2nd and 4th powers for later use:
      TAU2 = TAU * TAU
      TAU4 = TAU2 * TAU2

!    Calculate Henry's Coefficient (HP) in bars:
      if (TAU.GE.0.D0)THEN
         HP = 6.4D+03 - (2.07778D+03*TAU2) + (3.1111D+02*TAU4)
      ELSEif (TAU.LT.0.D0)THEN
         HP = 6.4D+03 - (2.14914D+03*TAU2) - (1.9543D+02*TAU4)
      ENDIF

!    Convert Henry's Coefficient to Pa:
      HP = HP * 1.0D+05

!    Convert temperature to K:
      T = TX + 2.7315D+02

!    Calculate mole fraction of CO2 (XMOLE):
      XCO2 = DEXP(DLOG(FC/HP)-(VID*(PCX-PS))/(R*T)) 

!    Calculate mass fraction of CO2 (XMASS):
      X1M = (XMWC*XCO2) / (((1.D0-XCO2)*XMWW)+(XCO2*XMWC))

      RETURN
end subroutine HENRY_sullivan

! ************************************************************************** !

subroutine SOLUT(PCX,TX,HSOL)
      implicit none
      
      PetscReal :: PCX,TX,HSOL
      PetscReal :: T,T2,T3,T4

!     This subroutine calculates the enthalpy of CO2 dissolution in
!     liquid water. The expression is from O'Sullivan et al. (1985).
!     The expression was created using a quadratic fit to data published
!     by Ellis and Goulding (1963).

      T  = 1.D-2 * TX
      T2 = T * T
      T3 = T * T2
      T4 = T * T3
      HSOL = -7.3696D-2 - 5.6405D-1*T + 7.0363D-1*T2 - &
             2.7882D-1*T3 + 4.2579D-2*T4
      HSOL = HSOL * 1.D6
      RETURN
end subroutine SOLUT

! ************************************************************************** !

subroutine DENMIX(TX,DW,X1M,D1M)
      implicit none
      
      PetscReal :: TX,DW,X1M,D1M
      PetscReal :: TX2,TX3,TX4,RHO,DC,XMWC,X2M

!     This subroutine returns density of CO2-H2O liquid mixture. The
!     expression is from Anderson et al. (1992).
!
!       Input:
!               TX   = Temperature in degrees C
!               DW   = Density of H2O in kg/m3
!               X1M  = Mass fraction of CO2 [-]
!
!       Output:
!               D1M  = Density of CO2-H2O mixture in kg/m3
!
!       Constants:
!               XMWC = Molecular weight of CO2 in kg/mol
!
!       Variables:
!               RHO  = Density of CO2 at saturation pressure in mol/cm3
!               DC   = Density of CO2 at saturation pressure in kg/m3
!               X2M  = Mass fraction H2O [-]

      PARAMETER(XMWC=4.40098D-02)

      if (X1M.LE.0.D0) THEN
        D1M = DW
        RETURN
      ENDIF

!    Calculate TX to the 2nd, 3rd and 4th powers for later use:
      TX2 = TX * TX
      TX3 = TX2 * TX
      TX4 = TX3 * TX

!    Calculate density of CO2 (RHO) at saturation pressure in mol/cm3:
      RHO = 1.D0/(3.736D+01 - (7.109D-02*TX) - (3.812D-05*TX2) + &
            (3.296D-06*TX3) - (3.702D-09*TX4))

!    Convert RHO to kg/m3 (DC):
      DC = RHO * 1.0D+06 * XMWC

!    Calculate mass fraction of H2O:
      X2M = 1.D0 - X1M

!    Calculate density of CO2-H2O mixture in kg/m3:
      D1M = (DW*DC) / (X1M*DW + X2M*DC)

      RETURN
end subroutine DENMIX

! ************************************************************************** !

subroutine VISCO2(TX,DC,VC)
      implicit none
!
!     This subroutine calculates the viscosity of pure CO2 as a function 
!     of temperature and density of CO2.  The expressions for calculating
!     the viscosity come from empirical equations provided in Vesovic et
!     al.(1990) and Fenghour et al. (1998).
!     The critical point enhancement for the viscosity of CO2
!     has been neglected since it is weak and restricted to a very small
!     region around the critical point.
!
!       Input:
!               TX     = Temperature in degrees C
!               DC     = Density of CO2 in kg/m3
!
!       Output:
!               VC     = Viscosity of CO2 in Pa-s       
!
!       Constants:
!               Ai     = Coefficients of the correlation of the
!                        zero-density viscosity 
!               ESCL   = Energy scaling parameter in K
!                      = epsilon/kappa
!               Dij    = Coefficients of the correlation of the
!                        excess viscosity 
!
!       Variables:
!               T      = Temperature in K
!               TSTAR  = (kappa*T)/epsilon = T/ESCL [-]
!               ETA0   = Zero-density viscosity in muPa-s
!               DETA   = Excess viscosity in muPa-s
!
      PetscReal :: TX,DC,VC
      PetscReal :: A0,A1,A2,A3,A4,ESCL,D11,D21,D64,D81,D82
      PetscReal :: T,DC2,DC6,DC8,TSTAR,TSTAR3,BETA1,BETA2,BETA3,BETA4
      PetscReal :: EXS,ETA0,DETA
      
      PARAMETER(A0   =  2.35156D-01)
      PARAMETER(A1   = -4.91266D-01)
      PARAMETER(A2   =  5.211155D-02)
      PARAMETER(A3   =  5.347906D-02)
      PARAMETER(A4   = -1.537102D-02)
      PARAMETER(ESCL =  2.51196D+02)
      PARAMETER(D11  =  0.4071119D-02)
      PARAMETER(D21  =  0.7198037D-04)
      PARAMETER(D64  =  0.2411697D-16)
      PARAMETER(D81  =  0.2971072D-22)
      PARAMETER(D82  = -0.1627888D-22)

!    Convert temperature from degrees C to K:
      T = TX + 2.7315D+02

!    Calculate DC to 2nd, 6th, and 8th powers:
      DC2 = DC*DC
      DC6 = DC2*DC2*DC2
      DC8 = DC6*DC2

!    Calculate TSTAR and 3rd power:
      TSTAR = T/ESCL
      TSTAR3=TSTAR*TSTAR*TSTAR

!    Calculate ln(TSTAR) and 2nd, 3rd, and 4th powers:
      BETA1 = DLOG(TSTAR)
      BETA2 = BETA1*BETA1
      BETA3 = BETA2*BETA1
      BETA4 = BETA3*BETA1

!    Calculate zero-density limit viscosity in muPa-s:
      EXS = DEXP(A0+(A1*BETA1)+(A2*BETA2)+(A3*BETA3)+(A4*BETA4))
      ETA0 = (1.00697D0 * DSQRT(T)) / EXS 

!    Calculate excess viscosity in muPa-s:
      DETA = (D11*DC)+(D21*DC2)+((D64*DC6)/TSTAR3)+(D81*DC8)+ &
             ((D82*DC8)/TSTAR) 

!    Calculate total viscosity in muPa-s:
      VC = ETA0 + DETA

!    Convert viscosity from muPa-s to Pa-s:
      VC = VC * 1.0D-06
      RETURN
end subroutine VISCO2

! ************************************************************************** !

subroutine SAT(T,P)
!--------- Fast SAT M.J.O'Sullivan - 17 SEPT 1990 ---------
!
      use PFLOTRAN_Constants_module, only : H2O_CRITICAL_PRESSURE, &
                                            H2O_CRITICAL_TEMPERATURE
      implicit none

      PetscReal :: T,P,A1,A2,A3,A4,A5,A6,A7,A8,A9
      PetscReal :: TC,X1,X2,SC,PC_
      
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9/ &
      -7.691234564,-2.608023696E1,-1.681706546E2,6.423285504E1, &
      -1.189646225E2,4.167117320,2.097506760E1,1.E9,6./

      if (T.LT.1..OR.T.GT.500.) GOTO 10
      TC=(T+273.15d0)/H2O_CRITICAL_TEMPERATURE
      X1=1.d0-TC
      X2=X1*X1
      SC=A5*X1+A4
      SC=SC*X1+A3
      SC=SC*X1+A2
      SC=SC*X1+A1
      SC=SC*X1
      PC_=EXP(SC/(TC*(1.d0+A6*X1+A7*X2))-X1/(A8*X2+A9))
      P=PC_*H2O_CRITICAL_PRESSURE
      RETURN
   10 continue
      WRITE(6,1) ' ',T
    1 FORMAT(A1,'TEMPERATURE = ',E12.6,'  OUT OF RANGE IN SAT ')
      RETURN
end subroutine SAT

! ************************************************************************** !

subroutine COWAT0(TF,PP,D,U)
!--------- Fast COWAT M.J.O'Sullivan - 17 SEPT 1990 ---------
      use PFLOTRAN_Constants_module, only : H2O_CRITICAL_PRESSURE, &
                                            H2O_CRITICAL_TEMPERATURE

      implicit none
      
      PetscReal :: TF,PP,D,U
      PetscReal :: A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12, &
                   A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23, &
                   SA1,SA2,SA3,SA4,SA5,SA6,SA7,SA8,SA9,SA10,SA11,SA12
      PetscReal :: TKR,TKR2,TKR3,TKR4,TKR5,TKR6,TKR7,TKR8,TKR9,TKR10,TKR11, &
        TKR18,TKR19,TKR20,PNMR,PNMR2,PNMR3,PNMR4,Y,YD,Z,ZP,V,VMKR, &
        CC1,CC2,CC4,CC8,CC10,CZ,SNUM,PRT1,PRT2,PRT3,PRT4,PRT5, &
        AA1,BB1,BB2,EE1,EE3,H,DD1,DD2,DD4,ENTR,PAR1,PAR2,PAR3,PAR4,PAR5

      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12 / &
      6.824687741E3,-5.422063673E2,-2.096666205E4,3.941286787E4, &
      -13.466555478E4,29.707143084E4,-4.375647096E5,42.954208335E4, &
      -27.067012452E4,9.926972482E4,-16.138168904E3,7.982692717E0/
      DATA A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23 / &
      -2.616571843E-2,1.522411790E-3,2.284279054E-2,2.421647003E2, &
      1.269716088E-10,2.074838328E-7,2.174020350E-8,1.105710498E-9, &
      1.293441934E1,1.308119072E-5,6.047626338E-14/
      DATA SA1,SA2,SA3,SA4,SA5,SA6,SA7,SA8,SA9,SA10,SA11,SA12 / &
      8.438375405E-1,5.362162162E-4,1.720000000E0,7.342278489E-2, &
      4.975858870E-2,6.537154300E-1,1.150E-6,1.51080E-5, &
      1.41880E-1,7.002753165E0,2.995284926E-4,2.040E-1   /

      TKR=(TF+273.15)/H2O_CRITICAL_TEMPERATURE
      TKR2=TKR*TKR
      TKR3=TKR*TKR2
      TKR4=TKR2*TKR2
      TKR5=TKR2*TKR3
      TKR6=TKR4*TKR2
      TKR7=TKR4*TKR3
      TKR8=TKR4*TKR4
      TKR9=TKR4*TKR5
      TKR10=TKR4*TKR6
      TKR11=TKR*TKR10
      TKR19=TKR8*TKR11
      TKR18=TKR8*TKR10
      TKR20=TKR10*TKR10
      PNMR=PP/H2O_CRITICAL_PRESSURE
      PNMR2=PNMR*PNMR
      PNMR3=PNMR*PNMR2
      PNMR4=PNMR*PNMR3
      Y=1.d0-SA1*TKR2-SA2/TKR6
      ZP=SA3*Y*Y-2.*SA4*TKR+2.*SA5*PNMR
      if (ZP.LT.0.) GOTO 1
!     19 September 1990.  double on VAX, single on CRAY
!      Z=Y+ DSQRT(ZP)
      Z=Y+SQRT(ZP)
      CZ=Z**(5./17.)
      PAR1=A12*SA5/CZ
      CC1=SA6-TKR
      CC2=CC1*CC1
      CC4=CC2*CC2
      CC8=CC4*CC4
      CC10=CC2*CC8
      AA1=SA7+TKR19
      PAR2=A13+A14*TKR+A15*TKR2+A16*CC10+A17/AA1
      PAR3=(A18+2.*A19*PNMR+3.*A20*PNMR2)/(SA8+TKR11)
      DD1=SA10+PNMR
      DD2=DD1*DD1
      DD4=DD2*DD2
      PAR4=A21*TKR18*(SA9+TKR2)*(-3./DD4+SA11)
      PAR5=3.*A22*(SA12-TKR)*PNMR2+4.*A23/TKR20*PNMR3
      VMKR=PAR1+PAR2-PAR3-PAR4+PAR5
      V=VMKR*3.17E-3
      D=1.d0/V
      YD=-2.d0*SA1*TKR+6.*SA2/TKR7
      SNUM= A10+A11*TKR
      SNUM=SNUM*TKR + A9
      SNUM=SNUM*TKR + A8
      SNUM=SNUM*TKR + A7
      SNUM=SNUM*TKR + A6
      SNUM=SNUM*TKR + A5
      SNUM=SNUM*TKR + A4
      SNUM=SNUM*TKR2 - A2
      PRT1=A12*(Z*(17.*(Z/29.-Y/12.)+5.*TKR*YD/12.)+SA4*TKR- &
      (SA3-1.)*TKR*Y*YD)/CZ
      PRT2=PNMR*(A13-A15*TKR2+A16*(9.*TKR+SA6)*CC8*CC1 &
      +A17*(19.*TKR19+AA1)/(AA1*AA1))
      BB1=SA8+TKR11
      BB2=BB1*BB1
      PRT3=(11.*TKR11+BB1)/BB2*(A18*PNMR+A19* &
      PNMR2+A20*PNMR3)
      EE1=SA10+PNMR
      EE3=EE1*EE1*EE1
      PRT4=A21*TKR18*(17.*SA9+19.*TKR2)*(1./EE3+SA11*PNMR)
      PRT5=A22*SA12*PNMR3+21.*A23/TKR20*PNMR4
      ENTR= A1*TKR - SNUM +PRT1+PRT2-PRT3+PRT4+PRT5
      H=ENTR*70120.4
      U=H-PP*V
      RETURN
    1 continue
      WRITE(6,2)TF
      WRITE(7,2)TF
    2 FORMAT(' TEMPERATURE = ',E12.6,'  OUT OF RANGE IN COWAT')
!  100 FORMAT(1H ,5X,A6,2X,E20.10)
!  102 FORMAT(1H ,5X,A6,5X,I2,2X,E20.10)
      RETURN
end subroutine COWAT0

! ************************************************************************** !

subroutine SUPST(T,P,D,U)
!--------- Fast SUPST M.J.O'Sullivan - 17 SEPT 1990 ---------
! SUPST    1.0 S     1 February  1991
! VAPOR DENSITY AND INTERNAL ENERGY AS FUNCTION OF TEMPERATURE AND 
! PRESSURE (M. OS.)
      use PFLOTRAN_Constants_module, only : H2O_CRITICAL_PRESSURE, &
                                            H2O_CRITICAL_TEMPERATURE
      
      implicit none
      
      PetscReal :: T,P,D,U
      PetscReal :: I1
      PetscReal :: B0,B01,B02,B03,B04,B05,B11,B12,B21,B22,B23,B31,B32,B41,B42, &
                   B51,B52,B53,B61,B62,B71,B72,B81,B82, &
                   B90,B91,B92,B93,B94,B95,B96,SB,SB61,SB71,SB81,SB82
      PetscReal :: THETA,THETA2,THETA3,THETA4
      PetscReal :: BETA,BETA2,BETA3,BETA4,BETA5,BETA6,BETA7,BETAL,DBETAL
      PetscReal :: X,X2,X3,X4,X5,X6,X8,X10,X11,X14,X18,X19,X24,X27
      PetscReal :: R,R2,R4,R6,R10,SD1,SD2,SD3,SC,SN,CHI1,CHI2,V
      PetscReal :: SN6,SN7,SN8,OS1,OS2,OS5,OS6,OS7,EPS2,H

      DATA B0,B01,B02,B03,B04,B05/ &
      16.83599274,28.56067796,-54.38923329,0.4330662834,-0.6547711697, &
      8.565182058E-2/
      DATA B11,B12,B21,B22,B23,B31,B32,B41,B42/ &
      6.670375918E-2,1.388983801,8.390104328E-2,2.614670893E-2, &
      -3.373439453E-2,4.520918904E-1,1.069036614E-1,-5.975336707E-1, &
      -8.847535804E-2/
      DATA B51,B52,B53,B61,B62,B71,B72,B81,B82/ &
      5.958051609E-1,-5.159303373E-1,2.075021122E-1,1.190610271E-1, &
      -9.867174132E-2,1.683998803E-1,-5.809438001E-2,6.552390126E-3, &
      5.710218649E-4/
      DATA B90,B91,B92,B93,B94,B95,B96/ &
      1.936587558E2,-1.388522425E3,4.126607219E3,-6.508211677E3, &
      5.745984054E3,-2.693088365E3,5.235718623E2/
      DATA SB,SB61,SB71,SB81,SB82/ &
      7.633333333E-1,4.006073948E-1,8.636081627E-2,-8.532322921E-1, &
      3.460208861E-1/

      THETA=(T+273.15)/H2O_CRITICAL_TEMPERATURE
      BETA=P/H2O_CRITICAL_PRESSURE
      I1=4.260321148
      X=EXP(SB*(1.d0-THETA))

      X2=X*X
      X3=X2*X
      X4=X3*X
      X5=X4*X
      X6=X5*X
      X8=X6*X2
      X10=X6*X4
      X11=X10*X
      X14=X10*X4
      X18=X14*X4
      X19=X18*X
      X24=X18*X6
      X27=X24*X3

      THETA2=THETA*THETA
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA

      BETA2=BETA*BETA
      BETA3=BETA2*BETA
      BETA4=BETA3*BETA
      BETA5=BETA4*BETA
      BETA6=BETA5*BETA
      BETA7=BETA6*BETA

      BETAL=15.74373327-34.17061978*THETA+19.31380707*THETA2
      DBETAL=-34.17061978+38.62761414*THETA
      R=BETA/BETAL
      R2=R*R
      R4=R2*R2
      R6=R4*R2
      R10=R6*R4

      CHI2=I1*THETA/BETA
      SC=(B11*X10+B12)*X3
      CHI2=CHI2-SC
      SC=B21*X18+B22*X2+B23*X
      CHI2=CHI2-2*BETA*SC
      SC=(B31*X8+B32)*X10
      CHI2=CHI2-3*BETA2*SC
      SC=(B41*X11+B42)*X14
      CHI2=CHI2-4*BETA3*SC
      SC=(B51*X8+B52*X4+B53)*X24
      CHI2=CHI2-5*BETA4*SC

      SD1=1./BETA4+SB61*X14
      SD2=1./BETA5+SB71*X19
      SD3=1./BETA6+(SB81*X27+SB82)*X27

      SN=(B61*X+B62)*X11
!over CHI2=CHI2-SN/SD12*4/BETA5
      chi2=chi2-(sn/sd1*4/beta5)/sd1
      SN=(B71*X6+B72)*X18
!over CHI2=CHI2-SN/SD22*5/BETA6
      chi2=chi2-(sn/sd2*5/beta6)/sd2
      SN=(B81*X10+B82)*X14
!over CHI2=CHI2-SN/SD32*6/BETA7
      chi2=chi2-(sn/sd3*6/beta7)/sd3
      SC=B96
      SC=SC*X+B95
      SC=SC*X+B94
      SC=SC*X+B93
      SC=SC*X+B92
      SC=SC*X+B91
      SC=SC*X+B90
      CHI2=CHI2+11.*R10*SC
      V=CHI2*0.00317
      D=1./V

      OS1=SB*THETA
      EPS2=0.0+B0*THETA-(-B01+B03*THETA2+2*B04*THETA3+3*B05*THETA4)
      SC=(B11*(1.+13.*OS1)*X10+B12*(1.+3.*OS1))*X3
      EPS2=EPS2-BETA*SC
      SC=B21*(1.+18.*OS1)*X18+B22*(1.+2.*OS1)*X2+B23*(1.+OS1)*X
      EPS2=EPS2-BETA2*SC
      SC=(B31*(1.+18.*OS1)*X8+B32*(1.+10.*OS1))*X10
      EPS2=EPS2-BETA3*SC
      SC=(B41*(1.+25.*OS1)*X11+B42*(1.+14.*OS1))*X14
      EPS2=EPS2-BETA4*SC
      SC=(B51*(1.+32.*OS1)*X8+B52*(1.+28.*OS1)*X4+ &
       B53*(1.+24.*OS1))*X24
      EPS2=EPS2-BETA5*SC

      SN6=14.*SB61*X14
      SN7=19.*SB71*X19
      SN8=(54.*SB81*X27+27.*SB82)*X27
      OS5= 1+11.*OS1-OS1*SN6/SD1
      SC=(B61*X*(OS1+OS5)+B62*OS5)*(X11/SD1)
      EPS2=EPS2-SC
      OS6= 1.+24.*OS1-OS1*SN7/SD2
      SC=(B71*X6*OS6+B72*(OS6-6.*OS1))*(X18/SD2)
      EPS2=EPS2-SC
      OS7= 1.+24.*OS1-OS1*SN8/SD3
      SC=(B81*X10*OS7+B82*(OS7-10.* OS1))*(X14/SD3)
      EPS2=EPS2-SC
      OS2=1+THETA*10.0*DBETAL/BETAL
      SC= (OS2+6*OS1)*B96
      SC=SC*X + (OS2+5*OS1)*B95
      SC=SC*X + (OS2+4*OS1)*B94
      SC=SC*X + (OS2+3*OS1)*B93
      SC=SC*X + (OS2+2*OS1)*B92
      SC=SC*X + (OS2+OS1)*B91
      SC=SC*X + OS2*B90
      EPS2=EPS2+BETA*R10*SC
      H=EPS2*70120.4
      U=H-P*V
      RETURN
end subroutine SUPST

! ************************************************************************** !

subroutine TSAT(PX,TX00,TS)
      implicit none
      
!     SATURATION TEMPERATURE TS AT PRESSURE PX.

      PetscReal :: PX,TX00,TX0,TS,PS,DT,TSD,PSD
      
      TX0=TX00
      if (TX0.NE.0.) GOTO 2
!
!-----COME HERE TO OBTAIN ROUGH STARTING VALUE FOR ITERATION.
      TX0=4606./(24.02-LOG(PX)) - 273.15
      TX0=MAX(TX0,5.d0)

    2 CONTINUE
      TS=TX0
      DT=TS*1.E-8
      TSD=TS+DT

    1 CONTINUE

      CALL SAT(TS,PS)

      if (ABS((PX-PS)/PX).LE.1.E-10) RETURN

      TSD=TS+DT
      CALL SAT(TSD,PSD)
      TS=TS+(PX-PS)*DT/(PSD-PS)

      goto 1
      
end subroutine TSAT

! ************************************************************************** !

subroutine SIGMA(T,ST)
      use PFLOTRAN_Constants_module, only : H2O_CRITICAL_TEMPERATURE
      implicit none
!
!-----COMPUTE SURFACE TENSION OF WATER, USING THE
!     "INTERNATIONAL REPRESENTATION OF THE SURFACE TENSION OF
!                                               WATER SUBSTANCE" (1975).

      PetscReal :: T,ST
      
      if (T.GE.374.15) GOTO 1
      ST=1.-0.625*(374.15-T)/H2O_CRITICAL_TEMPERATURE
      ST=ST*.2358*((374.15-T)/H2O_CRITICAL_TEMPERATURE)**1.256
      RETURN

    1 CONTINUE
      ST=0.
      RETURN
end subroutine SIGMA

! ************************************************************************** !

subroutine VIS(T,P,D,VW,VS,PS)
      implicit none
      
!     VISCOSITY OF LIQUID WATER AND VAPOR AS FUNCTION OF
!     TEMPERATURE AND PRESSURE

      PetscReal :: T,P,D,VW,VS,PS
      PetscReal :: EX,PHI,AM,V1
      
      EX=247.8/(T+133.15)
      PHI=1.0467*(T-31.85)
      AM=1.+PHI*(P-PS)*1.E-11
      VW=1.E-7*AM*241.4*10.**EX

      V1=.407*T+80.4
      if (T.LE.350.) VS=1.E-7*(V1-D*(1858.-5.9*T)*1.E-3)
      if (T.GT.350.) VS=1.E-7*(V1+.353*D+676.5E-6*D**2+102.1E-9*D**3)
      RETURN
end subroutine VIS

! ************************************************************************** !

subroutine VISW0(T,P,PS,VW)
      implicit none

!     VISCOSITY OF LIQUID WATER AS FUNCTION OF
!     TEMPERATURE AND PRESSURE

      PetscReal :: T,P,PS,VW
      PetscReal :: EX,PHI,AM

      EX=247.8/(T+133.15)
      PHI=1.0467*(T-31.85)
      AM=1.+PHI*(P-PS)*1.E-11
      VW=1.E-7*AM*241.4*10.**EX

      RETURN
end subroutine VISW0

! ************************************************************************** !

subroutine VISS(T,P,D,VS)
      implicit none

!     VISCOSITY OF VAPOR AS FUNCTION OF
!     TEMPERATURE AND PRESSURE

      PetscReal :: T,P,D,VS,V1
      
      V1=.407*T+80.4
      if (T.LE.350.) VS=1.E-7*(V1-D*(1858.-5.9*T)*1.E-3)
      if (T.GT.350.) VS=1.E-7*(V1+.353*D+676.5E-6*D**2+102.1E-9*D**3)

      RETURN
end subroutine VISS

! ************************************************************************** !

subroutine THERC(T,P,D,CONW,CONS,PS)
      implicit none

!     THERMAL CONDUCTIVITY OF WATER AND VAPOR AS FUNCTION OF
!     TEMPERATURE AND PRESSURE

      PetscReal :: T,P,D,CONW,CONS,PS
      PetscReal :: A0,A1,A2,A3,A4,B0,B1,B2,B3,C0,C1,C2,C3,T0,T1,T2,T3,T4
      PetscReal :: CON1,CON2,CON3,CONS1

      DATA A0,A1,A2,A3,A4/-922.47,2839.5,-1800.7,525.77,-73.440/
      DATA B0,B1,B2,B3/-.94730,2.5186,-2.0012,.51536/
      DATA C0,C1,C2,C3/1.6563E-3,-3.8929E-3,2.9323E-3,-7.1693E-4/
      DATA T0/273.15/

      T1=(T+273.15)/T0
      T2=T1*T1
      T3=T2*T1
      T4=T3*T1

!     if (P-PS.LT.0.) GOTO1
      CON1=A0+A1*T1+A2*T2+A3*T3+A4*T4
      CON2=(P-PS)*(B0+B1*T1+B2*T2+B3*T3)*1.E-5
      CON3=(P-PS)*(P-PS)*(C0+C1*T1+C2*T2+C3*T3)*1.E-10
      CONW=(CON1+CON2+CON3)*1.E-3
      CON1=17.6+5.87E-2*T
      CON2=1.04E-4*T*T
      CON3=4.51E-8*T**3
      CONS1=1.E-3*(CON1+CON2-CON3)
      CONS=CONS1+1.E-6*(103.51+.4198*T-2.771E-5*T*T)*D &
      +1.E-9*D*D*2.1482E14/T**4.2

!     PRINT 10,T,P,PS,CON
!   10 FORMAT(5H T = ,E12.6,5H P = ,E12.6,6H PS = ,E12.6,7H CON = ,E12.6)

      RETURN
!    1 CONTINUE
!     PRINT 2,T,P,PS
!    2 FORMAT(8H AT T = ,E12.6,5H P = ,E12.6,19H IS LESS THAN PS = ,E12.6)
      RETURN
end subroutine THERC
end module co2eos_module
