  Module MPHYSICS_PARAMETERS
! This module declares parameters used in microphysical process. 
! (Used in physics_interface.f90 & physics_v10d.f90)

  USE shr_kind_mod, only: dbl_kind => shr_kind_r8
  USE parmsld,      only: channel_l,nk1
  
  IMPLICIT NONE
     
  INTEGER, PARAMETER :: IM = channel_l, KM = nk1
     
  REAL (KIND=dbl_kind) :: RI50
  REAL (KIND=dbl_kind) :: &
         CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR,CRACI,CSACI,CGACW,CGACI,CRACW,   &
         CSSUB(5),CGSUB(5),CREVP(5),CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,  &
         ES0,CES0,C1BRG,C2BRG,RMI50,RMI40,CPLC,CPLS,CPLF,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),&
         DESWT(151),DESIT(151),TT,TTD 
                                                  
  REAL (KIND=dbl_kind) :: &
         PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ, &
         VK, Z0_ROUGHNESS            
  REAL (KIND=dbl_kind) :: VCONR,VCONS,VCONG
  REAL (KIND=dbl_kind) :: DTL,TDT,RHOSFC
     
  REAL (kind=dbl_kind), PARAMETER :: SFCRHO=1.1_8
    
  end module mphysics_parameters
  