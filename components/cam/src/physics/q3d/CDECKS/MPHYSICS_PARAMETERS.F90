  Module MPHYSICS_PARAMETERS
! This module declares parameters used in microphysical process. 
! (Used in physics_interface.f90 & physics_v10d.f90)

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE parmsld,      only: channel_l,nk1
  
  IMPLICIT NONE
     
  INTEGER, PARAMETER :: IM = channel_l, KM = nk1
  
  REAL (kind=r8), PARAMETER :: SFCRHO=1.225_r8
  ! standard density at 1013 mb
     
  REAL (KIND=r8) :: RI50
  REAL (KIND=r8) :: &
         CRACS,CSACR,CGACR,CGACS,ACCO(3,4),CSACW,CWACS,CIACR,CRACI,CSACI,CGACW,CGACI,CRACW,   &
         CSSUB(5),CGSUB(5),CREVP(5),CGFR(2),CSMLT(5),CGMLT(5),CGWET(4),CRAUT(2),QI0,QS0,QL0,  &
         ES0,CES0,C1BRG,C2BRG,RMI50,RMI40,CPLC,CPLS,CPLF,TDIF,ESW00,ESI00,ESWT(151),ESIT(151),&
         DESWT(151),DESIT(151),TT,TTD 
                                                  
  REAL (KIND=r8) :: &
         PIE,GRAV,VDIFU,TCOND,RVAPR,RDRYA,VISK,HLTS,HLTC,HLTF,CH2O,CICE,TICE,CP,EPS,TTFRZ
  REAL (KIND=r8) :: VCONR,VCONS,VCONG
  REAL (KIND=r8) :: DTL,TDT,RHOSFC
     
  end module mphysics_parameters
  