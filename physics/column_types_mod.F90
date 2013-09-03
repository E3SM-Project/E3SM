module column_types_mod

  use hybvcoord_mod,   only : hvcoord_t 
  use hybrid_mod,      only : hybrid_t
  use time_mod,        only : TimeLevel_t
  use dimensions_mod,  only : nlev, nlevp, np
  use kinds,           only : real_kind, int_kind

  implicit none

  public 


  type, public :: ColumnDataEmanuel_t
     ! in
     real (kind=real_kind)             :: T(nlev)
     real (kind=real_kind)             :: Q(nlev)
     real (kind=real_kind)             :: QS(nlev)
     real (kind=real_kind)             :: U(nlev)
     real (kind=real_kind)             :: V(nlev)
     real (kind=real_kind)             :: TRA(nlev,1)
     real (kind=real_kind)             :: P(nlev)
     real (kind=real_kind)             :: PH(nlevp)
     real (kind=real_kind)             :: R(nlev)
     
     ! out
     real (kind=real_kind)             :: FT(nlev)
     real (kind=real_kind)             :: FQ(nlev)
     real (kind=real_kind)             :: FU(nlev)
     real (kind=real_kind)             :: FV(nlev)
     real (kind=real_kind)             :: FTRA(nlev,1)
     real (kind=real_kind)             :: PRECIP   
     real (kind=real_kind)             :: WD
     real (kind=real_kind)             :: TPRIME
     real (kind=real_kind)             :: QPRIME

  end type ColumnDataEmanuel_t
  
  type, public :: ColumnModelEmanuel_t
     type (ColumnDataEmanuel_t) :: col(np,np)
     logical                    :: INIT = .FALSE.
  end type ColumnModelEmanuel_t
  
  type, public :: HeldSuarezForcing_t
     logical                    :: INIT = .FALSE.
  end type HeldSuarezForcing_t

  type, public :: AquaplanetForcing_t
     logical                    :: INIT = .FALSE.
  end type AquaplanetForcing_t

  type, public :: DataMulticloud_t
     ! Need to be computed from t=0 sounding
     real (kind=real_kind)         :: Qt0,Qt1,Qt2
     real (kind=real_kind)         :: phi(nlev,2),psi(nlev,2)
     real (kind=real_kind)         :: psitrunc(nlev,2)
     real (kind=real_kind)         :: cstmode(nlev)
     real (kind=real_kind)         :: Aphi(2,2)
     real (kind=real_kind)         :: Apsi(2,2)
     real (kind=real_kind)         :: psi1avg,psi2avg
     real (kind=real_kind)         :: psi1avgtrunc,psi2avgtrunc    
     real (kind=real_kind)         :: Q0c,Q0R2,Tau_evap_Inv,m0     
  end type DataMulticloud_t
  type, public :: ColumnModelMulticloud_t
     type (DataMulticloud_t)        :: D
     !real (kind=real_kind)         :: Qt0,Qt1,Qt2
     !real (kind=real_kind)         :: phi(nlev,2),psi(nlev,2)
     !real (kind=real_kind)         :: psitrunc(nlev,2)
     !real (kind=real_kind)         :: cstmode(nlev)
     !real (kind=real_kind)         :: Aphi(2,2)
     !real (kind=real_kind)         :: Apsi(2,2)
     !real (kind=real_kind)         :: psi1avg,psi2avg
     !real (kind=real_kind)         :: psi1avgtrunc,psi2avgtrunc    
     !real (kind=real_kind)         :: Q0c,Q0R2,Tau_evap_Inv,m0
     
     ! Either from namelist or hard coded in initialization
     real (kind=real_kind)         :: alpha,Q0R1,TstarMinTeb,Tminus,Tplus,Tau_Radiation,Tau_Damping
     real (kind=real_kind)         :: CD0,U0,h,TebMinTem,Tau_Convec, Tau_strat, Tau_congest
     real (kind=real_kind)         :: alpha_strat,alpha_congest,mu,gamma2,gamma2prime,alpha2
     real (kind=real_kind)         :: a0,a1,a2,a0prime
     real (kind=real_kind)         :: Htall,csr
     real (kind=real_kind)         :: lambdastar,A,B
     real (kind=real_kind)         :: Tau_r,relaxation
     real (kind=real_kind)         :: vertmask(nlev) 
     real (kind=real_kind)         :: sAlpha,sT,sBeta
     integer                       :: closure
     integer                       :: shutconvec

     logical                       :: INIT = .FALSE.

  end type ColumnModelMulticloud_t

  type, public :: ColumnModel_t
     type (ColumnModelEmanuel_t)       :: cm_em
     type (ColumnModelMulticloud_t)    :: cm_mc
     type (AquaplanetForcing_t)        :: cm_aq
     type (HeldSuarezForcing_t)        :: cm_hs
     type (hvcoord_t)                  :: hvcoord
     type (hybrid_t), pointer          :: hybrid
     type (TimeLevel_t), pointer       :: tl
     integer                           :: nets,nete
  end type ColumnModel_t
  

end module column_types_mod
