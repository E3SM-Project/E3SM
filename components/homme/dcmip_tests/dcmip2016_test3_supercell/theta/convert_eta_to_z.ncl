load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

  ; copy data file and append "_z"
  system("cp ./movies/dcmip2016_test31.nc ./movies/dcmip2016_test31_z.nc")

  ; open in read/write mode
  f   = addfile("./movies/dcmip2016_test31_z.nc", "w")

  ;_________________________________________________
  ; read in fields and invert vertical order (so k=1 is bottom level)
  geo   = f->geo(:,::-1,:,:)
  t     = f->T  (:,::-1,:,:)
  Th    = f->Th (:,::-1,:,:)
  u     = f->u  (:,::-1,:,:)
  v     = f->v  (:,::-1,:,:)
  w     = f->w  (:,::-1,:,:)
  p     = f->pnh(:,::-1,:,:)
  qv    = f->Q  (:,::-1,:,:)
  qv    = f->Q2 (:,::-1,:,:)
  qv    = f->Q3 (:,::-1,:,:)

  ;_________________________________________________
  ; Interpolate fields from eta levels to z levels

  z_eta = geo/9.81
  z     = f->lev
  z     = ispan(250,19750,500)

  print("z="+z)
  print("z_eta="+z_eta(:,0,0))
  print("interpolating from geopotential levels to evenly spaced z")

  ; interpolate to z and write to file
  f->Th    = linint1_n_Wrap(z_eta,Th   ,False,z,0,1)
  f->t     = linint1_n_Wrap(z_eta,t    ,False,z,0,1)
  f->u     = linint1_n_Wrap(z_eta,u    ,False,z,0,1)
  f->v     = linint1_n_Wrap(z_eta,v    ,False,z,0,1)
  f->w     = linint1_n_Wrap(z_eta,w    ,False,z,0,1)
  f->p     = linint1_n_Wrap(z_eta,p    ,False,z,0,1)
  f->qv    = linint1_n_Wrap(z_eta,qv   ,False,z,0,1)
  f->qc    = linint1_n_Wrap(z_eta,qc   ,False,z,0,1)
  f->qr    = linint1_n_Wrap(z_eta,qc   ,False,z,0,1)

end
