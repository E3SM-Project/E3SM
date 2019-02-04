load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

  ; copy data file and append "_z"
  print("cp ./movies/dcmip2016_test31.nc ./movies/dcmip2016_test31_z.nc")
  system("cp ./movies/dcmip2016_test31.nc ./movies/dcmip2016_test31_z.nc")

  ; open in read/write mode
  f   = addfile("./movies/dcmip2016_test31_z.nc", "w")

  ;_________________________________________________
  ; read in fields and invert vertical order (so k=1 is bottom level)
  ; output_varnames1  ='T','p','ps','pnh','geo','u','v','w','Th','Q','Q2','Q3','Q4'   ! variables to write to file


  time  = f->time
  nt    = dimsizes(time)
  print("nt="+nt)

  geo   = f->geo(:,::-1,:,:)
  z_eta = geo/9.81
  z     = f->lev
  z     = ispan(250,19750,500)

  print("z="+z)
  print("z_eta="+z_eta(0,:,0,0))
  print("interpolating from geopotential levels to evenly spaced z")

  ; interpolate to z and write to file
  print("converting T")
  T     = f->T  (:,::-1,:,:)
  f->T     = linint1_n_Wrap(z_eta,T    ,False,z,0,1)

  print("converting p")
  p     = f->pnh(:,::-1,:,:)
  f->p     = linint1_n_Wrap(z_eta,p    ,False,z,0,1)

  print("converting pnh")
  pnh   = f->pnh(:,::-1,:,:)
  f->pnh   = linint1_n_Wrap(z_eta,pnh  ,False,z,0,1)

  print("converting geo")
  geo   = f->geo(:,::-1,:,:)
  f->geo   = linint1_n_Wrap(z_eta,geo  ,False,z,0,1)

  print("converting u")
  u     = f->u  (:,::-1,:,:)
  f->u     = linint1_n_Wrap(z_eta,u    ,False,z,0,1)

  print("converting v")
  v     = f->v  (:,::-1,:,:)
  f->v     = linint1_n_Wrap(z_eta,v    ,False,z,0,1)

  print("converting w")
  w     = f->w  (:,::-1,:,:)
  f->w     = linint1_n_Wrap(z_eta,w    ,False,z,0,1)

  print("converting Th")
  Th    = f->Th (:,::-1,:,:)
  f->Th    = linint1_n_Wrap(z_eta,Th   ,False,z,0,1)

  print("converting qv")
  qv    = f->Q  (:,::-1,:,:)
  f->Q     = linint1_n_Wrap(z_eta,qv   ,False,z,0,1)

  print("converting qc")
  qc    = f->Q2 (:,::-1,:,:)
  f->Q2    = linint1_n_Wrap(z_eta,qc   ,False,z,0,1)

  print("converting qr")
  qr    = f->Q3 (:,::-1,:,:)
  f->Q3    = linint1_n_Wrap(z_eta,qr   ,False,z,0,1)

end
