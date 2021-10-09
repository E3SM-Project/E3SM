;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  CaseName     = "E3SMv1_FC5_NDG_RADCPL_OPT1_DT1800.ne30_ne30.cori-knl"
  ModelOutRoot = "/compyfs/zhan391/202010_SciDAC_simulation/regrid/climo"
  FileName     = ModelOutRoot+"/"+CaseName+"/"+CaseName+"_ANN_climo.nc"
 
  f         = addfile(FileName,"r")
  lat       = f->lat
  lon       = f->lon
  nlev      = dimsizes(f->ilev)
  nvars     = 8 
  data       = new((/nlev,nvars/),double)
  data(:,0)  = ispan(0,nlev-1,1)
  data(:,1)  = f->hyai
  data(:,2)  = f->hybi
  data(:,3)  = f->ilev
  data(1:nlev-1,4) = f->lev
  data(0,4)        = -9999
  tdz        = stdatmus_p2tdz(f->ilev)
  print(tdz)
  tdz(0,:) = tdz(0,:) + 273.15
  do i = 0,2,1
    data(:,i+5) = (/tdz(i,:)/)
  end do 
 
 ;data(:,5:)   = data(::-1,5:) 
  fmtf   = nlev + "f15.6"  
  opt  = True
  opt@title  = "       N             a[pa]           b           ph[hPa]        pf [hPa]      Temperature[K]     Density[kg/m^3]      Geometric Altidue [m]"
  opt@fout = "e3sm_L72_model_level_definitions.txt"
  write_matrix (data, fmtf, opt)
  
end

