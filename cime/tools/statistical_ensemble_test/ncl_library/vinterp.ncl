load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

undef("get_interpolation_levels")
function get_interpolation_levels(interp_type)
; For interpolating in vertical direction, plevs are the same as plvlM in AMWG
; diagnostics (functions_vertical_xz) "for comparison with another model 
; interpolate the model data to these 25 constant pressure levels" (they also
; use log interpolation)
begin
  plevs = (/30., 50., 70., 100., 150., 200., 250., 300., 350., 400., 450., \
            500., 550., 600., 650., 700., 750., 800., 850., 875., 900., 925., \
            950., 975., 1000./)
  plevs@units="mb"
  plevs@long_name = "pressure"
  
  interp_type = 2 ; 1 = LINEAR, 2 = LOG, 3 = LOG LOG
  return(plevs)
end

undef("vertical_interp")
function vertical_interp(f, varname, plevs, interp_type, extrap)
; Just a wrapper for vinth2p_nodes
begin
  return(vinth2p_nodes(f->$varname$(1,:,:), f->hyam, f->hybm, plevs, \
                       f->PS(1,:), interp_type, 1000., 1, extrap))
end
