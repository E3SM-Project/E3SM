load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
begin

whichVar="PS" ; PS or WIND

;configs=(/"reed_bryan","reed_reed","kessler_bryan","kessler_reed"/)
configs=(/"kessler_bryan","kessler_bryan2"/)
do ii = 0,dimsizes(configs)-1

all_files = systemfunc ("ls /glade/u/home/zarzycki/scratch/cam5_work/dcmip2/"+configs(ii)+"/dcmip2*.nc")
a      = addfiles (all_files, "r")   ; note the "s" of addfile

;========================
; choose how files are combined and read in variable across files
;========================  
ListSetType (a, "cat")             ; concatenate or "merge" (default)

if (whichVar .eq. "PS") then
  PS=a[:]->PS(:,{0.0:90.0},:)                   ; just take out NH since we know storm is there
  PS=PS/100.                                 ; convert from Pa to hPa/mb
  var_vtime = dim_min_n_Wrap(PS,(/1,2/))     ; for now, just find min over lat/lon domain at each time

else if (whichVar .eq. "WIND") then
  UBOT=a[:]->U(:,dimsizes(a->lev)-1,{0.0:90.0},:) ; dimsizes(a->lev)-1 subsets bottom level (CAM ordered top to bottom)
  VBOT=a[:]->V(:,dimsizes(a->lev)-1,{0.0:90.0},:)
  WIND=UBOT
  WIND=sqrt(UBOT^2+VBOT^2)
  WIND@long_name="Horizontal wind magnitude"
  var_vtime = dim_max_n_Wrap(WIND,(/1,2/))

else
  print("unsupported var choice... exiting...")
  exit
end if
end if

if (ii .eq. 0)
  var_vtime_all=new((/dimsizes(configs),dimsizes(var_vtime)/),typeof(var_vtime))
  time = var_vtime&time
end if

var_vtime_all(ii,:)=var_vtime

delete(var_vtime)
delete(a)
delete(all_files)

end do

wks   = gsn_open_wks ("x11","xy")               ; send graphics to PNG file

res                  = True                     ; plot mods desired
res@gsnDraw          = False
res@gsnFrame         = False
res@tiMainString     = "DCMIP2 CAM-SE "+whichVar+" v. time"          ; add title

; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .90                   ; move units right
 res@pmLegendOrthogonalPosF = -0.5                  ; more neg = down
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.25                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
 res@xyExplicitLegendLabels = configs

plot  = gsn_csm_xy (wks,time,var_vtime_all,res) ; create plot
draw(plot)
frame(wks)

end
