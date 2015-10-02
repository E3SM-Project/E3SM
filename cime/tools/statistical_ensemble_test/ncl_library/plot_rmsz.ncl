load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

undef("plot_RMSZscores")
procedure plot_RMSZscores(RMSZScore, fens, vars, machine, compiler, plot_type)
begin

  ps_plot = (plot_type.eq."ps").or.(plot_type.eq."eps")
  nvars = dimsizes(vars)

  plot_min = 0.77d
  plot_max = 1.49d
  nbins = 12
  score = (RMSZScore-plot_min)/(plot_max-plot_min)
  input_dims = getfiledimsizes(fens)
  nfile = input_dims(ind(getvardims(fens).eq."ens_size"))
  ;max_RMSE = new(nfile, double)
  max_RMSE = new(nfile, float)

  res = True
  res@trYMaxF = 40
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnHistogramClassIntervals = fspan(plot_min, plot_max, nbins+1)
  res@tiXAxisString = "Root Mean Square Z-score"
  res@tiMainFuncCode = "~"
  res@tmXBLabelAngleF = 60.
  res@gsFillColor = "gray"

  if (min(RMSZScore).lt.plot_min) then
    print((/"WARNING: your smallest raw score is smaller than the ensemble minimum"/))
    res@trXMinF = min(score)
  end if
  if (max(RMSZScore).gt.plot_max) then
    print((/"WARNING: your largest raw score is larger than the ensemble maximum"/))
    res@trXMaxF = max(score)
  end if

  txres=True
  txres@txFontHeightF = 0.0225

  if (.not.ps_plot) then
    wks = gsn_open_wks(plot_type, machine+"."+compiler+".validation_plots")
  end if

  all_vars = new(nvars, integer)
  do v = 0, nvars-1
    v_in_f = ind(chartostring(fens->vars).eq.vars(v))
    all_vars(v) = v_in_f
    val = round(1000.*RMSZScore(v),0)/1000.+""
    large_count = num(ind(fens->RMSZ(v_in_f,:).gt.RMSZScore(v)))
    res@tiMainString = "Raw score for "+vars(v)+" = "+val+"~C~Run is using the "+compiler+" compiler on "+machine+"~C~"+large_count+"/101 members produce a larger RMSZ"

    if (ps_plot) then
      wks = gsn_open_wks(plot_type, machine+"."+compiler+"."+vars(v)+".validation_plots")
    end if
    plot=gsn_histogram(wks,fens->RMSZ(v_in_f,:), res)
    X_on_plot=gsn_add_text(wks, plot, "X",score(v),1,txres)
    val = round(100.*RMSZScore(v),0)/100.+""
    text_on_plot=gsn_add_text(wks, plot, val,score(v),3,txres)
    draw(plot)
    frame(wks)
    if (ps_plot) then
      delete(wks)
    end if
  end do

  ; Max over all vars
; Commented out because of wide variance of RMSZ scores
if (False)
  do i=0,nfile-1
    max_RMSE(i) = max(fens->RMSZ(all_vars,i))
  end do
  large_count = num(ind(dim_max_n(fens->RMSZ(all_vars,:),0).gt.max(RMSZScore)))
  val = round(1000.*max(RMSZScore),0)/1000.+""
  res@tiMainString = "Max RMSZ = "+val+"~C~Run is using the "+compiler+" compiler on "+machine+"~C~"+large_count+"/101 members produce a larger max RMSZ"
  if (ps_plot) then
    wks = gsn_open_wks(plot_type, machine+"."+compiler+".max.validation_plots")
  end if
  plot=gsn_histogram(wks,max_RMSE, res)
  X_on_plot=gsn_add_text(wks, plot, "X",max(score),1,txres)
  val = round(100.*max(RMSZScore),0)/100.+""
  text_on_plot=gsn_add_text(wks, plot, val,max(score),3,txres)
  draw(plot)
  frame(wks)
end if

end
