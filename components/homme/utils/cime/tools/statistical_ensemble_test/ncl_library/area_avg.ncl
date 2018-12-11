undef("area_avg_SE")
function area_avg_SE(data, weight)
; Calculate global mean for single level of CAM-SE output
begin
  area2d = (/sum(data*weight(:))/)
  if (typeof(area2d).eq."double") then
    return(doubletofloat(area2d))
  else
    return(area2d)
  end if
end

undef("area_avg_FV")
function area_avg_FV(data, weight)
; Calculate global mean for single level of CAM-FV output
begin
  area2d = (/wgt_areaave(data, weight, 1.0, 1)/)
  if (typeof(area2d).eq."double") then
    return(doubletofloat(area2d))
  else
    return(area2d)
  end if
end

undef("calculate_global_mean")
function calculate_global_mean(runfile, ensemble_file, var_name3d, var_name2d, tind, vinterp, FailTest, VarInEns, verbose)
begin

; (1) Initial setup
;     (1a) Open run_file and ensemble_file

  frun = addfile(run_file,"r")
  fens = addfile(ensemble_file,"r")

;     (1b) Store dimensions of input fields
  input_dims = getfiledimsizes(frun)
  input_dim_names = getvardims(frun)
  nlev = -1
  ncol = -1
  nlat = -1
  nlon = -1
  ilev = -1
  icol = -1
  ilat = -1
  ilon = -1

  do i=0,dimsizes(input_dims)-1
    if (input_dim_names(i).eq."lev") then
      nlev = input_dims(i)
      ilev = i
    end if
    if (input_dim_names(i).eq."ncol") then
      ncol = input_dims(i)
      icol = i
    end if
    if (input_dim_names(i).eq."lat") then
      nlat = input_dims(i)
      ilat = i
    end if
    if (input_dim_names(i).eq."lon") then
      nlon = input_dims(i)
      ilon = i
    end if
  end do
  if (nlev.eq.-1) then
    print((/"ERROR: can not find nlev dimension in "+filelist(0)/))
    status_exit(1)
  end if
  if (((nlon.eq.-1).or.(nlat.eq.-1)).and.(ncol.eq.-1)) then
    print((/"ERROR: can not find ncol or nlat/nlon dimension in "+filelist(0)/))
    status_exit(1)
  end if
  if (ncol.ne.-1) then
    one_spatial_dim = True
  else
    one_spatial_dim = False
  end if

  nvars3d = dimsizes(var_name3d)
  nvars2d = dimsizes(var_name2d)

  ;    (1c) Compute Area weights
  if (one_spatial_dim) then
    area_wgt = frun->area/sum(frun->area)
  else
    area_wgt = frun->gw
  end if

  ; (2) Compute means
  GlobalMean2D = new(nvars2d, float)
  GlobalMean3D = new(nvars3d, float)
  GlobalMean3D_lev = new((/nvars3d, nlev/), float)

  do v=0,nvars3d-1
    if (isfilevar(frun, var_name3d(v))) then
      output = frun->$var_name3d(v)$
      do k=0,nlev-1
        if (one_spatial_dim) then
          GlobalMean3D_lev(v,k) = area_avg_SE(output(tind, k, :), area_wgt)
        else
          GlobalMean3D_lev(v,k) = area_avg_FV(output(tind, k, :, :), area_wgt)
        end if
        GlobalMean3D(v) = avg(GlobalMean3D_lev(v,:))
      end do
      delete(output)
    end if
  end do

  do v=0,nvars2d-1
    if (isfilevar(frun, var_name2d(v))) then
      output = frun->$var_name2d(v)$
      if (one_spatial_dim) then
        GlobalMean2D(v) = area_avg_SE(output(tind, :), area_wgt)
      else
        GlobalMean2D(v) = area_avg_FV(output(tind, :, :), area_wgt)
      end if
      delete(output)
    end if
  end do

  ; (3) Compare to ensemble
  v = 0
  do i=0, nvars3d-1
    v_in_f = ind(chartostring(fens->vars).eq.var_name3d(i))
    if (VarInEns(v)) then
      FailTest(v) = ((GlobalMean3D(i).gt.max(fens->global_mean(v_in_f,:)))  .or. \
                     (GlobalMean3D(i).lt.min(fens->global_mean(v_in_f,:))))
    end if
    v = v+1
  end do

  do i=0, nvars2d-1
    v_in_f = ind(chartostring(fens->vars).eq.var_name2d(i))
    if (VarInEns(v)) then
      FailTest(v) = ((GlobalMean2D(i).gt.max(fens->global_mean(v_in_f,:)))  .or. \
                     (GlobalMean2D(i).lt.min(fens->global_mean(v_in_f,:))))
    end if
    FailTest(v) = FailTest(v).and.VarInEns(v)
    v = v+1
  end do

  return(array_append_record(GlobalMean3D, GlobalMean2D,0))
end
