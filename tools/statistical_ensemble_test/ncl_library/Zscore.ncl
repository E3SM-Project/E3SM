load "$NCL_VALID_LIB_DIR/global.ncl"

undef("calc_Z")
function calc_Z(val, avg, stddev, count, flag)
; Calculate Z-score, print warning if standard deviation = 0
begin
  dims = dimsizes(val)
  return_val = new(dims, float)
  if (all(stddev.le.tol)) then
    if (flag) then
      print((/"WARNING: ALL standard dev = 0"/))
      flag = False
    end if
    count = count+num(stddev.le.tol)
    return_val = 0.
  else
    if (any(stddev.le.tol)) then
      if (flag) then
        print((/"WARNING: some standard dev = 0"/))
        flag = False
      end if
      count = count+num(stddev.le.tol)
      if (dimsizes(dims).eq.2) then
        ; lat,lon coords instead of ncol coords
        nlat = dims(0)
        nlon = dims(1)
        do i=0,nlat-1
          do j=0,nlon-1
            if (stddev(i,j).le.tol) then
              return_val(i,j) = 0.
            else
              return_val(i,j) =(val(i,j)-avg(i,j))/stddev(i,j)
            end if
          end do
        end do
      else
        return_val(ind(stddev.le.tol)) = 0.
        inds = ind(stddev.gt.tol)
        return_val(inds) = (val(inds)-avg(inds))/stddev(inds)
      end if
    else
      return_val = (val-avg)/stddev
    end if
  end if

  return(return_val)
end

