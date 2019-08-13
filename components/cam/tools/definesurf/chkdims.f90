subroutine chkdims (fileid, name, varid, londimid, latdimid, timdimid, verbose)

  implicit none

  include 'netcdf.inc'

  integer fileid, varid, londimid, latdimid
  integer timdimid
  logical verbose
  character*(*) name

  integer ret
  integer ndims, dimids(nf_max_dims)

  ret = nf_inq_varid (fileid, name, varid)

  if (ret.eq.NF_NOERR) then

    dimids(:) = -999
    ret = nf_inq_varndims (fileid, varid, ndims)
    ret = nf_inq_vardimid (fileid, varid, dimids)

    if (ret.ne.NF_NOERR) then
      write(6,*)'NF_INQ_VAR failed for ',name
      call handle_error (ret)
    end if

    if (ndims.eq.3 .and. dimids(3).ne.timdimid) then
      write(6,*)'3rd dim of ', name, ' must be time'
      call endrun
    end if
      
    if (dimids(1).ne.londimid .or. dimids(2).ne.latdimid) then
      write(6,*)'Dims of ', name,' must be lon by lat'
      call endrun
    end if

    if (verbose) write(6,*)'Overwriting existing ',name,' with hi-res topo'

  else

    dimids(1) = londimid
    dimids(2) = latdimid
    dimids(3) = timdimid
    if (verbose) write(6,*)name,' does not exist on netcdf file: Creating.'
    ret = nf_redef (fileid)
    ret = nf_def_var (fileid, name, NF_DOUBLE, 3, dimids, varid)
    if (ret.ne.NF_NOERR) call handle_error (ret)
    ret = nf_enddef (fileid)

  end if
end subroutine chkdims
