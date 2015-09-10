module netcdf_coord_file
  use netcdf
  implicit none
  private
  public write_netcdf_coord_file
contains
  subroutine write_netcdf_coord_file(filename, nlev, p0, hyai, hybi, hyam, hybm)
    character*(*), intent(in) :: filename
    integer, intent(in) :: nlev
    real*8 , intent(in) :: p0, hyai(:), hybi(:), hyam(:), hybm(:)
    integer :: levid, ilevid, lev_id, ilev_id, hyai_id, hybi_id, hyam_id, hybm_id, ncid, p0_id
    real*8 :: lev(nlev), ilev(nlev+1)
    integer :: k

    call check(nf90_create(path= trim(filename), cmode=nf90_clobber, ncid=ncid ))

    call check(nf90_def_dim( ncid, "lev", len=nlev, dimid=levid))

    call check(nf90_def_var( ncid, "lev", nf90_double, dimids=levid, varid=lev_id))
    call check(nf90_put_att( ncid, lev_id, "long_name", "hybrid level at midpoints (P0*(A+B))"))
    call check(nf90_put_att( ncid, lev_id, "units", "level"))
    call check(nf90_put_att( ncid, lev_id, "positive","down"))
    call check(nf90_put_att( ncid, lev_id, "standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
    call check(nf90_put_att( ncid, lev_id, "formula_terms","a: hyam b: hybm p0: P0 "))


    call check(nf90_def_dim( ncid, "ilev", len=nlev+1, dimid=ilevid))

    call check(nf90_def_var( ncid, "ilev", nf90_double, dimids=ilevid, varid=ilev_id))
    call check(nf90_put_att( ncid, ilev_id, "long_name", "hybrid level at interfaces (P0*(A+B))"))
    call check(nf90_put_att( ncid, ilev_id, "units", "level"))
    call check(nf90_put_att( ncid, ilev_id, "positive","down"))
    call check(nf90_put_att( ncid, ilev_id, "standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
    call check(nf90_put_att( ncid, ilev_id, "formula_terms","a: hyai b: hybi p0: P0 "))


    call check(nf90_def_var( ncid, "hyai", nf90_double, dimids=ilevid, varid=hyai_id))
    call check(nf90_put_att( ncid, hyai_id, "long_name", "hybrid A coefficient at layer interfaces"))

    call check(nf90_def_var( ncid, "hybi", nf90_double, dimids=ilevid, varid=hybi_id))
    call check(nf90_put_att( ncid, hybi_id, "long_name", "hybrid B coefficient at layer interfaces"))

    call check(nf90_def_var( ncid, "hyam", nf90_double, dimids=levid, varid=hyam_id))
    call check(nf90_put_att( ncid, hyam_id, "long_name", "hybrid A coefficient at layer midpoints"))

    call check(nf90_def_var( ncid, "hybm", nf90_double, dimids=levid, varid=hybm_id))
    call check(nf90_put_att( ncid, hybm_id, "long_name", "hybrid B coefficient at layer midpoints"))


    call check(nf90_def_var( ncid, "P0", nf90_double, varid=p0_id))
    call check(nf90_put_att( ncid, p0_id, "long_name", "reference pressure"))
    call check(nf90_put_att( ncid, p0_id, "units", "Pa"))

    call check(nf90_enddef(ncid))

    do k=1,nlev
       lev(k) = p0*(hyam(k)+hybm(k))
       ilev(k) = p0*(hyai(k)+hybi(k))
    end do
    ilev(nlev+1)=p0*(hyai(nlev+1)+hybi(nlev+1))

    call check(nf90_put_var(ncid,lev_id, lev))
    call check(nf90_put_var(ncid,ilev_id,ilev))

    call check(nf90_put_var(ncid,p0_id,p0*100.))
    call check(nf90_put_var(ncid,hyai_id,hyai))
    call check(nf90_put_var(ncid,hybi_id,hybi))
    call check(nf90_put_var(ncid,hyam_id,hyam))
    call check(nf90_put_var(ncid,hybm_id,hybm))

    call check(nf90_close(ncid))
    

    
  end subroutine write_netcdf_coord_file

 



  subroutine check(status, line)
    integer, intent(in) ::status
    integer, intent(in), optional :: line

    if(status /= nf90_noerr) then
       if(present(line)) then
          print *, trim(nf90_strerror(status)),' at line ',line, ' of file ',__FILE__
       else
          print *, trim(nf90_strerror(status))
       endif
       stop __FILE__
    end if
  end subroutine check
end module netcdf_coord_file
