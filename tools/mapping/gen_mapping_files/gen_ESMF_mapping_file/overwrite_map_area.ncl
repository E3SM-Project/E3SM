begin

  if (any(.not.isdefined((/"map_file","src_file","dst_file"/)))) then
    dq = str_get_dq()
    print((/"Usage: ncl overwrite_map_area.ncl map_file=\"+dq+"SCRIP_MAP\"+dq+" src_file=\"+dq+"SRC_SCRIP_GRID\"+dq+" dst_file=\"+dq+"DST_SCRIP_GRID\"+dq/))
    status_exit(1)
  end if
  fsrc=addfile(src_file,"r")
  fdst=addfile(dst_file,"r")
  fmap=addfile(map_file,"w")

  setfileoption(fmap,"DefineMode",False)
  if (isfilevar(fsrc,"grid_area")) then
    fmap->area_a = (/fsrc->grid_area/)
  end if
  if (isfilevar(fsrc,"grid_imask")) then
    fmap->frac_a = 1.d*(/fsrc->grid_imask/)
  end if
  if (isfilevar(fdst,"grid_area")) then
    fmap->area_b = (/fdst->grid_area/)
  end if
  if (isfilevar(fdst,"grid_imask")) then
    fmap->frac_b = 1.d*(/fdst->grid_imask/)
  end if

  ; Get area-weighting right!
  if (lupdate_S) then
    print((/"Changing S in "+map_file/))
    if (isfilevar(fsrc,"grid_area").and.isfilevar(fdst,"grid_area")) then
      ns = getfilevardimsizes(fmap,"S")
      do n=0,ns-1
        isrc=fmap->col(n)-1
        idst=fmap->row(n)-1
        wgt=fmap->area_a(isrc)/fmap->area_b(idst)
        fmap->S(n)=fmap->S(n)*wgt
      end do
    end if
  else
    print((/"Not changing S in "+map_file+" (set lupdate_S=True to edit)"/))
  end if
end
