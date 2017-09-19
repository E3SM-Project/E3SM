undef("print_stats")
procedure print_stats(cnt, zscores, varname)
begin
  ; Print stats to screen
  if (cnt.gt.0) then
    print("Total number of points with 0 variance in "+varname+" = "+cnt)
  end if
  print((/"RMSZ for "+varname+": mean = "+avg(zscores)+", stddev = "+stddev(zscores)/))
end
