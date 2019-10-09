;
; Compare two files without a time-axis and report if they are identical or not.
;
;  Erik Kluzek
;  May/04/2011
;  $Id: getregional_datasets.ncl 25177 2010-10-16 05:12:30Z erik $
;  $HeadURL;
;
begin
  ; ===========================================================================================================
  ;
  ; IMPORTANT NOTE: EDIT THE FOLLOWING TO CUSTOMIZE or use ENV VARIABLE SETTINGS
  ; Edit the following as needed
  ;
  ; Input files to compare
  ;
  file1      = getenv("MYFILE1")      ; Get the first file
  file2      = getenv("MYFILE2")      ; Get the second file
  debug_str  = getenv("DEBUG");       ; Don't run just -- debug
  print_str  = getenv("PRINT");       ; Do Extra printing for debugging
  break_str  = getenv("BREAKONDIFF"); ; Break on first difference

  if ( ismissing(file1) )then
     file1 = "../interpinic/clmi.BCN.1949-01-01_10x15_USGS_simyr1850_c100322.nc"
  end if
  if ( ismissing(file2) )then
     file2 = file1
  end if
  if ( ismissing(print_str) )then
     printn = False;
  else
     if ( print_str .eq. "TRUE" )then
        printn = True;
     else
        printn = False;
     end if
  end if
  if ( ismissing(break_str) )then
     breakondiff = False;
  else
     if ( break_str .eq. "TRUE" )then
        breakondiff = True;
     else
        breakondiff = False;
     end if
  end if
  if ( ismissing(debug_str) )then
     debug = False;
  else
     if ( debug_str .eq. "TRUE" )then
        print( "DEBUG is TRUE do extra printing AND do NOT execute -- just print what WOULD happen" );
        debug  = True;
        printn = True;
     else
        debug = False;
     end if
  end if
function commonlist( list1 [*]:string, list2 [*]:string )
;
; Figure out the list of common names in the two lists
;
local i, j;
begin
  common = new( (/ min( (/ dimsizes(list1), dimsizes(list2) /) ) /), "string" );
  n = 0;
  do i = 0, dimsizes(list1)-1
     found = False
     do j = 0, dimsizes(list2)-1
        if ( list2(j) .eq. list1(i) )then
           found = True;
           break;
        end if
     end do
     if ( found )then
        common(n) = list1(i);
        n = n + 1
     end if
  end do
  return( common(0:n-1) );
end
  ;
  ; Open files
  ;
  print( "Compare "+file1+" to "+file2 );
  if ( systemfunc("test -f "+file1+"; echo $?" ) .ne. 0 )then
     print( "Input file does not exist or not found: "+file1 );
     exit
  end if
  if ( printn .eq. True )then
     print( "file1:"+file1 );
  end if
  nc1 = addfile( file1, "r" );
  if ( systemfunc("test -f "+file2+"; echo $?" ) .ne. 0 )then
     print( "Input file does not exist or not found: "+file2 );
     exit
  end if
  if ( printn .eq. True )then
     print( "file2:"+file2 );
  end if
  nc2 = addfile( file2, "r" );

  ;
  ; Get the common list of variables
  ;
  varnames1 = getfilevarnames( nc1 );
  varnames2 = getfilevarnames( nc2 );

  vars = commonlist( varnames1, varnames2 );
  if ( printn .eq. True )then
     print( "var: "+vars )
  end if

  ;
  ; Go through variables and compare
  ;
  differ = False;
  do i = 0, dimsizes(vars)-1
     if ( printn .eq. True )then
        print( "compare var: "+vars(i) )
     end if
     if ( any( dimsizes(nc1->$vars(i)$) .ne. dimsizes(nc2->$vars(i)$) ) )then
        if ( printn .eq. True )then
           print( "var dimensions differ: "+vars(i) )
        end if
        differ = True;
        break;
     end if
     typ = typeof( nc1->$vars(i)$ );
     if ( typ .ne. "character" .and. typ .ne. "logical" )then
        if ( any( abs((/nc1->$vars(i)$/) - (/nc2->$vars(i)$/)) .gt. 0.0d00 ) )then
           differ = True;
           if ( printn .eq. True )then
              diff = abs((/nc1->$vars(i)$/) - (/nc2->$vars(i)$/));
              print( "  var differs:  "+vars(i) )
              print( "  max abs diff: "+max(diff)+" avg diff: "+avg(diff) )
              delete(diff)
           end if
        end if
     else 
        if ( any((/nc1->$vars(i)$/) .ne. (/nc2->$vars(i)$/)) )then
           differ = True;
           if ( printn .eq. True )then
              print( "var differs: "+vars(i) )
           end if
        end if
     end if
     if ( differ .and. breakondiff )then
        break
     end if
  end do
  if ( differ )then
     print( "  diff_test: the two files seem to be DIFFERENT " )
  else
     print( "  diff_test: the two files seem to be IDENTICAL " )
  end if

end
