#!/bin/tcsh -f

set HOMME = /ptmp/jedwards/homme_trunk

cd $HOMME/bld
rm -Rf Depends.pe.*
rm -Rf Depends.sw.*

cd $HOMME
cp -f Params.inc Params.inc.bak

foreach OMP  ( _NO_OMP_THREADS _OMP_THREADS )
foreach RESTART ( _PRESTART _RESTART )
foreach MOVIE ( _PIO _PIO_INTERP _NETCDF  )
foreach MPI ( _MPI _NO_MPI )

   echo $OMP $RESTART $MOVIE $MPI

   mv Params.inc temp.inc
   sed s/OMP.\*/"OMP = $OMP"/ temp.inc |\
   sed s/RESTART.\*/"RESTART = $RESTART"/  |\
   sed s/MPI.\*/"MPI = $MPI"/  |\
   sed s/GRID_STAG.\*/"GRID_STAG = _NONSTAGGER"/  |\
   sed s/MOVIE.\*/"MOVIE = $MOVIE"/  > Params.inc
   cd $HOMME/build.AIX
   gmake dep_pe
   cd $HOMME

end
end
end
end

#
#  shallow water depends. 
#
foreach OMP  ( _NO_OMP_THREADS _OMP_THREADS )
foreach RESTART ( _PRESTART _RESTART )
foreach MOVIE (  _PIO_INTERP _NETCDF _PIO )
foreach MPI ( _MPI _NO_MPI )
foreach STAG ( _NONSTAGGER _STAGGER )

   echo =====================================================================
   echo $OMP $RESTART $MOVIE $MPI $STAG
   echo =====================================================================

   mv Params.inc temp.inc
   sed s/OMP.\*/"OMP = $OMP"/ temp.inc |\
   sed s/RESTART.\*/"RESTART = $RESTART"/  |\
   sed s/MPI.\*/"MPI = $MPI"/  |\
   sed s/GRID_STAG.\*/"GRID_STAG = $STAG"/  |\
   sed s/MOVIE.\*/"MOVIE = $MOVIE"/  > Params.inc

   cd $HOMME/build.AIX
   gmake dep_sw
   cd $HOMME



end
end
end
end
end




cp -f Params.inc.bak Params.inc 
