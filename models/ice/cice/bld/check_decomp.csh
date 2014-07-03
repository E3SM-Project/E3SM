#!/bin/csh

set res = (gx3v7 gx1v6 tx0.1v2 ne240np4)
set pe1 = (1 1 1 1)
set pex = (500 4000 60000 60000)

foreach grid (1 2 3 4)

foreach thrd (1 4)

set n = $pe1[$grid]
echo " "
echo " res pes = task x thrd : nx ny bx by maxb type set"
echo " "

while ($n < $pex[$grid])

 foreach check (1 2 3)

   if ($check == 1) then
      @ tasks = $n
   else if ($check == 2) then
      @ tasks = $n + 1
   else if ($check == 3) then
      @ tasks = ($n * 3) / 2
   endif

   @ pes = $tasks * $thrd

   set config = `  ./generate_cice_decomp.pl -res $res[$grid] -nproc $tasks -thrds $thrd -output all`

#   if ($config[1] >= 0) then
      echo $res[$grid] ${pes} = ${tasks} x ${thrd} : $config 
#   endif

 end

 @ n = $n * 2

end

end

end
