#!/bin/csh

set res = gx1v5
set n = 0
set nmax = 2400
set file = $res.$nmax

echo "check_decomp on $res from $n to $nmax pes:" >! $file
echo "" >>& $file

while ($n < $nmax)

 @ n++
 set config = `./generate_pop_decomp.pl -res $res -nproc $n`

 if ($config[1] >= 0) then
    echo gx1v5 $n : $config  >>& $res.$nmax
 endif

end
