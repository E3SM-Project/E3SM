set term pslatex color solid rotate auxfile
set output "w_profile.pslatex"
set key below
set xlab "distance along profile [km]"
set ylab "vertical velocity [ma$^{-1}$]
set missing "NaN"
plot "velo/calc2.data" u ($1):($3) t "upper kinetic BC" w l, \
     "velo/e1-mm2.1.data" u ($1):($3) t "uncorrected" w l, \
     "velo/e1-mm.1.data" u ($1):($3) t "corrected" w l
