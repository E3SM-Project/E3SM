set term pslatex color solid rotate auxfile
set output "wt_sigma.pslatex"
set key below
set multiplot
set lmargin 10
set origin 0,0
set size square 0.5,1
set xtics 0.2
set xlab "vertical velocity [ma$^{-1}$]"
set ylab "normalised height"
plot "velo/velo_uncor.data" u ($2):(1-$1) t "uncorrected" w l 2, "velo/velo_cor.data" u ($2):(1-$1) t "" w l 3
set origin 0.5,0
set size square 0.5,1
set xtics 5
set xlab "temperature [$^{\\circ}$C]"
plot "velo/temp_uncor.data" u ($2):(1-$1) t "" w l 2, "velo/temp_cor.data" u ($2):(1-$1) t "corrected" w l 3
set nomultiplot