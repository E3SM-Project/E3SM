reset

set term pslatex color solid rotate aux
set output "stress-strain.pslatex"

load "fit.gp"

set size 1.,0.8
set lmargin 7
set bmargin 3
set multiplot
set nokey

set label "A" at 50.10, 106.50
set label "B" at 28.30, 77.50
set label "C" at 56.70, 82.00
set label "D" at 40.60, 64.60
set label "E" at 24.20, 50.80
set label "F" at 51.00, 39.40
set label "G" at 19.20, 26.60

set size 0.4,0.8
set ylabel "shear stress [kPa]"
set label 1 'A' at 2.5,110
plot [0:60][0:120] i1(x,30) t "" w l 0, i1(x,25) t "" w l 0, i1(x,20) t "" w l 0, i1(x,15) t "" w l 0, i1(x,10) t "" w l 0, i1(x,5) t "" w l 0,\
	"stress-strain.data" u ($3):($2) w p 1, t0(x) w l 2
set origin 0.3,0.
set format y ""
set ylabel ""

set xlabel "effective pressure [kPa]"
set label 1 'B' at 2.5,110
plot [0:60][0:120] i2(x,30) t "" w l 0, i2(x,25) t "" w l 0, i2(x,20) t "" w l 0, i2(x,15) t "" w l 0, i2(x,10) t "" w l 0, i2(x,5) t "" w l 0,\
	"stress-strain.data" u ($3):($2) w p 1, t0(x) w l 2

set origin 0.6,0
set xlabel ""
set label 1 'C' at 2.5,110
plot [0:60][0:120] i3(x,30) t "" w l 0, i3(x,25) t "" w l 0, i3(x,20) t "" w l 0, i3(x,15) t "" w l 0, i3(x,10) t "" w l 0, i3(x,5) t "" w l 0,\
	"stress-strain.data" u ($3):($2) w p 1, t0(x) w l 2

set nomultiplot
