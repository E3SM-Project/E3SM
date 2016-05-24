reset

set term pslatex color solid rotate aux
set output "stress-strain-fit.pslatex"

load "fit.gp"

set size square 
set key bottom right 
set xlabel 'observed $\dot\epsilon$ [a$^{-1}$]'
set ylabel 'calculated $\dot\epsilon$ [a$^{-1}$]'
plot [0:32][0:32] x t "" w l 0, \
	"stress-strain.data" u ($4):(f1($2,$3)) t "Model A" w p 1, \
	"stress-strain.data" u ($4):(f2($2,$3)) t "Model B" w p 2, \
	"stress-strain.data" u ($4):(f3($2,$3)) t "Model C" w p 3
