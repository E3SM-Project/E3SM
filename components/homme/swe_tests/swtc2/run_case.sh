
#!/bin/bash


mkdir -p run-$2
cd run-$2
rm -rf *

mpirun -np $1 ../../build/sweqx-$2/sweqx.cpu < ../swtc2.nl > swtc2.log
