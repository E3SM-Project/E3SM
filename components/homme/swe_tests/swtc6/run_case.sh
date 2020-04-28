
#!/bin/bash


mkdir -p run-$2
cd run-$2
rm -rf *

mpirun.mpich -np $1 ../../build/sweqx-$2/sweqx.cpu < ../swtc6.nl > swtc6.log
