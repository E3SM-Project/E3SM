#!/usr/bin/csh
#./testpio_bench.pl --iofmt pnc --pecount 128 --bench POPC --numIO 8
#./testpio_bench.pl --maxiter 2 --iofmt pnc --pecount 64 --bench POPD --numIO 12

./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 1000 --bench POPC  --numIO 160
./testpio_bench.pl --maxiter 5  --iofmt pnc --pecount 1000 --bench POPD  --numIO 160
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 832 --bench CAM05 --numIO 140
#./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 1000 --bench WRFB  --numIO 160
