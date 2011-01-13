#!/usr/bin/csh
#./testpio_bench.pl --iofmt pnc --pecount 128 --bench POPC --numIO 8
#./testpio_bench.pl --maxiter 2 --iofmt pnc --pecount 64 --bench POPD --numIO 12

./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 512 --bench POPC  --numIO 80
./testpio_bench.pl --maxiter 2  --iofmt pnc --pecount 512 --bench POPD  --numIO 80
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 416 --bench CAM05 --numIO 70
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 500 --bench WRFB  --numIO 80
