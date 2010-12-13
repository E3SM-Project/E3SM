#!/usr/bin/csh
#./testpio_bench.pl --iofmt pnc --pecount 128 --bench POPC --numIO 8
#./testpio_bench.pl --maxiter 2 --iofmt pnc --pecount 64 --bench POPD --numIO 12

./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 64 --bench POPC --numIO 6
./testpio_bench.pl --maxiter 2 --iofmt pnc --pecount 64 --bench POPD --numIO 6
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 64 --bench CAM05 --numIO 6
#./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 256 --bench WRFB  --numIO 40

