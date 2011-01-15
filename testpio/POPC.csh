#!/usr/bin/csh
set id = `date "+%m%d%y-%H%M"`
#set host = 'kraken'
set host = 'hopper'
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 64  --bench POPC  --numIO 12  --log ${host}.0064.pnc.iotask_12.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 128 --bench POPC  --numIO 20  --log ${host}.0128.pnc.iotask_20.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 256 --bench POPC  --numIO 40  --log ${host}.0256.pnc.iotask_40.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 512 --bench POPC  --numIO 80  --log ${host}.0512.pnc.iotask_80.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 1000 --bench POPC --numIO 160 --log ${host}.1000.pnc.iotask_160.log.${id}
