#!/usr/bin/csh
set id = `date "+%m%d%y-%H%M"`
set host = 'kraken'
#set host = 'hopper'
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 125  --bench WRFB --numIO 32  --log ${host}.0125.pnc.iotask_32.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 250  --bench WRFB --numIO 40  --log ${host}.0250.pnc.iotask_40.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 500  --bench WRFB --numIO 80  --log ${host}.0500.pnc.iotask_80.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 1000 --bench WRFB --numIO 160 --log ${host}.1000.pnc.iotask_160.log.${id}
