#!/usr/bin/csh
set id = `date "+%m%d%y-%H%M"`
#set host = 'kraken'
set host = 'hopper'
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 80  --bench POPB --numIO 30  --log ${host}.0080.pnc.iotask_30.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 160 --bench POPB --numIO 30  --log ${host}.0160.pnc.iotask_30.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 320 --bench POPB --numIO 50  --log ${host}.0320.pnc.iotask_50.log.${id}
./testpio_bench.pl --maxiter 10 --iofmt pnc --pecount 640 --bench POPB --numIO 100 --log ${host}.0640.pnc.iotask_100.log.${id}

