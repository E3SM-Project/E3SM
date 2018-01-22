import os, time
from mpi4py import MPI


comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()


#output_dir = '/lustre/or-hydra/cades-ccsi/scratch/xyk/TESTRD4'
output_dir = '/lustre/or-hydra/cades-ccsi/scratch/dmricciuto/'

#sites=['AU-Tum','BR-Sa1','FI-Hyy','FR-Pue','GF-Guy','RU-Cok','RU-Fyo','RU-Sam','US-Blo','US-Ha1','US-PFa','US-SRM','US-UMB','US-Wkg','ZA-Kru']
#sites=['BR-Sa1','FI-Hyy','GF-Guy','US-Ha1','US-UMB']
sites=['AT-Neu','AU-Tum','BE-Bra','BE-Vie','BR-Sa1','CA-Gro','CZ-BK1','DE-Hai','DE-Tha','DK-Sor','DK-ZaH','FI-Hyy','FR-Pue','GF-Guy','IT-Lav','IT-Ren','IT-Ro2','NL-Loo','RU-Cok','RU-Fyo','RU-Sam','US-Blo','US-Ha1','US-Los','US-Oho','US-PFa','US-SRM','US-Syv','US-Ton','US-UMB','US-Var','US-WCr','US-Wkg', 'ZA-Kru']
#sites = ['AT-Neu','AU-Tum']

time_shift=[1, 10, 1, 1, -4, -6, 1, 1, 1, 1, 1, -2, 2, 1, -3, 1, 1, 1, 1, 11, 3, 3, -8, -5, -6, -5, -6, -8, -6, -8, -5, -8, -6, -7, 2]

cases=['171206_ECAv1_CNP']
compsets=['ICB20TRCNPRDCTCBC']
endyear=2014

var_list='NEE,GPP,ER,EFLX_LH_TOT,FSH,TBOT,FSDS,WIND,RAIN'

n_sites = len(sites)

if (rank == 0):  #master
   #os.system("rm "+cases[0]+"*.nc")
   n_done = 0
   for n_job in range(1,size):
       comm.send(n_job,dest=n_job,tag=1)
       comm.send(0,    dest=n_job,tag=2)
       time.sleep(1)
   for n_job in range(size,n_sites+1):
       process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
       thisjob = comm.recv(source=process, tag=4)
       n_done = n_done+1
       comm.send(n_job, dest=process, tag=1)
       comm.send(0,     dest=process, tag=2)
   while(n_done < n_sites):
       process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
       thisjob = comm.recv(source=process, tag=4)
       n_done = n_done+1
       comm.send(-1, dest=process, tag=1)
       comm.send(-1, dest=process, tag=2)
   times = ['hourly','daily','monthly','annual']
   #times=['monthly','annual']
   for t in times:
       os.system("ncecat ./plots/"+cases[0]+'/'+t+"/"+cases[0]+"_??-???_model.nc "+cases[0]+"_"+t+"_model.nc")
       os.system("ncecat ./plots/"+cases[0]+"/"+t+"/"+cases[0]+"_??-???_obs.nc "+cases[0]+"_"+t+"_obs.nc")
       os.system("ncwa -O -a gridcell -d gridcell,0,0 "+cases[0]+"_"+t+"_model.nc "+cases[0]+"_"+t+"_model.nc")
       os.system("ncwa -O -a gridcell -d gridcell,0,0 "+cases[0]+"_"+t+"_obs.nc "+cases[0]+"_"+t+"_obs.nc")
       os.system("ncrename -h -O -d record,gridcell "+cases[0]+"_"+t+"_model.nc "+cases[0]+"_"+t+"_model.nc")
       os.system("ncrename -h -O -d record,gridcell "+cases[0]+"_"+t+"_obs.nc "+cases[0]+"_"+t+"_obs.nc")
       #os.system("rm "+cases[0]+"_??-???_"+t+"_*.nc")
   MPI.Finalize()
else:
    status=0
    while status == 0:
         myjob  = comm.recv(source=0, tag=1)
         status = comm.recv(source=0, tag=2)
         if (status == 0):
             
             #monthly
             os.system('python plotcase.py --case '+cases[0]+' --compset '+compsets[0]+' --site ' \
                +sites[myjob-1]+' --obs --vars '+var_list+' --ystart 1991 --yend '+str(endyear)+' --csmdir ' \
                +output_dir+' --noplot')     
             #hourly
             os.system('python plotcase.py --case '+cases[0]+' --compset '+compsets[0]+' --site ' \
                +sites[myjob-1]+' --obs --vars '+var_list+' --ystart 1991 --yend '+str(endyear)+' --h1 --csmdir ' \
                +output_dir+' --noplot --timezone '+str(time_shift[myjob-1]))
             #daily
             os.system('python plotcase.py --case '+cases[0]+' --compset '+compsets[0]+' --site ' \
                +sites[myjob-1]+' --obs --vars '+var_list+' --ystart 1991 --yend '+str(endyear)+' --h2 --csmdir ' \
                +output_dir+' --noplot')
             #annual
             os.system('python plotcase.py --case '+cases[0]+' --compset '+compsets[0]+' --site ' \
                +sites[myjob-1]+' --obs --vars '+var_list+' --ystart 1991 --yend '+str(endyear)+' --h4 --csmdir ' \
                +output_dir+' --noplot')

             comm.send(rank, dest=0,tag=3)
             comm.send(myjob,dest=0,tag=4) 
    MPI.Finalize() 
