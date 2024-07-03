# luca's verbatim example works from any dir now!
from pyeamxx import pyscream as ps
import numpy as np

ps.init()

v = np.zeros(shape=(10,2))

for i in range(0,10):
    for j in range(0,2):
        v[i,j] = i*2 + j + 1 

print(f"v:{v}")

f = ps.Field("T_mid",v)

f.print()

f.cleanup()
ps.finalize()
