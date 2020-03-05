import numpy
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

fig = plt.gcf()
nRow = 6
nCol = 2
iTime = [8, 16]
nu = ['0.01', '0.1', '1', '10', '100', '200']
time = ['hour 8', 'hour 16']

fig, axs = plt.subplots(nRow, nCol, figsize=(
    5.3 * nCol, 2.0 * nRow), constrained_layout=True)

for iRow in range(nRow):
    ncfile = Dataset('output_' + str(iRow + 1) + '.nc', 'r')
    var = ncfile.variables['temperature']
    xtime = ncfile.variables['xtime']
    for iCol in range(nCol):
        ax = axs[iRow, iCol]
        dis = ax.imshow(var[iTime[iCol], 0:512:4, :].T, extent=[
                        0, 120, 20, 0], aspect=2, cmap='jet', vmin=5, vmax=30)
        if iRow == nRow - 1:
            ax.set_xlabel('x, km')
        if iCol == 0:
            ax.set_ylabel('depth, m')
        if iCol == nCol - 1:
            fig.colorbar(dis, ax=axs[iRow, iCol], aspect=5)
        ax.set_title(time[iCol] + ", " + r"$\nu_h=$" + nu[iRow])
    ncfile.close()

plt.savefig('sections_lock_exchange.png', bbox_inches='tight')
