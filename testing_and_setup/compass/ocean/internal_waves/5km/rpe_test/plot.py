import numpy
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

fig = plt.gcf()
nRow = 4
nCol = 2
nu = ['0.01', '1', '15', '150']
iTime = [1, 2]
time = ['day 10', 'day 20']

fig, axs = plt.subplots(nRow, nCol, figsize=(
    4.0 * nCol, 3.7 * nRow), constrained_layout=True)

for iRow in range(nRow):
    ncfile = Dataset('output_' + str(iRow + 1) + '.nc', 'r')
    var = ncfile.variables['temperature']
    xtime = ncfile.variables['xtime']
    for iCol in range(nCol):
        ax = axs[iRow, iCol]
        dis = ax.imshow(var[iTime[iCol], 0::4, :].T, extent=[
                        0, 250, 500, 0], aspect='0.5', cmap='jet', vmin=10, vmax=20)
        if iRow == nRow - 1:
            ax.set_xlabel('x, km')
        if iCol == 0:
            ax.set_ylabel('depth, m')
        if iCol == nCol - 1:
            fig.colorbar(dis, ax=axs[iRow, iCol], aspect=10)
        ax.set_title(time[iCol] + ", " + r"$\nu_h=$" + nu[iRow])
    ncfile.close()

plt.savefig('sections_internal_waves.png')
