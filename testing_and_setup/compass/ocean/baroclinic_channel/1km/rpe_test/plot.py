import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

fig = plt.gcf()
nRow = 1  # 2
nCol = 5
nu = ['1', '5', '10', '100', '200']
iTime = [0]
time = ['20']

# ---nx,ny for 10 km
#nx = 16
#ny = 50

# ---nx,ny for 4 km
#nx = 40
#ny = 126

# ---nx,ny for 1 km
nx = 160
ny = 500

fig, axs = plt.subplots(nRow, nCol, figsize=(
    2.1 * nCol, 5.0 * nRow), constrained_layout=True)

for iCol in range(nCol):
    for iRow in range(nRow):
        ncfile = Dataset('output_' + str(iCol + 1) + '.nc', 'r')
        var = ncfile.variables['temperature']
        var1 = np.reshape(var[iTime[iRow], :, 0], [ny, nx])
        # --- flip in y-dir
        var = np.flipud(var1)

        # --- Every other row in y needs to average two neighbors in x on planar hex mesh
        var_avg = var
        for j in range(0, ny, 2):
            for i in range(1, nx - 2):
                var_avg[j, i] = (var[j, i + 1] + var[j, i]) / 2.0

        if nRow == 1:
            ax = axs[iCol]
        if nRow > 1:
            ax = axs[iRow, iCol]
        dis = ax.imshow(
            var_avg,
            extent=[
                0,
                160,
                0,
                500],
            cmap='jet',
            vmin=11.8,
            vmax=13.0)
        ax.set_title("day " + time[iRow] + ", " + r"$\nu_h=$" + nu[iCol])
        ax.set_xticks(np.arange(0, 161, step=40))
        ax.set_yticks(np.arange(0, 501, step=50))

        if iRow == nRow - 1:
            ax.set_xlabel('x, km')
        if iCol == 0:
            ax.set_ylabel('y, km')
        if iCol == nCol - 1:
            if nRow == 1:
                fig.colorbar(dis, ax=axs[nCol - 1], aspect=40)
            if nRow > 1:
                fig.colorbar(dis, ax=axs[iRow, nCol - 1], aspect=40)
        ncfile.close()

if nx == 16:
    res = '10'
if nx == 40:
    res = '4'
if nx == 160:
    res = '1'

plt.savefig("sections_baroclinic_channel_" + res + "km.png")
