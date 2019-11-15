import numpy
import xarray
import os
import copy
import progressbar
import subprocess
import glob

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors


class TimeSeriesPlotter(object):
    '''
    A plotter object to hold on to some info needed for plotting time series
    from ISOMIP+ simulation results

    Attributes
    ----------
    inFolder : str
        The folder with simulation results

    outFolder : str
        The folder where images will be written

    expt : int
        The number of the experiment (0 for Ocean0, 1 for Ocean1, etc.)
    '''
    def __init__(self, inFolder='.', outFolder='plots', expt=0):
        '''
        Create a plotter object to hold on to some info needed for plotting
        time series from ISOMIP+ simulation results

        Parameters
        ----------
        inFolder : str, optional
            The folder with simulation results

        outFolder : str, optional
            The folder where images will be written

        expt : int, optional
            The number of the experiment (0 for Ocean0, 1 for Ocean1, etc.)

        '''

        self.inFolder = inFolder
        self.outFolder = outFolder
        self.expt = expt

        self.dsMesh = xarray.open_dataset('{}/init.nc'.format(self.inFolder))

        self.ds = xarray.open_mfdataset(
            '{}/timeSeriesStatsMonthly*.nc'.format(self.inFolder),
            concat_dim='Time')

        try:
            os.makedirs(self.outFolder)
        except OSError:
            pass

    def plot_melt_time_series(self, sshMax=None):
        '''
        Plot a series of image for each of several variables related to melt
        at the ice shelf-ocean interface: mean melt rate, total melt flux,
        mean thermal driving, mean friction velocity
        '''

        rho_fw = 1000.
        secPerYear = 365*24*60*60

        areaCell = self.dsMesh.areaCell
        iceMask = self.ds.timeMonthly_avg_landIceFraction
        meltFlux = self.ds.timeMonthly_avg_landIceFreshwaterFlux
        if sshMax is not None:
            ssh = self.ds.timeMonthly_avg_ssh
            iceMask = iceMask.where(ssh < sshMax)

        totalMeltFlux = (meltFlux*areaCell*iceMask).sum(dim='nCells')
        totalArea = (areaCell*iceMask).sum(dim='nCells')
        meanMeltRate = totalMeltFlux/totalArea/rho_fw*secPerYear
        self.plot_time_series(meanMeltRate, 'mean melt rate', 'meanMeltRate',
                              'm/yr')

        self.plot_time_series(1e-6*totalMeltFlux, 'total melt flux',
                              'totalMeltFlux', 'kT/yr')

        da = (self.ds.timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature -
              self.ds.timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature)
        da = (da*areaCell*iceMask).sum(dim='nCells')/totalArea

        self.plot_time_series(da, 'mean thermal driving',
                              'meanThermalDriving', 'deg C')

        da = self.ds.timeMonthly_avg_landIceFrictionVelocity
        da = (da*areaCell*iceMask).sum(dim='nCells')/totalArea

        self.plot_time_series(da, 'mean friction velocity',
                              'meanFrictionVelocity', 'm/s')

    def plot_time_series(self, da, nameInTitle, prefix, units=None,
                         figsize=[12, 6], color=None, overwrite=True):

        fileName = '{}/{}.png'.format(self.outFolder, prefix)
        if(not overwrite and os.path.exists(fileName)):
            return

        nTime = da.sizes['Time']
        time = numpy.arange(nTime)/12.

        plt.figure(1, figsize=figsize)
        plt.plot(time, da.values, color=color)

        if units is None:
            ylabel = nameInTitle
        else:
            ylabel = '{} ({})'.format(nameInTitle, units)

        plt.ylabel(ylabel)
        plt.xlabel('time (yrs)')

        plt.savefig(fileName)
        plt.close()


class MoviePlotter(object):
    '''
    A plotter object to hold on to some info needed for plotting images from
    ISOMIP+ simulation results

    Attributes
    ----------
    inFolder : str
        The folder with simulation results

    outFolder : str
        The folder where images will be written

    expt : int
        The number of the experiment (0 for Ocean0, 1 for Ocean1, etc.)

    sectionY : float
        The location along the y axis of a transect in the x-z plane to plot

    cmap : Colormap or str
        A color map for plots

    dsMesh : ``xarray.Dataset``
        A data set with mesh data

    ds : ``xarray.Dataset``
        A data set with the montly-mean simulation results

    oceanMask : ``numpy.ndarray``
        A mask of cells that are in the ocean domain (probably all ones)

    cavityMask : ``numpy.ndarray``
        A mask of cells that are in the sub-ice-shelf cavity

    oceanPatches : ``PatchCollection``
        A set of polygons covering ocean cells

    cavityPatches : ``PatchCollection``
        A set of polygons covering only cells in the cavity

    X, Z : ``numpy.ndarray``
        The horiz. and vert. coordinates of the x-z cross section

    sectionMask : ``numpy.ndarray``
        A mask for the cross section indicating where values are valid (i.e.
        above the bathymetry)
    '''

    def __init__(self, inFolder='.', outFolder='plots', expt=0,
                 sectionY=40e3, cmap=None):
        '''
        Create a plotter object to hold on to some info needed for plotting
        images from ISOMIP+ simulation results

        Parameters
        ----------
        inFolder : str, optional
            The folder with simulation results

        outFolder : str, optional
            The folder where images will be written

        expt : int, optional
            The number of the experiment (0 for Ocean0, 1 for Ocean1, etc.)

        sectionY : float, optional
            The location along the y axis of a transect in the x-z plane to
            plot

        cmap : Colormap or str
            A color map to plot, default is based on the rainbow map from
            Ferret
        '''
        self.inFolder = inFolder
        self.outFolder = outFolder
        self.expt = expt
        self.sectionY = sectionY

        if cmap is None:
            self.cmap = _make_ferret_colormap()
        else:
            self.cmap = cmap

        self.dsMesh = xarray.open_dataset('{}/init.nc'.format(self.inFolder))

        landIceMask = self.dsMesh.landIceFraction.isel(Time=0) > 0.01
        self.oceanMask = self.dsMesh.maxLevelCell-1 >= 0
        self.cavityMask = numpy.logical_and(self.oceanMask, landIceMask)

        self.oceanPatches = _compute_cell_patches(
            self.dsMesh, self.oceanMask, self.cmap)
        self.cavityPatches = _compute_cell_patches(
            self.dsMesh, self.cavityMask, self.cmap)

        self.sectionCellIndices = _compute_section_cell_indices(self.sectionY,
                                                                self.dsMesh)

        self.ds = xarray.open_mfdataset(
            '{}/timeSeriesStatsMonthly*.nc'.format(self.inFolder),
            concat_dim='Time')

        self._compute_section_x_z()

    def plot_barotropic_streamfunction(self, vmin=None, vmax=None):
        '''
        Plot a series of image of the barotropic streamfunction

        Parameters
        ----------
        vmin, vmax : float, optional
            The minimum and maximum values for the colorbar, defaults are
            chosen depending on ``expt``
        '''

        ds = xarray.open_dataset('{}/barotropicStreamfunction.nc'.format(
            self.inFolder))

        if vmin is None or vmax is None:
            if self.expt in [0, 1]:
                vmin = -1
                vmax = 1
            else:
                vmin = -0.5
                vmax = 0.5

        nTime = ds.sizes['Time']
        widgets = ['plotting barotropic streamfunction: ',
                   progressbar.Percentage(), ' ',
                   progressbar.Bar(), ' ', progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets,
                                      maxval=nTime).start()

        for tIndex in range(nTime):
            self.update_date(tIndex)
            bsf = ds.bsfCell.isel(Time=tIndex)
            outFileName = '{}/bsf/bsf_{:04d}.png'.format(
                self.outFolder, tIndex+1)
            self._plot_horiz_field(bsf, title='barotropic streamfunction (Sv)',
                                   outFileName=outFileName, oceanDomain=True,
                                   vmin=vmin, vmax=vmax)
            bar.update(tIndex+1)
        bar.finish()

    def plot_overturning_streamfunction(self, vmin=-0.3, vmax=0.3):
        '''
        Plot a series of image of the overturning streamfunction

        Parameters
        ----------
        vmin, vmax : float, optional
            The minimum and maximum values for the colorbar
        '''

        ds = xarray.open_dataset('{}/overturningStreamfunction.nc'.format(
            self.inFolder))

        nTime = ds.sizes['Time']
        widgets = ['plotting overturning streamfunction: ',
                   progressbar.Percentage(), ' ',
                   progressbar.Bar(), ' ', progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets,
                                      maxval=nTime).start()

        for tIndex in range(nTime):
            self.update_date(tIndex)
            osf = ds.osf.isel(Time=tIndex)
            outFileName = '{}/osf/osf_{:04d}.png'.format(self.outFolder,
                                                         tIndex+1)
            x = _interp_extrap_corner(ds.x.values)
            z = _interp_extrap_corner(ds.z.values)
            self._plot_vert_field(
                x, z, osf, title='overturning streamfunction (Sv)',
                outFileName=outFileName, vmin=vmin, vmax=vmax)
            bar.update(tIndex+1)
        bar.finish()

    def plot_melt_rates(self, vmin=-100., vmax=100.):
        '''
        Plot a series of image of the melt rate

        Parameters
        ----------
        vmin, vmax : float, optional
            The minimum and maximum values for the colorbar
        '''
        rho_fw = 1000.
        secPerYear = 365*24*60*60

        da = secPerYear/rho_fw*self.ds.timeMonthly_avg_landIceFreshwaterFlux

        self.plot_horiz_series(da, 'melt rate', prefix='meltRate',
                               oceanDomain=False, units='m/yr', vmin=vmin,
                               vmax=vmax)

    def plot_ice_shelf_boundary_variables(self):
        '''
        Plot a series of image for each of several variables related to the
        ice shelf-ocean interface: heat flux from the ocean, heat flux into the
        ice, thermal driving, haline driving, and the friction velocity under
        ice
        '''

        self.plot_horiz_series(self.ds.timeMonthly_avg_landIceHeatFlux,
                               'heat flux from ocean to ice-ocean interface',
                               prefix='oceanHeatFlux',
                               oceanDomain=False, units='W/s',
                               vmin=-1e3, vmax=1e3)

        self.plot_horiz_series(self.ds.timeMonthly_avg_heatFluxToLandIce,
                               'heat flux into ice at ice-ocean interface',
                               prefix='iceHeatFlux',
                               oceanDomain=False, units='W/s',
                               vmin=-1e1, vmax=1e1)

        da = (self.ds.timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerTemperature -
              self.ds.timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceTemperature)
        self.plot_horiz_series(da, 'thermal driving',
                               prefix='thermalDriving',
                               oceanDomain=False, units='deg C',
                               vmin=-2, vmax=2)

        da = (self.ds.timeMonthly_avg_landIceBoundaryLayerTracers_landIceBoundaryLayerSalinity -
              self.ds.timeMonthly_avg_landIceInterfaceTracers_landIceInterfaceSalinity)
        self.plot_horiz_series(da, 'haline driving',
                               prefix='halineDriving',
                               oceanDomain=False, units='PSU',
                               vmin=-10, vmax=10)

        self.plot_horiz_series(self.ds.timeMonthly_avg_landIceFrictionVelocity,
                               'friction velocity',
                               prefix='frictionVelocity',
                               oceanDomain=True, units='m/s',
                               vmin=0, vmax=0.05)

    def plot_temperature(self):
        '''
        Plot a series of images of temperature at the sea surface or
        ice-ocean interface, sea floor and in an x-z section
        '''

        da = self.ds.timeMonthly_avg_activeTracers_temperature
        self.plot_3d_field_top_bot_section(da,
                                           nameInTitle='temperature',
                                           prefix='Temp', units='deg C',
                                           vmin=-2.5, vmax=1.0)

    def plot_salinity(self):
        '''
        Plot a series of images of salinity at the sea surface or
        ice-ocean interface, sea floor and in an x-z section
        '''

        da = self.ds.timeMonthly_avg_activeTracers_salinity
        self.plot_3d_field_top_bot_section(da,
                                           nameInTitle='salinity',
                                           prefix='Salinity', units='PSU',
                                           vmin=33.8, vmax=34.7)

    def plot_potential_density(self):
        '''
        Plot a series of images of salinity at the sea surface or
        ice-ocean interface, sea floor and in an x-z section
        '''

        da = self.ds.timeMonthly_avg_potentialDensity
        self.plot_3d_field_top_bot_section(da,
                                           nameInTitle='potential density',
                                           prefix='PotRho', units='kg/m^3',
                                           vmin=1027., vmax=1028.)

    def plot_haney_number(self):
        '''
        Plot a series of images of the Haney number rx1 at the sea surface or
        ice-ocean interface, sea floor and in an x-z section
        '''

        ds = xarray.open_dataset('{}/haney.nc'.format(self.inFolder))

        self.plot_3d_field_top_bot_section(ds.haneyCell,
                                           nameInTitle='Haney Number (rx1)',
                                           prefix='Haney', units=None,
                                           vmin=0., vmax=8.)

    def plot_horiz_series(self, da, nameInTitle, prefix, oceanDomain,
                          units=None, vmin=None, vmax=None):
        '''
        Plot a series of image of a given variable

        Parameters
        ----------
        da : ``xarray.DataArray``
            The data array of time series to plot

        nameInTitle : str
            The name of the variable to use in the title and the progress bar

        prefix : str
            The nae of the variable to use in the subfolder and file prefix

        oceanDomain : bool
            True if the variable is for the full ocean, False if only for the
            cavity

        units : str, optional
            The units of the variable to be included in the title

        vmin, vmax : float, optional
            The minimum and maximum values for the colorbar
        '''

        nTime = self.ds.sizes['Time']
        widgets = ['plotting {}: '.format(nameInTitle),
                   progressbar.Percentage(), ' ',
                   progressbar.Bar(), ' ', progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets,
                                      maxval=nTime).start()

        for tIndex in range(nTime):
            self.update_date(tIndex)
            field = da.isel(Time=tIndex).values
            outFileName = '{}/{}/{}_{:04d}.png'.format(
                self.outFolder, prefix, prefix, tIndex+1)
            if units is None:
                title = nameInTitle
            else:
                title = '{} ({})'.format(nameInTitle, units)
            self._plot_horiz_field(field, title=title, outFileName=outFileName,
                                   oceanDomain=oceanDomain, vmin=vmin,
                                   vmax=vmax)
            bar.update(tIndex+1)
        bar.finish()

    def plot_3d_field_top_bot_section(self, da, nameInTitle, prefix,
                                      units=None, vmin=None, vmax=None):
        '''
        Plot a series of images of a given 3D variable showing the value
        at the top (sea surface or ice-ocean interface), sea floor and in an
        x-z section

        Parameters
        ----------
        da : ``xarray.DataArray``
            The data array of time series to plot

        nameInTitle : str
            The name of the variable to use in the title and the progress bar

        prefix : str
            The nae of the variable to use in the subfolder and file prefix

        units : str, optional
            The units of the variable to be included in the title

        vmin, vmax : float, optional
            The minimum and maximum values for the colorbar
        '''

        if vmin is None:
            vmin = da.min()
        if vmax is None:
            vmax = da.max()

        self.plot_horiz_series(da.isel(nVertLevels=0),
                               'top {}'.format(nameInTitle),
                               'top{}'.format(prefix), oceanDomain=True,
                               vmin=vmin, vmax=vmax)

        maxLevelCell = self.dsMesh.maxLevelCell-1

        daBot = xarray.DataArray(da)

        daBot.coords['verticalIndex'] = \
            ('nVertLevels',
             numpy.arange(daBot.sizes['nVertLevels']))

        # mask only the values with the right vertical index
        daBot = daBot.where(daBot.verticalIndex == maxLevelCell)

        # Each vertical layer has at most one non-NaN value so the "sum"
        # over the vertical is used to collapse the array in the vertical
        # dimension
        daBot = daBot.sum(dim='nVertLevels').where(maxLevelCell >= 0)

        self.plot_horiz_series(daBot,
                               'bot {}'.format(nameInTitle),
                               'bot{}'.format(prefix), oceanDomain=True,
                               vmin=vmin, vmax=vmax)

        daSection = da.isel(nCells=self.sectionCellIndices)

        nTime = self.ds.sizes['Time']
        widgets = ['plotting {} section: '.format(nameInTitle),
                   progressbar.Percentage(), ' ',
                   progressbar.Bar(), ' ', progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets,
                                      maxval=nTime).start()

        for tIndex in range(nTime):
            self.update_date(tIndex)
            mask = numpy.logical_not(self.sectionMask)
            field = numpy.ma.masked_array(daSection.isel(Time=tIndex).values.T,
                                          mask=mask)
            outFileName = '{}/section{}/section{}_{:04d}.png'.format(
                self.outFolder, prefix, prefix, tIndex+1)
            if units is None:
                title = nameInTitle
            else:
                title = '{} ({}) along section at y={:g} km'.format(
                    nameInTitle, units, 1e-3*self.sectionY)
            self._plot_vert_field(self.X, self.Z[tIndex, :, :],
                                  field, title=title,
                                  outFileName=outFileName,
                                  vmin=vmin, vmax=vmax)
            bar.update(tIndex+1)
        bar.finish()

    def images_to_movies(self, outFolder, framesPerSecond=30, extension='mp4',
                         overwrite=True):
        '''
        Convert all the image sequences into movies with ffmpeg
        '''
        try:
            os.makedirs('{}/logs'.format(outFolder))
        except OSError:
            pass

        framesPerSecond = '{}'.format(framesPerSecond)

        for fileName in sorted(glob.glob(
                '{}/*/*0001.png'.format(self.outFolder))):
            prefix = os.path.basename(fileName)[:-9]
            outFileName = '{}/{}.{}'.format(outFolder, prefix, extension)
            if not overwrite and os.path.exists(outFileName):
                continue

            imageFileTemplate = '{}/{}/{}_%04d.png'.format(self.outFolder,
                                                           prefix, prefix)
            logFileName = '{}/logs/{}.log'.format(outFolder, prefix)
            with open(logFileName, 'w') as logFile:
                args = ['ffmpeg', '-y', '-r', framesPerSecond,
                        '-i', imageFileTemplate, '-b:v', '32000k',
                        '-r', framesPerSecond, '-pix_fmt', 'yuv420p',
                        outFileName]
                print('running {}'.format(' '.join(args)))
                subprocess.check_call(args, stdout=logFile, stderr=logFile)

    def update_date(self, tIndex):
        xtime = self.ds.xtime_startMonthly.isel(Time=tIndex).values
        xtime = ''.join(str(xtime.astype('U'))).strip()
        year = xtime[0:4]
        month = xtime[5:7]
        self.date = '{}-{}'.format(year, month)

    def _plot_horiz_field(self, field, title, outFileName, oceanDomain=True,
                          vmin=None, vmax=None, figsize=[9, 3]):

        try:
            os.makedirs(os.path.dirname(outFileName))
        except OSError:
            pass

        if os.path.exists(outFileName):
            return

        if oceanDomain:
            localPatches = copy.copy(self.oceanPatches)
            localPatches.set_array(field[self.oceanMask])
        else:
            localPatches = copy.copy(self.cavityPatches)
            localPatches.set_array(field[self.cavityMask])

        localPatches.set_edgecolor('face')
        localPatches.set_clim(vmin=vmin, vmax=vmax)

        plt.figure(figsize=figsize)
        ax = plt.subplot('111')
        ax.add_collection(localPatches)
        plt.colorbar(localPatches)
        plt.axis([0, 500, 0, 1000])
        ax.set_aspect('equal')
        ax.autoscale(tight=True)
        plt.title('{} {}'.format(title, self.date))
        plt.tight_layout(pad=0.5)
        plt.savefig(outFileName)
        plt.close()

    def _plot_vert_field(self, inX, inZ, field, title, outFileName, vmin=None,
                         vmax=None, figsize=[9, 5]):
        try:
            os.makedirs(os.path.dirname(outFileName))
        except OSError:
            pass

        if os.path.exists(outFileName):
            return

        plt.figure(figsize=figsize)
        ax = plt.subplot('111')
        plt.pcolormesh(1e-3*inX, inZ, field, vmin=vmin, vmax=vmax,
                       cmap=self.cmap)
        plt.colorbar()
        ax.autoscale(tight=True)
        plt.ylim([numpy.amin(inZ), 20])
        plt.title('{} {}'.format(title, self.date))
        plt.tight_layout(pad=0.5)
        plt.savefig(outFileName)
        plt.close()

    def _compute_section_x_z(self):
        x = _interp_extrap_corner(self.dsMesh.xCell[self.sectionCellIndices])
        nx = len(x)
        nVertLevels = self.dsMesh.sizes['nVertLevels']
        nTime = self.ds.sizes['Time']
        self.X = numpy.zeros((nVertLevels+1, nx))
        for zIndex in range(nVertLevels+1):
            self.X[zIndex, :] = x

        self.sectionMask = numpy.zeros((nVertLevels, nx-1))
        for zIndex in range(nVertLevels):
            maxLevelCell = self.dsMesh.maxLevelCell.isel(
                nCells=self.sectionCellIndices) - 1
            self.sectionMask[zIndex, :] = maxLevelCell >= zIndex

        self.Z = numpy.zeros((nTime, nVertLevels+1, nx))
        for tIndex in range(nTime):
            layerThickness = self.ds.timeMonthly_avg_layerThickness.isel(
                Time=tIndex, nCells=self.sectionCellIndices).values
            bottomDepth = self.dsMesh.bottomDepth.isel(
                nCells=self.sectionCellIndices).values
            self.Z[tIndex, -1, :] = -_interp_extrap_corner(bottomDepth)
            for zIndex in range(nVertLevels-1, -1, -1):
                layerThicknessSection = _interp_extrap_corner(
                    layerThickness[:, zIndex])
                self.Z[tIndex, zIndex, :] = self.Z[tIndex, zIndex+1, :] + \
                    layerThicknessSection


def _compute_cell_patches(dsMesh, mask, cmap):
    patches = []
    nVerticesOnCell = dsMesh.nEdgesOnCell.values
    verticesOnCell = dsMesh.verticesOnCell.values - 1
    xVertex = dsMesh.xVertex.values
    yVertex = dsMesh.yVertex.values
    for iCell in range(dsMesh.sizes['nCells']):
        if(not mask[iCell]):
            continue
        nVert = nVerticesOnCell[iCell]
        vertexIndices = verticesOnCell[iCell, :nVert]
        vertices = numpy.zeros((nVert, 2))
        vertices[:, 0] = 1e-3*xVertex[vertexIndices]
        vertices[:, 1] = 1e-3*yVertex[vertexIndices]

        polygon = Polygon(vertices, True)
        patches.append(polygon)

    p = PatchCollection(patches, cmap=cmap, alpha=1.)

    return p


def _compute_section_cell_indices(y, dsMesh):
    xCell = dsMesh.xCell.values
    yCell = dsMesh.yCell.values
    xMin = numpy.amin(xCell)
    xMax = numpy.amax(xCell)
    xs = numpy.linspace(xMin, xMax, 10000)
    cellIndices = []
    for x in xs:
        distanceSquared = (x - xCell)**2 + (y-yCell)**2
        index = numpy.argmin(distanceSquared)
        if(len(cellIndices) == 0 or cellIndices[-1] != index):
            cellIndices.append(index)

    return numpy.array(cellIndices)


def _make_ferret_colormap():
    red = numpy.array([[0, 0.6],
                       [0.15, 1],
                       [0.35, 1],
                       [0.65, 0],
                       [0.8, 0],
                       [1, 0.75]])

    green = numpy.array([[0, 0],
                         [0.1, 0],
                         [0.35, 1],
                         [1, 0]])

    blue = numpy.array([[0, 0],
                        [0.5, 0],
                        [0.9, 0.9],
                        [1, 0.9]])

    colorCount = 21
    ferretColorList = numpy.ones((colorCount, 4), float)
    ferretColorList[:, 0] = numpy.interp(numpy.linspace(0, 1, colorCount),
                                         red[:, 0], red[:, 1])
    ferretColorList[:, 1] = numpy.interp(numpy.linspace(0, 1, colorCount),
                                         green[:, 0], green[:, 1])
    ferretColorList[:, 2] = numpy.interp(numpy.linspace(0, 1, colorCount),
                                         blue[:, 0], blue[:, 1])
    ferretColorList = ferretColorList[::-1, :]

    cmap = colors.LinearSegmentedColormap.from_list('ferret', ferretColorList,
                                                    N=255)
    return cmap


def _interp_extrap_corner(inField):
    '''Interpolate/extrapolate a 1D field from grid centers to grid corners'''

    outField = numpy.zeros(len(inField) + 1)
    outField[1:-1] = 0.5 * (inField[0:-1] + inField[1:])
    # extrapolate the ends
    outField[0] = 1.5 * inField[0] - 0.5 * inField[1]
    outField[-1] = 1.5 * inField[-1] - 0.5 * inField[-2]
    return outField
