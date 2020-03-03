import numpy as np

class GeoGrid:
    def __init__(self, lon: np.ndarray, lat: np.ndarray):
        """
        Constructor.
        :param lon: longitude of the grid in radians, as numpy array
        :param lat: latitude of the grid in radians, as numpy array
        """
        self.lon = lon
        self.lat = lat
        self.ncells = len(lon)




    # '''
    # A class that defines the structure, location, extent, and resolution of a geographic grid.
    # The grid is not the same as a geospatial raster, but is related in that, while a raster numbers vertical cells
    # starting from the top of the raster, the grid cells are numbered from the bottom. That is, a raster is oriented
    # like a raster of pixels, while the geographic grid is oriented like a regular Cartesian grid of cells. The
    # data in the grid is contained in a two-dimensional NumPy array. Because of this, the grid cell is indexed like
    # a Fortran array (column major indexing, i.e. i=column, j=row).
    # '''
    # def __init__(self, lon, lat, nlon, nlat, cellsize, defaultValue=0.0):
    #     '''
    #     Constructor.
    #     :param lon: Lower-left longitude of the grid in decimal degrees.
    #     :param lat: Lower-left latitude of the grid in decimal degrees.
    #     :param nlon: The number of cells in longitude.
    #     :param nlat: The number of cells in latitude.
    #     :param cellsize: The size of a cell in the grid.
    #     '''
    #     self.lon = lon
    #     self.lat = lat
    #     self.nlon = nlon
    #     self.nlat = nlat
    #     self.cellsize = cellsize
    #     self.defaultValue = defaultValue
    #     self.grid = np.zeros([nlat,nlon],dtype=np.float64)
    #     self.bounds = [self.lon, self.lon + self.nlon*self.cellsize,
    #         self.lat, self.lat + self.nlat*self.cellsize]
    #
    #
    # def put(self,i,j,v):
    #     if self.indexInside(i,j):
    #         self.grid[self.nlat-j-1,i]=v
    #
    # def getByIndex(self,i,j):
    #     if self.indexInside(i,j):
    #         return self.grid[self.nlat-j-1,i]
    #     else:
    #         return self.defaultValue
    #
    # def getByCoordinate(self,lon,lat):
    #     if self.coordinateInside(lon,lat):
    #         index = self.getIndex(lon,lat)
    #         return self.getByIndex(index[0],index[1])
    #     else:
    #         return self.defaultValue
    #
    # def clear(self):
    #     self.grid.fill(0.0)
    #
    # def indexInside(self,i,j):
    #     if i>=0 and i<self.nlon and j>=0 and j<self.nlat:
    #         return True
    #     else:
    #         return False
    #
    # def coordinateInside(self,lon,lat):
    #     if lon>=self.bounds[0] and lon<=self.bounds[1] and lat>=self.bounds[2] and lat<=self.bounds[3]:
    #         return True
    #     else:
    #         return False
    #
    # def getOrigin(self):
    #     return [self.lon,self.lat]
    #
    # def getCenter(self,i,j):
    #     clon = self.lon + (i+0.5)*self.cellsize
    #     clat = self.lat + (j+0.5)*self.cellsize
    #     return [clon,clat]
    #
    # def getIndex(self,lon,lat):
    #     i = int((lon-self.lon)/self.cellsize)
    #     j = int((lat-self.lat)/self.cellsize)
    #     return [i,j]
    #
