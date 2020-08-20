from structures.geogrid import GeoGrid

def test_geogrid():
    lon = -106.0
    lat = 35
    nlon = 8
    nlat = 4
    cellsize = 1.0
    defaultValue = -1.0
    grid = GeoGrid(lon,lat,nlon,nlat,cellsize,defaultValue = defaultValue)
    assert grid.lon == lon
    assert grid.lat == lat
    assert grid.nlon == nlon
    assert grid.nlat == nlat
    assert grid.cellsize == cellsize
    assert defaultValue == defaultValue

    l = int(nlat/2)
    k = int(nlon/2)
    for j in range(0,l):
        for i in range(0,k):
            grid.put(i,j,1.0)
        for i in range(k,nlon):
            grid.put(i,j,2.0)
    for j in range(l,nlat):
        for i in range(0,k):
            grid.put(i,j,3.0)
        for i in range(k,nlon):
            grid.put(i,j,4.0)

    for j in range(0,l):
        for i in range(0,k):
            assert grid.getByIndex(i,j) == 1.0
        for i in range(k,nlon):
            assert grid.getByIndex(i,j) == 2.0
    for j in range(l,nlat):
        for i in range(0,k):
            assert grid.getByIndex(i,j) == 3.0
        for i in range(k,nlon):
            assert grid.getByIndex(i,j) == 4.0

    testcell = [3,3]
    center = grid.getCenter(testcell[0],testcell[1])
    centerx = lon + (testcell[0]+0.5)*cellsize
    centery = lat + (testcell[1]+0.5)*cellsize
    assert center[0] == centerx
    assert center[1] == centery

    index = grid.getIndex(centerx,centery)
    assert index[0] == testcell[0]
    assert index[1] == testcell[1]

    value = grid.getByIndex(testcell[0],testcell[1])
    testcoords = grid.getCenter(testcell[0],testcell[1])
    valuec = grid.getByCoordinate(testcoords[0],testcoords[1])
    assert value == valuec

    origin = grid.getOrigin()
    assert origin[0] == lon
    assert origin[1] == lat

    bounds = grid.bounds
    assert bounds[0] == lon
    assert bounds[1] == lon + nlon*cellsize
    assert bounds[2] == lat
    assert bounds[3] == lat + nlat*cellsize

    assert grid.indexInside(-1,l) == False
    assert grid.indexInside(k,l) == True
    assert grid.indexInside(nlon,l) == False
    assert grid.indexInside(k,-1) == False
    assert grid.indexInside(k,l) == True
    assert grid.indexInside(k,nlat) == False

    assert grid.coordinateInside(bounds[0]+cellsize,bounds[2]+cellsize) == True
    assert grid.coordinateInside(bounds[0]-cellsize,bounds[2]+cellsize) == False
    assert grid.coordinateInside(bounds[0]+cellsize,bounds[2]-cellsize) == False

    assert grid.coordinateInside(bounds[1]-cellsize,bounds[2]+cellsize) == True
    assert grid.coordinateInside(bounds[1]-cellsize,bounds[2]-cellsize) == False
    assert grid.coordinateInside(bounds[1]+cellsize,bounds[2]+cellsize) == False

    assert grid.coordinateInside(bounds[0]+cellsize,bounds[3]-cellsize) == True
    assert grid.coordinateInside(bounds[0]+cellsize,bounds[3]+cellsize) == False
    assert grid.coordinateInside(bounds[0]-cellsize,bounds[3]+cellsize) == False

    assert grid.coordinateInside(bounds[1]-cellsize,bounds[3]-cellsize) == True
    assert grid.coordinateInside(bounds[1]-cellsize,bounds[3]+cellsize) == False
    assert grid.coordinateInside(bounds[1]+cellsize,bounds[3]-cellsize) == False

    grid.clear()
    for j in range(0,nlat):
        for i in range(0,nlon):
            assert grid.getByIndex(i,j) == 0.0

