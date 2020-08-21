from geopy.distance import geodesic
from utils.gis import geodistkm

def test_gis():
    albuquerque = [35.0844, -106.6504] #(lat,lon)
    los_alamos = [35.8800, -106.3031] #(lat,lon)

    result1 = geodesic(albuquerque,los_alamos).km
    result2 = geodistkm(albuquerque[1],albuquerque[0],los_alamos[1],los_alamos[0])

    assert result1 == result2
