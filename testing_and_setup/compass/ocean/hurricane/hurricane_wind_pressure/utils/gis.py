from geopy.distance import geodesic

def geodistkm(x1,y1,x2,y2):
    '''
    Returns the geodesic distance in km given two pairs of (lon, lat) coordinates.
    Note: Because it uses geopy, the coordinate order is reversed to (lat,lon)
    before calling the geopy function.
    :param x1: lon of the first point.
    :param y1: lat of the first point.
    :param x2: lon of the second point.
    :param y2: lat of the second point.
    :return: Geodesic distance between the two points in km.
    '''
    return geodesic((y1,x1),(y2,x2)).km
