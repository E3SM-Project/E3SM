import pytest
from hurricane_model.hurricane import Hurricane

def test_hurricane():
    center = [1.0,2.0]  # Position of the eye (lon,lat) in decimal degrees.
    extent = 100.0  # The maximum extent of the hurricane in kilometers.
    vforward = [3.0, 4.0]  # Forward velocity [ve, vn] in km/hr.
    pcentral = 200.0  # Central pressure in millibars.
    deltap = 50.0  # Pressure difference in millibars.
    vmax = 15.0  # The maximum gradient wind speed in km/hr.
    b = 1.2  # The Holland parameter, conventionally in the range [0.5,2.5].

    hurricane = Hurricane(center,extent)
    hurricane.setVForward(vforward[0],vforward[1])
    hurricane.setPCentral(pcentral)
    hurricane.setDeltaP(deltap)
    hurricane.setVMax(vmax)
    hurricane.setB(b)

    assert hurricane.center == center
    assert hurricane.extent == extent
    assert hurricane.vforward == vforward
    assert hurricane.pcentral == pcentral
    assert hurricane.deltap == deltap
    assert hurricane.vmax == vmax
    assert hurricane.b == b
