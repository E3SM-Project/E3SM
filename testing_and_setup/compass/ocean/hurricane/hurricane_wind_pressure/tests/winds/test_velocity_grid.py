import numpy as np
from structures.geogrid import GeoGrid
from profile_model.radialprofiles import HollandWindSpeedProfile
from winds.parameters import Parameters
from winds.velocities import Velocities
import matplotlib.pyplot as plt

def test_velocity_grid():
    # Grid of x, y points
    n = 50
    nr = 200
    rmax = 40
    cmin, cmax = -200 , 200
    cellsize = (cmax-cmin)/n
    x = np.linspace(cmin, cmax, n)
    y = np.linspace(cmin, cmax, n)
    U = GeoGrid(cmin,cmin,n,n,cellsize)
    V = GeoGrid(cmin,cmin,n,n,cellsize)

    params = Parameters()
    b = 1.4
    hc = [0,0]
    vf = [0,10]
    deltap = 100
    coriol = False
    profile = HollandWindSpeedProfile(nr,2*cmax,rmax,deltap,params.rho,params.getCoriolisMid(),b,coriolis=coriol)
    vels = Velocities(vf[0],vf[1],profile.getVmax())
    for j in range(0,n):
        for i in range(0,n):
            pt = U.getCenter(i,j)
            r = np.sqrt(pow(pt[0]-hc[0],2)+pow(pt[1]-hc[1],2))
            vg = profile.getValue(r)
            vv = vels.compute_wind_vector(vg,pt[0],pt[1])
            U.put(i,j,vv[0])
            V.put(i,j,vv[1])

    assert True # If we made it to here.

    fig = plt.figure()

    ax = fig.add_subplot(131)
    ax.plot(profile.rvals, profile.profile)

    ax.set(xlabel='r (km)', ylabel='wind speed (km/hr)',
           title='Radial Wind')
    ax1 = fig.add_subplot(133)

    # Plot the streamlines.
    # Matplotlib assume an ordinary row ordering, so the rows must be reversed before plotting.
    Ug = U.grid
    Vg = V.grid
    Uplt = np.zeros([n,n])
    Vplt = np.zeros([n,n])
    for j in range(0,n):
        jp = n-j-1
        for i in range(0,n):
            Uplt[jp,i]=Ug[j,i]
            Vplt[jp,i]=Vg[j,i]

    Vmag = np.sqrt(Ug*Ug+Vg*Vg)
    ax1.streamplot(x, y, Uplt, Vplt, color=Vmag, cmap='Spectral')
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')
    ax1.set_xlim(cmin,cmax)
    ax1.set_ylim(cmin,cmax)
    ax1.set_aspect('equal')
    plt.title('Wind Vectors')

    plt.show()

