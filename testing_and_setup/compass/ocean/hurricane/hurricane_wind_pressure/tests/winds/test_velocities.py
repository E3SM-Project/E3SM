from winds.velocities import Velocities
import math

def test_velocities():
    # Forward velocity in km/hr.
    vfe = -1.0 # Eastward .
    vfn = 0.0 # Northward.
    vg = 1.0 # Tangential gradient wind speed in km/hr.

    veloc = Velocities(vfe,vfn)
    r = 1.0 # Unit circle about the origin.
    np = 360
    dtheta = 2*math.pi/np
    with open('test_velocities_out.csv','wt') as out:
        out.write('x,y,vx,vy,r,theta_degrees\n')
        for i in range(0,np):
            theta = i*dtheta
            degrees = 180.0*theta/math.pi
            x = r*math.cos(theta)
            y = r*math.sin(theta)
            v = veloc.compute_wind_vector(vg,x,y)
            out.write(str(x)+','+str(y)+','+str(v[0])+','+str(v[1])+','+str(r)+','+str(degrees)+'\n')

