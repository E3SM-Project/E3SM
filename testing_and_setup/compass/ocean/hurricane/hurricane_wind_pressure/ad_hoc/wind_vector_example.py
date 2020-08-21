import sys
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.patches import Circle
import math

def W(x, y):
    """Return the wind vector given a wind speed."""
    r = np.sqrt(x*x+y*y)
    v = V(r)
    if r>0:
        costheta = x/r
        sintheta = y/r
        return [-sintheta*v,costheta*v]
    else:
        return [0,0]

def V(r):
    return 2*r*r*np.exp(-r)

def example(n):
    # Grid of x, y points
    nx, ny = n, n
    x = np.linspace(-2, 2, nx)
    y = np.linspace(-2, 2, ny)

    # Wind field vector components U,V
    U, V = np.zeros((ny, nx)), np.zeros((ny, nx))
    for j in range(ny-1,-1,-1):
        for i in range(0,nx):
            vv = W(x[i],y[j])
            U[j,i]=vv[0]
            V[j,i]=vv[1]

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)

    # Plot the streamlines.
    ax1.streamplot(x, y, U, V, color=np.sqrt(U*U+V*V), cmap='Spectral')
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$')
    ax1.set_xlim(-2,2)
    ax1.set_ylim(-2,2)
    ax1.set_aspect('equal')
    plt.title('Tangential Wind Vectors')
    plt.show()

if __name__=='__main__':
    example(8)

