import numpy as np
import matplotlib.pyplot as plt

def example():
    x,y = np.linspace(-1,1,2), np.linspace(-1,1,2)
    A, B = np.zeros((2,2)), np.zeros((2,2))
    A[0,0]=1
    B[0,0]=-1
    A[0,1]=1
    B[0,1]=1
    A[1,0]=-1
    B[1,0]=-1
    A[1,1]=-1
    B[1,1]=1

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot the streamlines.
    ax.streamplot(x,y,A,B)

    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    ax.set_aspect('equal')
    plt.show()

if __name__=='__main__':
    example()

