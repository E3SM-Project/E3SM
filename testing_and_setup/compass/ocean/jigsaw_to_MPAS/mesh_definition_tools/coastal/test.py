import numpy as np

x = np.linspace(-1,1,5)
y = np.linspace(-1,1,3)

nx = len(x)
ny = len(y)

X,Y = np.meshgrid(x,y)
XY = np.vstack([X.ravel(), Y.ravel()]).T
print X
print ""
print Y
print ""
print XY
print ""

X2,Y2 = np.vsplit(XY.T,2)
print X2
print ""
print Y2
print ""

X3 = np.reshape(X2,(ny,nx))
Y3 = np.reshape(Y2,(ny,nx))
print X3
print ""
print Y3
print ""

