# Redi_verification
# test functions for Redi tendancy terms
# Mark Petersen, LANL, Nov 2019

# followed example at:
# https://pythonhosted.org/algopy/symbolic_differentiation.html

import sympy as sp
x,y,z = sp.symbols('x y z')
fkx,fky,fkz = sp.symbols('fkx fky fkz')
Tkx,Tky,Tkz = sp.symbols('Tkx Tky Tkz')
Lx,Ly,Lz = sp.symbols('Lx Ly Lz')
A,rho0,alpha = sp.symbols('A rho0 alpha')
T0,T1,T0eos= sp.symbols('T0 T1 T0eos')

# define the passive debug tracer
def phiDef(x,y,z):
    return A \
       *sp.sin(x*fkx*2*sp.pi/Lx) \
       *sp.sin(y*fky*2*sp.pi/Ly) \
       *sp.sin(z*fkz*2*sp.pi/Lz)
phi = phiDef(x,y,z)
    
# Define Temperature, with sines
#def TDef(x,y,z):
#    return A \
#       *sp.sin(x*Tkx*2*sp.pi/Lx) \
#       *sp.sin(y*Tky*2*sp.pi/Ly) \
#       *sp.sin(z*Tkz*2*sp.pi/Lz)

# Define Temperature, linear
def TDef(x,y,z):
    return T0 + T1 \
       *x*Tkx/Lx \
       *y*Tky/Ly \
       *z*Tkz/Lz
T = TDef(x,y,z)

# density using the linear equation of state
def rhoDef(x,y,z):
    return rho0 - alpha * (T - T0eos)
rho = rhoDef(x,y,z)

# isopycnal slope in the x direction
def Sx(x,y,z):
    return rho.diff(x)/rho.diff(z)

# isopycnal slope in the y direction
def Sy(x,y,z):
    return rho.diff(y)/rho.diff(z)

# Term 1 div(grad(phi))
def term1(x,y,z):
    return \
    phi.diff(x).diff(x) + phi.diff(y).diff(y)

# Term 2: div( S dphi/dz)
def term2(x,y,z):
    return \
    (Sx(x,y,z)*phi.diff(z)).diff(x) + \
    (Sy(x,y,z)*phi.diff(z)).diff(y)

# Term 3: d/dz( S.grad(phi))
def term3(x,y,z):
    return \
    (Sx(x,y,z)*phi.diff(x)).diff(z) + \
    (Sy(x,y,z)*phi.diff(y)).diff(z)

print('term 1 = ',term1(x,y,z))
print('term 2 = ',term2(x,y,z))
print('term 3 = ',term3(x,y,z))

