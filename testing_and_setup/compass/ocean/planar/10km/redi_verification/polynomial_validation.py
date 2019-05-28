# Redi_verification
# test functions for Redi tendancy terms
# Mark Petersen, LANL, Nov 2019

# followed example at:
# https://pythonhosted.org/algopy/symbolic_differentiation.html

import sympy as sp
import numpy as np
x,y,z = sp.symbols('x y z')
fkx,fky,fkz = sp.symbols('fkx fky fkz')
Tkx,Tky,Tkz = sp.symbols('Tkx Tky Tkz')
Lx,Ly,Lz = sp.symbols('Lx Ly Lz')
Tx,Ty,Tz = sp.symbols('Tx Ty Tz')
A,rho0,alpha = sp.symbols('A rho0 alpha')
T0,T1,T0eos= sp.symbols('T0 T1 T0eos')
P0,T1,T0eos= sp.symbols('P0 T1 T0eos')

# define the passive debug tracer
def PDef(x,y,z):
    return P0*(x**px * y**py * z**pz)

# define the passive debug tracer
def TDef(x,y,z):
    return Tx*x**pTx + Ty*y**pTy + Tz*z**pTz

# density using the linear equation of state
def rhoDef(x,y,z):
    return rho0 - alpha * (T - T0eos)

# isopycnal slope in the x direction
def SxDef(x,y,z):
    return rho.diff(x)/rho.diff(z)

# isopycnal slope in the y direction
def SyDef(x,y,z):
    return rho.diff(y)/rho.diff(z)

# Term 1 div(grad(P))
def term1Def(x,y,z):
    return \
    P.diff(x).diff(x) + P.diff(y).diff(y)

# Term 2: div( S dP/dz)
def term2Def(x,y,z):
    return \
    (Sx*P.diff(z)).diff(x) + \
    (Sy*P.diff(z)).diff(y)

# Term 3: d/dz( S.grad(P))
def term3Def(x,y,z):
    return \
    (Sx*P.diff(x)).diff(z) + \
    (Sy*P.diff(y)).diff(z)

def printColumns(a):
    itemLen = np.zeros(20,'int')
    for row in a:
        for i in range(len(row)):
            itemLen[i] = max(itemLen[i],len(row[i]))
    for row in a:
        for i in range(len(row)):
            print(row[i].ljust(itemLen[i]),' ',end='')
        print()

def sympyList2Str(a):
    for i in range(len(a)):
       a[i] = str(a[i])
    return a

#print('term 1 = ',term1(x,y,z))
#print('term 2 = ',term2(x,y,z))
#print('term 3 = ',term3(x,y,z))

py=0
pTy=0
a=[['powers','tracer P','Temperature T','Slope Sx','term 1','term 2','term 3']]
for pz in [1, 2]:
    for px in [1, 2]:
        for pTz in [1, 2]:
            for pTx in [1, 2]:
                P = PDef(x,y,z)
                T = TDef(x,y,z)
                rho = rhoDef(x,y,z)
                Sx = SxDef(x,y,z)
                Sy = SyDef(x,y,z)
                term1 = term1Def(x,y,z)
                term2 = term2Def(x,y,z)
                term3 = term3Def(x,y,z)
                const = 'px='+str(px)+';py='+str(py)+';pz='+str(pz)+';'+ \
                        'pTx='+str(pTx)+';pTy='+str(pTy)+';pTz='+str(pTz)+';'
                a.append(sympyList2Str([const,P,T,Sx,term1,term2,term3]))

printColumns(a)
