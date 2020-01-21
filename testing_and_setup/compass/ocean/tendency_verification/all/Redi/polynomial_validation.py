#!/usr/bin/env python
'''
This script creates test functions for verification of Redi tendancy terms.
Mark Petersen, LANL, Nov 2019

followed example at:
https://pythonhosted.org/algopy/symbolic_differentiation.html
'''
import sympy as sp
import numpy as np

# define algabraic variables:
x, y, z = sp.symbols('x y z')
fkx, fky, fkz = sp.symbols('fkx fky fkz')
x0, y0, z0 = sp.symbols('x0 y0 z0')
Tkx, Tky, Tkz = sp.symbols('Tkx Tky Tkz')
Lx, Ly, Lz = sp.symbols('Lx Ly Lz')
Tx, Ty, Tz = sp.symbols('Tx Ty Tz')
A, rho0, alpha = sp.symbols('A rho0 alpha')
T0, T1, T0eos = sp.symbols('T0 T1 T0eos')
P0, T1, T0eos = sp.symbols('P0 T1 T0eos')


def PDef(x, y, z):
    # define the passive debug tracer
    return P0 * ((x + x0)**px * y**py * z**pz)


def TDef(x, y, z):
    # define the passive debug tracer
    return Tx * (x + x0)**pTx + Ty * y**pTy + Tz * z**pTz


def rhoDef(x, y, z):
    # density using the linear equation of state
    return rho0 - alpha * (T - T0eos)


def SxDef(x, y, z):
    # isopycnal slope in the x direction
    return -rho.diff(x) / rho.diff(z)


def SyDef(x, y, z):
    # isopycnal slope in the y direction
    return -rho.diff(y) / rho.diff(z)


def term1Def(x, y, z):
    # Term 1 div(grad(P))
    return \
        P.diff(x).diff(x) + P.diff(y).diff(y)


def term2Def(x, y, z):
    # Term 2: div( S dP/dz)
    return \
        (Sx * P.diff(z)).diff(x) + \
        (Sy * P.diff(z)).diff(y)


def term3Def(x, y, z):
    # Term 3: d/dz( S.grad(P))
    return \
        (Sx * P.diff(x)).diff(z) + \
        (Sy * P.diff(y)).diff(z)


def printColumns(a):
    # Print columns neatly, with enough white space to fit widest entry in each
    # column.
    itemLen = np.zeros(20, 'int')
    for row in a:
        for i in range(len(row)):
            itemLen[i] = max(itemLen[i], len(row[i]))
    for row in a:
        for i in range(len(row)):
            print(row[i].ljust(itemLen[i]), ' ', end='')
        print()


def sympyList2Str(a):
    for i in range(len(a)):
        a[i] = str(a[i])
    return a


j = 1
py = 0
pTy = 0
a = [['#', 'Temperature T', 'tracer P', 'Slope Sx', 'term 1', 'term 2', 'term 3']]
for pTz in [1, 2]:
    for pTx in [1, 2]:
        for pz in [1, 2]:
            for px in [1, 2]:
                P = PDef(x, y, z)
                T = TDef(x, y, z)
                rho = rhoDef(x, y, z)
                Sx = SxDef(x, y, z)
                Sy = SyDef(x, y, z)
                term1 = term1Def(x, y, z)
                term2 = term2Def(x, y, z)
                term3 = term3Def(x, y, z)
                const = 'px=' + str(px) + ';py=' + str(py) + ';pz=' + str(pz) + ';' + \
                        'pTx=' + str(pTx) + ';pTy=' + str(pTy) + ';pTz=' + str(pTz) + ';'
                a.append(sympyList2Str([j, T, P, Sx, term1, term2, term3]))
                j += 1

printColumns(a)
