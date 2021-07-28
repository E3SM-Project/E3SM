from netCDF4 import Dataset
import numpy as np
from scipy.special import sph_harm
import math
from math import sin, cos, tan, pi, fabs, pow, sqrt, factorial

#-------------------------------------------------------------------------------

def legendre_polynomial(z, m, l):

    c0   =   0.0
    p5   =   0.5
    c1   =   1.0
    c1p5 =   1.5
    c2   =   2.0
    c3   =   3.0
    c4   =   4.0
    c5   =   5.0
    c6   =   6.0
    c7   =   7.0
    c8   =   8.0
    c9   =   9.0
    c12  =  12.0
    c14  =  14.0
    c15  =  15.0
    c18  =  18.0
    c20  =  20.0
    c21  =  21.0
    c24  =  24.0
    c28  =  28.0
    c30  =  30.0
    c35  =  35.0
    c42  =  42.0
    c45  =  45.0
    c56  =  56.0
    c60  =  60.0
    c63  =  63.0
    c70  =  70.0
    c84  =  84.0
    c90  =  90.0
    c105 = 105.0
    c108 = 108.0
    c140 = 140.0
    c168 = 168.0
    c252 = 252.0
    c280 = 280.0
    c420 = 420.0
    c504 = 504.0
    c756 = 756.0
    c945 = 945.0

    s = sqrt(c1 - z*z)

    if (l == 0 and m == 0):

        poly = c1

    elif (l == 1 and m == 0):

        poly = z

    elif (l == 1 and m == 1):

        poly = s

    elif (l == 2 and m == 0):

        poly = p5 * (c3 * pow(z,2) - c1)

    elif (l == 2 and m == 1):

        poly = c3 * z * s

    elif (l == 2 and m == 2):

        poly = c3 * (c1 - pow(z,2))

    elif (l == 3 and m == 0):

        poly = p5 * z * (c5 * pow(z,2) - c3)

    elif (l == 3 and m == 1):

        poly = c1p5 * (c5 * pow(z,2) - c1) * s

    elif (l == 3 and m == 2):

        poly = c15 * z * (c1 - pow(z,2))

    elif (l == 3 and m == 3):

        poly = c15 * pow(s,3)

    elif (l == 4 and m == 0):

        poly = (c1/c8) * (c35 * pow(z,4) - c30 * pow(z,2) + c3)

    elif (l == 4 and m == 1):

        poly = (c5/c2) * (c7 * pow(z,3) - c3 * z) * s

    elif (l == 4 and m == 2):

        poly = (c15/c2) * (c7 * pow(z,2) - c1) * (c1 - pow(z,2))

    elif (l == 4 and m == 3):

        poly = c105 * z * pow(s,3)

    elif (l == 4 and m == 4):

        poly = c105 * pow(s,4)

    elif (l == 5 and m == 0):

        poly = (c1/c8) * z * (c63 * pow(z,4) - c70 * pow(z,2) + c15)

    elif (l == 5 and m == 1):

        poly = (c15/c8) * s * (c21 * pow(z,4) - c14 * pow(z,2) + c1)

    elif (l == 5 and m == 2):

        poly = (c105/c2) * z * (c1 - pow(z,2)) * (c3 * pow(z,2) - c1)

    elif (l == 5 and m == 3):

        poly = (c105/c2) * pow(s,3) * (c9 * pow(z,2) - c1)

    elif (l == 5 and m == 4):

        poly = c945 * z * pow(s,4)

    elif (l == 5 and m == 5):

        poly = c945 * pow(s,5)

    else:

        raise RuntimeError("legendre_polynomial: Unknown input")

    return poly

#-------------------------------------------------------------------------------

def dlegendre_polynomial_dz(z, m, l):

    c0   =   0.0
    p5   =   0.5
    c1   =   1.0
    c1p5 =   1.5
    c2   =   2.0
    c3   =   3.0
    c4   =   4.0
    c5   =   5.0
    c6   =   6.0
    c7   =   7.0
    c8   =   8.0
    c9   =   9.0
    c12  =  12.0
    c14  =  14.0
    c15  =  15.0
    c18  =  18.0
    c20  =  20.0
    c21  =  21.0
    c24  =  24.0
    c28  =  28.0
    c30  =  30.0
    c35  =  35.0
    c42  =  42.0
    c45  =  45.0
    c56  =  56.0
    c60  =  60.0
    c63  =  63.0
    c70  =  70.0
    c84  =  84.0
    c90  =  90.0
    c105 = 105.0
    c108 = 108.0
    c140 = 140.0
    c168 = 168.0
    c252 = 252.0
    c280 = 280.0
    c420 = 420.0
    c504 = 504.0
    c756 = 756.0
    c945 = 945.0

    s = sqrt(c1 - pow(z,2))
    if (s == 0):
        return 0.0
    dsdz = -z / s

    if (l == 0 and m == 0):

        #poly = c1
        poly = c0

    elif (l == 1 and m == 0):

        #poly = z
        poly = c1

    elif (l == 1 and m == 1):

        #poly = s
        poly = dsdz

    elif (l == 2 and m == 0):

        #poly = p5 * (c3 * pow(z,2) - c1)
        poly = c3 * z

    elif (l == 2 and m == 1):

        #poly = c3 * z * s
        poly = c3 * z * dsdz + c3 * s

    elif (l == 2 and m == 2):

        #poly = c3 * (c1 - pow(z,2))
        poly = -c6 * z

    elif (l == 3 and m == 0):

        #poly = p5 * z * (c5 * pow(z,2) - c3)
        poly = p5 * (c15 * z - c3)

    elif (l == 3 and m == 1):

        #poly = c1p5 * (c5 * pow(z,2) - c1) * s
        poly = c15 * z * s + c1p5 * (c5 * pow(z,2) - c1) * dsdz

    elif (l == 3 and m == 2):

        #poly = c15 * z * (c1 - pow(z,2))
        poly = c15 * (c1 - c3 * pow(z,2))

    elif (l == 3 and m == 3):

        #poly = c15 * pow(s,3)
        poly = c45 * pow(s,2) * dsdz

    elif (l == 4 and m == 0):

        #poly = (c1/c8) * (c35 * pow(z,4) - c30 * pow(z,2) + c3)
        poly = (c1/c8) * (c140 * pow(z,3) - c60 * z)

    elif (l == 4 and m == 1):

        #poly = (c5/c2) * (c7 * pow(z,3) - c3 * z) * s
        poly = (c5/c2) * ((c21 * pow(z,2) - c3) * s + (c7 * pow(z,3) - 3 * z) * dsdz)

    elif (l == 4 and m == 2):

        #poly = (c15/c2) * (c7 * pow(z,2) - c1) * (c1 - pow(z,2))
        poly = (c15/c2) * (c14 * z * (c1 - pow(z,2)) - c2 * z * (c7 * pow(z,2) - c1))

    elif (l == 4 and m == 3):

        #poly = c105 * z * pow(s,3)
        poly = c105 * (pow(s,3) + c3 * z * pow(s,2) * dsdz)

    elif (l == 4 and m == 4):

        #poly = c105 * pow(s,4)
        poly = c420 * pow(s,3) * dsdz

    elif (l == 5 and m == 0):

        #poly = (c1/c8) * z * (c63 * pow(z,4) - c70 * pow(z,2) + c15)
        poly = (c1/c8) * ((c63 * pow(z,4) - c70 * pow(z,2) + c15) + z * (c252 * pow(z,3) - c140 * z))

    elif (l == 5 and m == 1):

        #poly = (c15/c8) * s * (c21 * pow(z,4) - c14 * pow(z,2) + c1)
        poly = (c15/c8) * (s * (c84 * pow(z,3) - c28 * z) + dsdz * (c21 * pow(z,4) - c14 * pow(z,2) + 1))

    elif (l == 5 and m == 2):

        #poly = (c105/c2) * z * (c1 - pow(z,2)) * (c3 * pow(z,2) - c1)
        poly = (c105/c2) * ((c1 - pow(z,2)) * (c3 * pow(z,2) - c1) - c2 * pow(z,2) * (c3 * pow(z,2) - c1) + c6 * pow(z,2) * (c1 - pow(z,2)))

    elif (l == 5 and m == 3):

        #poly = (c105/c2) * pow(s,3) * (c9 * pow(z,2) - c1)
        poly = (c105/c2) * (c3 * pow(s,2) * dsdz * (c9 * pow(z,2) - c1) + c18 * pow(s,3) * z)

    elif (l == 5 and m == 4):

        #poly = c945 * z * pow(s,4)
        poly = c945 * (pow(s,4) + c4 * z * pow(s,3) * dsdz)

    elif (l == 5 and m == 5):

        #poly = c945 * pow(s,5)
        poly = c945 * (c5 * pow(s,4) * dsdz)

    else:

        raise RuntimeError("legendre_polynomial: Unknown input")

    return poly

#-------------------------------------------------------------------------------

def d2legendre_polynomial_dz2(z, m, l):

    c0   =   0.0
    p5   =   0.5
    c1   =   1.0
    c1p5 =   1.5
    c2   =   2.0
    c3   =   3.0
    c4   =   4.0
    c5   =   5.0
    c6   =   6.0
    c7   =   7.0
    c8   =   8.0
    c9   =   9.0
    c12  =  12.0
    c14  =  14.0
    c15  =  15.0
    c18  =  18.0
    c20  =  20.0
    c21  =  21.0
    c24  =  24.0
    c28  =  28.0
    c30  =  30.0
    c35  =  35.0
    c42  =  42.0
    c45  =  45.0
    c56  =  56.0
    c60  =  60.0
    c63  =  63.0
    c70  =  70.0
    c84  =  84.0
    c90  =  90.0
    c105 = 105.0
    c108 = 108.0
    c140 = 140.0
    c168 = 168.0
    c252 = 252.0
    c280 = 280.0
    c420 = 420.0
    c504 = 504.0
    c756 = 756.0
    c945 = 945.0

    s = sqrt(c1 - pow(z,2))
    if (s == 0):
        return 0.0

    dsdz = -z / s
    d2sdz2 = -c1 / s + (z * dsdz) / pow(s,2)

    if (l == 0 and m == 0):

        #poly = c1
        poly = c0

    elif (l == 1 and m == 0):

        #poly = z
        poly = c0

    elif (l == 1 and m == 1):

        #poly = s
        poly = d2sdz2

    elif (l == 2 and m == 0):

        #poly = p5 * (c3 * pow(z,2) - c1)
        poly = c3

    elif (l == 2 and m == 1):

        #poly = c3 * z * s
        poly = c3 * z * d2sdz2 + c6 * dsdz

    elif (l == 2 and m == 2):

        #poly = c3 * (c1 - pow(z,2))
        poly = -c6

    elif (l == 3 and m == 0):

        #poly = p5 * z * (c5 * pow(z,2) - c3)
        poly = c15 / 2

    elif (l == 3 and m == 1):

        #poly = c1p5 * (c5 * pow(z,2) - c1) * s
        poly = c15 * s + c30 * z * dsdz + c1p5 * (c5 * pow(z,2) - c1) * d2sdz2

    elif (l == 3 and m == 2):

        #poly = c15 * z * (c1 - pow(z,2))
        poly = -c90 * z

    elif (l == 3 and m == 3):

        #poly = c15 * pow(s,3)
        poly = c45 * (c2 * s * dsdz + pow(s,2) * d2sdz2)

    elif (l == 4 and m == 0):

        #poly = (c1/c8) * (c35 * pow(z,4) - c30 * pow(z,2) + c3)
        poly = (c1/c8) * (c420 * pow(z,2) - c60)

    elif (l == 4 and m == 1):

        #poly = (c5/c2) * (c7 * pow(z,3) - c3 * z) * s
        poly = (c5/c2) * (c42 * z * s + c2 * (c21 * pow(z,2) - c3) * dsdz + (c7 * pow(z,3) - c3 * z) * d2sdz2)

    elif (l == 4 and m == 2):

        #poly = (c15/c2) * (c7 * pow(z,2) - c1) * (c1 - pow(z,2))
        poly = (c15/c2) * (c14 * (c1 - pow(z,2)) - c56 * pow(z,2) - c2 * (c7 * pow(z,2) - c1))

    elif (l == 4 and m == 3):

        #poly = c105 * z * pow(s,3)
        poly = c105 * (c6 * pow(s,2) * dsdz + c6 * z * s * dsdpow(z,2) + c3 * z * pow(s,2) * d2sdz2)

    elif (l == 4 and m == 4):

        #poly = c105 * pow(s,4)
        poly = c420 * (c3 * pow(s,2) * dsdz + pow(s,3) * d2sdz2)

    elif (l == 5 and m == 0):

        #poly = (c1/c8) * z * (c63 * pow(z,4) - c70 * pow(z,2) + c15)
        poly = (c1/c8) * (c504 * pow(z,3) - c280 * z + z * (c756 * pow(z,2) - c140))

    elif (l == 5 and m == 1):

        #poly = (c15/c8) * s * (c21 * pow(z,4) - c14 * pow(z,2) + c1)
        poly = (c15/c8) * (s * (c252 * pow(z,2) - c28) + dsdz * (c168 * pow(z,3) - c56 * z) + d2sdz2 * (c21 * pow(z,4) - c14 * pow(z,2) + 1))

    elif (l == 5 and m == 2):

        #poly = (c105/c2) * z * (c1 - pow(z,2)) * (c3 * pow(z,2) - c1)
        poly = (c105/c2) * (c18 * z * (c1 - pow(z,2)) - c6 * z * (c3 * pow(z,2) - 1) - c24 * pow(z,3))

    elif (l == 5 and m == 3):

        #poly = (c105/c2) * pow(s,3) * (c9 * pow(z,2) - 1)
        poly = (c105/c2) * (c6 * s * pow(dsdz,2) * (c9 * pow(z,2) - c1) + c3 * pow(s,2) * d2sdz2 * (c9 * pow(z,2) - c1) + c108 * pow(s,2) * dsdz * z + c18 * pow(s,3))

    elif (l == 5 and m == 4):

        #poly = c945 * z * pow(s,4)
        poly = c945 * (c8 * pow(s,3) * dsdz + c12 * z * pow(s,2) * dsdpow(z,2) + c4 * z * pow(s,3) * d2sdz2)

    elif (l == 5 and m == 5):

        #poly = c945 * pow(s,5)
        poly = c945 * (c20 * pow(s,3) * pow(dsdz,2) + c5 * pow(s,4) * d2sdz2)

    else:

        raise RuntimeError("legendre_polynomial: Unknown input")

    return poly

#-------------------------------------------------------------------------------

def K(m, l):

    Kr = sqrt((float(2 * l + 1) * factorial(l-fabs(m))) / (4.0 * pi * factorial(l+abs(m))))

    return Kr

#-------------------------------------------------------------------------------

def harmonic_derivatives(theta, phi, m, l):

    if (m > 0):

       #sphericalHarmonic = sqrt(2.0) * K(m,l) * cos( real(m,RKIND) * phi) * legendre_polynomial(cos(theta),  m, l)

       dsh_dtheta        =                    cos( float(m) * phi) * dlegendre_polynomial_dz(cos(theta),  m, l) * (-sin(theta))
       dsh_dphi          = -float(m)        * sin( float(m) * phi) * legendre_polynomial(cos(theta),  m, l)
       d2sh_dthetadtheta =                    cos( float(m) * phi) * ((-cos(theta)) * dlegendre_polynomial_dz(cos(theta),  m, l) + pow(sin(theta),2) * d2legendre_polynomial_dz2(cos(theta),  m, l))
       d2sh_dthetadphi   = -float(m)        * sin( float(m) * phi) * dlegendre_polynomial_dz(cos(theta),  m, l) * (-sin(theta))
       d2sh_dphidtheta   = -float(m)        * sin( float(m) * phi) * dlegendre_polynomial_dz(cos(theta),  m, l) * (-sin(theta))
       d2sh_dphidphi     = -pow(float(m),2) * cos( float(m) * phi) * legendre_polynomial(cos(theta),  m, l)

       dsh_dtheta        = dsh_dtheta        * sqrt(2.0) * K(m,l)
       dsh_dphi          = dsh_dphi          * sqrt(2.0) * K(m,l)
       d2sh_dthetadtheta = d2sh_dthetadtheta * sqrt(2.0) * K(m,l)
       d2sh_dthetadphi   = d2sh_dthetadphi   * sqrt(2.0) * K(m,l)
       d2sh_dphidtheta   = d2sh_dphidtheta   * sqrt(2.0) * K(m,l)
       d2sh_dphidphi     = d2sh_dphidphi     * sqrt(2.0) * K(m,l)

    elif (m == 0):

       #sphericalHarmonic =                   K(m,l) *                             legendre_polynomial(cos(theta),  m, l)

       dsh_dtheta        = K(m,l) * dlegendre_polynomial_dz(cos(theta),  m, l) * (-sin(theta))
       dsh_dphi          = 0.0
       d2sh_dthetadtheta = K(m,l) * ((-cos(theta)) * dlegendre_polynomial_dz(cos(theta),  m, l) + pow(sin(theta),2) * d2legendre_polynomial_dz2(cos(theta),  m, l))
       d2sh_dthetadphi   = 0.0
       d2sh_dphidtheta   = 0.0
       d2sh_dphidphi     = 0.0

    elif (m < 0):

       #sphericalHarmonic = sqrt(2.0) * K(m,l) * sin(-real(m,RKIND) * phi) * legendre_polynomial(cos(theta), -m, l)

       dsh_dtheta        =                    sin(-float(m) * phi) * dlegendre_polynomial_dz(cos(theta),  -m, l) * (-sin(theta))
       dsh_dphi          = -float(m)        * cos(-float(m) * phi) * legendre_polynomial(cos(theta),  -m, l)
       d2sh_dthetadtheta =                    sin(-float(m) * phi) * ((-cos(theta)) * dlegendre_polynomial_dz(cos(theta),  -m, l) + pow(sin(theta),2) * d2legendre_polynomial_dz2(cos(theta),  -m, l))
       d2sh_dthetadphi   = -float(m)        * cos(-float(m) * phi) * dlegendre_polynomial_dz(cos(theta),  -m, l) * (-sin(theta))
       d2sh_dphidtheta   = -float(m)        * cos(-float(m) * phi) * dlegendre_polynomial_dz(cos(theta),  -m, l) * (-sin(theta))
       d2sh_dphidphi     = -pow(float(m),2) * sin(-float(m) * phi) * legendre_polynomial(cos(theta),  -m, l)

       dsh_dtheta        = dsh_dtheta        * sqrt(2.0) * K(m,l)
       dsh_dphi          = dsh_dphi          * sqrt(2.0) * K(m,l)
       d2sh_dthetadtheta = d2sh_dthetadtheta * sqrt(2.0) * K(m,l)
       d2sh_dthetadphi   = d2sh_dthetadphi   * sqrt(2.0) * K(m,l)
       d2sh_dphidtheta   = d2sh_dphidtheta   * sqrt(2.0) * K(m,l)
       d2sh_dphidphi     = d2sh_dphidphi     * sqrt(2.0) * K(m,l)

    return dsh_dtheta, dsh_dphi, d2sh_dthetadtheta, d2sh_dthetadphi, d2sh_dphidtheta, d2sh_dphidphi

#-------------------------------------------------------------------------------

def spherical_harmonic_geographical_derivatives(lon, lat, m, l):

    theta = pi / 2.0 - lat
    phi = lon

    dsh_dtheta, dsh_dphi, d2sh_dthetatheta, d2sh_dthetaphi, d2sh_dphitheta, d2sh_dphiphi = harmonic_derivatives(theta, phi, m, l)

    dphi_dlon = 1.0
    dtheta_dlat = -1.0

    dsh_dlon     = dsh_dphi         * dphi_dlon
    dsh_dlat     = dsh_dtheta       * dtheta_dlat
    d2sh_dlonlon = d2sh_dphiphi     * dphi_dlon   * dphi_dlon
    d2sh_dlonlat = d2sh_dphitheta   * dphi_dlon   * dtheta_dlat
    d2sh_dlatlon = d2sh_dthetaphi   * dphi_dlon   * dtheta_dlat
    d2sh_dlatlat = d2sh_dthetatheta * dtheta_dlat * dtheta_dlat

    return dsh_dlon, dsh_dlat, d2sh_dlonlon, d2sh_dlonlat, d2sh_dlatlon, d2sh_dlatlat

#-------------------------------------------------------------------------------

def spherical_harmonic(theta, phi, m, l):

    if (m > 0):
       sphericalHarmonic = sqrt(2.0) * K(m,l) * cos( float(m) * phi) * legendre_polynomial(cos(theta),  m, l)
    elif (m == 0):
       sphericalHarmonic =             K(m,l) *                        legendre_polynomial(cos(theta),  m, l)
    elif (m < 0):
       sphericalHarmonic = sqrt(2.0) * K(m,l) * sin(-float(m) * phi) * legendre_polynomial(cos(theta), -m, l)

    return sphericalHarmonic

#-------------------------------------------------------------------------------

def spherical_harmonics_geographical(lon, lat, m, l):

    theta = pi * 0.5 - lat
    phi = lon

    return spherical_harmonic(theta, phi, m, l)

#-------------------------------------------------------------------------------

def velocities_strains_analytical(lat, lon, mu, lu, mv, lv):

    r = 1.0

    u = spherical_harmonics_geographical(lon, lat, mu, lu)
    v = spherical_harmonics_geographical(lon, lat, mv, lv)

    du_dlon, du_dlat, d2u_dlonlon, d2u_dlonlat, d2u_dlatlon, d2u_dlatlat = spherical_harmonic_geographical_derivatives(lon, lat, mu, lu)
    dv_dlon, dv_dlat, d2v_dlonlon, d2v_dlonlat, d2v_dlatlon, d2v_dlatlat = spherical_harmonic_geographical_derivatives(lon, lat, mv, lv)

    strain11 = (1.0 / (r * cos(lat))) * (du_dlon - v * sin(lat))
    strain22 = dv_dlat / r
    strain12 = (0.5 / r) * (du_dlat + u * tan(lat) + dv_dlon / cos(lat))

    dstrain11_dlon = (1.0 / (r * cos(lat))) * (d2u_dlonlon - dv_dlon * sin(lat))
    dstrain22_dlat = d2v_dlatlat / r

    dstrain12_dlon = (0.5 / r) * (d2u_dlatlon + du_dlon * tan(lat) + d2v_dlonlon / cos(lat))
    dstrain12_dlat = (0.5 / r) * (d2u_dlatlat + du_dlat * tan(lat) + u / pow(cos(lat),2) + d2v_dlatlon / cos(lat) + dv_dlon * (tan(lat) / cos(lat)))

    divu = (1.0 / (r * cos(lat))) * dstrain11_dlon + \
           (1.0 / r)              * dstrain12_dlat - \
           (2.0 / r) * tan(lat)   * strain12

    divv = (1.0 / (r * cos(lat))) * dstrain12_dlon + \
           (1.0 / r)              * dstrain22_dlat + \
           (1.0 / r) * tan(lat)   * strain11       - \
           (1.0 / r) * tan(lat)   * strain22

    return u, v, strain11, strain22, strain12, divu, divv

#-------------------------------------------------------------------------------

def grid_rotation_forward(x, y, z, rotateCartesianGrid):

    # rotate xyz coordinates from geographical grid to rotated grid with poles on real equator

    if (rotateCartesianGrid):

       xp = -z
       yp = y
       zp = x

    else:

       xp = x
       yp = y
       zp = z

    return xp, yp, zp

#-------------------------------------------------------------------------------

def latlon_from_xyz(x, y, z, r):

    # given xyz coordinates determine the latitude and longitude

    lon = math.atan2(y, x)
    lat = math.asin(z/r)

    return lat, lon

#-------------------------------------------------------------------------------

def create_ic():

    mu = 3
    lu = 5

    mv = 2
    lv = 4

    gridSizes = [2562, 10242, 40962, 163842]

    rotateCartesianGrid = True
    r = 1.0

    for gridSize in gridSizes:

        print("  Gridsize: ", gridSize)

        # input
        filenameIn = "x1.%i.grid.nc" %(gridSize)

        fileIn = Dataset(filenameIn,"r")

        nCells = len(fileIn.dimensions["nCells"])
        nVertices = len(fileIn.dimensions["nVertices"])

        xCell = fileIn.variables["xCell"][:]
        yCell = fileIn.variables["yCell"][:]
        zCell = fileIn.variables["zCell"][:]

        xVertex = fileIn.variables["xVertex"][:]
        yVertex = fileIn.variables["yVertex"][:]
        zVertex = fileIn.variables["zVertex"][:]

        latCell = fileIn.variables["latCell"][:]
        lonCell = fileIn.variables["lonCell"][:]

        latVertex = fileIn.variables["latVertex"][:]
        lonVertex = fileIn.variables["lonVertex"][:]

        fileIn.close()

        # velocities
        uVelocity = np.zeros(nVertices)
        vVelocity = np.zeros(nVertices)

        strain11VertexAnalytical = np.zeros(nVertices)
        strain22VertexAnalytical = np.zeros(nVertices)
        strain12VertexAnalytical = np.zeros(nVertices)

        stressDivergenceUAnalytical = np.zeros(nVertices)
        stressDivergenceVAnalytical = np.zeros(nVertices)

        for iVertex in range(0, nVertices):

            #print("iVertex: ", iVertex, latVertex[iVertex], lonVertex[iVertex])

            xp, yp, zp = grid_rotation_forward(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], rotateCartesianGrid)
            lat, lon = latlon_from_xyz(xp, yp, zp, r)

            u, v, strain11, strain22, strain12, divu, divv = velocities_strains_analytical(lat, lon, mu, lu, mv, lv)

            uVelocity[iVertex] = u
            vVelocity[iVertex] = v

            strain11VertexAnalytical[iVertex] = strain11
            strain22VertexAnalytical[iVertex] = strain22
            strain12VertexAnalytical[iVertex] = strain12

            stressDivergenceUAnalytical[iVertex] = divu
            stressDivergenceVAnalytical[iVertex] = divv

        strain11CellAnalytical = np.zeros(nCells)
        strain22CellAnalytical = np.zeros(nCells)
        strain12CellAnalytical = np.zeros(nCells)

        for iCell in range(0, nCells):

            #print("iCell: ", iCell, latCell[iCell], lonCell[iCell])

            xp, yp, zp = grid_rotation_forward(xCell[iCell], yCell[iCell], zCell[iCell], rotateCartesianGrid)
            lat, lon = latlon_from_xyz(xp, yp, zp, r)

            u, v, strain11, strain22, strain12, divu, divv = velocities_strains_analytical(lat, lon, mu, lu, mv, lv)

            strain11CellAnalytical[iCell] = strain11
            strain22CellAnalytical[iCell] = strain22
            strain12CellAnalytical[iCell] = strain12

        solveVelocityPrevious = np.ones(nVertices,dtype="i")


        # output
        filenameOut = "ic_%i.nc" %(gridSize)

        fileOut = Dataset(filenameOut, "w", format="NETCDF3_CLASSIC")

        fileOut.createDimension("nVertices", nVertices)
        fileOut.createDimension("nCells", nCells)

        var = fileOut.createVariable("uVelocity","d",dimensions=["nVertices"])
        var[:] = uVelocity[:]

        var = fileOut.createVariable("vVelocity","d",dimensions=["nVertices"])
        var[:] = vVelocity[:]

        var = fileOut.createVariable("solveVelocityPrevious","i",dimensions=["nVertices"])
        var[:] = solveVelocityPrevious[:]

        var = fileOut.createVariable("strain11VertexAnalytical","d",dimensions=["nVertices"])
        var[:] = strain11VertexAnalytical[:]

        var = fileOut.createVariable("strain22VertexAnalytical","d",dimensions=["nVertices"])
        var[:] = strain22VertexAnalytical[:]

        var = fileOut.createVariable("strain12VertexAnalytical","d",dimensions=["nVertices"])
        var[:] = strain12VertexAnalytical[:]

        var = fileOut.createVariable("strain11CellAnalytical","d",dimensions=["nCells"])
        var[:] = strain11CellAnalytical[:]

        var = fileOut.createVariable("strain22CellAnalytical","d",dimensions=["nCells"])
        var[:] = strain22CellAnalytical[:]

        var = fileOut.createVariable("strain12CellAnalytical","d",dimensions=["nCells"])
        var[:] = strain12CellAnalytical[:]

        var = fileOut.createVariable("stressDivergenceUAnalytical","d",dimensions=["nVertices"])
        var[:] = stressDivergenceUAnalytical[:]

        var = fileOut.createVariable("stressDivergenceVAnalytical","d",dimensions=["nVertices"])
        var[:] = stressDivergenceVAnalytical[:]

        fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ic()
