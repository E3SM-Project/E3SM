/*
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_okubo_weiss_eigenvalues
!
!> \brief MPAS ocean analysis core member: okubo_weiss
!> \author Andre Schmeisser
!> \date   August 2014
!> \details
!>  MPAS ocean analysis core member: okubo_weiss
!
!-----------------------------------------------------------------------
*/

#include <math.h>

#ifdef UNDERSCORE
#define compute_ev_2 compute_ev_2_
#define compute_ev_3 compute_ev_3_
#else
#ifdef DOUBLEUNDERSCORE
#define compute_ev_2 compute_ev_2__
#define compute_ev_3 compute_ev_3__
#endif
#endif

#ifdef SINGLE_PRECISION
    typedef float real;
#else
    typedef double real;
#endif

#define swap(A, B) tmp = A; A = B; B = tmp;

static inline void sort_descending_complex_3(real* wr, real* wi)
{
    real tmp;
    if (wr[0] < wr[1])
    {
        swap(wr[0],  wr[1]);
        swap(wi[0],  wi[1]);
    }
    if (wr[1] < wr[2])
    {
        swap(wr[1],  wr[2]);
        swap(wi[1],  wi[2]);
    }
    if (wr[0] < wr[1])
    {
        swap(wr[0],  wr[1]);
        swap(wi[0],  wi[1]);
    }
}

/*
!***********************************************************************
!
!  routine compute_ev_2
!
!> \brief   Compute the eigenvalues of real 2x2 matrix
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  Compute the eigenvalues of real 2x2 matrix
!
!-----------------------------------------------------------------------
*/
void compute_ev_2(real A[4], real wr[2], real wi[2])
{
    real a = A[0];
    real b = A[1];
    real c = A[2];
    real d = A[3];
    real trA = a+d;
    real detA = a*d - b*c;
    real discr = trA*trA - 4*detA;
    real root;
    if (discr >= 0)
    {
        root = sqrt(discr);
        wr[0] = (trA + root) / 2;
        wr[1] = (trA - root) / 2;
        wi[0] = wi[1] = 0;
    }
    else
    {
        wr[0] = wr[1] = trA / 2;
        if (fabs(discr) < 1e-10)
        {
            wi[0] = wi[1] = 0;
        }
        else
        {
            root = sqrt(-discr);
            wi[0] =   root / 2;
            wi[1] = - root / 2;
        }
    }

}

/*
!***********************************************************************
!
!  routine compute_ev_3
!
!> \brief   Compute the eigenvalues of real 3x3 matrix
!> \author  Andre Schmeisser
!> \date    August 2014
!> \details 
!>  Compute the eigenvalues of real 3x3 matrix
!>    characteristic polynomial: det(A-x*I) = 0 for eigenvalues x
!>    x^3 + x^2 (- Tr(A)) + x (-1/2*(Tr(A^2)-Tr(A)^2)) - det(A) = 0
!>    x^3 + b*x^2 + c*x + d = 0
!> 
!>    reduce do depressed cubic t^3 + pt + qt = 0
!>    x = t - b/3
!>    p = (3c - b^2)/3 = c - b^2/3
!>    q = (2b^3 - 9bc + 27d)/27 = (2b^3 - 9bc)/27 + d
!> 
!
!-----------------------------------------------------------------------
*/
void compute_ev_3(real* mat, real* wr, real* wi )
{
    /* find value to normalize with, to protect against over-/underflow */
    double maxAbsVal = 0.;
    int i;

    for (i = 0; i < 9; i++)
        if (fabs(mat[i]) > maxAbsVal)
            maxAbsVal = fabs(mat[i]);

    if (maxAbsVal == 0.)
    {
        wr[0] = wr[1] = wr[2] = 0.;
        wi[0] = wi[1] = wi[2] = 0.;
        return;
    }

    double normVal = maxAbsVal;
    double A[9];
    for (i = 0; i < 9; i++)
        A[i] = mat[i] / normVal;

    /*
    characteristic polynomial: det(A-x*I) = 0 for eigenvalues x
    x^3 + x^2 (- Tr(A)) + x (-1/2*(Tr(A^2)-Tr(A)^2)) - det(A) = 0
    x^3 + b*x^2 + c*x + d = 0
    */

    double b = -A[0] -A[4] -A[8];
    double c = A[0]*A[4] + A[0]*A[8] + A[4]*A[8] - A[2]*A[6] - A[1]*A[3] - A[5]*A[7];
    double det =   A[0] * A[4] * A[8]
                 + A[1] * A[5] * A[6]
                 + A[2] * A[3] * A[7]
                 - A[2] * A[4] * A[6]
                 - A[1] * A[3] * A[8]
                 - A[0] * A[5] * A[7];
    double d = -det;

    /*
    reduce do depressed cubic t^3 + pt + qt = 0
    x = t - b/3
    p = (3c - b^2)/3 = c - b^2/3
    q = (2b^3 - 9bc + 27d)/27 = (2b^3 - 9bc)/27 + d
    */

    double Q = (b*b - 3*c)/9;
    double R = (b*(2*b*b - 9*c) + 27*d) / 54;
    double RR = R*R;
    double QQQ = Q*Q*Q;
    double b3 = b / 3;
    if (RR < QQQ)
    {
        // three real roots
        double theta = acos(R/sqrt(QQQ));
        double f = -2 * sqrt(Q);
        wr[0] = f * cos(theta/3) -b3;
        wr[1] = f * cos((theta + 2*M_PI)/3) -b3;
        wr[2] = f * cos((theta - 2*M_PI)/3) -b3;
        wi[0] = wi[1] = wi[2] = 0;
    }
    else
    {
        // one real root, two complex conjugates
        double sign = R >= 0 ? 1 : -1;
        double rad = fabs(R) - sqrt(RR - QQQ);
        double root3;
        if (rad >= 0)
            root3 = pow(rad, 1./3.);
        else
            root3 = -pow(-rad, 1./3.);
        double A = -sign*root3;
        double B = A == 0 ? 0 : Q / A;
        wr[0] = (A+B) - b3;
        wi[0] = 0;
        
        wr[1] = -0.5*(A+B)-b3;
        wi[1] = 0.5*sqrt(3)*(A-B);

        wr[2] =  wr[1];
        wi[2] = -wi[1];
    }
    sort_descending_complex_3(wr, wi);

    for (i = 0; i < 3; i++)
    {
        wr[i] *= normVal;
        wi[i] *= normVal;
    }
}

