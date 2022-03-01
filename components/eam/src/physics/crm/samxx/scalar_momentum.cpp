/*!-----------------------------------------------------------------------
! Purpose:
!
! Explicit Scalar Momentum Transport (ESMT)
! Transport large-scale momentum as non-conserved scalars,
! including estimate of 3D pressure gradient force. This was
! implemented to help avoid using any physics modules outside
! of CRM and radiation routines. These non-SP routines can be
! problematic when changing the vertical grid.
! See Tulich (2015) for further details on ESMT
!
! Revision history:
! Nov, 2017 - Walter Hannah - Lawrence Livermore National Lab
!             initial version based on crmtracers.F90
!             Possoin solver and fft routines provided by Stefan Tulich
!---------------------------------------------------------------------------*/
#include "scalar_momentum.h"

#ifdef MMF_ESMT

void scalar_momentum_pgf( real4d& scalar_wind, real4d& tend ) {
   /*!------------------------------------------------------------------
   ! Purpose: calculate pgf for scalar momentum transport
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   ! adapted from SP-WRF code by Stefan Tulich
   !------------------------------------------------------------------*/
   auto &ncrms = :: ncrms;
   auto &z     = :: z;
   auto &zi    = :: zi;
   auto &w     = :: w;
   auto &rho   = :: rho;
   auto &dx    = ::dx;

   real constexpr pi = 3.14159;
   
   real1d k_arr("k_arr",nx);
   real2d dz_loc("dz_loc",nzm+1,ncrms);
   real3d scalar_wind_avg("scalar_wind_avg",nzm,ny,ncrms);
   real3d shear("shear",nzm,ny,ncrms);
   real4d a("a",nzm,ny,nx,ncrms);
   real4d b("b",nzm,ny,nx,ncrms);
   real4d c("c",nzm,ny,nx,ncrms);
   real4d w_i("w_i",nzm,ny,nx,ncrms);
   real4d pgf("pgf",nzm,ny,nx,ncrms);
   int nx2 = nx+2;
   real4d w_hat("w_hat",nzm,ny,nx2,ncrms);
   real4d pgf_hat("pgf_hat",nzm,ny,nx2,ncrms);

   // The loop over "y" points is mostly unessary, since ESMT
   // is for 2D CRMs, but it is useful for directly comparing
   // ESMT tendencies to fully resolved 3D momentum transport

   // Calculate layer thickness
   //do k = 1,nzm
   //  for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<2>(nzm+1,ncrms) , YAKL_LAMBDA (int k, int icrm) {
      if (k < nzm) {
         dz_loc(k,icrm) = zi(k+1,icrm)-zi(k,icrm);
      } else{
         dz_loc(k,icrm) = zi(k,icrm)-zi(k-1,icrm); 
      }
   });


   //-----------------------------------------
   // Initialize stuff for averaging
   //-----------------------------------------
   //do k=1,nzm
   //  for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
      scalar_wind_avg(k,j,icrm) = 0.0;
      shear(k,j,icrm) = 0.0;
   });

   //-----------------------------------------
   // Calculate shear of wind profile
   // averaged over the "x" dimension
   // defined on scalar levels
   //-----------------------------------------
   //do k=1,nzm
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
      // real tmp = scalar_wind(k,j,i,icrm) / real(nx);
      yakl::atomicAdd( scalar_wind_avg(k,j,icrm) , scalar_wind(k,j,i,icrm) / real(nx) );
      // note that w is on interface levels - need to interpolate to mid-levels
      w_i(k,j,i,icrm) = ( w(k,j,i,icrm) + w(k+1,j,i,icrm) )/2.0 ;
   });

   //do k = 1,nzm
   //  for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<3>(nzm-1,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
      if ( k>0) {
        shear(k,j,icrm) = ( scalar_wind_avg(k+1,j,icrm) - scalar_wind_avg(k-1,j,icrm) )/(z(k+1,icrm)-z(k-1,icrm));
      }
      shear(0,j,icrm) = ( scalar_wind_avg(1,j,icrm) - scalar_wind_avg(0,j,icrm) ) / ( z(1,icrm) - z(0,icrm) );
   });

   /*!------------------------------------------------------------------------
   ! Use Poisson solver to calculate pressure gradient force (PGF)
   ! pgf is diagnosed from w * d/dz(scalar_wind) using the poisson equation
   ! (see Wu and Yanai 1994)
   !------------------------------------------------------------------------*/

   //-----------------------------------------
   // compute forward fft of w
   //-----------------------------------------
   yakl::RealFFT1D<nx> fftx;
   fftx.init();

   // for (int k=0; k<nzm; k++) {
   //   for (int j=0; j<ny; j++) {
   //     for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_LAMBDA (int k, int j, int icrm) {
      SArray<real,1,nx+2> ftmp;
      for (int i=0; i<nx ; i++) { ftmp(i) = w_i(k,j,i,icrm); }
      fftx.forward(ftmp, fftx.trig, yakl::FFT_SCALE_ECMWF);
      for (int i=0; i<nx2; i++) { w_hat(k,j,i,icrm) = ftmp(i); }
   });

   //-----------------------------------------
   //-----------------------------------------
   int nh;

   if((nx%2) == 0) nh = nx/2 - 1;
   if((nx%2) != 0) nh = (nx-1)/2;

   parallel_for( nh, YAKL_LAMBDA (int j) {
      k_arr(2*j+1) = 2.*pi*real(j+1)/(real(nx)*dx);   //cos
      k_arr(2*j+2) = 2.*pi*real(j+1)/(real(nx)*dx);   //sin
      if (j==0) {
         k_arr(j) = 0.0;
      }
      if ( j==(nh-1) && (nx%2)==0) { 
         k_arr(nx-1) = 2.*pi/(2.*dx); //nyquist wavelength for even n
      }  
   });

   //-----------------------------------------
   // solve vertical structure equation
   // for each zonal wavelength (k_arr)
   // solution method involves constructing
   // a tridiagonal matrix
   //-----------------------------------------
   //do k=1,nzm
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm,ny,nx2,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      pgf_hat(k,j,i,icrm) = 0.;
   });

   //-----------------------------------------
   // Loop through wavelengths
   //-----------------------------------------
   //do k=1,nzm
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (i>0) {
         a(k,j,i,icrm) = dz_loc(k+1,icrm) / ( dz_loc(k+1,icrm) + dz_loc(k,icrm) );
         // the factor of 1.25 crudely accounts for difference between 2D and 3D updraft geometry
         b(k,j,i,icrm) = -0.5 * 1.25 * pow(k_arr(i), 2.0) * dz_loc(k,icrm) * dz_loc(k+1,icrm) - 1.0;
         c(k,j,i,icrm) = dz_loc(k,icrm) / ( dz_loc(k+1,icrm) + dz_loc(k,icrm) );
         pgf_hat(k,j,i,icrm) = pow(k_arr(i), 2.) * w_hat(k,j,i,icrm) * shear(k,j,icrm) * dz_loc(k,icrm) * dz_loc(k+1,icrm);

         //lower boundary condition (symmetric)
         if (k == 0) {
            b(0,j,i,icrm) =  b(0,j,i,icrm) + a(0,j,i,icrm);
            a(0,j,i,icrm) = 0.0;
         }

         //upper boundary condition (symmetric)
         if (k == nzm-1) {
            b(nzm-1,j,i,icrm) = b(nzm-1,j,i,icrm) + a(nzm-1,j,i,icrm);
            c(nzm-1,j,i,icrm) = 0.0;
         }
      }
   });

   //-----------------------------------------
   // gaussian elimination with no pivoting
   //-----------------------------------------
   //do k=1,nzm-1
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int i, int icrm) {
      if (i>0) {
         b(k+1,j,i,icrm) = b(k+1,j,i,icrm) - a(k+1,j,i,icrm) / b(k,j,i,icrm) * c(k,j,i,icrm);
         pgf_hat(k+1,j,i,icrm) = pgf_hat(k+1,j,i,icrm) - a(k+1,j,i,icrm) / b(k,j,i,icrm) * pgf_hat(k,j,i,icrm);
      }
   });

   //-----------------------------------------
   // backward substitution
   //-----------------------------------------
   //do k=1,nzm-1
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm-1,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (i>0) {
         // f90: do k=nzm-1,1,-1
         // cpp: for (auto k = nzm-2; k > -1; --k) {
         int kk = nzm-1 - k;
         if (kk == nzm-1) {
            pgf_hat(kk,j,i,icrm) = pgf_hat(kk,j,i,icrm) / b(kk,j,i,icrm);
         } else {
            pgf_hat(kk,j,i,icrm) = ( pgf_hat(kk,j,i,icrm) - c(kk,j,i,icrm) * pgf_hat(kk+1,j,i,icrm) ) / b(kk,j,i,icrm);
         }
      }
   });

   // Note sure what this part does... 
   // something about the Nyquist freq and whether nx is odd or even
   if (nx%2 == 0) {
      //do k=1,nzm
      //  for (int j=0; j<ny; j++) {
      //    for (int i=0; i<nx; i++) {
      //      for (int icrm=0; icrm<ncrms; icrm++) {
      parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
         pgf_hat(k,j,nx-1,icrm) = pgf_hat(k,j,nx-1,icrm) / 2.0;
      });
   }

   //-----------------------------------------
   // invert fft of pgf_hat to get pgf
   //-----------------------------------------
   // for (int k=0; k<nzm; k++) {
   //   for (int j=0; j<ny; j++) {
   //     for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<3>(nzm,ny,ncrms) , YAKL_DEVICE_LAMBDA (int k, int j, int icrm) {
      SArray<real,1,nx+2> ftmp;
      for (int i=0; i<nx2; i++) { ftmp(i) = pgf_hat(k,j,i,icrm); }
      fftx.inverse(ftmp, fftx.trig, yakl::FFT_SCALE_ECMWF);
      for (int i=0; i<nx ; i++) { pgf(k,j,i,icrm) = ftmp(i); }
   });

   //-----------------------------------------
   // Compute final tendency
   //-----------------------------------------
   //do k=1,nzm
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {t
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      if (k == 0) {
         tend(k,j,i,icrm) = 0.0;
      } else {
         tend(k,j,i,icrm) = -1.0 * pgf(k,j,i,icrm) * rho(k,icrm);
      }
   });

}

void scalar_momentum_tend() {
   /*!------------------------------------------------------------------
   * Purpose: Calculate pressure gradient effects on scalar momentum
   * Author: Walter Hannah - Lawrence Livermore National Lab
   *------------------------------------------------------------------*/
   auto &dtn       = :: dtn;
   auto &ncrms     = :: ncrms;
   auto &u_esmt    = :: u_esmt;
   auto &v_esmt    = :: v_esmt;
   
   real4d u_esmt_pgf_3D("u_esmt_pgf_3D",nzm,ny,nx,ncrms);
   real4d v_esmt_pgf_3D("v_esmt_pgf_3d",nzm,ny,nx,ncrms); 

   // Calculate pressure gradient force tendency
   scalar_momentum_pgf(u_esmt,u_esmt_pgf_3D);
   scalar_momentum_pgf(v_esmt,v_esmt_pgf_3D);

   // Add PGF tendency
   //do k=1,nzm
   //  for (int j=0; j<ny; j++) {
   //    for (int i=0; i<nx; i++) {
   //      for (int icrm=0; icrm<ncrms; icrm++) {
   parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
     u_esmt(k,j+offy_s,i+offx_s,icrm) = u_esmt(k,j+offy_s,i+offx_s,icrm) + u_esmt_pgf_3D(k,j,i,icrm)*dtn;
     v_esmt(k,j+offy_s,i+offx_s,icrm) = v_esmt(k,j+offy_s,i+offx_s,icrm) + v_esmt_pgf_3D(k,j,i,icrm)*dtn;
   });

}
#endif
