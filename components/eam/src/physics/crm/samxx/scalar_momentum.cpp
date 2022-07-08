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
void esmt_fft_forward(real2d& arr_in, real1d& k_out, real2d& arr_out) {
/*!------------------------------------------------------------------
 * Purpose: calculate forward FFT transform
 * Author: Walter Hannah - Lawrence Livermore National Lab
 * adapted from SP-WRF code provided by Stefan Tulich
 *------------------------------------------------------------------*/
   auto& dx  = ::dx;

   int lensave, ier, nh, nl;
   int lot, jump, n, inc, lenr, lensav, lenwrk;
   real wsave[nx+15];

   // naming convention follows fftpack5 routines
   n = nx;
   lot = 1;
   lensav = n+15;
   inc = 1;
   lenr = nx;
   jump = nx;
   lenwrk = lenr;

#ifndef USE_ORIG_FFT
   realHost2d arrinHost  = arr_in.createHostCopy();
   realHost2d arroutHost = arr_out.createHostCopy();

   // initialization for FFT
   rfft1i(n,wsave,lensav,ier);

   //  do the forward transform
   //do k = 1,nzm
   for (auto k = 0; k < nzm; ++k) {
      real atmp[nx];
      real tmp[nx];
      for (auto i = 0; i < nx; ++i) { atmp[i] = arrinHost(k,i); }
      rfft1f( n, inc, atmp, lenr, wsave, lensav, tmp, lenwrk, ier );
      for (int i=0; i<nx; i++) { arroutHost(k,i) = atmp[i]; }
   }

   arroutHost.deep_copy_to(arr_out);
#else
    yakl::FFT<nx> fft;
    parallel_for( nzm , YAKL_LAMBDA (int k) {
      real atmp[nx+2];
      real tmp[nx];

      for (int i=0; i<nx ; i++) { atmp[i+1] = arr_in(k,i); }
      fft.forward(atmp, tmp);
      for (int i=0; i<nx; i++) { arr_out(k,i) = atmp[i+1];  }
    });
#endif

   if((n%2) == 0) {
      nh = n/2 - 1;
   } else {
      nh = (n-1)/2;
   }

   parallel_for( nh, YAKL_LAMBDA (int j) {
      k_out(2*j+1) = 2.*pi*real(j+1)/(real(n)*dx);   //cos
      k_out(2*j+2) = 2.*pi*real(j+1)/(real(n)*dx);   //sin
      if ((n%2) == 0) {
        k_out(n-1) =  2.*pi/(2.*dx);   //nyquist wavelength for even n
      }
      if (j==0) k_out(0) = 0.0;
   });
}

void esmt_fft_backward(real2d& arr_in, real2d& arr_out) {
/*!------------------------------------------------------------------
 *  Purpose: calculate backward FFT transform
 *  Author: Walter Hannah - Lawrence Livermore National Lab
 *  adapted from SP-WRF code provided by Stefan Tulich
 *------------------------------------------------------------------*/
   int lensave, ier, nh, nl;
   int lot, jump, n, inc, lenr, lensav, lenwrk;
   real wsave[nx+15];

   // naming convention follows fftpack5 routines
   n = nx;
   lot = 1;
   lensav = n+15;
   inc = 1;
   lenr = nx;
   jump = nx;
   lenwrk = lenr;

#ifndef USE_ORIG_FFT
   realHost2d arrinHost  = arr_in.createHostCopy();
   realHost2d arroutHost = arr_out.createHostCopy();

   // initialization for FFT
   rfft1i(n,wsave,lensav,ier);

   //  do the backward transform
   //do k = 1,nzm
   for (auto k = 0; k < nzm; ++k) {
      real atmp[nx];
      real tmp[nx];
      for (auto i = 0; i < nx; ++i) atmp[i] = arrinHost(k,i);
      rfft1b( n, inc, atmp, lenr, wsave, lensav, tmp, lenwrk, ier );
      for (auto i = 0; i < nx; ++i) arroutHost(k,i) = atmp[i];
   }

   arroutHost.deep_copy_to(arr_out);
#else
    yakl::FFT<nx> fft;
    parallel_for( nzm , YAKL_LAMBDA (int k) {
      real atmp[nx+2];
      real tmp[nx];

      for (int i=0; i<nx ; i++) { atmp[i+1] = arr_in(k,i); }
      fft.forward(atmp, tmp);
      for (int i=0; i<nx; i++) { arr_out(k,i) = atmp[i+1];  }
    });
#endif
}


void scalar_momentum_pgf( int icrm, real3d& u_s, real3d& tend ) {
 /*!------------------------------------------------------------------
   ! Purpose: calculate pgf for scalar momentum transport
   ! Author: Walter Hannah - Lawrence Livermore National Lab
   ! adapted from SP-WRF code by Stefan Tulich
   !------------------------------------------------------------------*/
   auto &ncrms = :: ncrms;
   auto &z     = :: z;
   auto &pres  = :: pres;
   auto &zi    = :: zi;
   auto &w     = :: w;
   auto &rho   = :: rho;
   
   real dampwt;
   real1d k_arr("k_arr",nx);
   real1d u_s_avg("u_s_avg",nzm);
   real1d shr("shr",nzm);
   real2d a("a",nzm,nx);
   real2d b("b",nzm,nx);
   real2d c("c",nzm,nx);
   real2d rhs("rhs",nzm,nx);
   real2d w_i("w_i",nzm,nx);
   real2d w_hat("w_hat",nzm,nx);
   real2d pgf_hat("pgf_hat",nzm,nx);
   real2d pgf("pgf",nzm,nx);
   real1d dz("dz",nzm+1);

   // The loop over "y" points is mostly unessary, since ESMT
   // is for 2D CRMs, but it is useful for directly comparing
   // ESMT tendencies to fully resolved 3D momentum transport

   // Calculate layer thickness
   //do k = 1,nzm
   parallel_for( nzm+1, YAKL_LAMBDA (int k) {
     if (k < nzm) {
      dz(k) = zi(k+1,icrm)-zi(k,icrm);
     }
     dz(nzm) = zi(nzm,icrm)-zi(nzm-1,icrm);
   });

   //do j=1,ny
   for (auto j = 0; j < ny; ++j) {
      //-----------------------------------------
      // Initialize stuff for averaging
      //-----------------------------------------
      parallel_for( nzm, YAKL_LAMBDA (int k) {
        u_s_avg(k) = 0.0;
        shr(k)     = 0.0;
      });

      //-----------------------------------------
      // Calculate shear of domain average profile
      // defined on scalar levels
      //-----------------------------------------
      //do k = 1,nzm
      parallel_for( nzm, YAKL_LAMBDA (int k) {
         for (auto i = 0; i < nx; ++i) {
            u_s_avg(k) = u_s_avg(k) + u_s(k,j,i);
           // yakl::atomicAdd( u_s_avg(k) , u_s(k,j,i) );

            // note that w is on interface levels
            w_i(k,i) = ( w(k,j,i,icrm) + w(k+1,j,i,icrm) )/2.0;
         }
         u_s_avg(k) = u_s_avg(k) / real(nx);
      });

      //do k = 2,nzm-1
      parallel_for( nzm-1, YAKL_LAMBDA (int k) {
         if ( k>0) {
           shr(k) = ( u_s_avg(k+1) - u_s_avg(k-1) )/(z(k+1,icrm)-z(k-1,icrm));
         }
        shr(0) = ( u_s_avg(1) - u_s_avg(0) ) / ( z(1,icrm) - z(0,icrm) );
      });

      /*!------------------------------------------------------------------------
      ! Use Poisson solver to calculate pressure gradient force (PGF)
      ! pgf is diagnosed from w * du_si/dz using the poisson equation
      ! (see Wu and Yanai 1994)
      !------------------------------------------------------------------------*/

      //-----------------------------------------
      // compute forward fft of w
      //-----------------------------------------
      esmt_fft_forward(w_i, k_arr, w_hat);

      yakl::fence();
      //-----------------------------------------
      // solve vertical structure equation
      // for each zonal wavelength (k_arr)
      // solution method involves constructing
      // a tridiagonal matrix
      //-----------------------------------------
      parallel_for( SimpleBounds<2>(nzm,nx) , YAKL_LAMBDA (int k, int i) {
         pgf_hat(k,i) = 0.;
       });

      // Loop through wavelengths
      //do i = 2,nx
      parallel_for( SimpleBounds<2>(nzm,nx), YAKL_LAMBDA (int k, int i) {
        if (i>0) {
           a(k,i) = dz(k+1) / ( dz(k+1) + dz(k) );
           // the factor of 1.25 crudely accounts for difference between 2D and 3D updraft geometry
           b(k,i) = -0.5 * 1.25 * pow(k_arr(i), 2.0) * dz(k) * dz(k+1) - 1.0;
           c(k,i) = dz(k) / ( dz(k+1) + dz(k) );
           rhs(k,i) = pow(k_arr(i), 2.) * w_hat(k,i) * shr(k) * dz(k) * dz(k+1);

           //lower boundary condition (symmetric)
           if (k == 0) {
             b(0,i) =  b(0,i) + a(0,i);
             a(0,i) = 0.0;
           }

           //upper boundary condition (symmetric)
           if (k == nzm-1) {
             b(nzm-1,i) = b(nzm-1,i) + a(nzm-1,i);
             c(nzm-1,i) = 0.0;
           }
        }
     });

      // gaussian elimination with no pivoting
      //do k = 1,nzm-1
      parallel_for( nx, YAKL_LAMBDA (int i) {
       if (i>0) {
         for (auto k = 0; k < nzm-1; ++k) {
           b(k+1,i) = b(k+1,i) - a(k+1,i) / b(k,i) * c(k,i);
           rhs(k+1,i) = rhs(k+1,i) - a(k+1,i) / b(k,i) * rhs(k,i);
         }
       }
      });


      parallel_for( nx, YAKL_LAMBDA (int i) {
       if (i>0) {
         // backward substitution
         rhs(nzm-1,i) = rhs(nzm-1,i) / b(nzm-1,i);
         //do k=nzm-1,1,-1
         for (auto k = nzm-2; k > -1; --k) {
           rhs(k,i) = ( rhs(k,i) - c(k,i) * rhs(k+1,i) ) / b(k,i);
         }
       }
      });

         //do k = 1,nzm
       parallel_for( SimpleBounds<2>(nzm,nx), YAKL_LAMBDA (int k, int i) {
         pgf_hat(k,i) = rhs(k,i);
       });


      // Note sure what this part does... 
      // something about the Nyquist freq and whether nx is odd or even
      if (nx%2 == 0) {
        // do k=1,nzm
         parallel_for( nzm, YAKL_LAMBDA (int k) {
           pgf_hat(k,nx-1) = pgf_hat(k,nx-1) / 2.0;
         }); // k
      }

      //-----------------------------------------
      // invert fft of pgf_hat to get pgf
      //-----------------------------------------
      esmt_fft_backward(pgf_hat, pgf);

      //-----------------------------------------
      // Compute final tendency
      //----------------------------------------
      //do k = 1,nzm
      //   do i = 1,nx
      parallel_for( SimpleBounds<2>(nzm,nx) , YAKL_LAMBDA (int k, int i) {
         if (k == 0) {
            tend(k,j,i) = 0.0;
         } else {
            tend(k,j,i) = -1.0 * pgf(k,i) * rho(k,icrm);
         }
      });
   } // j
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
   
   real3d u_esmt_pgf_3D("u_esmt_pgf_3D",nzm,ny,nx);
   real3d v_esmt_pgf_3D("v_esmt_pgf_3d",nzm,ny,nx);

   real3d u_esmt_tmp("u_esmt_tmp", nzm,ny,nx);
   real3d v_esmt_tmp("v_esmt_tmp", nzm,ny,nx);

   //do icrm = 1 , ncrms
   for (auto icrm = 0; icrm < ncrms; ++icrm) {
     parallel_for( SimpleBounds<3>(nzm,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        u_esmt_tmp(k,j,i) = u_esmt(k,j+offy_s,i+offx_s,icrm);
        v_esmt_tmp(k,j,i) = v_esmt(k,j+offy_s,i+offx_s,icrm);
     });

      scalar_momentum_pgf(icrm,u_esmt_tmp,u_esmt_pgf_3D);
      scalar_momentum_pgf(icrm,v_esmt_tmp,v_esmt_pgf_3D);
   
     // Add PGF tendency
     //do k=1,nzm
     //   do j=1,ny
     //      do i=1,nx
     parallel_for( SimpleBounds<3>(nzm,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        u_esmt(k,j+offy_s,i+offx_s,icrm) = u_esmt(k,j+offy_s,i+offx_s,icrm) + u_esmt_pgf_3D(k,j,i)*dtn;
        v_esmt(k,j+offy_s,i+offx_s,icrm) = v_esmt(k,j+offy_s,i+offx_s,icrm) + v_esmt_pgf_3D(k,j,i)*dtn;
     });
  }
}
#endif
