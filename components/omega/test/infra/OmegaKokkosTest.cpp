//===-- Test driver for OMEGA Kokkos -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Kokkos
///
/// This driver tests the OmegaKokkos capabilities.
/// The code is based on Kokkos Tutorial
///
//
//===-----------------------------------------------------------------------===/

#include <iostream>

#include "OmegaKokkos.h"

using namespace OMEGA;

int main(int argc, char **argv) {

   int RetVal  = 0;
   int N       = 4096;
   int M       = 1024;
   int S       = 4194304;
   int nrepeat = 100;

   try {

      Kokkos::initialize(argc, argv);
      {

         // Allocate y, x vectors and Matrix A on device.
         Array1DR8 y("y", N);
         Array1DR8 x("x", M);
         Array2DR8 A("A", N, M);

         deepCopy(y, 1);
         deepCopy(x, 1);
         deepCopy(A, 1);
         /*
         #ifdef OMEGA_TARGET_DEVICE
               // Create host mirrors of device views.
               Array1DR8::HostMirror h_y = createHostMirror( d_y );
               Array1DR8::HostMirror h_x = createHostMirror( d_x );
               Array2DR8::HostMirror h_A = createHostMirror( d_A );
         #endif

               // Initialize y vector on host.
               for ( int i = 0; i < N; ++i ) {
                 y( i ) = 1;
               }

               // Initialize x vector on host.
               for ( int i = 0; i < M; ++i ) {
                 x( i ) = 1;
               }

               // Initialize A matrix on host.
               for ( int j = 0; j < N; ++j ) {
                 for ( int i = 0; i < M; ++i ) {
                   A( j, i ) = 1;
                 }
               }

         #ifdef OMEGA_TARGET_DEVICE
               // Deep copy host views to device views.
               Kokkos::deep_copy( d_y, y );
               Kokkos::deep_copy( d_x, x );
               Kokkos::deep_copy( d_A, A );
         #endif
         */

         // Timer products.
         Kokkos::Timer timer;

         for (int repeat = 0; repeat < nrepeat; repeat++) {

            // Application: <y,Ax> = y^T*A*x
            double result = 0;

            parallelReduce(
                "yAx", {N},
                KOKKOS_LAMBDA(int j, double &update) {
                   double temp2 = 0;

                   for (int i = 0; i < M; ++i) {
                      temp2 += A(j, i) * x(i);
                   }

                   update += y(j) * temp2;
                },
                result);

            // Output result.
            if (repeat == (nrepeat - 1)) {
               std::cout << "  Computed result for " << N << " x " << M
                         << " is " << result << std::endl;
            }

            const double solution = (double)N * (double)M;

            if (result != solution) {
               std::cout << "  FAIL: result( " << result << " ) != solution( "
                         << solution << " )" << std::endl;
               RetVal -= -1;
            }

            // Calculate time.
            double time = timer.seconds();

            // Calculate bandwidth.
            // Each matrix A row (each of length M) is read once.
            // The x vector (of length M) is read N times.
            // The y vector (of length N) is read once.
            // double Gbytes = 1.0e-9 * double( sizeof(double) * ( 2 * M * N + N
            // ) );
            double Gbytes = 1.0e-9 * double(sizeof(double) * (M + M * N + N));

            // Print results (problem size, time and bandwidth in GB/s).
            std::cout << "  N( " << N << " ) M( " << M << " ) nrepeat ( "
                      << nrepeat << " ) problem( " << Gbytes * 1000
                      << " MB ) time( " << time << " s ) bandwidth( "
                      << Gbytes * nrepeat / time << " GB/s )" << std::endl;
         }

         std::cout << "OmegaKokkos test: PASS" << std::endl;
      }
      Kokkos::finalize();

   } catch (const std::exception &Ex) {
      std::cout << Ex.what() << ": FAIL" << std::endl;
      RetVal -= -1;
   } catch (...) {
      std::cout << "Unknown: FAIL" << std::endl;
      RetVal -= -1;
   }

   return RetVal;
}
