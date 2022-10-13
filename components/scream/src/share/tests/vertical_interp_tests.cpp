#include <catch2/catch.hpp>

#include "share/util/scream_vertical_interpolation.hpp"

using namespace scream;
using namespace vinterp;

template <typename S>
using SmallPack = Pack<S,SCREAM_SMALL_PACK_SIZE>;
using Spack = SmallPack<Real>;

using Smask = ekat::Mask<Spack::n>;

template<int N>
void run(){
  //This function first tests cases where the output levels are the same
  //as the input levels to make sure you get the same output as input
  //In this case only two cases are investigated:
  //1) n_layers_src = n_layers_tgt = 2*N
  //2) n_layers_src = n_layers_tgt = 2*N+1
  //where N is the pack size
  //Then the function tests 4 different scenarios where the 
  //the output levels are different from the input levels:
  //1) n_layers_src= 2*N,   n_layers_tgt= 2*N
  //2) n_layers_src= 2*N+1, n_layers_tgt= 2*N
  //3) n_layers_src= 2*N,   n_layers_tgt= 2*N-1
  //4) n_layers_src= 2*N+1, n_layers_tgt= 2*N-1
  //For each scenario target levels are at the midpoint and so should be the average
  //of the source layers. 

  int n_layers_src[4] = {2*N,2*N+1,2*N,2*N+1};
  int n_layers_tgt[4] = {2*N,2*N,2*N-1,2*N-1};

  printf (" -- Testing vertical interpolation with Pack size %d --\n",N);

  for (double perturb : {0, 1}){
    printf ("    -> Target pressure levels: p_tgt = p_src + %f\n",perturb);
    for (int i=0; i<4; i++){
      if (perturb == 0){
        //The 3rd and 4th test are redundant in the case of the
        //same target levels as source levels
        if(i > 1){ break;}
        //In this case the target levels are the same as the source levels
        n_layers_tgt[i] = n_layers_src[i];

        printf ("      -> Testing %d source layers, %d target layers\n",
                n_layers_src[i],n_layers_tgt[i]);
      }
      else{
        printf ("      -> Testing %d source layers, %d target layers\n",
                n_layers_src[i],n_layers_tgt[i]);
      }
      auto npacks_src = ekat::PackInfo<N>::num_packs(n_layers_src[i]);
      auto npacks_tgt = ekat::PackInfo<N>::num_packs(n_layers_tgt[i]);
      auto p_tgt = view_1d<Pack<Real,N>>("",npacks_tgt);
      auto p_tgt_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt));
      auto tmp_src = view_2d<Pack<Real,N>>("",2,npacks_src);
      auto tmp_src_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src));
      auto p_src = view_2d<Pack<Real,N>>("",2,npacks_src);
      auto p_src_s = Kokkos::create_mirror_view(ekat::scalarize(p_src));
      auto out = view_2d<Pack<Real,N>>("",2,npacks_tgt);
      auto out_s = Kokkos::create_mirror_view(ekat::scalarize(out));
      auto mask = view_2d<Mask<N>>("",2,npacks_tgt);
  
      //Set target levels    
      for (int lev=0; lev<(n_layers_tgt[i]-1); lev++){
        p_tgt_s(lev) = 2*lev+perturb;
      } 
      //Set source levels and source input (tmp_src)
      for (int col=0; col<2; col++){
        for (int lev=0; lev<n_layers_src[i]; lev++){
          p_src_s(col,lev) = 2.*lev;
          tmp_src_s(col,lev) = 100.*(col+1) + 2.*lev;
        }//end of looping over levels
      }//end of looping over columns
  
      perform_vertical_interpolation(p_src,
                                     p_tgt,
                                     tmp_src,
                                     out,
                                     mask,
                                     n_layers_src[i],
                                     n_layers_tgt[i]);

      //Check that output of interpolation is as expected
      for(int col=0; col<2; col++){
        for(int lev=0; lev<(n_layers_tgt[i]-1); lev++){
          if (perturb == 0){
            REQUIRE(out_s(col,lev) == tmp_src_s(col,lev));
          }
          else{
            auto avg_val = (tmp_src_s(col,lev+1)+tmp_src_s(col,lev))/2.;
            REQUIRE(out_s(col,lev) == avg_val);
  	  }
        }//end of looping over levels
      }//end of looping over columns
  
      //Take subview and run through the 1d interpolator function 
      //and make sure get same thing back
      ekat::LinInterp<Real,N> vert_interp(2,n_layers_src[i],n_layers_tgt[i]);
      const int num_vert_packs = p_tgt.extent(0);
      const auto policy = ESU::get_default_team_policy(1, num_vert_packs);
      Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
         	       KOKKOS_LAMBDA(MemberType const& team) {
          const int icol = team.league_rank();
          view_1d<Pack<Real,N>> x1=ekat::subview(p_src, icol);
          view_1d<Pack<Real,N>> in=ekat::subview(tmp_src, icol);
          view_1d<Pack<Real,N>> out_1d=ekat::subview(out, icol);
          view_1d<Mask<N>> msk=ekat::subview(mask, icol);
          perform_vertical_interpolation_impl_1d(x1,
                                                 p_tgt,
                                                 in,
                                                 out_1d,
                                                 msk,
                                                 n_layers_src[i],
                                                 n_layers_tgt[i],
                                                 icol,
                                                 masked_val,
                                                 team,
                                                 vert_interp); 
      });
      Kokkos::fence();

      //Check that 1d interpolator output is consistent with what is expected
      for(int col=0; col<2; col++){
        for(int lev=0; lev<(n_layers_tgt[i]-1); lev++){
          if (perturb == 0){
            REQUIRE(out_s(col,lev) == tmp_src_s(col,lev));
          }
  	  else{
            auto avg_val = (tmp_src_s(col,lev+1)+tmp_src_s(col,lev))/2.;
            REQUIRE(out_s(col,lev) == avg_val);
  	  }
        }//end of looping over levels
      }//end of looping over columns
    }//end of looping over pack size tests
  }//end of looping over perturbed or not
}

TEST_CASE("main_vertical_interpolation_test"){
  //This test first checks that you get the same answer if source and target levels 
  //are the same
  //It then tests that you get the expected answer for nlevs % packsize==0 and nlevs%packsize!=0
  //and for different pack sizes. 
  run<1> ();
  if (SCREAM_PACK_SIZE>1) {
    run<SCREAM_PACK_SIZE> ();
  }
}

TEST_CASE("testing_masking"){
  //This test performs 3 tests:
  //1) That the interpolation is working properly using 2d views, 
  //   including the masking of out-of-bounds values
  //2) It checks the interpolation is working properly with 
  //   a user defined masking value
  //3) It checks that the interpolation is working properly when 
  //   using the 1d interpolation function
  const int n_layers_src = 9;
  const int n_layers_tgt = 17;
  const int N = SCREAM_PACK_SIZE;

  auto npacks_src = ekat::PackInfo<N>::num_packs(n_layers_src);
  auto npacks_tgt = ekat::PackInfo<N>::num_packs(n_layers_tgt);
  auto p_tgt = view_1d<Pack<Real,N>>("",npacks_tgt);
  auto p_tgt_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt));
  auto tmp_src = view_2d<Pack<Real,N>>("",2,npacks_src);
  auto tmp_src_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src));
  auto p_src = view_2d<Pack<Real,N>>("",2,npacks_src);
  auto p_src_s = Kokkos::create_mirror_view(ekat::scalarize(p_src));
  auto out = view_2d<Pack<Real,N>>("",2,npacks_tgt);
  auto out_s = Kokkos::create_mirror_view(ekat::scalarize(out));
  auto mask = view_2d<Mask<N>>("",2,npacks_tgt);

  //Fist test to see if interpolate properly using 2d views
  //Also test that when out-of-bounds returns masked values

  //Set target levels from 25-105
  for (int i=0; i<17; i++){
    p_tgt_s(i) = 25.0 + i*5.0;
  }
 
  //Set source levels:
  //For column 1 from 20-100
  //For column 2 from 30-110
  //Set input variable, in this case
  //set tmp_src from 200-280 by 10, with final level at 240 so 
  //can test interpolation with increasing/decreasing values
  for (int i=0; i<9; i++){
    p_src_s(0,i) = 20.0 + i*10.0;
    p_src_s(1,i) = 30.0 + i*10.0;
    if (i<8){
      tmp_src_s(0,i) = 200.+i*10.;
      tmp_src_s(1,i) = 200.+i*10.;
    }
    else{
      tmp_src_s(0,i) = 240.;
      tmp_src_s(1,i) = 240.;
    }
  }

  perform_vertical_interpolation(p_src,
				 p_tgt,
				 tmp_src,
				 out,
				 mask,
				 n_layers_src,
				 n_layers_tgt);

  Real correct_val[2][17];
  for (int i=0; i<15; i++){
    correct_val[0][i] = 205+i*5.;
    correct_val[1][i] = 195+i*5.;
  }

  correct_val[0][14] = 255.;
  correct_val[0][15] = 240.;
  correct_val[0][16] = masked_val;
  correct_val[1][0]  = masked_val;
  correct_val[1][15] = 270.;
  correct_val[1][16] = 255.;

  for(int col=0; col<2; col++){
    for(int lev=0; lev<17; lev++){
      REQUIRE(out_s(col,lev) == correct_val[col][lev]);
    }
  }

  //Test to see if get same answer when call 1D interpolation function 
  //instead of 2D interpolation function
  auto out_1d_test = view_2d<Pack<Real,N>>("",2,npacks_tgt);
  auto out_1d_test_s = Kokkos::create_mirror_view(ekat::scalarize(out_1d_test));
  auto mask_1d_test = view_2d<Mask<N>>("",2,npacks_tgt);

  ekat::LinInterp<Real,N> vert_interp(2,n_layers_src,n_layers_tgt);
  const int num_vert_packs = p_tgt.extent(0);
  const auto policy = ESU::get_default_team_policy(2, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      view_1d<Pack<Real,N>> x1=ekat::subview(p_src, icol);
      view_1d<Pack<Real,N>> in=ekat::subview(tmp_src, icol);
      view_1d<Pack<Real,N>> out_1d=ekat::subview(out_1d_test, icol);
      view_1d<Mask<N>> msk=ekat::subview(mask_1d_test, icol);
      perform_vertical_interpolation_impl_1d(x1,
                                             p_tgt,
                                             in,
                                             out_1d,
                                             msk,
                                             n_layers_src,
                                             n_layers_tgt,
                                             icol,
                                             masked_val,
                                             team,
                                             vert_interp);
  });
  Kokkos::fence();
  
  for(int col=0; col<2; col++){
    for(int lev=0; lev<17; lev++){
      REQUIRE(out_1d_test_s(col,lev) == correct_val[col][lev]);
    }
  }

  //Check to see if choose different masked value than default that it returns as expected
  auto out_same = view_2d<Pack<Real,N>>("",2,npacks_tgt);
  auto out_same_s = Kokkos::create_mirror_view(ekat::scalarize(out_same));
  auto mask_same = view_2d<Mask<N>>("",2,npacks_tgt);
  Real mod_mask_val = -999.;
  perform_vertical_interpolation(p_src,
				 p_tgt,
				 tmp_src,
				 out_same,
				 mask_same,
				 n_layers_src,
				 n_layers_tgt,
                                 mod_mask_val);
  correct_val[0][16] = -999.;
  correct_val[1][0] = -999.;

  for(int col=0; col<2; col++){
    for(int lev=0; lev<17; lev++){
      REQUIRE(out_same_s(col,lev) == correct_val[col][lev]);
    }
  }

}

