#include <catch2/catch.hpp>

#include "share/util/scream_vertical_interpolation.hpp"

using namespace scream;
using namespace vinterp;

template <typename S>
using SmallPack = Pack<S,SCREAM_SMALL_PACK_SIZE>;
using Spack = SmallPack<Real>;

using Smask = ekat::Mask<Spack::n>;

template<int P>
void run(){
  //This function first tests cases where the output levels are the same
  //as the input levels to make sure you get the same output as input
  //In this case only two cases are investigated:
  //1) n_layers_src = n_layers_tgt = 2*P
  //2) n_layers_src = n_layers_tgt = 2*P+1
  //where P is the pack size
  //Then the function tests 4 different scenarios where the 
  //the output levels are different from the input levels:
  //1) n_layers_src= 2*P,   n_layers_tgt= 2*P
  //2) n_layers_src= 2*P+1, n_layers_tgt= 2*P
  //3) n_layers_src= 2*P,   n_layers_tgt= 2*P-1
  //4) n_layers_src= 2*P+1, n_layers_tgt= 2*P-1
  //For each scenario target levels are at the midpoint and so should be the average
  //of the source layers. 

  int n_layers_src[4] = {2*P,2*P+1,2*P,2*P+1};
  int n_layers_tgt[4] = {2*P,2*P,2*P-1,2*P-1};

  printf (" -- Testing vertical interpolation with Pack size %d --\n",P);

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
      auto npacks_src = ekat::PackInfo<P>::num_packs(n_layers_src[i]);
      auto npacks_tgt = ekat::PackInfo<P>::num_packs(n_layers_tgt[i]);
      auto p_tgt = view_1d<Pack<Real,P>>("",npacks_tgt);
      auto p_tgt_h = Kokkos::create_mirror_view(p_tgt);
      auto p_tgt_h_s = ekat::scalarize(p_tgt_h);
      auto tmp_src = view_Nd<Pack<Real,P>,2>("",2,npacks_src);
      auto tmp_src_h = Kokkos::create_mirror_view(tmp_src);
      auto tmp_src_h_s = ekat::scalarize(tmp_src_h);
      auto p_src = view_Nd<Pack<Real,P>,2>("",2,npacks_src);
      auto p_src_h = Kokkos::create_mirror_view(p_src);
      auto p_src_h_s = ekat::scalarize(p_src_h);
      auto out = view_Nd<Pack<Real,P>,2>("",2,npacks_tgt);
      auto out_h = Kokkos::create_mirror_view(out);
      auto out_h_s = ekat::scalarize(out_h);
      auto mask = view_Nd<Mask<P>,2>("",2,npacks_tgt);
      auto mask_h = Kokkos::create_mirror_view(mask);
  
      //Set target levels    
      for (int lev=0; lev<(n_layers_tgt[i]); lev++){
        p_tgt_h_s(lev) = 2*lev+perturb;
      } 
      //Set source levels and source input (tmp_src)
      for (int col=0; col<2; col++){
        for (int lev=0; lev<n_layers_src[i]; lev++){
          p_src_h_s(col,lev) = 2.*lev;
          tmp_src_h_s(col,lev) = 100.*(col+1) + 2.*lev;
        }//end of looping over levels
      }//end of looping over columns
  
      Kokkos::deep_copy(p_src,p_src_h);
      Kokkos::deep_copy(p_tgt,p_tgt_h);
      Kokkos::deep_copy(tmp_src,tmp_src_h);
      perform_vertical_interpolation<Real,P,2>(p_src,
                                     p_tgt,
                                     tmp_src,
                                     out,
                                     mask,
                                     n_layers_src[i],
                                     n_layers_tgt[i]);

      Kokkos::deep_copy(out_h,out);
      Kokkos::deep_copy(mask_h,mask);

      //Check that output of interpolation is as expected
      for(int col=0; col<2; col++){
        for(int lev=0; lev<(n_layers_tgt[i]-1); lev++){
          const int ipack = lev / P;
          const int jpack  = lev % P;
          REQUIRE(mask_h(col,ipack)[jpack] == false);

          if (perturb == 0){
            REQUIRE(out_h_s(col,lev) == tmp_src_h_s(col,lev));
          }
          else{
            auto avg_val = (tmp_src_h_s(col,lev+1)+tmp_src_h_s(col,lev))/2.;
            REQUIRE(out_h_s(col,lev) == avg_val);
  	  }
        }//end of looping over levels
      }//end of looping over columns

      auto out_1d_test = view_Nd<Pack<Real,P>,2>("",2,npacks_tgt);
      auto out_1d_test_h = Kokkos::create_mirror_view(out_1d_test);
      auto out_1d_test_h_s = ekat::scalarize(out_1d_test_h);
      auto mask_1d_test = view_Nd<Mask<P>,2>("",2,npacks_tgt);
      auto mask_1d_test_h = Kokkos::create_mirror_view(mask_1d_test);
  
      //Take subview and run through the 1d interpolator function 
      //and make sure get same thing back
      ekat::LinInterp<Real,P> vert_interp(2,n_layers_src[i],n_layers_tgt[i]);
      const int num_vert_packs = p_tgt.extent(0);
      const auto policy = ESU::get_default_team_policy(2, num_vert_packs);
      auto loc_layers_src = n_layers_src[i];
      auto loc_layers_tgt = n_layers_tgt[i];
      Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
         	       KOKKOS_LAMBDA(MemberType const& team) {
          const int icol = team.league_rank();
          const auto x1=ekat::subview(p_src, icol);
          const auto in=ekat::subview(tmp_src, icol);
          const auto out_1d=ekat::subview(out_1d_test, icol);
          const auto msk=ekat::subview(mask_1d_test, icol);
          apply_interpolation_impl_1d<Real,P>(x1,
                                              p_tgt,
                                              in,
                                              out_1d,
                                              msk,
                                              loc_layers_src,
                                              loc_layers_tgt,
                                              icol,
                                              masked_val,
                                              team,
                                              vert_interp); 
      });
      Kokkos::fence();

      Kokkos::deep_copy(out_1d_test_h,out_1d_test);
      Kokkos::deep_copy(mask_1d_test_h,mask_1d_test);

      //Check that 1d interpolator output is consistent with what is expected
      for(int col=0; col<2; col++){
        for(int lev=0; lev<(n_layers_tgt[i]-1); lev++){
          const int ipack = lev / P;
          const int jpack  = lev % P;
          REQUIRE(mask_1d_test_h(col,ipack)[jpack] == false);

          if (perturb == 0){
            REQUIRE(out_1d_test_h_s(col,lev) == tmp_src_h_s(col,lev));
          }
  	  else{
            auto avg_val = (tmp_src_h_s(col,lev+1)+tmp_src_h_s(col,lev))/2.;
            REQUIRE(out_1d_test_h_s(col,lev) == avg_val);
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

template<int P>
void check_mask(const view_Nd_host<Mask<P>,2>& mask, int col, int lev){
  const int ipack = lev / P;
  const int jpack  = lev % P;
  if ((col == 0 && lev == 16) || (col == 1 && lev == 0)){
    REQUIRE(mask(col,ipack)[jpack] == true);
  }
  else{
    REQUIRE(mask(col,ipack)[jpack] == false);
  }
}

TEST_CASE("testing_masking"){
  printf (" -- Testing Masking --\n");
  //This test performs 3 tests:
  //1) That the interpolation is working properly using 2d views, 
  //   including the masking of out-of-bounds values
  //2) It checks the interpolation is working properly with 
  //   a user defined masking value
  //3) It checks that the interpolation is working properly when 
  //   using the 1d interpolation function
  const int n_layers_src = 9;
  const int n_layers_tgt = 17;
  const int P = SCREAM_PACK_SIZE;

  auto npacks_src  = ekat::PackInfo<P>::num_packs(n_layers_src);
  auto npacks_tgt  = ekat::PackInfo<P>::num_packs(n_layers_tgt);
  auto p_tgt       = view_1d<Pack<Real,P>>("",npacks_tgt);
  auto p_tgt_h     = Kokkos::create_mirror_view(p_tgt);
  auto p_tgt_h_s   = ekat::scalarize(p_tgt_h);
  auto tmp_src     = view_Nd<Pack<Real,P>,2>("",2,npacks_src);
  auto tmp_src_h   = Kokkos::create_mirror_view(tmp_src);
  auto tmp_src_h_s = ekat::scalarize(tmp_src_h);
  auto p_src       = view_Nd<Pack<Real,P>,2>("",2,npacks_src);
  auto p_src_h     = Kokkos::create_mirror_view(p_src);
  auto p_src_h_s   = ekat::scalarize(p_src_h);
  auto out         = view_Nd<Pack<Real,P>,2>("",2,npacks_tgt);
  auto out_h       = Kokkos::create_mirror_view(out);
  auto mask        = view_Nd<Mask<P>,2>("",2,npacks_tgt);
  auto mask_h      = Kokkos::create_mirror_view(mask);

  //Fist test to see if interpolate properly using 2d views
  //Also test that when out-of-bounds returns masked values

  //Set target levels from 25-105
  for (int i=0; i<17; i++){
    p_tgt_h_s(i) = 25.0 + i*5.0;
  }
 
  //Set source levels:
  //For column 1 from 20-100
  //For column 2 from 30-110
  //Set input variable, in this case
  //set tmp_src from 200-280 by 10, with final level at 240 so 
  //can test interpolation with increasing/decreasing values
  for (int i=0; i<9; i++){
    p_src_h_s(0,i) = 20.0 + i*10.0;
    p_src_h_s(1,i) = 30.0 + i*10.0;
    if (i<8){
      tmp_src_h_s(0,i) = 200.+i*10.;
      tmp_src_h_s(1,i) = 200.+i*10.;
    }
    else{
      tmp_src_h_s(0,i) = 240.;
      tmp_src_h_s(1,i) = 240.;
    }
  }

  Kokkos::deep_copy(p_src,p_src_h);
  Kokkos::deep_copy(p_tgt,p_tgt_h);
  Kokkos::deep_copy(tmp_src,tmp_src_h);

  perform_vertical_interpolation<Real,P,2>(p_src,
				 p_tgt,
				 tmp_src,
				 out,
				 mask,
				 n_layers_src,
				 n_layers_tgt);

  Kokkos::deep_copy(out_h,out);
  Kokkos::deep_copy(mask_h,mask);

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

  auto out_h_s     = ekat::scalarize(out_h);
  for(int col=0; col<2; col++){
    for(int lev=0; lev<17; lev++){
      REQUIRE(out_h_s(col,lev) == correct_val[col][lev]);
      check_mask<P>(mask_h,col,lev);
    }
  }

  //Test to see if get same answer when call 1D interpolation function 
  //instead of 2D interpolation function
  auto out_1d_test = view_Nd<Pack<Real,P>,2>("",2,npacks_tgt);
  auto out_1d_test_h = Kokkos::create_mirror_view(out_1d_test);
  auto out_1d_test_h_s = ekat::scalarize(out_1d_test_h);
  auto mask_1d_test = view_Nd<Mask<P>,2>("",2,npacks_tgt);
  auto mask_1d_test_h = Kokkos::create_mirror_view(mask_1d_test);

  ekat::LinInterp<Real,P> vert_interp(2,n_layers_src,n_layers_tgt);
  const int num_vert_packs = p_tgt.extent(0);
  const auto policy = ESU::get_default_team_policy(2, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      const auto x1=ekat::subview(p_src, icol);
      const auto in=ekat::subview(tmp_src, icol);
      const auto out_1d=ekat::subview(out_1d_test, icol);
      const auto msk=ekat::subview(mask_1d_test, icol);
      apply_interpolation_impl_1d<Real,P>(x1,
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

  Kokkos::deep_copy(out_1d_test_h,out_1d_test);
  Kokkos::deep_copy(mask_1d_test_h,mask_1d_test);
  
  for(int col=0; col<2; col++){
    for(int lev=0; lev<17; lev++){
      REQUIRE(out_1d_test_h_s(col,lev) == correct_val[col][lev]);
      check_mask<P>(mask_h,col,lev);
    }
  }

  //Check to see if choose different masked value than default that it returns as expected
  auto out_usr_msk = view_Nd<Pack<Real,P>,2>("",2,npacks_tgt);
  auto out_usr_msk_h = Kokkos::create_mirror_view(out_usr_msk);
  auto out_usr_msk_h_s = ekat::scalarize(out_usr_msk_h);
  auto mask_usr_msk = view_Nd<Mask<P>,2>("",2,npacks_tgt);
  auto mask_usr_msk_h = Kokkos::create_mirror_view(mask_usr_msk);
  Real mod_mask_val = -999.;
  perform_vertical_interpolation<Real,P,2>(p_src,
				 p_tgt,
				 tmp_src,
				 out_usr_msk,
				 mask_usr_msk,
				 n_layers_src,
				 n_layers_tgt,
                                 mod_mask_val);
  correct_val[0][16] = -999.;
  correct_val[1][0] = -999.;

  Kokkos::deep_copy(out_usr_msk_h,out_usr_msk);
  Kokkos::deep_copy(mask_usr_msk_h,mask_usr_msk);

  for(int col=0; col<2; col++){
    for(int lev=0; lev<17; lev++){
      REQUIRE(out_usr_msk_h_s(col,lev) == correct_val[col][lev]);
      check_mask<P>(mask_h,col,lev);
    }
  }

}

