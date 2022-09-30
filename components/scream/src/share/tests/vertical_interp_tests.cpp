#include <catch2/catch.hpp>

#include "share/util/scream_vertical_interpolation.hpp"

TEST_CASE("main_vertical_interpolation_test"){

  using namespace scream;
  using namespace vinterp;
  const int n_layers_src = 3;
  const int n_layers_tgt = 4;
  
  auto npacks_src = ekat::PackInfo<Spack::n>::num_packs(n_layers_src);
  auto npacks_tgt = ekat::PackInfo<Spack::n>::num_packs(n_layers_tgt);
  auto p_tgt = view_1d<Spack>("",npacks_tgt);
  auto p_tgt_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt));
  auto tmp_src = view_2d<Spack>("",3,npacks_src);
  auto tmp_src_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src));
  auto p_src = view_2d<Spack>("",3,npacks_src);
  auto p_src_s = Kokkos::create_mirror_view(ekat::scalarize(p_src));
  auto out = view_2d<Spack>("",3,npacks_tgt);
  auto out_s = Kokkos::create_mirror_view(ekat::scalarize(out));
  auto mask = view_2d<Smask>("",3,npacks_tgt);

  p_tgt_s(0) = 20.0;
  p_tgt_s(1) = 25.0;
  p_tgt_s(2) = 35.0;
  p_tgt_s(3) = 45.0;

  p_src_s(0,0) = 20.0;
  p_src_s(0,1) = 30.0;
  p_src_s(0,2) = 40.0;

  p_src_s(1,0) = 30.0;
  p_src_s(1,1) = 40.0;
  p_src_s(1,2) = 50.0;

  p_src_s(2,0) = 10.0;
  p_src_s(2,1) = 20.0;
  p_src_s(2,2) = 30.0;

  for (int i=0; i<3; i++){
    tmp_src_s(i,0)=260.;
    tmp_src_s(i,1)=270.;
    tmp_src_s(i,2)=240.;
  }

  perform_vertical_interpolation(p_src,
				 p_tgt,
				 tmp_src,
				 out,
				 mask,
				 n_layers_src,
				 n_layers_tgt);
  
  Real correct_val[3][4];
  correct_val[0][0] = 260.0;
  correct_val[0][1] = 265.0;
  correct_val[0][2] = 255.0;
  correct_val[0][3] = masked_val;
  correct_val[1][0] = masked_val;
  correct_val[1][1] = masked_val;
  correct_val[1][2] = 265.0;
  correct_val[1][3] = 255.0;
  correct_val[2][0] = 270.0;
  correct_val[2][1] = 255.0;
  correct_val[2][2] = masked_val;
  correct_val[2][3] = masked_val;
  
  for(int col=0; col<3; col++){
    for(int lev=0; lev<4; lev++){
      REQUIRE(out_s(col,lev) == correct_val[col][lev]);
    }
  }


  //Test to see if 1-D works
  auto out_1d_test = view_2d<Spack>("",3,npacks_tgt);
  auto out_1d_test_s = Kokkos::create_mirror_view(ekat::scalarize(out_1d_test));
  auto mask_1d_test = view_2d<Smask>("",3,npacks_tgt);


  ekat::LinInterp<Real,Spack::n> vert_interp(3,n_layers_src,n_layers_tgt);
  const int num_vert_packs = p_tgt.extent(0);
  const auto policy = ESU::get_default_team_policy(3, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      view_1d<Spack> x1=ekat::subview(p_src, icol);
      view_1d<Spack> in=ekat::subview(tmp_src, icol);
      view_1d<Spack> out_1d=ekat::subview(out_1d_test, icol);
      view_1d<Smask> msk=ekat::subview(mask_1d_test, icol);
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
  
  for(int col=0; col<3; col++){
    for(int lev=0; lev<4; lev++){
      REQUIRE(out_1d_test_s(col,lev) == correct_val[col][lev]);
    }
  }



  //Check to see if choose different masked value than default that it returns as expected
  auto out_same = view_2d<Spack>("",3,npacks_tgt);
  auto out_same_s = Kokkos::create_mirror_view(ekat::scalarize(out_same));
  auto mask_same = view_2d<Smask>("",3,npacks_tgt);
  Real mod_mask_val = -999.;
  perform_vertical_interpolation(p_src,
				 p_tgt,
				 tmp_src,
				 out_same,
				 mask_same,
				 n_layers_src,
				 n_layers_tgt,
                                 mod_mask_val);
  correct_val[0][3] = mod_mask_val;
  correct_val[1][0] = mod_mask_val;
  correct_val[1][1] = mod_mask_val;
  correct_val[2][2] = mod_mask_val;
  correct_val[2][3] = mod_mask_val;


  for(int col=0; col<3; col++){
    for(int lev=0; lev<4; lev++){
      REQUIRE(out_same_s(col,lev) == correct_val[col][lev]);
    }
  }
  

  
}


TEST_CASE("pack_size_test"){

  //Testing that get expected answer for nlevs % packsize==0 and nlevs%packsize!=0
  //Can I use an arbitrary pack size like SCREAM SMALL Pack size?
  //Update this to take arbitrary pack size (i.e. pack size from SCREAM Small), 
  //check that it can run with 1 or 16
  using namespace scream;
  using namespace vinterp;

  //First test looks at case when nlevs % packsize==0 for both target and source
  const int N = SCREAM_SMALL_PACK_SIZE;
  //const int N = 16;
  //const int n_layers_src = 8;
  //const int n_layers_tgt = 8;
  const int n_layers_src = 2*N;
  const int n_layers_tgt = 2*N;


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
  auto mask = view_2d<Pmask<N>>("",2,npacks_tgt);

  for (int lev=0; lev<(n_layers_tgt-1); lev++){
    p_tgt_s(lev) = lev*2+1;
  }
  //p_tgt_s(2*N-1) = 14.;
  p_tgt_s(2*N-1) = 2*(n_layers_tgt-1);
  for (int col=0; col<2; col++){
    for (int lev=0; lev<n_layers_src; lev++){
      p_src_s(col,lev) = 2.*lev;
      tmp_src_s(col,lev) = 100.*(col+1) + 2.*lev;
    }
  }
  perform_vertical_interpolation(p_src,
				 p_tgt,
				 tmp_src,
				 out,
				 mask,
				 n_layers_src,
				 n_layers_tgt);
  for(int col=0; col<2; col++){
    for(int lev=0; lev<(n_layers_tgt-1); lev++){
      auto avg_val = (tmp_src_s(col,lev+1)+tmp_src_s(col,lev))/2.;
      if (lev < (n_layers_tgt-1)){
        REQUIRE(out_s(col,lev) == avg_val);
      }
      else{
        REQUIRE(out_s(col,lev) == tmp_src_s(col,lev));
      }
    }
  }

  //Second test looks at case when nlevs % packsize!=0 for source 
  //and nlevs % packsize==0 for target
  //const int n_layers_src_2 = 9;
  //const int n_layers_tgt_2 = 8;
  const int n_layers_src_2 = 2*N+1;
  const int n_layers_tgt_2 = 2*N;

  auto npacks_src_2 = ekat::PackInfo<N>::num_packs(n_layers_src_2);
  auto npacks_tgt_2 = ekat::PackInfo<N>::num_packs(n_layers_tgt_2);
  auto p_tgt_2 = view_1d<Pack<Real,N>>("",npacks_tgt_2);
  auto p_tgt_2_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt_2));
  auto tmp_src_2 = view_2d<Pack<Real,N>>("",2,npacks_src_2);
  auto tmp_src_2_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src_2));
  auto p_src_2 = view_2d<Pack<Real,N>>("",2,npacks_src_2);
  auto p_src_2_s = Kokkos::create_mirror_view(ekat::scalarize(p_src_2));
  auto out_2 = view_2d<Pack<Real,N>>("",2,npacks_tgt_2);
  auto out_2_s = Kokkos::create_mirror_view(ekat::scalarize(out_2));
  auto mask_2 = view_2d<Pmask<N>>("",2,npacks_tgt_2);

  for (int lev=0; lev<(n_layers_tgt_2-1); lev++){
    p_tgt_2_s(lev) = lev*2+1;
  }
  p_tgt_2_s(n_layers_tgt_2-1) = 2*n_layers_tgt_2;
  for (int col=0; col<2; col++){
    for (int lev=0; lev<n_layers_src_2; lev++){
      p_src_2_s(col,lev) = 2.*lev;
      tmp_src_2_s(col,lev) = 100.*(col+1) + 2.*lev;
    }//end of looping over levels
  }//end of looping over columns
  //std::cout<<"n_layers_src_2: "<<n_layers_src_2<<std::endl;
  perform_vertical_interpolation(p_src_2,
				 p_tgt_2,
				 tmp_src_2,
				 out_2,
				 mask_2,
				 n_layers_src_2,
				 n_layers_tgt_2);
  for(int col=0; col<2; col++){
    for(int lev=0; lev<n_layers_tgt_2; lev++){
      auto avg_val = (tmp_src_2_s(col,lev+1)+tmp_src_2_s(col,lev))/2.;
      if (lev < (n_layers_tgt_2-1)){
        REQUIRE(out_2_s(col,lev) == avg_val);
	//std::cout<<"out_2_s(col,lev)"<<out_2_s(col,lev)<<std::endl;
      }
      else if(lev == (n_layers_tgt_2-1)){
        REQUIRE(out_2_s(col,lev) == tmp_src_2_s(col,lev+1));
	//std::cout<<"out_2_s(col,lev) lev 7: "<<out_2_s(col,lev)<<std::endl;
      }
    }//end of looping over levels
  }//end of looping over columns

  
  const int n_layers_src_3 = 2*N;
  const int n_layers_tgt_3 = 2*N-1;

  auto npacks_src_3 = ekat::PackInfo<N>::num_packs(n_layers_src_3);
  auto npacks_tgt_3 = ekat::PackInfo<N>::num_packs(n_layers_tgt_3);
  auto p_tgt_3 = view_1d<Pack<Real,N>>("",npacks_tgt_3);
  auto p_tgt_3_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt_3));
  auto tmp_src_3 = view_2d<Pack<Real,N>>("",2,npacks_src_3);
  auto tmp_src_3_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src_3));
  auto p_src_3 = view_2d<Pack<Real,N>>("",2,npacks_src_3);
  auto p_src_3_s = Kokkos::create_mirror_view(ekat::scalarize(p_src_3));
  auto out_3 = view_2d<Pack<Real,N>>("",2,npacks_tgt_3);
  auto out_3_s = Kokkos::create_mirror_view(ekat::scalarize(out_3));
  auto mask_3 = view_2d<Pmask<N>>("",2,npacks_tgt_3);

  for (int lev=0; lev<(n_layers_tgt_3-1); lev++){
    p_tgt_3_s(lev) = lev*2+1;
  }
  p_tgt_3_s(n_layers_tgt_3-1) = 2*n_layers_tgt_3;
  for (int col=0; col<2; col++){
    for (int lev=0; lev<n_layers_src_3; lev++){
      p_src_3_s(col,lev) = 2.*lev;
      tmp_src_3_s(col,lev) = 100.*(col+1) + 2.*lev;
    }//end of looping over levels
  }//end of looping over columns
  //std::cout<<"n_layers_src_3: "<<n_layers_src_3<<std::endl;
  perform_vertical_interpolation(p_src_3,
				 p_tgt_3,
				 tmp_src_3,
				 out_3,
				 mask_3,
				 n_layers_src_3,
				 n_layers_tgt_3);
  for(int col=0; col<2; col++){
    for(int lev=0; lev<n_layers_tgt_3; lev++){
      auto avg_val = (tmp_src_3_s(col,lev+1)+tmp_src_3_s(col,lev))/2.;
      if (lev < (n_layers_tgt_3-1)){
        REQUIRE(out_3_s(col,lev) == avg_val);
	//std::cout<<"out_3_s(col,lev)"<<out_3_s(col,lev)<<std::endl;
      }
      else if(lev == (n_layers_tgt_3-1)){
        REQUIRE(out_3_s(col,lev) == tmp_src_3_s(col,lev+1));
	//std::cout<<"out_3_s(col,lev) lev 7: "<<out_3_s(col,lev)<<std::endl;
      }
    }//end of looping over levels
  }//end of looping over columns
  
  const int n_layers_src_4 = 2*N+1;
  const int n_layers_tgt_4 = 2*N-1;

  auto npacks_src_4 = ekat::PackInfo<N>::num_packs(n_layers_src_4);
  auto npacks_tgt_4 = ekat::PackInfo<N>::num_packs(n_layers_tgt_4);
  auto p_tgt_4 = view_1d<Pack<Real,N>>("",npacks_tgt_4);
  auto p_tgt_4_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt_4));
  auto tmp_src_4 = view_2d<Pack<Real,N>>("",2,npacks_src_4);
  auto tmp_src_4_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src_4));
  auto p_src_4 = view_2d<Pack<Real,N>>("",2,npacks_src_4);
  auto p_src_4_s = Kokkos::create_mirror_view(ekat::scalarize(p_src_4));
  auto out_4 = view_2d<Pack<Real,N>>("",2,npacks_tgt_4);
  auto out_4_s = Kokkos::create_mirror_view(ekat::scalarize(out_4));
  auto mask_4 = view_2d<Pmask<N>>("",2,npacks_tgt_4);

  for (int lev=0; lev<(n_layers_tgt_4-1); lev++){
    p_tgt_4_s(lev) = lev*2+1;
  }
  p_tgt_4_s(n_layers_tgt_4-1) = 2*n_layers_tgt_4;
  for (int col=0; col<2; col++){
    for (int lev=0; lev<n_layers_src_4; lev++){
      p_src_4_s(col,lev) = 2.*lev;
      tmp_src_4_s(col,lev) = 100.*(col+1) + 2.*lev;
    }//end of looping over levels
  }//end of looping over columns
  //std::cout<<"n_layers_src_4: "<<n_layers_src_4<<std::endl;
  perform_vertical_interpolation(p_src_4,
				 p_tgt_4,
				 tmp_src_4,
				 out_4,
				 mask_4,
				 n_layers_src_4,
				 n_layers_tgt_4);
  for(int col=0; col<2; col++){
    for(int lev=0; lev<n_layers_tgt_4; lev++){
      auto avg_val = (tmp_src_4_s(col,lev+1)+tmp_src_4_s(col,lev))/2.;
      if (lev < (n_layers_tgt_4-1)){
        REQUIRE(out_4_s(col,lev) == avg_val);
	//std::cout<<"out_4_s(col,lev)"<<out_4_s(col,lev)<<std::endl;
      }
      else if(lev == (n_layers_tgt_4-1)){
        REQUIRE(out_4_s(col,lev) == tmp_src_4_s(col,lev+1));
	//std::cout<<"out_4_s(col,lev) lev 7: "<<out_4_s(col,lev)<<std::endl;
      }
    }//end of looping over levels
  }//end of looping over columns

}//end of pack size test



TEST_CASE("same_output_grid_test"){

  //Check to see if same levels for target and source returns same data
  using namespace scream;
  using namespace vinterp;
  const int n_layers_src = 4;
  const int n_layers_tgt = 4;

  auto npacks_src = ekat::PackInfo<Spack::n>::num_packs(n_layers_src);
  auto npacks_tgt = ekat::PackInfo<Spack::n>::num_packs(n_layers_tgt);
  auto p_tgt = view_1d<Spack>("",npacks_tgt);
  auto p_tgt_s = Kokkos::create_mirror_view(ekat::scalarize(p_tgt));
  auto tmp_src = view_2d<Spack>("",2,npacks_src);
  auto tmp_src_s = Kokkos::create_mirror_view(ekat::scalarize(tmp_src));
  auto p_src = view_2d<Spack>("",2,npacks_src);
  auto p_src_s = Kokkos::create_mirror_view(ekat::scalarize(p_src));
  auto out = view_2d<Spack>("",2,npacks_tgt);
  auto out_s = Kokkos::create_mirror_view(ekat::scalarize(out));
  auto mask = view_2d<Smask>("",2,npacks_tgt);


  p_tgt_s(0) = 20.0;
  p_tgt_s(1) = 25.0;
  p_tgt_s(2) = 35.0;
  p_tgt_s(3) = 45.0;

  p_src_s(0,0) = 20.0;
  p_src_s(0,1) = 25.0;
  p_src_s(0,2) = 35.0;
  p_src_s(0,3) = 45.0;

  p_src_s(1,0) = 20.0;
  p_src_s(1,1) = 25.0;
  p_src_s(1,2) = 35.0;
  p_src_s(1,3) = 45.0;

  for (int i=0; i<2; i++){
    tmp_src_s(i,0)=260.;
    tmp_src_s(i,1)=270.;
    tmp_src_s(i,2)=240.;
    tmp_src_s(i,3)=250.;
  }

  perform_vertical_interpolation(p_src,
				 p_tgt,
				 tmp_src,
				 out,
				 n_layers_src,
				 n_layers_tgt);
  
  for(int col=0; col<2; col++){
    for(int lev=0; lev<4; lev++){
      REQUIRE(out_s(col,lev) == tmp_src_s(col,lev));
    }
  }


  //Take subview and run through the 1d interpolator function and make sure get same thing back
  ekat::LinInterp<Real,Spack::n> vert_interp(1,n_layers_src,n_layers_tgt);
  const int num_vert_packs = p_tgt.extent(0);
  const auto policy = ESU::get_default_team_policy(1, num_vert_packs);
  Kokkos::parallel_for("scream_vert_interp_setup_loop", policy,
     	       KOKKOS_LAMBDA(MemberType const& team) {
      const int icol = team.league_rank();
      view_1d<Spack> x1=ekat::subview(p_src, icol);
      view_1d<Spack> in=ekat::subview(tmp_src, icol);
      view_1d<Spack> out_1d=ekat::subview(out, icol);
      view_1d<Smask> msk=ekat::subview(mask, icol);
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
    for(int lev=0; lev<4; lev++){
      REQUIRE(out_s(col,lev) == tmp_src_s(col,lev));
    }
  }

}
