#ifndef LAPIS_MODULE_H
#define LAPIS_MODULE_H
#include <Kokkos_Core.hpp>
extern Kokkos::View<float[72][64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_72x64xf32;
extern Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_72xf32_0;
extern Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_64xf32_0;
extern Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_72xf32;
extern Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_64xf32;
extern Kokkos::View<float[64][72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_64x72xf32_0;
extern Kokkos::View<float[144][64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_144x64xf32;
extern Kokkos::View<float[64][72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> __constant_64x72xf32;
struct GlobalViews_forward {
  GlobalViews_forward() {
    m__constant_72x64xf32 = __constant_72x64xf32;
    m__constant_72xf32_0 = __constant_72xf32_0;
    m__constant_64xf32_0 = __constant_64xf32_0;
    m__constant_72xf32 = __constant_72xf32;
    m__constant_64xf32 = __constant_64xf32;
    m__constant_64x72xf32_0 = __constant_64x72xf32_0;
    m__constant_144x64xf32 = __constant_144x64xf32;
    m__constant_64x72xf32 = __constant_64x72xf32;
  }
  Kokkos::View<float[72][64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_72x64xf32;
  Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_72xf32_0;
  Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_64xf32_0;
  Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_72xf32;
  Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_64xf32;
  Kokkos::View<float[64][72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_64x72xf32_0;
  Kokkos::View<float[144][64], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_144x64xf32;
  Kokkos::View<float[64][72], Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> m__constant_64x72xf32;
};
// Return the total amount of scratch that the function uses
// This is the upper bound to what forward_L0_scratch_required can return
KOKKOS_INLINE_FUNCTION constexpr int forward_total_scratch_required() {
  return 2024;
}

// Find L1 address shift: for allocations spilling into level 1 scratch,
// this is subtracted from the allocation address to find its relative address within L1.
KOKKOS_INLINE_FUNCTION constexpr int forward_L1_shift(int L0_scratch_max) {
  int tmp = 2024;
  if(328 > L0_scratch_max && 72 < tmp)
    tmp = 72;
  if(584 > L0_scratch_max && 328 < tmp)
    tmp = 328;
  if(840 > L0_scratch_max && 584 < tmp)
    tmp = 584;
  if(872 > L0_scratch_max && 584 < tmp)
    tmp = 584;
  if(1160 > L0_scratch_max && 872 < tmp)
    tmp = 872;
  if(1160 > L0_scratch_max && 872 < tmp)
    tmp = 872;
  if(72 > L0_scratch_max && 0 < tmp)
    tmp = 0;
  if(1448 > L0_scratch_max && 872 < tmp)
    tmp = 872;
  if(2024 > L0_scratch_max && 1448 < tmp)
    tmp = 1448;
  if(1448 > L0_scratch_max && 872 < tmp)
    tmp = 872;
  return tmp;
}

// Find the actual level 0 scratch required by the function, assuming this limit.
// The answer will be at most L0_scratch_max but never larger.
KOKKOS_INLINE_FUNCTION constexpr int forward_L0_scratch_required(int L0_scratch_max) {
  int tmp = 0;
  if(328 <= L0_scratch_max && 328 > tmp)
    tmp = 328;
  if(584 <= L0_scratch_max && 584 > tmp)
    tmp = 584;
  if(840 <= L0_scratch_max && 840 > tmp)
    tmp = 840;
  if(872 <= L0_scratch_max && 872 > tmp)
    tmp = 872;
  if(1160 <= L0_scratch_max && 1160 > tmp)
    tmp = 1160;
  if(1160 <= L0_scratch_max && 1160 > tmp)
    tmp = 1160;
  if(72 <= L0_scratch_max && 72 > tmp)
    tmp = 72;
  if(1448 <= L0_scratch_max && 1448 > tmp)
    tmp = 1448;
  if(2024 <= L0_scratch_max && 2024 > tmp)
    tmp = 2024;
  if(1448 <= L0_scratch_max && 1448 > tmp)
    tmp = 1448;
  return tmp;
}

// Find the actual level 1 scratch required by the function, assuming this limit for L0 scratch.
// This has no strict upper bound.
KOKKOS_INLINE_FUNCTION constexpr int forward_L1_scratch_required(int L0_scratch_max) {
  return 2024 - forward_L1_shift(L0_scratch_max);
}

template<typename ExecSpace, int L0_scratch_max, typename ViewArg0, typename ViewArg1, typename ViewArg2, typename ViewArg3>
KOKKOS_INLINE_FUNCTION void forward(const typename Kokkos::TeamPolicy<ExecSpace>::member_type& team, const GlobalViews_forward& globals, const ViewArg0& v1, const ViewArg1& v2, const ViewArg2& v3, const ViewArg3& v4, char* scratch0, char* scratch1) {
  constexpr int l1_cutoff = forward_L1_shift(L0_scratch_max);
  const auto& v5 = globals.m__constant_64x72xf32;
  const auto& v6 = globals.m__constant_144x64xf32;
  const auto& v7 = globals.m__constant_64x72xf32_0;
  const auto& v8 = globals.m__constant_64xf32;
  const auto& v9 = globals.m__constant_72xf32;
  const auto& v10 = globals.m__constant_64xf32_0;
  const auto& v11 = globals.m__constant_72xf32_0;
  const auto& v12 = globals.m__constant_72x64xf32;
  constexpr bool v13_spill = 328 > L0_scratch_max;
  Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v13((float*) (v13_spill ? (scratch1 + 72 - l1_cutoff) : (scratch0 + 72)));
  ;
  constexpr bool v14_spill = 584 > L0_scratch_max;
  Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v14((float*) (v14_spill ? (scratch1 + 328 - l1_cutoff) : (scratch0 + 328)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 64),
  [=](size_t v15) {
    v14(v15) = 0.0e+00f;
  });
  team.team_barrier();
  constexpr bool v16_spill = 840 > L0_scratch_max;
  Kokkos::View<float[64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v16((float*) (v16_spill ? (scratch1 + 584 - l1_cutoff) : (scratch0 + 584)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 64),
  [=](size_t v17) {
    float v18 = v14(v17);
    v16(v17) = v18;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, 64),
  [=](size_t v19) {
    float v20 = v16(v19);
    float v21;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 72),
    [=](size_t v22, float& lreduce1) {
      float v23 = v3(v22);
      float v24 = v12(v22, v19);
      float v25 = v23 * v24;
      float v26 = lreduce1 + v25;
      lreduce1 = v26;
      ;
    }, (v21));
    Kokkos::Sum<float> v21_joiner(v21);
    v21_joiner.join(v21, v20);
    ;
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v16(v19) = v21;
    });
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 64),
  [=](size_t v27) {
    float v28 = v16(v27);
    float v29 = v8(v27);
    float v30 = v28 + v29;
    v13(v27) = v30;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 64),
  [=](size_t v31) {
    float v32 = v13(v31);
    bool v33 = (Kokkos::isnan(v32) || Kokkos::isnan(0.0e+00f)) || (v32 > 0.0e+00f);
    float v34 = v33? v32 : 0.0e+00f;
    v13(v31) = v34;
  });
  team.team_barrier();
  constexpr bool v35_spill = 872 > L0_scratch_max;
  Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v35((float*) (v35_spill ? (scratch1 + 584 - l1_cutoff) : (scratch0 + 584)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v36) {
    v35(v36) = 0.0e+00f;
  });
  team.team_barrier();
  constexpr bool v37_spill = 1160 > L0_scratch_max;
  Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v37((float*) (v37_spill ? (scratch1 + 872 - l1_cutoff) : (scratch0 + 872)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v38) {
    float v39 = v35(v38);
    v37(v38) = v39;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, 72),
  [=](size_t v40) {
    float v41 = v37(v40);
    float v42;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 64),
    [=](size_t v43, float& lreduce1) {
      float v44 = v13(v43);
      float v45 = v7(v43, v40);
      float v46 = v44 * v45;
      float v47 = lreduce1 + v46;
      lreduce1 = v47;
      ;
    }, (v42));
    Kokkos::Sum<float> v42_joiner(v42);
    v42_joiner.join(v42, v41);
    ;
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v37(v40) = v42;
    });
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v48) {
    float v49 = v37(v48);
    float v50 = v9(v48);
    float v51 = v49 + v50;
    v2(v48) = v51;
  });
  team.team_barrier();
  constexpr bool v52_spill = 1160 > L0_scratch_max;
  Kokkos::View<float[72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v52((float*) (v52_spill ? (scratch1 + 872 - l1_cutoff) : (scratch0 + 872)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v53) {
    float v54 = v2(v53);
    float v55 = -v54;
    float v56 = Kokkos::exp(v55);
    float v57 = v56 + 1.000000000e+00f;
    float v58 = 1.000000000e+00f / v57;
    v52(v53) = v58;
  });
  team.team_barrier();
  constexpr bool v59_spill = 72 > L0_scratch_max;
  Kokkos::View<bool[72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v59((bool*) (v59_spill ? (scratch1 + 0 - l1_cutoff) : (scratch0 + 0)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v60) {
    float v61 = v52(v60);
    double v62 = (double) v61;
    bool v63 = (v62 > 5.00000000000000000e-01);
    v59(v60) = v63;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v64) {
    bool v65 = v59(v64);
    float v66 = (float) v65;
    v2(v64) = v66;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v67) {
    float v68 = v2(v67);
    float v69 = v52(v67);
    float v70 = v68 - v69;
    v2(v67) = v70;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v71) {
    float v72 = v2(v71);
    float v73 = v52(v71);
    float v74 = v72 + v73;
    v1(v71) = v74;
  });
  team.team_barrier();
  constexpr bool v75_spill = 1448 > L0_scratch_max;
  Kokkos::View<float[144], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v75((float*) (v75_spill ? (scratch1 + 872 - l1_cutoff) : (scratch0 + 872)));
  ;
  constexpr bool v76_spill = 2024 > L0_scratch_max;
  Kokkos::View<float[144], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v76((float*) (v76_spill ? (scratch1 + 1448 - l1_cutoff) : (scratch0 + 1448)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 144),
  [=](size_t v77) {
    float v78 = v75(v77);
    v76(v77) = v78;
  });
  team.team_barrier();
  Kokkos::LayoutStride v79_layout(72, 1);
  Kokkos::View<float[72], Kokkos::LayoutStride, Kokkos::AnonymousSpace> v79(v76.data() + 0, v79_layout);
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v80) {
    float v81 = v4(v80);
    v79(v80) = v81;
  });
  team.team_barrier();
  constexpr bool v82_spill = 1448 > L0_scratch_max;
  Kokkos::View<float[144], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v82((float*) (v82_spill ? (scratch1 + 872 - l1_cutoff) : (scratch0 + 872)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 144),
  [=](size_t v83) {
    float v84 = v76(v83);
    v82(v83) = v84;
  });
  team.team_barrier();
  Kokkos::LayoutStride v85_layout(72, 1);
  Kokkos::View<float[72], Kokkos::LayoutStride, Kokkos::AnonymousSpace> v85(v82.data() + 72, v85_layout);
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v86) {
    float v87 = v1(v86);
    v85(v86) = v87;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, 64),
  [=](size_t v88) {
    float v89 = v14(v88);
    float v90;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 144),
    [=](size_t v91, float& lreduce1) {
      float v92 = v82(v91);
      float v93 = v6(v91, v88);
      float v94 = v92 * v93;
      float v95 = lreduce1 + v94;
      lreduce1 = v95;
      ;
    }, (v90));
    Kokkos::Sum<float> v90_joiner(v90);
    v90_joiner.join(v90, v89);
    ;
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v14(v88) = v90;
    });
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 64),
  [=](size_t v96) {
    float v97 = v14(v96);
    float v98 = v10(v96);
    float v99 = v97 + v98;
    v13(v96) = v99;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 64),
  [=](size_t v100) {
    float v101 = v13(v100);
    bool v102 = (Kokkos::isnan(v101) || Kokkos::isnan(0.0e+00f)) || (v101 > 0.0e+00f);
    float v103 = v102? v101 : 0.0e+00f;
    v13(v100) = v103;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, 72),
  [=](size_t v104) {
    float v105 = v35(v104);
    float v106;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 64),
    [=](size_t v107, float& lreduce1) {
      float v108 = v13(v107);
      float v109 = v5(v107, v104);
      float v110 = v108 * v109;
      float v111 = lreduce1 + v110;
      lreduce1 = v111;
      ;
    }, (v106));
    Kokkos::Sum<float> v106_joiner(v106);
    v106_joiner.join(v106, v105);
    ;
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v35(v104) = v106;
    });
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, 72),
  [=](size_t v112) {
    float v113 = v35(v112);
    float v114 = v11(v112);
    float v115 = v113 + v114;
    v2(v112) = v115;
  });
  team.team_barrier();
  return;
}



extern "C" void lapis_initialize();
extern "C" void lapis_finalize();
#endif
