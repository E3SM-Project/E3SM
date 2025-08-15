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
  return 2280;
}

// Find L1 address shift: for allocations spilling into level 1 scratch,
// this is subtracted from the allocation address to find its relative address within L1.
KOKKOS_INLINE_FUNCTION constexpr int forward_L1_shift(int L0_scratch_max) {
  int tmp = 2280;
  if(328 > L0_scratch_max && 72 < tmp)
    tmp = 72;
  if(584 > L0_scratch_max && 328 < tmp)
    tmp = 328;
  if(840 > L0_scratch_max && 584 < tmp)
    tmp = 584;
  if(1096 > L0_scratch_max && 840 < tmp)
    tmp = 840;
  if(1128 > L0_scratch_max && 840 < tmp)
    tmp = 840;
  if(1416 > L0_scratch_max && 1128 < tmp)
    tmp = 1128;
  if(1416 > L0_scratch_max && 1128 < tmp)
    tmp = 1128;
  if(72 > L0_scratch_max && 0 < tmp)
    tmp = 0;
  if(1704 > L0_scratch_max && 1416 < tmp)
    tmp = 1416;
  if(1704 > L0_scratch_max && 1128 < tmp)
    tmp = 1128;
  if(2280 > L0_scratch_max && 1704 < tmp)
    tmp = 1704;
  if(1704 > L0_scratch_max && 1128 < tmp)
    tmp = 1128;
  if(840 > L0_scratch_max && 584 < tmp)
    tmp = 584;
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
  if(1096 <= L0_scratch_max && 1096 > tmp)
    tmp = 1096;
  if(1128 <= L0_scratch_max && 1128 > tmp)
    tmp = 1128;
  if(1416 <= L0_scratch_max && 1416 > tmp)
    tmp = 1416;
  if(1416 <= L0_scratch_max && 1416 > tmp)
    tmp = 1416;
  if(72 <= L0_scratch_max && 72 > tmp)
    tmp = 72;
  if(1704 <= L0_scratch_max && 1704 > tmp)
    tmp = 1704;
  if(1704 <= L0_scratch_max && 1704 > tmp)
    tmp = 1704;
  if(2280 <= L0_scratch_max && 2280 > tmp)
    tmp = 2280;
  if(1704 <= L0_scratch_max && 1704 > tmp)
    tmp = 1704;
  if(840 <= L0_scratch_max && 840 > tmp)
    tmp = 840;
  return tmp;
}

// Find the actual level 1 scratch required by the function, assuming this limit for L0 scratch.
// This has no strict upper bound.
KOKKOS_INLINE_FUNCTION constexpr int forward_L1_scratch_required(int L0_scratch_max) {
  return 2280 - forward_L1_shift(L0_scratch_max);
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
  Kokkos::View<float[1][64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v13((float*) (v13_spill ? (scratch1 + 72 - l1_cutoff) : (scratch0 + 72)));
  ;
  constexpr bool v14_spill = 584 > L0_scratch_max;
  Kokkos::View<float[1][64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v14((float*) (v14_spill ? (scratch1 + 328 - l1_cutoff) : (scratch0 + 328)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 64),
  [=](size_t v15, size_t v16) {
    v14(v15, v16) = 0.0e+00f;
  });
  team.team_barrier();
  constexpr bool v17_spill = 840 > L0_scratch_max;
  Kokkos::View<float[1][64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v17((float*) (v17_spill ? (scratch1 + 584 - l1_cutoff) : (scratch0 + 584)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 64),
  [=](size_t v18, size_t v19) {
    float v20 = v14(v18, v19);
    v17(v18, v19) = v20;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadMDRange(team, 1, 64),
  [=](size_t v21, size_t v22) {
    float v23 = v17(v21, v22);
    float v24;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 72),
    [=](size_t v25, float& lreduce1) {
      float v26 = v3(v21, v25);
      float v27 = v12(v25, v22);
      float v28 = v26 * v27;
      float v29 = lreduce1 + v28;
      lreduce1 = v29;
      ;
    }, (v24));
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v17(v21, v22) = v24;
    });
  });
  team.team_barrier();
  constexpr bool v30_spill = 1096 > L0_scratch_max;
  Kokkos::View<float[1][64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v30((float*) (v30_spill ? (scratch1 + 840 - l1_cutoff) : (scratch0 + 840)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 64),
  [=](size_t v31, size_t v32) {
    float v33 = v17(0, v32);
    float v34 = v8(v32);
    float v35 = v33 + v34;
    v30(v31, v32) = v35;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 64),
  [=](size_t v36, size_t v37) {
    float v38 = v30(0, v37);
    bool v39 = (Kokkos::isnan(v38) || Kokkos::isnan(0.0e+00f)) || (v38 > 0.0e+00f);
    float v40 = v39? v38 : 0.0e+00f;
    v13(v36, v37) = v40;
  });
  team.team_barrier();
  constexpr bool v41_spill = 1128 > L0_scratch_max;
  Kokkos::View<float[1][72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v41((float*) (v41_spill ? (scratch1 + 840 - l1_cutoff) : (scratch0 + 840)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v42, size_t v43) {
    v41(v42, v43) = 0.0e+00f;
  });
  team.team_barrier();
  constexpr bool v44_spill = 1416 > L0_scratch_max;
  Kokkos::View<float[1][72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v44((float*) (v44_spill ? (scratch1 + 1128 - l1_cutoff) : (scratch0 + 1128)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v45, size_t v46) {
    float v47 = v41(v45, v46);
    v44(v45, v46) = v47;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadMDRange(team, 1, 72),
  [=](size_t v48, size_t v49) {
    float v50 = v44(v48, v49);
    float v51;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 64),
    [=](size_t v52, float& lreduce1) {
      float v53 = v13(v48, v52);
      float v54 = v7(v52, v49);
      float v55 = v53 * v54;
      float v56 = lreduce1 + v55;
      lreduce1 = v56;
      ;
    }, (v51));
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v44(v48, v49) = v51;
    });
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v57, size_t v58) {
    float v59 = v44(0, v58);
    float v60 = v9(v58);
    float v61 = v59 + v60;
    v2(v57, v58) = v61;
  });
  team.team_barrier();
  constexpr bool v62_spill = 1416 > L0_scratch_max;
  Kokkos::View<float[1][72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v62((float*) (v62_spill ? (scratch1 + 1128 - l1_cutoff) : (scratch0 + 1128)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v63, size_t v64) {
    float v65 = v2(0, v64);
    float v66 = -v65;
    float v67 = Kokkos::exp(v66);
    float v68 = v67 + 1.000000000e+00f;
    float v69 = 1.000000000e+00f / v68;
    v62(v63, v64) = v69;
  });
  team.team_barrier();
  constexpr bool v70_spill = 72 > L0_scratch_max;
  Kokkos::View<bool[1][72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v70((bool*) (v70_spill ? (scratch1 + 0 - l1_cutoff) : (scratch0 + 0)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v71, size_t v72) {
    float v73 = v62(0, v72);
    double v74 = (double) v73;
    bool v75 = (v74 > 5.00000000000000000e-01);
    v70(v71, v72) = v75;
  });
  team.team_barrier();
  constexpr bool v76_spill = 1704 > L0_scratch_max;
  Kokkos::View<float[1][72], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v76((float*) (v76_spill ? (scratch1 + 1416 - l1_cutoff) : (scratch0 + 1416)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v77, size_t v78) {
    bool v79 = v70(0, v78);
    float v80 = (float) v79;
    v76(v77, v78) = v80;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v81, size_t v82) {
    float v83 = v76(0, v82);
    float v84 = v62(0, v82);
    float v85 = v83 - v84;
    v2(v81, v82) = v85;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v86, size_t v87) {
    float v88 = v2(0, v87);
    float v89 = v62(0, v87);
    float v90 = v88 + v89;
    v1(v86, v87) = v90;
  });
  team.team_barrier();
  constexpr bool v91_spill = 1704 > L0_scratch_max;
  Kokkos::View<float[1][144], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v91((float*) (v91_spill ? (scratch1 + 1128 - l1_cutoff) : (scratch0 + 1128)));
  ;
  constexpr bool v92_spill = 2280 > L0_scratch_max;
  Kokkos::View<float[1][144], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v92((float*) (v92_spill ? (scratch1 + 1704 - l1_cutoff) : (scratch0 + 1704)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 144),
  [=](size_t v93, size_t v94) {
    float v95 = v91(v93, v94);
    v92(v93, v94) = v95;
  });
  team.team_barrier();
  Kokkos::LayoutStride v96_layout(1, 144, 72, 1);
  Kokkos::View<float[1][72], Kokkos::LayoutStride, Kokkos::AnonymousSpace> v96(v92.data() + 0, v96_layout);
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v97, size_t v98) {
    float v99 = v4(v97, v98);
    v96(v97, v98) = v99;
  });
  team.team_barrier();
  constexpr bool v100_spill = 1704 > L0_scratch_max;
  Kokkos::View<float[1][144], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v100((float*) (v100_spill ? (scratch1 + 1128 - l1_cutoff) : (scratch0 + 1128)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 144),
  [=](size_t v101, size_t v102) {
    float v103 = v92(v101, v102);
    v100(v101, v102) = v103;
  });
  team.team_barrier();
  Kokkos::LayoutStride v104_layout(1, 144, 72, 1);
  Kokkos::View<float[1][72], Kokkos::LayoutStride, Kokkos::AnonymousSpace> v104(v100.data() + 72, v104_layout);
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v105, size_t v106) {
    float v107 = v1(v105, v106);
    v104(v105, v106) = v107;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadMDRange(team, 1, 64),
  [=](size_t v108, size_t v109) {
    float v110 = v14(v108, v109);
    float v111;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 144),
    [=](size_t v112, float& lreduce1) {
      float v113 = v100(v108, v112);
      float v114 = v6(v112, v109);
      float v115 = v113 * v114;
      float v116 = lreduce1 + v115;
      lreduce1 = v116;
      ;
    }, (v111));
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v14(v108, v109) = v111;
    });
  });
  team.team_barrier();
  constexpr bool v117_spill = 840 > L0_scratch_max;
  Kokkos::View<float[1][64], Kokkos::LayoutRight, Kokkos::AnonymousSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v117((float*) (v117_spill ? (scratch1 + 584 - l1_cutoff) : (scratch0 + 584)));
  ;
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 64),
  [=](size_t v118, size_t v119) {
    float v120 = v14(0, v119);
    float v121 = v10(v119);
    float v122 = v120 + v121;
    v117(v118, v119) = v122;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 64),
  [=](size_t v123, size_t v124) {
    float v125 = v117(0, v124);
    bool v126 = (Kokkos::isnan(v125) || Kokkos::isnan(0.0e+00f)) || (v125 > 0.0e+00f);
    float v127 = v126? v125 : 0.0e+00f;
    v13(v123, v124) = v127;
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadMDRange(team, 1, 72),
  [=](size_t v128, size_t v129) {
    float v130 = v41(v128, v129);
    float v131;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, 0, 64),
    [=](size_t v132, float& lreduce1) {
      float v133 = v13(v128, v132);
      float v134 = v5(v132, v129);
      float v135 = v133 * v134;
      float v136 = lreduce1 + v135;
      lreduce1 = v136;
      ;
    }, (v131));
    Kokkos::single(Kokkos::PerThread(team), [&]() {
      v41(v128, v129) = v131;
    });
  });
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorMDRange(team, 1, 72),
  [=](size_t v137, size_t v138) {
    float v139 = v41(0, v138);
    float v140 = v11(v138);
    float v141 = v139 + v140;
    v2(v137, v138) = v141;
  });
  team.team_barrier();
  return;
}



extern "C" void lapis_initialize();
extern "C" void lapis_finalize();
#endif
