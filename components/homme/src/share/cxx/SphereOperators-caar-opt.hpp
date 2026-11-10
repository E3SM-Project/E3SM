#ifndef HOMMEXX_SPHERE_OPERATORS_CAAR_OPT_HPP
#define HOMMEXX_SPHERE_OPERATORS_CAAR_OPT_HPP

#include "SphereOperators.hpp"

namespace Homme {

static constexpr int NPNP = NP * NP;

#ifdef KOKKOS_ENABLE_CUDA
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace,Kokkos::LaunchBounds<1024,1> >;
#else
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
#endif

#ifndef WARP_SIZE

#if defined(KOKKOS_ENABLE_HIP)
#define WARP_SIZE 64
#elif defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_SYCL)
#define WARP_SIZE 32
#else
#define WARP_SIZE 1
#endif

#endif

#if (WARP_SIZE == 1)

#define SPHERE_BLOCK_START3(B,X,Y,Z) \
  Scalar X; Scalar sbo##X[NP][NP][NUM_LEV]; \
  Scalar Y; Scalar sbo##Y[NP][NP][NUM_LEV]; \
  Scalar Z; Scalar sbo##Z[NP][NP][NUM_LEV]; \
  for (int ix = 0; ix < NP; ix++) for(int iy = 0; iy < NP; iy++) { \
    B.update(ix, iy); \
    for (int iz = 0; iz < NUM_LEV; iz++) { \
      B.z = iz;

#define SPHERE_BLOCK_START0(B) \
  for (int ix = 0; ix < NP; ix++) for(int iy = 0; iy < NP; iy++) { \
    B.update(ix,iy); \
    for (int iz = 0; iz < NUM_LEV; iz++) { \
      B.z = iz;

#define SPHERE_BLOCK_MIDDLE3(B,X,Y,Z) \
      sbo##X[B.x][B.y][B.z] = X; \
      sbo##Y[B.x][B.y][B.z] = Y; \
      sbo##Z[B.x][B.y][B.z] = Z; \
    }\
  } \
  for (int ix = 0; ix < NP; ix++) for(int iy = 0; iy < NP; iy++) { \
    B.update(ix,iy); \
    for (int iz = 0; iz < NUM_LEV; iz++) { \
      B.z = iz; \
      X = sbo##X[B.x][B.y][B.z]; \
      Y = sbo##Y[B.x][B.y][B.z]; \
      Z = sbo##Z[B.x][B.y][B.z];

#define SPHERE_BLOCK_MIDDLE0(B) \
    } \
  } \
  for (int ix = 0; ix < NP; ix++) for(int iy = 0; iy < NP; iy++) { \
    B.update(ix,iy); \
    for (int iz = 0; iz < NUM_LEV; iz++) { \
      B.z = iz;

#define SPHERE_BLOCK_END() } }

#else

static_assert(VECTOR_SIZE == 1, "VECTOR_SIZE != 1");

// SPHERE_BLOCK macros for GPU (WARP_SIZE > 1).
//
// When NUM_PHYSICAL_LEV is not a multiple of WARP_SIZE (e.g. nlev=30, warp=32),
// threads with B.z >= NUM_LEV are padding threads.  The original formulation
// used an early `return` before the team_barrier(), which violates the CUDA
// requirement that all threads in a CTA must reach __syncthreads().  This
// caused non-deterministic results between otherwise-identical runs.
//
// Fix: guard the body with `if (B.z < NUM_LEV)` blocks so that all 512
// threads always reach the barrier.  Local variables X/Y/Z are declared in
// outer scope (default-initialized) so they remain accessible after the
// barrier for the active (z < NUM_LEV) threads.

#define SPHERE_BLOCK_START3(B,X,Y,Z) \
      Scalar X = {}, Y = {}, Z = {}; \
      if (B.z < NUM_LEV) {
#define SPHERE_BLOCK_START0(B) if (B.z < NUM_LEV) {

#define SPHERE_BLOCK_MIDDLE3(B,X,Y,Z) } B.t.team_barrier(); if (B.z < NUM_LEV) {
#define SPHERE_BLOCK_MIDDLE0(B)       } B.t.team_barrier(); if (B.z < NUM_LEV) {

#define SPHERE_BLOCK_END() }

static constexpr int SPHERE_BLOCK_LEV = WARP_SIZE;
static constexpr int SPHERE_BLOCK = NPNP * SPHERE_BLOCK_LEV;
static constexpr int SPHERE_BLOCKS_PER_COL = (NUM_PHYSICAL_LEV - 1) / SPHERE_BLOCK_LEV + 1;

#endif

namespace SphereOuter {
#if WARP_SIZE == 1
  using Team = TeamPolicy::member_type;
  template <typename F>
  void parallel_for(const int num_elems, F f) {
    Kokkos::parallel_for(
      __PRETTY_FUNCTION__,
      TeamPolicy(num_elems, 1),
      f);
  }
#else
  using Team = int;
  template <typename F>
  void parallel_for(const int num_elems, F f) {
    f(num_elems);
  }
#endif
}

struct SphereGlobal {

  const ExecViewManaged<const Real *[2][2][NP][NP]> d;
  const ExecViewManaged<const Real *[2][2][NP][NP]> dinv;
  const ExecViewManaged<const Real [NP][NP]> dvv;
  const ExecViewManaged<const Real *[NP][NP]> metdet;
  const Real scale_factor_inv;

  SphereGlobal(const SphereOperators &op):
    d(op.m_d),
    dinv(op.m_dinv),
    dvv(op.dvv),
    metdet(op.m_metdet),
    scale_factor_inv(op.m_scale_factor_inv)
  {}
};

struct SphereBlockOps;

#if (WARP_SIZE == 1)

struct SphereBlockScratch {
  const SphereBlockOps &b;
  Scalar v[NP][NP][NUM_LEV];
  KOKKOS_INLINE_FUNCTION SphereBlockScratch(const SphereBlockOps &);
  KOKKOS_INLINE_FUNCTION const Scalar &sv(int x, int y) const;
  KOKKOS_INLINE_FUNCTION Scalar &sv(int x, int y);
};

#else

using SphereBlockScratchView = Kokkos::View<
  Scalar[SPHERE_BLOCK_LEV][NP][NP],
  ExecSpace::scratch_memory_space,
  Kokkos::MemoryTraits<Kokkos::Unmanaged>
  >;

using SphereBlockScratchSubview = Kokkos::Subview<
  SphereBlockScratchView, int,
  std::remove_const_t<decltype(Kokkos::ALL)>,
  std::remove_const_t<decltype(Kokkos::ALL)>
  >;

struct SphereBlockScratch {
  SphereBlockScratchView v;
  SphereBlockScratchSubview sv;
  KOKKOS_INLINE_FUNCTION SphereBlockScratch(const SphereBlockOps &b);
};

#endif

struct SphereBlockOps {
  using Team = TeamPolicy::member_type;

  const SphereGlobal &g;
  const Team &t;
  Real scale_factor_inv;
  Real metdet;
  Real rrdmd;
  Real dinv[2][2];
  Real dvvx[NP];
  Real dvvy[NP];
  int e,x,y,z;

  KOKKOS_INLINE_FUNCTION SphereBlockOps(const SphereGlobal &sg, const Team &team):
    g(sg),
    t(team),
    scale_factor_inv(g.scale_factor_inv)
  {
#if (WARP_SIZE == 1)
    e = t.league_rank();
    x = y = z = 0;
#else
    const int lr = t.league_rank();
    e = lr / SPHERE_BLOCKS_PER_COL;
    const int iw = lr % SPHERE_BLOCKS_PER_COL;
    const int tr = t.team_rank();
    const int ixy = tr / SPHERE_BLOCK_LEV;
    const int dz = tr % SPHERE_BLOCK_LEV;
    z = dz + iw * SPHERE_BLOCK_LEV;
    update(ixy / NP, ixy % NP);
#endif
  }

  KOKKOS_INLINE_FUNCTION void update(const int ix, const int iy)
  {
    rrdmd = 0;
    x = ix;
    y = iy;
    metdet = g.metdet(e,x,y);
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) dinv[i][j] = g.dinv(e,i,j,x,y);
    for (int j = 0; j < NP; j++) {
      dvvx[j] = g.dvv(x,j);
      dvvy[j] = g.dvv(y,j);
    }
  }

  KOKKOS_INLINE_FUNCTION Scalar div(const SphereBlockScratch &t0, const SphereBlockScratch &t1)
  {
    // Separately compute du and dv to match Fortran bfb
    Scalar du = 0;
    Scalar dv = 0;
    for (int j = 0; j < NP; j++) {
      du += dvvy[j] * t0.sv(x,j);
      dv += dvvx[j] * t1.sv(j,y);
    }
    if (rrdmd == 0) rrdmd = (1.0 / metdet) * scale_factor_inv;
    return (du + dv) * rrdmd;
  }

  KOKKOS_INLINE_FUNCTION void divInit(SphereBlockScratch &t0, SphereBlockScratch &t1, const Scalar v0, const Scalar v1) const
  {
    t0.sv(x,y) = (dinv[0][0] * v0 + dinv[1][0] * v1) * metdet;
    t1.sv(x,y) = (dinv[0][1] * v0 + dinv[1][1] * v1) * metdet;
  }

  KOKKOS_INLINE_FUNCTION void grad(Scalar &g0, Scalar &g1, const SphereBlockScratch &t) const
  {
    Scalar s0 = 0;
    Scalar s1 = 0;
    for (int j = 0; j < NP; j++) {
      s0 += dvvy[j] * t.sv(x,j);
      s1 += dvvx[j] * t.sv(j,y);
    }
    s0 *= scale_factor_inv;
    s1 *= scale_factor_inv;
    g0 = dinv[0][0] * s0 + dinv[0][1] * s1;
    g1 = dinv[1][0] * s0 + dinv[1][1] * s1;
  }

  KOKKOS_INLINE_FUNCTION void gradInit(SphereBlockScratch &t0, const Scalar v0)
  {
    t0.sv(x,y) = v0;
  }

  KOKKOS_INLINE_FUNCTION Scalar vort(const SphereBlockScratch &t0, const SphereBlockScratch &t1)
  {
    Scalar du = 0;
    Scalar dv = 0;
    for (int j = 0; j < NP; j++) {
      du += dvvx[j] * t0.sv(j,y);
      dv += dvvy[j] * t1.sv(x,j);
    }
    if (rrdmd == 0) rrdmd = (1.0 / metdet) * scale_factor_inv;
    return (dv - du) * rrdmd;
  }

  KOKKOS_INLINE_FUNCTION void vortInit(SphereBlockScratch &t0, SphereBlockScratch &t1, const Scalar v0, const Scalar v1) const
  {
    t0.sv(x,y) = g.d(e,0,0,x,y) * v0 + g.d(e,0,1,x,y) * v1;
    t1.sv(x,y) = g.d(e,1,0,x,y) * v0 + g.d(e,1,1,x,y) * v1;
  }

  KOKKOS_INLINE_FUNCTION Scalar zabove(const Scalar *const here, const Scalar *const above) const
  {
    if (VECTOR_SIZE == 1) return above[0];
    Scalar out;
    for (int i = 0; i < VECTOR_END; i++) out[i] = here[0][i+1];
    out[VECTOR_END] = (z < NUM_LEV_P-1) ? above[0][0] : 0;
    return out;
  }

  KOKKOS_INLINE_FUNCTION Scalar zabove(const Scalar top, const Scalar *const here, const Scalar *const above) const
  {
    if (VECTOR_SIZE == 1) return (z < NUM_LEV-1) ? above[0] : top;
    Scalar out;
    for (int i = 0; i < VECTOR_END; i++) out[i] = here[0][i+1];
    if (z < NUM_LEV-1) {
      out[VECTOR_END] = above[0][0];
    } else {
      for (int i = (NUM_PHYSICAL_LEV-1)%VECTOR_SIZE; i < VECTOR_SIZE; i++) out[i] = top[i];
    }
    return out;
  }

  KOKKOS_INLINE_FUNCTION Scalar zbelow(const Scalar bottom, const Scalar *const below, const Scalar *const here) const
  {
    if (VECTOR_SIZE == 1) return (z > 0) ? below[0] : bottom;
    Scalar out;
    for (int i = 1; i < VECTOR_SIZE; i++) out[i] = here[0][i-1];
    if (z > 0) out[0] = below[0][VECTOR_END];
    else out[0] = bottom[0];
    return out;
  }

#if (WARP_SIZE == 1)
  template <typename F>
  static void parallel_for(const Team &team, int /*num_scratch*/, F f)
  {
    f(team);
  }
#else
  template <typename F>
  static void parallel_for(const int num_elems, const int num_scratch, F f)
  {
    Kokkos::parallel_for(
      __PRETTY_FUNCTION__,
      TeamPolicy(num_elems * SPHERE_BLOCKS_PER_COL, SPHERE_BLOCK).
      set_scratch_size(0, Kokkos::PerTeam(num_scratch * SphereBlockScratchView::shmem_size())),
      f);
  }
#endif

};

#if (WARP_SIZE == 1)

KOKKOS_INLINE_FUNCTION SphereBlockScratch::SphereBlockScratch(const SphereBlockOps &b):
  b(b)
{}

KOKKOS_INLINE_FUNCTION const Scalar &SphereBlockScratch::sv(const int x, const int y) const
{ return v[x][y][b.z]; }

KOKKOS_INLINE_FUNCTION Scalar &SphereBlockScratch::sv(const int x, const int y)
{ return v[x][y][b.z]; }

#else

KOKKOS_INLINE_FUNCTION SphereBlockScratch::SphereBlockScratch(const SphereBlockOps &b):
  v(b.t.team_scratch(0)),
  sv(Kokkos::subview(v, b.z % SPHERE_BLOCK_LEV, Kokkos::ALL, Kokkos::ALL))
{}

#endif

struct SphereCol {
  int e,x,y,z;

#if (WARP_SIZE == 1)
  struct Team { int e,x,y,z; };
  KOKKOS_INLINE_FUNCTION SphereCol(const Team &team):
    e(team.e),
    x(team.x),
    y(team.y),
    z(team.z)
  {}

  template <typename F>
  static void parallel_for(const SphereOuter::Team &outer, const int num_lev, F f)
  {
    const int e = outer.league_rank();
    for (int x = 0; x < NP; x++) for (int y = 0; y < NP; y++) {
      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(outer, num_lev),
        [=] (const int z) { f(Team{e,x,y,z}); });
    }
  }
#else
  using Team = TeamPolicy::member_type;

  KOKKOS_INLINE_FUNCTION SphereCol(const Team &team)
  {
    const int lr = team.league_rank();
    e = lr / NPNP;
    const int xy = lr % NPNP;
    x = xy / NP;
    y = xy % NP;
    z = team.team_rank();
  }

  template <typename F>
  static void parallel_for(const int num_elems, const int num_lev, F f)
  {
    Kokkos::parallel_for(
      __PRETTY_FUNCTION__,
      TeamPolicy(num_elems * NPNP, num_lev),
      f);
  }
#endif

  KOKKOS_INLINE_FUNCTION Scalar zplus(const Scalar *const here, const Scalar *const above) const
  {
    if (VECTOR_SIZE == 1) return above[0];
    Scalar out;
    for (int i = 0; i < VECTOR_END; i++) out[i] = here[0][i+1];
    out[VECTOR_END] = (z < NUM_LEV_P-1) ? above[0][0] : 0;
    return out;
  }
};

struct SphereColOps: SphereCol {
  Real scale_factor_inv;
  Real dinv[2][2];
  Real dvvx[NP];
  Real dvvy[NP];

  KOKKOS_INLINE_FUNCTION SphereColOps(const SphereGlobal &sg, const Team &team):
    SphereCol(team)
  {
    scale_factor_inv = sg.scale_factor_inv;
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) dinv[i][j] = sg.dinv(e,i,j,x,y);
    for (int j = 0; j < NP; j++) {
      dvvx[j] = sg.dvv(x,j);
      dvvy[j] = sg.dvv(y,j);
    }
  }

  template <typename OutView, typename InView>
  KOKKOS_INLINE_FUNCTION void grad(OutView &out, const InView &in, const int n) const
  {
    Scalar s0 = 0;
    Scalar s1 = 0;
    for (int j = 0; j < NP; j++) {
      s0 += dvvy[j] * in(e,n,x,j,z);
      s1 += dvvx[j] * in(e,n,j,y,z);
    }
    s0 *= scale_factor_inv;
    s1 *= scale_factor_inv;
    out(e,0,x,y,z) = dinv[0][0] * s0 + dinv[0][1] * s1;
    out(e,1,x,y,z) = dinv[1][0] * s0 + dinv[1][1] * s1;
  }

  KOKKOS_INLINE_FUNCTION Scalar zabove(const Scalar top, const Scalar *const in) const
  {
    // When z == NUM_LEV (only possible when NUM_LEV_P > NUM_LEV), `in` is nullptr
    // because z is one pack past a NUM_LEV-sized midpoint view.
    // When NUM_LEV == NUM_LEV_P and z == NUM_LEV_P-1, `in` points to the last
    // midpoint pack (valid for slots 0..NUM_PHYSICAL_LEV%VECTOR_SIZE-1), while
    // slot NUM_PHYSICAL_LEV%VECTOR_SIZE is the surface interface and takes `top`.
    if (z == NUM_LEV) return top;
    Scalar out = *in;
    if (z == NUM_LEV_P-1) out[NUM_PHYSICAL_LEV%VECTOR_SIZE] = top[NUM_PHYSICAL_LEV%VECTOR_SIZE];
    return out;
  }

  KOKKOS_INLINE_FUNCTION Scalar zbelow(const Real bottom, const Scalar *const below, const Scalar *const here) const
  {
    Scalar out;
    out[0] = (z == 0) ? bottom : below[0][VECTOR_END];
    for (int i = 1; i < VECTOR_SIZE; i++) out[i] = here ? here[0][i-1] : 0;
    return out;
  }

  KOKKOS_INLINE_FUNCTION Scalar ztop(const Real top, const Real in) const
  {
    Scalar out = in;
    if (z == NUM_LEV_P-1) out[NUM_PHYSICAL_LEV%VECTOR_SIZE] = top;
    return out;
  }

  KOKKOS_INLINE_FUNCTION Scalar ztrim(const Scalar bottom, const Scalar top, const Scalar in) const
  {
    Scalar out = in;
    if (z == 0) out[0] = bottom[0];
    if (z == NUM_LEV_P-1) out[NUM_PHYSICAL_LEV%VECTOR_SIZE] = top[NUM_PHYSICAL_LEV%VECTOR_SIZE];
    return out;
  }
};

struct SphereScanOps {
  TeamPolicy::member_type t;
  int e,x,y;

#if (WARP_SIZE == 1)
  struct Team {
    TeamPolicy::member_type t;
    int e,x,y;
  };

  KOKKOS_INLINE_FUNCTION SphereScanOps(const Team &team):
    t(team.t),
    e(team.e),
    x(team.x),
    y(team.y)
  {}

  template <typename F>
  static void parallel_for(const SphereOuter::Team &outer, F f)
  {
    const int e = outer.league_rank();
    for (int x = 0; x < NP; x++) for (int y = 0; y < NP; y++) {
      f(Team{outer,e,x,y});
    }
  }
#else
  using Team = TeamPolicy::member_type;

  KOKKOS_INLINE_FUNCTION SphereScanOps(const Team &team):
    t(team)
  {
    const int lr = t.league_rank();
    const int tr = t.team_rank();
    e = lr;
    x = tr / NP;
    y = tr % NP;
  }

  template <typename F>
  static void parallel_for(const int num_elems, F f)
  {
    Kokkos::parallel_for(
      __PRETTY_FUNCTION__,
      TeamPolicy(num_elems, NPNP, WARP_SIZE),
      f);
  }
#endif

  template <typename OutView, typename InView, typename EndView>
  KOKKOS_INLINE_FUNCTION void bwd_scan(OutView &out, const int n, const InView &in, const EndView &end) const
  {
    Dispatch<>::parallel_scan(
      t, NUM_PHYSICAL_LEV,
      [&](const int a, Real &sum, const bool last) {
        if (a == 0) out(e,n,x,y,NUM_PHYSICAL_LEV) = sum = end(e,x,y);
        const int z = NUM_PHYSICAL_LEV-a-1;
        sum += in(e,x,y,z);
        if (last) out(e,n,x,y,z) = sum;
      });
  }

  template <typename OutView, typename InView>
  KOKKOS_INLINE_FUNCTION void scan(OutView &out, const InView &in, const Real zero) const
  {
    Dispatch<>::parallel_scan(
      t, NUM_PHYSICAL_LEV,
      [&](const int z, Real &sum, const bool last) {
        if (z == 0) out(e,x,y,0) = sum = zero;
        sum += in(e,x,y,z);
        if (last) out(e,x,y,z+1) = sum;
      });
  }

  template <typename OutView, typename InView>
  KOKKOS_INLINE_FUNCTION void scan(OutView &out, const InView &in, const int n, const Real zero) const
  {
    Dispatch<>::parallel_scan(
      t, NUM_PHYSICAL_LEV,
      [&](const int z, Real &sum, const bool last) {
        if (z == 0) out(e,x,y,0) = sum = zero;
        sum += in(e,n,x,y,z);
        if (last) out(e,x,y,z+1) = sum;
      });
  }
};

} // namespace Homme

#endif // HOMMEXX_SPHERE_OPERATORS_CAAR_OPT_HPP
