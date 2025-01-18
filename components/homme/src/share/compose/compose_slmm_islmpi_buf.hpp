#ifndef INCLUDE_COMPOSE_SLMM_ISLMPI_BUF_HPP
#define INCLUDE_COMPOSE_SLMM_ISLMPI_BUF_HPP

namespace homme {
namespace islmpi {

const int nreal_per_2int = (2*sizeof(Int) + sizeof(Real) - 1) / sizeof(Real);

template <typename Buffer> SLMM_KIF
Int setbuf (Buffer& buf, const Int& os, const Int& i1, const Int& i2) {
  Int* const b = reinterpret_cast<Int*>(&buf(os));
  b[0] = i1;
  b[1] = i2;
  return nreal_per_2int;
}

template <typename Buffer> SLMM_KIF
Int setbuf (Buffer& buf, const Int& os, const Int& i1, const short& i2, const short& i3) {
  static_assert(sizeof(Int) >= 2*sizeof(short), "Need >= 2 shorts per Int");
  Int* const b = reinterpret_cast<Int*>(&buf(os));
  b[0] = i1;
  short* const b2 = reinterpret_cast<short*>(b+1);
  b2[0] = i2;
  b2[1] = i3;
  return nreal_per_2int;
}

template <typename Buffer> SLMM_KIF
Int setbuf (Buffer& buf, const Int& os, const Int& i1, const Int& i2,
            const bool final) {
  if (final) setbuf(buf, os, i1, i2);
  return nreal_per_2int;
}

template <typename Buffer> SLMM_KIF
Int setbuf (Buffer& buf, const Int& os, const Int& i1, const short& i2, const short& i3,
            const bool final) {
  if (final) setbuf(buf, os, i1, i2, i3);
  return nreal_per_2int;
}

template <typename Buffer> SLMM_KIF
Int getbuf (Buffer& buf, const Int& os, Int& i1, Int& i2) {
  const Int* const b = reinterpret_cast<const Int*>(&buf(os));
  i1 = b[0];
  i2 = b[1];
  return nreal_per_2int;
}

} // namespace islmpi
} // namespace homme

#endif
