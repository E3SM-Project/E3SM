#ifndef HOMMEXX_PACK_TRAITS_HPP
#define HOMMEXX_PACK_TRAITS_HPP

namespace Homme
{

template<typename PackType>
struct PackTraits {
private:
  // Note: I want any use of the general, non-specialized, PackTratis to generate
  //       a compiler error. However, simply static_assert(false,""); does not work,
  //       since false does not depend on the template parameter, so it can be
  //       immediately evaluated. On the other hand, any check that requires knowledge
  //       of PackType triggers SFINAE-like behavior, so it's safe.
  static constexpr bool dummy_false = (sizeof(PackType)<0);
  static_assert(dummy_false, "Error! Specialization of PackTraits missing for this PackType.\n");

public:
  // These are only to show what your specialization should expose.

  // The length of the pack
  static constexpr int pack_length = -1;

  // The underlying type
  using value_type = void;
};

} // namespace Homme

#endif // HOMMEXX_PACK_TRAITS_HPP
