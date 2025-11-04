#ifndef EAMXX_EVALUATE_EXPRESSION_HPP
#define EAMXX_EVALUATE_EXPRESSION_HPP

#include <share/expressions/base.hpp>

namespace scream {

namespace impl {

template<int N, typename Derived>
void evaluate (Expression<Derived>& e, Field& result)
{
  static_assert(N>=1 and N<=3, "Unsuppoerted rank.\n");
  using dev_t = Field::device_t;
  using kk_t  = KokkosTypes<dev_t>;
  using Policy1D = typename kk_t::RangePolicy;
  using PolicyMD = Kokkos::MDRangePolicy<kk_t::ExeSpace,Kokkos::Rank<N>>;
  using DataT  = typename ekat::DataND<Real,N>::type;

  const auto& fl = result.get_header().get_identifier().get_layout();

  int beg[N] = {0};
  int end[N];
  for (int i=0; i<N; ++i) end[i] = fl.dim(i);

  // Cast now, and capture the derived obj in the lambda, to make sure we get the
  // correct default copy constructor behavior
  auto impl = e.cast();
  impl.set_eval_layout(result.get_header().get_identifier().get_layout());
  auto res_v = result.get_view<DataT>();
  if constexpr (N==1) {
    Policy1D p(0,end[0]);
    auto eval = KOKKOS_LAMBDA (int i) {
      EvalData d {i};
      impl.set_eval_data(d);
      res_v(i) = impl.eval();
    };
    Kokkos::parallel_for(p,eval);
  } else if constexpr (N==2) {

    PolicyMD p(beg,end);
    auto eval = KOKKOS_LAMBDA (int i,int j) {
      EvalData d {i,j};
      impl.set_eval_data(d);
      res_v(i,j) = impl.eval();
    };
    Kokkos::parallel_for(p,eval);
  } else {
    PolicyMD p(beg,end);
    auto eval = KOKKOS_LAMBDA (int i,int j,int k) {
      EvalData d {i,j,k};
      impl.set_eval_data(d);
      res_v(i,j,k) = impl.eval();
    };
    Kokkos::parallel_for(p,eval);
  }
}

} // namespace impl

template<typename Derived>
void evaluate (Expression<Derived>& e, Field& result)
{
  EKAT_REQUIRE_MSG (result.rank()==e.num_indices(),
    "[evaluate] Error! Input expression and result field have different ranks.\n"
    " - field name: " + result.name() + "\n"
    " - field rank: " << result.rank() << "\n"
    " - expression rank: " << e.num_indices() + "\n");

  EKAT_REQUIRE_MSG (result.data_type()==DataType::RealType,
    "[evaluate_1d] Error! We currently only support expression templates for real-valued fields.\n"
    " - field name: " + result.name() + "\n"
    " - data type : " + e2str(result.data_type()) + "\n");

  if (e.num_indices()==1) {
    impl::evaluate<1>(e,result);
  } else if (e.num_indices()==2) {
    impl::evaluate<2>(e,result);
  } else {
    impl::evaluate<3>(e,result);
  }
}

} // namespace scream

#endif // EAMXX_EVALUATE_EXPRESSION_HPP
