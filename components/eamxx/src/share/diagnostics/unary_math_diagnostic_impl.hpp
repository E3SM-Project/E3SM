#include "unary_mth_diagnostic.hpp"

namespace scream {

namespace {
  template<int N>
  using MDRange = Kokkos::MDRangePolicy<typename KokkosTypes<DefaultDevice>::ExeSpace,Kokkos::Rank<N>>;
} // anonymous namespace

template <typename Op>
UnaryMathDiagnostic::
UnaryMathDiagnostic(const ekat::Comm& comm,
                    const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
  : AbstractDiagnostic(comm, params, grid)
{
  // Retrieve input field name from params (e.g. "field_name")
  m_field_in_names.push_back(m_params.get<std::string>("field_name"));
}

template <typename Op>
std::string UnaryMathDiagnostic::
name() const override {
  return std::string(Op::name()) + "(" + m_fields_in_names.front() + ")";
}

template <typename Op>
void UnaryMathDiagnostic::
initialize_impl() override {
  const auto& f_in = m_fields_in.at(m_fields_in_names.front());

  // Setup diagnostic output metadata (same layout/grid as input)
  auto fid = f_in.clone(name());
  m_diagnostic_output = Field(fid,true);
}

template <typename Op>
void UnaryMathDiagnostic::
compute_impl() override {
  const auto& f_in = m_fields_in.at(m_input_name);
  auto do_compute = [&](auto fv, auto dv) {
    if constexpr (fv.rank()==0) {
    } else if constexpr (fv.rank()==
  };

  auto dims = f_in.get_header().get_identifier().get_layout().dims();
  if (f_in.rank()==0) {
    Kokkos::RangePolicy<exec_space> policy(0,1);
    auto lhs = m_diagnostic_output.get_view<Real>();
    auto rhs = m_diagnostic_output.get_view<const Real>();
    auto lambda = KOKKOS_LAMBDA (int) { lhs() = op(rhs()); };
    Kokkos::parallel_for(policy,lambda);
  } else if (f_in.rank()==1) {
    Kokkos::RangePolicy<exec_space> policy(0,dims[0]);
    auto lhs = m_diagnostic_output.get_view<Real*>();
    auto rhs = m_diagnostic_output.get_view<const Real*>();
    auto lambda = KOKKOS_LAMBDA (int i) { lhs(i) = op(rhs(i)); };
    Kokkos::parallel_for(policy,lambda);
  } else if (f_in.rank()==2) {
    MDRange<2> policy({0,0},{dims[0],dims[1]});
    auto lhs = m_diagnostic_output.get_view<Real**>();
    auto rhs = m_diagnostic_output.get_view<const Real**>();
    auto lambda = KOKKOS_LAMBDA (int i, int j) { lhs(i,j) = op(rhs(i,j)); };
    Kokkos::parallel_for(policy,lambda);
  } else if (f_in.rank()==3) {
    MDRange<3> policy({0,0},{dims[0],dims[1],dims[2]});
    auto lhs = m_diagnostic_output.get_view<Real***>();
    auto rhs = m_diagnostic_output.get_view<const Real***>();
    auto lambda = KOKKOS_LAMBDA (int i, int j, int k) { lhs(i,j,k) = op(rhs(i,j,k)); };
    Kokkos::parallel_for(policy,lambda);
  } else {
    EKAT_ERROR_MSG ("Error! Unsupported rank in UnaryMathDiagnostic.\n");
  }
}

} // namespace scream
