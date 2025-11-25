#include "field_mask.hpp"

namespace scream {

namespace impl {

template<bool AND, typename LhsView, typename RhsView>
struct Combine
{
  using exec_space = typename Field::device_t::execution_space;
  template<int N>
  using MDRange = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<N>>;

  LhsView lhs;
  RhsView rhs;
  int     not_lhs;
  int     not_rhs;

  Combine(const LhsView& lhs_in, const RhsView& rhs_in, bool not_lhs_in, bool not_rhs_in)
  {
    lhs = lhs_in;
    rhs = rhs_in;
    not_lhs = static_cast<int>(not_lhs_in);
    not_rhs = static_cast<int>(not_rhs_in);
  }

  template<typename... Indices>
  KOKKOS_INLINE_FUNCTION
  void operator (Indices... idx) const {
    auto a = rhs(idx...) ^ not_lhs;
    auto b = lhs(idx...) ^ not_rhs;
    lhs(idx...) = AND ? a & b : a | b;
  }

  void run () {
    constexpr int N = LhsView::rank;
    if constexpr (N==0 or N==1) {
      Kokkos::RangePolicy<exec_space> p(0,lhs.size());
      Kokkos::parallel_for(p,*this);
    } else {
      int beg[N] = {};
      int end[N] = {};
      for (int i=0; i<N; ++i) end[i] = lhs.extent_int(i);

      MDRange p(beg,end);
      Kokkos::parallel_for(p,*this);
    }
  }
}

template<bool AND, typename LhsView, typename RhsView>
void combine(const LhsView& lhs, const RhsView& rhs, bool not_rhs)
{
  Combine<AND,LhsView,RhsView> functor(lhs,rhs,not_rhs);
  functor.run();
}

} // namespace detail

FieldMask::
FieldMask(const std::string& name, const FieldLayout& layout, const std::string& grid_name, bool allocate)
 : Field(FieldIdentifier(name,layout,ekat::units::Units::nondimensional(),grid_name,DataType::IntType))
{
  if (allocate)
    allocate_view();
}

FieldMask (const Field& f, const std::string& mname, bool allocate)
 : FieldMask(mname,
             f.get_header().get_identifier().get_layout(),
             f.get_header().get_identifier().get_grid_name(),
             allocate)
{
  // Nothing to do here
}

FieldMask (const Field& f)
 : Field(f)
{
  EKAT_REQUIRE_MSG (f.data_type()==DataType::IntType,
      "[FieldMask] Error! Cannot create a FieldMask from input Field. Bad data type.\n"
      " - field name: " + f.name() + "\n"
      " - field data type: " + e2str(f.data_type()) + "\n");
}

void FieldMask::
combine(const FieldMask& x, bool and_op, bool not_rhs)
{
  EKAT_REQUIRE_MSG (x.rank()==rank(),
      "[FieldMask::land] Error! Incompatible mask ranks.\n"
      " - lhs name: " + name() + "\n"
      " - rhs name: " + x.name() + "\n"
      " - lhs rank: " + std::to_string(rank()) + "\n"
      " - rhs rank: " + std::to_string(x.rank()) + "\n");

  EKAT_REQUIRE_MSG (not m_is_read_only,
      "[FieldMask::land] Error! Cannot update read-only mask.\n"
      " - mask name: " + name() + "\n");

  switch (rank()) {
    case 1:
    {
      auto xv = x.get_view<int*,Device>();
      auto yv = get_view<const int*,Device>();
      if (and_op)
        combine<true>(yv,xv,not_rhs);
      else
        combine<false>(yv,xv,not_rhs);
    } break;
    case 2:
    {
      auto xv = x.get_view<int**,Device>();
      auto yv = get_view<const int**,Device>();
      if (and_op)
        combine<true>(yv,xv,not_rhs);
      else
        combine<false>(yv,xv,not_rhs);
    } break;
    case 3:
    {
      auto xv = x.get_view<int***,Device>();
      auto yv = get_view<const int***,Device>();
      if (and_op)
        combine<true>(yv,xv,not_rhs);
      else
        combine<false>(yv,xv,not_rhs);
    } break;
    case 4:
    {
      auto xv = x.get_view<int****,Device>();
      auto yv = get_view<const int****,Device>();
      if (and_op)
        combine<true>(yv,xv,not_rhs);
      else
        combine<false>(yv,xv,not_rhs);
    } break;
    case 5:
    {
      auto xv = x.get_view<int*****,Device>();
      auto yv = get_view<const int*****,Device>();
      if (and_op)
        combine<true>(yv,xv,not_rhs);
      else
        combine<false>(yv,xv,not_rhs);
    } break;
    case 6:
    {
      auto xv = x.get_view<int******,Device>();
      auto yv = get_view<const int******,Device>();
      if (and_op)
        combine<true>(yv,xv,not_rhs);
      else
        combine<false>(yv,xv,not_rhs);
    } break;
    default:
      EKAT_ERROR_MSG (
          "[FieldMask::land] Error! Unsupported rank.\n"
          " - mask name: " + name() + "\n"
          " - mask rank: " + std::to_string(rank()) + "\n");
  }
  return *this;
}

FieldMask& FieldMask::operator&= (const FieldMask& x)
{
  return combine(x,true,false,false);
}


FieldMask& FieldMask::operator|= (const FieldMask& x)
{
  return combine(x,false,false,false);
}

FieldMask& FieldMask::flip ()
{
  return combine(*this,false,true); 
}

FieldMask FieldMask::operator!() const
{
  auto m = clone();
  return m.flip();
}

FieldMask FieldMask::clone(const std::string& name)
{
  return FieldMask(Field::clone(name));
}

template<int N>
ekat::ViewBroadcast<get_view_type<typename ekat::DataND<int,N>::type,Device>>
FieldMask::get_view (const FieldLayout& lt) const {
  using data_t = typename ekat::DataND<int,N>::type;
  using ret_t = ekat::ViewBroadcast<get_view_type<data_t,Device>>;

  const auto& my_lt = get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (lt.rank()>=my_lt.rank(),
      "[FieldMask::get_view] Error! Input layout has a rank lower than this mask.\n"
      " - mask name: " + name() + "\n"
      " - mask layout: " + my_lt.to_string() + "\n"
      " - input layout: " + lt.to_string() + "\n");

  const auto& dims = my_lt.dims();
  const auto& tags = my_lt.tags();
  auto bdims = lt.dims();
  auto btags = lt.tags();

  // We need to guard against ambiguous broadcast. E.g., bcasting a field with
  // layout ("col","cmp2","lev") to layout ("col","cmp2","cmp2","lev") is ambiguos,
  // as we don't know which of the two "cmp2" in the tgt layout is the original
  // vs bcasted dimension. We *could* use some info from the extent, but that's fragile,
  // as things may change going fwd in the src/tgt layout. Hence, we simply don't allow this scenario.
  // This implies that, if this field has N tags with name "foo" and N>1, then the target layout
  // must also have N tags with name "foo".
  // NOTE: we use dim NAMES rather than TAGS to give a tad more flexibility.
  auto lt_diff = lt;
  for (int i=0; i<my_lt.rank(); ++i) {
    const auto& n = my_lt.name(i);
    int my_count = ekat::count(my_lt.names(),n);
    int tgt_count = ekat::count(lt.names(),n);
    EKAT_REQUIRE_MSG(my_count>=tgt_count,
        "Error! Target layout is missing one of this field's layout dimensions.\n"
        " - field name: " + name() + "\n"
        " - field layout: " + my_lt.to_string() + "\n"
        " - target layout: " + lt.to_string() + "\n");
    EKAT_REQUIRE_MSG(my_count==1 or my_count==tgt_count,
        "Error! Target layout has more copies of one of this field's layout dimensions.\n"
        " - field name: " + name() + "\n"
        " - field layout: " + my_lt.to_string() + "\n"
        " - target layout: " + lt.to_string() + "\n");
    lt_diff.strip_dim(lt_diff.dim_idx(n));
  }

  int rank_diff = lt.rank()-my_lt.rank();


  for (int i=0,j=0; i<my_lt.rank() && j<lt.rank();) {
    if (btags[i]==tags[i] and ) {

    }
  }
  switch (m_src_field.rank()) {
    case 1:
    {
      auto v = m_src_field.get_view<int*,HD>();
      return broadcast(v,m_extents);
    }
    case 2:
    {
      auto v = m_src_field.get_view<int**,HD>();
      return broadcast(v,m_extents);
    }
  }
}

} // namespace scream
