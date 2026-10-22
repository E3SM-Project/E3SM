#include "share/remap/field_views_repo.hpp"

#include <ekat_assert.hpp>

namespace scream
{

void FieldViewsRepo::
init (const std::vector<Field>& fields)
{
  const int nf = fields.size();

  views  = views_t  ("views", nf);
  cviews = cviews_t ("cviews",nf);

  h_views  = Kokkos::create_mirror_view(views);
  h_cviews = Kokkos::create_mirror_view(cviews);

  m_subfield_info.clear();
  for (int i=0; i<nf; ++i) {
    const auto& h = fields[i].get_header();

    // A dynamic subfield will need some special treatment at runtime
    // Namely, we'll need to re-extract the view every time,
    // since the subview info may have changed
    if (h.get_parent()) {
      EKAT_REQUIRE_MSG (h.get_parent()->get_parent()==nullptr,
          "Error! We do not support remapping of subfields of other subfields.\n");
      m_subfield_info[i] = h.get_alloc_properties().get_subview_info();
    }

    set_views(i,fields[i]);
  }

  Kokkos::deep_copy(views,  h_views);
  Kokkos::deep_copy(cviews, h_cviews);
}

bool FieldViewsRepo::
subfields_info_have_changed (const std::vector<Field>& fields) const
{
  for (const auto& it : m_subfield_info) {
    const auto& f = fields[it.first];
    const auto& info_old = it.second;
    const auto& info_new = f.get_header().get_alloc_properties().get_subview_info();
    if ( not(info_old==info_new) ) {
      return true;
    }
  }
  return false;
}

void FieldViewsRepo::
update_subfields (const std::vector<Field>& fields)
{
  bool changed = false;
  for (auto& [i, svi] : m_subfield_info) {
    const auto& f = fields[i];
    const auto& new_svi = f.get_header().get_alloc_properties().get_subview_info();
    if ( not(svi==new_svi) ) {
      set_views(i,f);
      svi = new_svi;
      changed = true;
    }
  }
  if (changed) {
    Kokkos::deep_copy(views,  h_views);
    Kokkos::deep_copy(cviews, h_cviews);
  }
}

void FieldViewsRepo::
set_views (const int i, const Field& f)
{
  const auto rank = f.get_header().get_identifier().get_layout().rank();
  switch (rank) {
    case 1: h_cviews[i].v1d = f.get_view<const Real*>();      break;
    case 2: h_cviews[i].v2d = f.get_view<const Real**>();     break;
    case 3: h_cviews[i].v3d = f.get_view<const Real***>();    break;
    case 4: h_cviews[i].v4d = f.get_view<const Real****>();   break;
    case 5: h_cviews[i].v5d = f.get_view<const Real*****>();  break;
  }
  if (not f.is_read_only()) {
    switch (rank) {
      case 1: h_views[i].v1d = f.get_view<Real*>();      break;
      case 2: h_views[i].v2d = f.get_view<Real**>();     break;
      case 3: h_views[i].v3d = f.get_view<Real***>();    break;
      case 4: h_views[i].v4d = f.get_view<Real****>();   break;
      case 5: h_views[i].v5d = f.get_view<Real*****>();  break;
    }
  }
}

} // namespace scream
