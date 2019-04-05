#include "share/field/remap/remap_utils.hpp"
#include "share/mpi/scream_comm.hpp"

#include <algorithm>

namespace scream
{

// Shortcuts to avoid using long iterators syntax (like v.begin() and v.end())
template<typename T>
int count (const std::vector<T>& vector, const T& value) {
  return std::count(vector.begin(), vector.end(), value);
}

template<typename T>
std::vector<T> sort (const std::vector<T>& vector) {
  std::vector<T> copy(vector);
  std::sort(copy.begin(), copy.end());
  return copy;
}

template<typename T>
typename std::vector<T>::const_iterator
find (const std::vector<T>& vector, const T& value) {
  return std::find(vector.begin(), vector.end(), value);
}

template<typename T>
bool erase (std::vector<T>& vector, const T& value) {
  auto it = std::find(vector.begin(),vector.end(),value);
  if (it!=vector.end()) {
    vector.erase(it);
    return true;
  }
  return false;
}

// =================== Remap utilities ================== //

std::pair<PhysDyn,LayoutType> get_layout_specs (const FieldLayout& layout) {
  auto tags = layout.tags();

  const int n_element = count(tags,FieldTag::Element);
  const int n_column  = count(tags,FieldTag::Column);
  const int ngp = count(tags,FieldTag::GaussPoint);

  // Start from undefined/invalid
  auto result = std::make_pair(PhysDyn::Undefined,LayoutType::Invalid);

  if ( (n_element==0 && n_column==0) ||
        n_element>1 || n_column>1) {
    // This is not a valid layout.
    return result;
  }

  if (n_element>0) {
    // Remove the element tag
    erase(tags,FieldTag::Element);

    if (ngp!=2) {
      // Currently, we only support 2d and 3d fields, i.e., fields must have two
      // 'GaussPoint' tags. If they don't, we return invalid layout type.
      return result;
    }

    // A Dyn layout (potentially)
    result.first = PhysDyn::Dyn;

    // Remove the two 'GaussPoint' tags
    erase(tags,FieldTag::GaussPoint);
    erase(tags,FieldTag::GaussPoint);
  } else {
    if (ngp!=0) {
      // We don't allow to mix Column and GaussPoint tags.
      return result;
    }

    // A Phys layout (potentially)
    result.first = PhysDyn::Phys;

    // Remove the column tag
    erase(tags,FieldTag::Column);

  }

  // Get the size of what's left
  const auto size = tags.size();
  switch (size) {
    case 0:
      result.second = LayoutType::Scalar2D;
      break;
    case 1:
      // The only tag left should be 'Component', 'TimeLevel', 'Variable', or 'VerticalLevel
      if (tags[0]==FieldTag::Component || tags[0]==FieldTag::TimeLevel || tags[0]==FieldTag::Variable) {
        result.second = LayoutType::Vector2D;
      } else if (tags[0]==FieldTag::VerticalLevel) {
        result.second = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Possible scenarios:
      //  1) <Component|TimeLevel|Variable,VerticalLevel>
      //  2) <ComponentX,ComponentY>
      //  3) <Component,TimeLevel|Variable>
      //  4) <TimeLevel|Variable,TimeLevel|Variable>
      if (erase(tags,FieldTag::VerticalLevel)) {
        if (tags[0]==FieldTag::Component || tags[0]==FieldTag::TimeLevel || tags[0]==FieldTag::Variable) {
          result.second = LayoutType::Vector3D;
        }
      } else if (erase(tags,FieldTag::ComponentX)) {
        if (tags[0]==FieldTag::ComponentY) {
          result.second = LayoutType::Tensor2D;
        }
      } else if (erase(tags,FieldTag::Component)) {
        if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
          result.second = LayoutType::Tensor2D;
        }
      } else if (erase(tags,FieldTag::Variable)) {
        if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
          result.second = LayoutType::Tensor2D;
        }
      } else if (erase(tags,FieldTag::TimeLevel)) {
        if (tags[0]==FieldTag::TimeLevel) {
          result.second = LayoutType::Tensor2D;
        }
      }
      break;
    case 3:
      if (erase(tags,FieldTag::VerticalLevel)) {
        if (erase(tags,FieldTag::ComponentX)) {
          if (tags[0]==FieldTag::ComponentY) {
            result.second = LayoutType::Tensor3D;
          }
        } else if (erase(tags,FieldTag::Component)) {
          if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
            result.second = LayoutType::Tensor3D;
          }
        } else if (erase(tags,FieldTag::Variable)) {
          if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
            result.second = LayoutType::Tensor3D;
          }
        } else if (erase(tags,FieldTag::TimeLevel)) {
          if (tags[0]==FieldTag::TimeLevel) {
            result.second = LayoutType::Tensor3D;
          }
        }
      }
  }
  
  return result;
}

bool is_valid_layout (const FieldLayout& layout) {
  const auto& specs = get_layout_specs(layout);
  return specs.first!=PhysDyn::Undefined &&
         specs.second!=LayoutType::Invalid;
}

bool compatible_layout_types (const LayoutType& src, const LayoutType& tgt) {
  return src==tgt                 &&  // Same type
         src!=LayoutType::Invalid &&  // No invalid src type
         tgt!=LayoutType::Invalid;    // No invalid tgt type
}

bool is_tags_permutation (const FieldLayout& src, const FieldLayout& tgt) {
  return sort(src.tags())==sort(tgt.tags());
}

bool is_phys_dyn_remap (const FieldLayout& src, const FieldLayout& tgt) {
  const auto& src_pd = get_layout_specs(src).first;
  const auto& tgt_pd = get_layout_specs(tgt).first;
  return (src_pd==PhysDyn::Phys && tgt_pd==PhysDyn::Dyn) ||
         (src_pd==PhysDyn::Dyn  && tgt_pd==PhysDyn::Phys);
}

bool is_mpi_remap (const FieldLayout& src, const FieldLayout& tgt, const Comm& comm) {
  const auto& src_lt = get_layout_specs(src);
  const auto& tgt_lt = get_layout_specs(tgt);

  // Note: we do not allow to perform mpi redistribution unless the layouts are exactly the same.
  //       We could in principle do it together with, say, a phys-dyn remap, but it is much much
  //       easier to do the two remap in sequence rather than in one shot
  const bool same_specs = src_lt.first==tgt_lt.first && src.tags()==tgt.tags();
  if (!same_specs) {
    return false;
  }

  // Furthermore, we only allow mpi re-distribution over the Element tag for PhysDyn::Dyn, and
  // over Column for PhysDyn::Phys
  int pos;
  if (src_lt.first==PhysDyn::Dyn) {
    auto it = find(src.tags(),FieldTag::Element);
    pos = it - src.tags().begin();
  } else if (src_lt.first==PhysDyn::Phys) {
    auto it = find(src.tags(),FieldTag::Column);
    pos = it - src.tags().begin();
  } else {
    // No luck, this is a bad layout
    return false;
  }

  const int my_src_dim = src.dim(pos);
  const int my_tgt_dim = tgt.dim(pos);

  const int my_tgt_src_match = my_src_dim==my_tgt_dim ? 1 : 0;
  int global_tgt_src_match;
  MPI_Allreduce(&my_tgt_src_match,&global_tgt_src_match,1,MPI_INT,MPI_MIN,comm.mpi_comm());

  // For now, we assume it is an mpi re-distribution if at least one rank has two different
  // values in src and tgt for the distributed dimension (Element or Column).
  if (global_tgt_src_match) {
    return false;
  }

  int global_tgt_dim, global_src_dim;
  MPI_Allreduce(&my_src_dim,&global_src_dim,1,MPI_INT,MPI_SUM,comm.mpi_comm());
  MPI_Allreduce(&my_tgt_dim,&global_tgt_dim,1,MPI_INT,MPI_SUM,comm.mpi_comm());

  return global_src_dim==global_tgt_dim;
}

} // namespace scream
