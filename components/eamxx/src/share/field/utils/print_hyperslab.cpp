#include "share/field/field_utils.hpp"

namespace scream {

namespace impl {

template<typename T>
void print_field_hyperslab (const Field& f,
                            std::vector<FieldTag> tags,
                            std::vector<int> indices,
                            std::ostream& out,
                            const int orig_rank,
                            const size_t curr_idx)
{
  // General idea: call f.subfield with the proper index, and recurse
  // until all indices are exausted, then print the field that is left.
  //
  // We keep all the tags/indices we slice away, since we need them at
  // the end of recursion when we print the info of the field location.
  // E.g., if f has tags/dims <COL,CMP,LEV>/(2,3,4) and we subview at
  // tags/indices <COL,LEV>/(0,1), we want to print something like
  //   f(0,:,1):
  //     0.123, 0.456, 0.789

  EKAT_REQUIRE_MSG (tags.size()==indices.size(),
      "Error! Tags vector size differs from indices vector size.\n");

  const auto& layout = f.get_header().get_identifier().get_layout();

  // Get layout of original field (before all the slicing happened)
  auto get_orig_header = [&]() -> std::shared_ptr<const FieldHeader> {
    auto fh = f.get_header_ptr();
    while (fh->get_identifier().get_layout().rank()<orig_rank) {
      fh = fh->get_parent();
    }
    return fh;
  };

  // Get a vector of strings, which we'll use to print the info of the field
  // slice/entry currently printed. We can already fill the entries corresponding
  // to dimensions we sliced away.
  auto get_dims_str = [&](const FieldLayout& orig_layout) -> std::vector<std::string> {
    const int orig_rank = orig_layout.rank();
    const int num_tags  = tags.size();
    std::vector<std::string> dims_str (orig_rank,"");
    for (int ii=0, jj=0; ii<orig_rank && jj<num_tags; ++ii) {
      if (orig_layout.tag(ii)==tags[jj]) {
        // Was sliced. Store slice idx as a stirng
        dims_str[ii] = std::to_string(indices[jj]);
        ++jj;
      }
    }
    return dims_str;
  };

  // Get the dimensions (in the original field) that were *NOT* sliced away,
  // since we'll have to loop over them. We'll use these indices to add the
  // missing info in the dims_str vector<string> as we loop and print.
  auto get_dims_left = [&](const FieldLayout& orig_layout) -> std::vector<int> {
    const int orig_rank = orig_layout.rank();
    const int num_tags  = tags.size();
    std::vector<int> dims_left;
    for (int ii=0, jj=0; ii<orig_rank; ++ii) {
      if (jj<num_tags && orig_layout.tag(ii)==tags[jj]) {
        // Was sliced. Skip
        ++jj;
      } else {
        // Was not sliced.
        dims_left.push_back(ii);
      }
    }
    return dims_left;
  };

  constexpr int max_per_line = 5;
  if (curr_idx==tags.size()) {
    // All slices have been taken. Print the whole input field.
    const auto& orig_layout = get_orig_header()->get_identifier().get_layout();
    const auto dims_left = get_dims_left(orig_layout);
    auto dims_str = get_dims_str(orig_layout);

    // NOTE: because print_field_hyperslab() is only used for testing, we
    // generalize by always calling get_strided_view(), even if it could be
    // contiguous.
    // This allows us to call this function on any type of field,
    // including multi-slice subfields
    f.sync_to_host();
    const int rank = layout.rank();
    out << "     " << f.name() << orig_layout.to_string() << "\n\n";
    switch (rank) {
      case 0:
      {
        out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
        // NOTE: add ", " at the end, to make rank0 behave the same as other ranks,
        //       for the sake of any script trying to manipulate output
        out << "\n    " << f.get_strided_view<const T,Host>()() << ", \n";
        break;
      }
      case 1:
      {
        dims_str[dims_left[0]] = ":";
        out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
        auto v = f.get_strided_view<const T*,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          if (i%max_per_line==0) {
            out << "\n    ";
          }
          out << v(i) << ", ";
        }
        out << "\n";
        break;
      }
      case 2:
      {
        dims_str[dims_left[1]] = ":";
        auto v = f.get_strided_view<const T**,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
          for (int j=0; j<layout.dim(1); ++j) {
            if (j%max_per_line==0) {
              out << "\n    ";
            }
            out << v(i,j) << ", ";
          }
          out << "\n";
        }
        break;
      }
      case 3:
      {
        dims_str[dims_left[2]] = ":";
        auto v = f.get_strided_view<const T***,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          for (int j=0; j<layout.dim(1); ++j) {
            dims_str[dims_left[1]] = std::to_string(j);
            out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
            for (int k=0; k<layout.dim(2); ++k) {
              if (k%max_per_line==0) {
                out << "\n    ";
              }
              out << v(i,j,k) << ", ";
            }
            out << "\n";
          }
        }
        break;
      }
      case 4:
      {
        dims_str[dims_left[3]] = ":";
        auto v = f.get_strided_view<const T****,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          for (int j=0; j<layout.dim(1); ++j) {
            dims_str[dims_left[1]] = std::to_string(j);
            for (int k=0; k<layout.dim(2); ++k) {
              dims_str[dims_left[2]] = std::to_string(k);
              out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
              for (int l=0; l<layout.dim(3); ++l) {
                if (l%max_per_line==0) {
                  out << "\n    ";
                }
                out << v(i,j,k,l) << ", ";
              }
              out << "\n";
            }
          }
        }
        break;
      }
      case 5:
      {
        dims_str[dims_left[4]] = ":";
        auto v = f.get_strided_view<const T*****,Host>();
        for (int i=0; i<layout.dim(0); ++i) {
          dims_str[dims_left[0]] = std::to_string(i);
          for (int j=0; j<layout.dim(1); ++j) {
            dims_str[dims_left[1]] = std::to_string(j);
            for (int k=0; k<layout.dim(2); ++k) {
              dims_str[dims_left[2]] = std::to_string(k);
              for (int l=0; l<layout.dim(3); ++l) {
                dims_str[dims_left[3]] = std::to_string(l);
                out << "  " << f.name() << "(" << ekat::join(dims_str,",") << ")";
                for (int m=0; m<layout.dim(3); ++m) {
                  if (m%max_per_line==0) {
                    out << "\n    ";
                  }
                  out << v(i,j,k,l,m) << ", ";
                }
                out << "\n";
              }
            }
          }
        }
        break;
      }
      default:
        EKAT_ERROR_MSG (
            "Unsupported rank in print_field_hyperslab.\n"
            "  - field name  : " + f.name() + "\n"
            "  - field layout (upon slicing): " + layout.to_string() + "\n");
    }
  } else {
    auto tag = tags[curr_idx];
    auto idx = indices[curr_idx];

    auto it = ekat::find(layout.tags(),tag);
    EKAT_REQUIRE_MSG (it!=layout.tags().end(),
        "Error! Something went wrong while slicing field.\n"
        "  - field name  : " + f.name() + "\n"
        "  - field layout: " + layout.to_string() + "\n"
        "  - curr tag    : " + e2str(tag) + "\n");
    auto idim = std::distance(layout.tags().begin(),it);

    EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Cannot subview field for printing.\n"
        "  - field name  : " + f.name() + "\n"
        "  - field layout: " + layout.to_string() + "\n"
        "  - loc tags    : <" + ekat::join(tags,",") + ">\n"
        "  - loc indices : (" + ekat::join(indices,",") + ")\n");

    auto sub_f = f.subfield(idim,idx);
    return print_field_hyperslab<T>(sub_f,tags,indices,out,orig_rank,curr_idx+1);
  }
}

} // namespace impl

void print_field_hyperslab (const Field& f,
                            const std::vector<FieldTag>& tags,
                            const std::vector<int>& indices,
                            std::ostream& out)
{
  const auto dt = f.data_type();
  const auto rank = f.rank();
  const auto& fl = f.get_header().get_identifier().get_layout();

  EKAT_REQUIRE_MSG (rank>=static_cast<int>(tags.size()),
      "Error! Requested location incompatible with field rank.\n"
      "  - field name: " + f.name() + "\n"
      "  - field rank: " + std::to_string(rank) + "\n"
      "  - requested indices: (" + ekat::join(indices,",") + ")\n");

  const int num_indices = indices.size();
  for (int i=0; i<num_indices; ++i) {
    EKAT_REQUIRE_MSG ( indices[i]>=0 && indices[i]<fl.dim(tags[i],false),
      "Error! Requested slice index is out of bound.\n"
      "  - field name: " + f.name() + "\n"
      "  - field layout: " + fl.to_string() + "\n"
      "  - hyperslab tags: (" + ekat::join(tags2str(tags),",") + ")\n"
      "  - hyperslab indices: (" + ekat::join(indices,",") + ")\n");
  }

  switch (dt) {
    case DataType::IntType:
      impl::print_field_hyperslab<int>(f,tags,indices,out,rank,0);
      break;
    case DataType::FloatType:
      impl::print_field_hyperslab<float>(f,tags,indices,out,rank,0);
      break;
    case DataType::DoubleType:
      impl::print_field_hyperslab<double>(f,tags,indices,out,rank,0);
      break;
    default:
      EKAT_ERROR_MSG ("[print_field_hyperslab] Error! Invalid/unsupported data type.\n"
          " - field name: " + f.name() + "\n");
  }
}

} // namespace scream
