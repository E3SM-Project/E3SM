#include "field_identifier.hpp"

namespace scream
{

FieldIdentifier::FieldIdentifier (const std::string& name,
                                  const std::vector<FieldTag>& tags)
 : m_name (name)
 , m_tags (tags)
 , m_rank (tags.size())
 , m_dims (tags.size(),-1)
{
  update_identifier ();
}

void FieldIdentifier::set_dimension (const int idim, const int dimension) {
  error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  error::runtime_check(m_dims[idim] == -1,
                       "Error! You cannot reset field dimensions once set.\n");
  error::runtime_check(dimension>0, "Error! Dimensions must be positive.");
  m_dims[idim] = dimension;

  update_identifier ();
}


void FieldIdentifier::set_dimensions (const std::vector<int>& dims) {
  // Check, then set dims
  error::runtime_check(dims.size()==static_cast<size_t>(m_rank),
                       "Error! Input dimensions vector not properly sized.");
  for (int idim=0; idim<m_rank; ++idim) {
    set_dimension(idim,dims[idim]);
  }
}

void FieldIdentifier::update_identifier () {
  // Create a verbose identifier string.
  m_identifier = m_name + "<" + tag2string(m_tags[0]);
  for (int dim=1; dim<m_rank; ++dim) {
    m_identifier += "," + tag2string(m_tags[dim]);
  }
  m_identifier += ">(" + std::to_string(m_dims[0]);
  for (int dim=1; dim<m_rank; ++dim) {
    m_identifier += "," + std::to_string(m_dims[dim]);
  }
  m_identifier += ")";
}

bool operator== (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // First, compare names
  if (fid1.m_name != fid2.m_name) { return false; }

  // Then, compare ranks
  if (fid1.m_rank != fid2.m_rank) { return false; }

  // Check tags along each dimension
  for (int dim=0; dim<fid1.m_rank; ++dim) {
    if (fid1.m_tags[dim]!=fid2.m_tags[dim]) {
      return false;
    }
  }

  // Check each dimension
  for (int dim=0; dim<fid1.m_rank; ++dim) {
    if (fid1.m_dims[dim]!=fid2.m_dims[dim]) {
      return false;
    }
  }

  // Field identifiers are exactly the same. Return true.
  return true;
}

bool operator< (const FieldIdentifier& fid1, const FieldIdentifier& fid2) {
  // First, compare names
  if (fid1.m_name < fid2.m_name) {
    return true;
  } else if (fid1.m_name == fid2.m_name) {
    // If equal, compare ranks
    if (fid1.m_rank < fid2.m_rank) {
      return true;
    } else if (fid1.m_rank == fid2.m_rank) {
      // If still equal, check tags along each dimension
      for (int idim=0; idim<fid1.m_rank; ++idim) {
        if (fid1.m_tags[idim] < fid2.m_tags[idim]) {
          return true;
        } else if (fid2.m_tags[idim] > fid2.m_tags[idim]) {
          return false;
        }
      }

      // If still equal, check dimensions
      for (int idim=0; idim<fid1.m_rank; ++idim) {
        if (fid1.m_dims[idim] < fid2.m_dims[idim]) {
          return true;
        } else if (fid2.m_dims[idim] > fid2.m_dims[idim]) {
          return false;
        }
      }
    }
  }

  // Everything compared equal. Return false.
  return false;
}

} // namespace scream
