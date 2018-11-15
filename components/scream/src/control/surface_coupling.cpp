#include "surface_coupling.hpp"

namespace scream {

SurfaceCoupling::SurfaceCoupling (const ParameterList& /*params*/)
{
  // Grab what you need from the parameter list
}

void SurfaceCoupling::initialize (const Comm& comm) {
  m_comm = comm;
  // Initialize the FieldRepository (FR) for host fields
  // (i.e., fields that are I/O w.r.t the coupler)
}

void SurfaceCoupling::run ( /* inputs ? */ ) {
  // Recall that the surface coupling can (and usually does) happen
  // in the middle of an atm time step. Therefore, we first export
  // atm output fields to the coupler, then import atm input fields
  // The order is in a way pointless (the same fields should not be
  // both exported and imported), but keeping this order makes it
  // more coherent with what is happening.

  // Copy (device->host) data to coupler views

  // Copy (host->device) data from coupler views
}

void SurfaceCoupling::finalize ( /* inputs? */ ) {
  // clean up
}

}  // namespace scream
