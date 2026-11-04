#ifndef SCREAM_FIELD_PCH_HPP
#define SCREAM_FIELD_PCH_HPP

// Precompiled header for field library ETI files.
// This header bundles all the heavy template headers needed by the
// 36 ETI .cpp files, allowing them to reuse the precompiled content
// instead of independently parsing these headers.

#include "share/field/field.hpp"
#include "share/field/field_get_view_impl.hpp"
#include "share/field/field_update_impl.hpp"
#include "share/field/field_update_masked.hpp"
#include "share/field/field_update_fill_aware.hpp"

#endif // SCREAM_FIELD_PCH_HPP
