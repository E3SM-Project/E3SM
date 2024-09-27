#ifndef OMEGA_LOG_FORMATTERS_H
#define OMEGA_LOG_FORMATTERS_H
//===-- infra/LogFormatters.h - Omega specific log formatters --*- C++ -*-===//
//
/// \file
/// \brief Defines spdlog custom formatters
///
/// This header defines custom formatters for Omega logging.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include <spdlog/spdlog.h>

#ifdef OMEGA_DEBUG
#define GENERATE_FORMATTER(D, T)                                               \
   template <>                                                                 \
   struct fmt::formatter<OMEGA::Array##D##T> : fmt::formatter<std::string> {   \
      auto format(OMEGA::Array##D##T my,                                       \
                  format_context &ctx) -> decltype(ctx.out()) {                \
         return fmt::format_to(ctx.out(), "{}({}D:{})", my.label(), my.rank(), \
                               my.size());                                     \
      }                                                                        \
   };
#else
#define GENERATE_FORMATTER(D, T)                                             \
   template <>                                                               \
   struct fmt::formatter<OMEGA::Array##D##T> : fmt::formatter<std::string> { \
      auto format(OMEGA::Array##D##T my,                                     \
                  format_context &ctx) -> decltype(ctx.out()) {              \
         return fmt::format_to(ctx.out(), "{}", my.label());                 \
      }                                                                      \
   };
#endif

#define GENERATE_FORMATTER_DIM(D) \
   GENERATE_FORMATTER(D, I4)      \
   GENERATE_FORMATTER(D, I8)      \
   GENERATE_FORMATTER(D, R4)      \
   GENERATE_FORMATTER(D, R8)

GENERATE_FORMATTER_DIM(1D)
GENERATE_FORMATTER_DIM(2D)
GENERATE_FORMATTER_DIM(3D)
GENERATE_FORMATTER_DIM(4D)
GENERATE_FORMATTER_DIM(5D)

#undef GENERATE_FORMATTER_DIM
#undef GENERATE_FORMATTER

#endif // OMEGA_LOG_FORMATTERS_H
