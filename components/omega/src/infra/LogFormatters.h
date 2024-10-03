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
#define GENERATE_FORMATTER(ARR, DIM, TYPE)                                     \
   template <>                                                                 \
   struct fmt::formatter<OMEGA::ARR##DIM##TYPE>                                \
       : fmt::formatter<std::string> {                                         \
      auto format(OMEGA::ARR##DIM##TYPE my,                                    \
                  format_context &ctx) -> decltype(ctx.out()) {                \
         return fmt::format_to(ctx.out(), "{}({}D:{})", my.label(), my.rank(), \
                               my.size());                                     \
      }                                                                        \
   };
#else
#define GENERATE_FORMATTER(ARR, DIM, TYPE)                      \
   template <>                                                  \
   struct fmt::formatter<OMEGA::ARR##DIM##TYPE>                 \
       : fmt::formatter<std::string> {                          \
      auto format(OMEGA::ARR##DIM##TYPE my,                     \
                  format_context &ctx) -> decltype(ctx.out()) { \
         return fmt::format_to(ctx.out(), "{}", my.label());    \
      }                                                         \
   };
#endif

#define GENERATE_FORMATTER_DIM(ARR, DIM) \
   GENERATE_FORMATTER(ARR, DIM, I4)      \
   GENERATE_FORMATTER(ARR, DIM, I8)      \
   GENERATE_FORMATTER(ARR, DIM, R4)      \
   GENERATE_FORMATTER(ARR, DIM, R8)

#define GENERATE_FORMATTER_ARR(ARR) \
   GENERATE_FORMATTER_DIM(ARR, 1D)  \
   GENERATE_FORMATTER_DIM(ARR, 2D)  \
   GENERATE_FORMATTER_DIM(ARR, 3D)  \
   GENERATE_FORMATTER_DIM(ARR, 4D)  \
   GENERATE_FORMATTER_DIM(ARR, 5D)

GENERATE_FORMATTER_ARR(HostArray)

#ifdef OMEGA_TARGET_DEVICE
GENERATE_FORMATTER_ARR(Array)
#endif

#undef GENERATE_FORMATTER_ARR
#undef GENERATE_FORMATTER_DIM
#undef GENERATE_FORMATTER

#endif // OMEGA_LOG_FORMATTERS_H
