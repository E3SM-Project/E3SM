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

// TODO:
// 1. Use template to create formatter for various array types
// 2. Consider using some of the following for formatting
// View.rank()
// View.rank_dynamic()
// View.stride_(0, 1,2,3...)()
// View.span()
// View.size()
// View.span_is_contiguous()
// View.use_count()
// View.label()
// View.is_allocated()
// ExecSpace.name()
// ExecSpace.print_configuration(ostr);
// ExecSpace.print_configuration(ostr, detail);
// MemSpace.name()

template <>
struct fmt::formatter<OMEGA::HostArray1DReal> : fmt::formatter<std::string> {
   auto format(OMEGA::HostArray1DReal my, format_context &ctx)
       -> decltype(ctx.out()) {
#ifdef OMEGA_DEBUG
      return fmt::format_to(
          ctx.out(), "[data type of '{}' is HostArray1DReal.]", my.label());
#else
      return fmt::format_to(ctx.out(), "[data type of '' is HostArray1DReal.]");
#endif
   }
};

template <>
struct fmt::formatter<OMEGA::HostArray2DReal> : fmt::formatter<std::string> {
   auto format(OMEGA::HostArray2DReal my, format_context &ctx)
       -> decltype(ctx.out()) {
#ifdef OMEGA_DEBUG
      return fmt::format_to(
          ctx.out(), "[data type of '{}' is HostArray2DReal.]", my.label());
#else
      return fmt::format_to(ctx.out(), "[data type of '' is HostArray2DReal.]");
#endif
   }
};

#endif // OMEGA_LOG_FORMATTERS_H
