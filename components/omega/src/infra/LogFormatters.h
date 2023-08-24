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

template <>
struct fmt::formatter<OMEGA::ArrayHost1DReal> : fmt::formatter<std::string> {
   auto format(OMEGA::ArrayHost1DReal my, format_context &ctx)
       -> decltype(ctx.out()) {
#ifdef OMEGA_DEBUG
      return fmt::format_to(
          ctx.out(), "[data type of '{}' is ArrayHost1DReal.]", my.label());
#else
      return fmt::format_to(ctx.out(), "[data type of '' is ArrayHost1DReal.]");
#endif
   }
};

template <>
struct fmt::formatter<OMEGA::ArrayHost2DReal> : fmt::formatter<std::string> {
   auto format(OMEGA::ArrayHost2DReal my, format_context &ctx)
       -> decltype(ctx.out()) {
#ifdef OMEGA_DEBUG
      return fmt::format_to(
          ctx.out(), "[data type of '{}' is ArrayHost2DReal.]", my.label());
#else
      return fmt::format_to(ctx.out(), "[data type of '' is ArrayHost2DReal.]");
#endif
   }
};

#endif // OMEGA_LOG_FORMATTERS_H
