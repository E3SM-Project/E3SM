#ifndef OMEGA_ANALYSISOPFACTORY_H
#define OMEGA_ANALYSISOPFACTORY_H

//===-- analysis/AnalysisOpFactory.h - Operator factory ---------*- C++ -*-===//
//
/// \file
/// \brief Defines the AnalysisOpFactory for runtime operator creation
///
/// The AnalysisOpFactory provides a runtime registry that maps operator type
/// names to constructor functions, enabling dynamic instantiation of templated
/// operators without hard-coded switch statements. The factory:
/// - Supports decentralized registration via template helpers
/// - Performs type-safe dispatch based on Field metadata (scalar type, rank,
///   memory location)
/// - Enables extensibility by allowing new operators to integrate without
///   modifying orchestration code
///
/// Each operator template (e.g., SpatialMaxOp<ArrayT>) is registered for all
/// supported array type combinations via registerAllArrayVariants(). At
/// operator creation time, the factory inspects the primary upstream Field's
/// metadata and selects the matching templated specialization.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Config.h"

#include <string>

namespace OMEGA {

/// Macro defining all supported array type combinations for analysis operators.
/// Used by registerAllArrayVariants() to register operator variants for all
/// combinations of scalar type, rank, and memory location. The macro takes a
/// single argument X which should be a macro that processes each
/// (DType, Rank, MemLoc, ArrayT) tuple.
#define OMEGA_ANALYSIS_ARRAY_TYPES(X)                      \
   /* 1D Arrays */                                         \
   X(ArrayDataType::I4, 1, ArrayMemLoc::Both, Array1DI4)   \
   X(ArrayDataType::I4, 1, ArrayMemLoc::Device, Array1DI4) \
   X(ArrayDataType::I8, 1, ArrayMemLoc::Both, Array1DI8)   \
   X(ArrayDataType::I8, 1, ArrayMemLoc::Device, Array1DI8) \
   X(ArrayDataType::R4, 1, ArrayMemLoc::Both, Array1DR4)   \
   X(ArrayDataType::R4, 1, ArrayMemLoc::Device, Array1DR4) \
   X(ArrayDataType::R8, 1, ArrayMemLoc::Both, Array1DR8)   \
   X(ArrayDataType::R8, 1, ArrayMemLoc::Device, Array1DR8) \
   /* 2D Arrays */                                         \
   X(ArrayDataType::I4, 2, ArrayMemLoc::Both, Array2DI4)   \
   X(ArrayDataType::I4, 2, ArrayMemLoc::Device, Array2DI4) \
   X(ArrayDataType::I8, 2, ArrayMemLoc::Both, Array2DI8)   \
   X(ArrayDataType::I8, 2, ArrayMemLoc::Device, Array2DI8) \
   X(ArrayDataType::R4, 2, ArrayMemLoc::Both, Array2DR4)   \
   X(ArrayDataType::R4, 2, ArrayMemLoc::Device, Array2DR4) \
   X(ArrayDataType::R8, 2, ArrayMemLoc::Both, Array2DR8)   \
   X(ArrayDataType::R8, 2, ArrayMemLoc::Device, Array2DR8) \
   /* 3D Arrays */                                         \
   X(ArrayDataType::I4, 3, ArrayMemLoc::Both, Array3DI4)   \
   X(ArrayDataType::I4, 3, ArrayMemLoc::Device, Array3DI4) \
   X(ArrayDataType::I8, 3, ArrayMemLoc::Both, Array3DI8)   \
   X(ArrayDataType::I8, 3, ArrayMemLoc::Device, Array3DI8) \
   X(ArrayDataType::R4, 3, ArrayMemLoc::Both, Array3DR4)   \
   X(ArrayDataType::R4, 3, ArrayMemLoc::Device, Array3DR4) \
   X(ArrayDataType::R8, 3, ArrayMemLoc::Both, Array3DR8)   \
   X(ArrayDataType::R8, 3, ArrayMemLoc::Device, Array3DR8) //               \
//   /* 4D Arrays */                                                       \
//   X(ArrayDataType::I4, 4, ArrayMemLoc::Both, Array4DI4)                 \
//   X(ArrayDataType::I4, 4, ArrayMemLoc::Device, Array4DI4)               \
//   X(ArrayDataType::I8, 4, ArrayMemLoc::Both, Array4DI8)                 \
//   X(ArrayDataType::I8, 4, ArrayMemLoc::Device, Array4DI8)               \
//   X(ArrayDataType::R4, 4, ArrayMemLoc::Both, Array4DR4)                 \
//   X(ArrayDataType::R4, 4, ArrayMemLoc::Device, Array4DR4)               \
//   X(ArrayDataType::R8, 4, ArrayMemLoc::Both, Array4DR8)                 \
//   X(ArrayDataType::R8, 4, ArrayMemLoc::Device, Array4DR8)               \
//   /* 5D Arrays */                                                       \
//   X(ArrayDataType::I4, 5, ArrayMemLoc::Both, Array5DI4)                 \
//   X(ArrayDataType::I4, 5, ArrayMemLoc::Device, Array5DI4)               \
//   X(ArrayDataType::I8, 5, ArrayMemLoc::Both, Array5DI8)                 \
//   X(ArrayDataType::I8, 5, ArrayMemLoc::Device, Array5DI8)               \
//   X(ArrayDataType::R4, 5, ArrayMemLoc::Both, Array5DR4)                 \
//   X(ArrayDataType::R4, 5, ArrayMemLoc::Device, Array5DR4)               \
//   X(ArrayDataType::R8, 5, ArrayMemLoc::Both, Array5DR8)                 \
//   X(ArrayDataType::R8, 5, ArrayMemLoc::Device, Array5DR8)               \

/// AnalysisOpFactory is a singleton factory class that manages runtime
/// registration and creation of analysis operators. The factory maintains a
/// registry mapping operator type keys to constructor functions. Keys are
/// formed from base operator name, array type, and memory location (e.g.,
/// "SpatialMax_Array2DR8_Device"). This enables type-safe dispatch: the
/// factory inspects the primary upstream Field's metadata at creation time
/// and selects the correct templated specialization.
///
/// All methods are static; the underlying registry is a Meyer's singleton
/// guaranteed to be initialized before first use. Operators self-register
/// during program initialization via registerAllArrayVariants().
class AnalysisOpFactory {
 public:
   /// Function signature for operator constructors. Takes upstream field names
   /// and configuration options, returns a unique_ptr to a new operator
   /// instance.
   using CreatorFunc = std::function<std::unique_ptr<AnalysisOperator>(
       const std::vector<std::string> &UpstreamNames, Config Options)>;

   /// Registers a single operator variant in the factory by associating a
   /// string label (key) with a constructor function. The label typically
   /// includes the operator type, array type, and memory location. Aborts if
   /// the label is already registered (duplicate key).
   static void registerOperator(
       const std::string &Label, ///< [in] unique key for this operator variant
       CreatorFunc Creator       ///< [in] constructor function for this variant
   );

   /// Creates an operator instance by inspecting the primary upstream Field's
   /// metadata (scalar type, rank, memory location) to build a fully-qualified
   /// type key, looking up the constructor in the registry, and invoking it
   /// with the provided arguments. Aborts if the operator type or array variant
   /// is not registered.
   static std::unique_ptr<AnalysisOperator>
   createOp(const std::string &OpType, ///< [in] operator type name
            const std::vector<std::string>
                &UpstreamNames, ///< [in] upstream field names
            Config Options      ///< [in] operator configuration
   );

   /// Returns a list of all registered operator variant keys. Useful for
   /// validation and generating informative error messages when an operator
   /// type is not found.
   static std::vector<std::string> availableOperators();

   /// Checks whether an operator variant with the given key is registered.
   /// Returns true if found, false otherwise.
   static bool
   hasOperator(const std::string &Type ///< [in] operator variant key
   );

   /// Registers all array type variants of a templated operator class by
   /// expanding OMEGA_ANALYSIS_ARRAY_TYPES over all supported (DType, Rank,
   /// MemLoc, ArrayT) combinations. For each combination, constructs a key
   /// from baseName + "_" + ArrayT + "_" + MemLoc and registers a lambda
   /// that instantiates OperatorTemplate<ArrayT>. This enables a single call
   /// to register dozens of operator variants.
   template <template <typename> class OperatorTemplate>
   static void registerAllArrayVariants(
       const std::string &baseName ///< [in] base operator name
   ) {
// Define a macro that registers one variant with registerOperator
#define REGISTER_VARIANT(dtype, rank, memloc, ArrayT)                     \
   registerOperator(baseName + "_" #ArrayT + "_" #memloc,                 \
                    [](const std::vector<std::string> &names, Config c) { \
                       return std::make_unique<OperatorTemplate<ArrayT>>( \
                           names, c);                                     \
                    });

      // Expand the macro over all array type combinations
      OMEGA_ANALYSIS_ARRAY_TYPES(REGISTER_VARIANT)
#undef REGISTER_VARIANT
   } // end registerAllArrayVariants

   /// Registers only 1D array type variants of a templated operator class.
   /// Used for operators that only support 1D arrays. Uses SFINAE to prevent
   /// template instantiation for non-1D array types at compile time.

   // Helper: enabled only for 1D arrays - registers the operator variant
   template <template <typename> class OperatorTemplate, typename ArrayT>
   static std::enable_if_t<ArrayRank<ArrayT>::Is1D, void>
   register1DVariantHelper(const std::string &baseName,
                           const std::string &typeName,
                           const std::string &memLocName) {
      registerOperator(baseName + "_" + typeName + "_" + memLocName,
                       [](const std::vector<std::string> &names, Config c) {
                          return std::make_unique<OperatorTemplate<ArrayT>>(
                              names, c);
                       });
   }

   // Helper: enabled only for non-1D arrays - does nothing (no-op)
   template <template <typename> class OperatorTemplate, typename ArrayT>
   static std::enable_if_t<!ArrayRank<ArrayT>::Is1D, void>
   register1DVariantHelper(const std::string &, const std::string &,
                           const std::string &) {
      // No-op for 2D, 3D, etc. - prevents template instantiation
   }

   // Public method that expands over all array types but only registers 1D
   template <template <typename> class OperatorTemplate>
   static void
   register1DVariants(const std::string &baseName ///< [in] base operator name
   ) {
// Define a macro that calls the SFINAE helper for each array type
#define REGISTER_1D_VARIANT(dtype, rank, memloc, ArrayT)                \
   register1DVariantHelper<OperatorTemplate, ArrayT>(baseName, #ArrayT, \
                                                     #memloc);

      // Expand over all array types; SFINAE filters out non-1D at compile time
      OMEGA_ANALYSIS_ARRAY_TYPES(REGISTER_1D_VARIANT)
#undef REGISTER_1D_VARIANT
   } // end register1DVariants

   /// Registers only 2D array type variants of a templated operator class.
   /// Used for operators that only support 2D arrays. Uses SFINAE to prevent
   /// template instantiation for non-2D array types at compile time.

   // Helper: enabled only for 2D arrays - registers the operator variant
   template <template <typename> class OperatorTemplate, typename ArrayT>
   static std::enable_if_t<ArrayRank<ArrayT>::Is2D, void>
   register2DVariantHelper(const std::string &baseName,
                           const std::string &typeName,
                           const std::string &memLocName) {
      registerOperator(baseName + "_" + typeName + "_" + memLocName,
                       [](const std::vector<std::string> &names, Config c) {
                          return std::make_unique<OperatorTemplate<ArrayT>>(
                              names, c);
                       });
   }

   // Helper: enabled only for non-2D arrays - does nothing (no-op)
   template <template <typename> class OperatorTemplate, typename ArrayT>
   static std::enable_if_t<!ArrayRank<ArrayT>::Is2D, void>
   register2DVariantHelper(const std::string &, const std::string &,
                           const std::string &) {
      // No-op for 1D, 3D, etc. - prevents template instantiation
   }

   // Public method that expands over all array types but only registers 2D
   template <template <typename> class OperatorTemplate>
   static void
   register2DVariants(const std::string &baseName ///< [in] base operator name
   ) {
// Define a macro that calls the SFINAE helper for each array type
#define REGISTER_2D_VARIANT(dtype, rank, memloc, ArrayT)                \
   register2DVariantHelper<OperatorTemplate, ArrayT>(baseName, #ArrayT, \
                                                     #memloc);

      // Expand over all array types; SFINAE filters out non-2D at compile time
      OMEGA_ANALYSIS_ARRAY_TYPES(REGISTER_2D_VARIANT)
#undef REGISTER_2D_VARIANT
   } // end register2DVariants

   /// Registers only 2D Real array type variants of a templated operator class.
   /// Used for operators that only support 2D Real arrays. Uses SFINAE to
   /// prevent template instantiation for non-2D or non-Real array types at
   /// compile time.

   // Helper: enabled only for 2D Real arrays - registers the operator variant
   template <template <typename> class OperatorTemplate, typename ArrayT>
   static std::enable_if_t<
       ArrayRank<ArrayT>::Is2D &&
           std::is_same_v<typename ArrayT::non_const_value_type, Real>,
       void>
   register2DRealVariantHelper(const std::string &baseName,
                               const std::string &typeName,
                               const std::string &memLocName) {
      registerOperator(baseName + "_" + typeName + "_" + memLocName,
                       [](const std::vector<std::string> &names, Config c) {
                          return std::make_unique<OperatorTemplate<ArrayT>>(
                              names, c);
                       });
   }

   // Helper: enabled only for non-2D or non-Real arrays - does nothing (no-op)
   template <template <typename> class OperatorTemplate, typename ArrayT>
   static std::enable_if_t<
       !ArrayRank<ArrayT>::Is2D ||
           !std::is_same_v<typename ArrayT::non_const_value_type, Real>,
       void>
   register2DRealVariantHelper(const std::string &, const std::string &,
                               const std::string &) {
      // No-op for 1D, 3D, or non-Real types - prevents template instantiation
   }

   // Public method that expands over all array types but only registers 2D Real
   template <template <typename> class OperatorTemplate>
   static void register2DRealVariants(
       const std::string &baseName ///< [in] base operator name
   ) {
// Define a macro that calls the SFINAE helper for each array type
#define REGISTER_2D_REAL_VARIANT(dtype, rank, memloc, ArrayT)               \
   register2DRealVariantHelper<OperatorTemplate, ArrayT>(baseName, #ArrayT, \
                                                         #memloc);

      // Expand over all array types; SFINAE filters out non-2D-Real at compile
      // time
      OMEGA_ANALYSIS_ARRAY_TYPES(REGISTER_2D_REAL_VARIANT)
#undef REGISTER_2D_REAL_VARIANT
   } // end register2DRealVariants

 private:
   /// Returns a reference to the static registry map (Meyer's singleton).
   /// The map associates operator variant keys with constructor functions.
   /// Static local variable ensures the registry is initialized exactly once
   /// before first use, avoiding static initialization order problems.
   static std::map<std::string, CreatorFunc> &registry() {
      static std::map<std::string, CreatorFunc> Registry;
      return Registry;
   }

   /// Helper function to construct an array type name string from Field
   /// metadata (data type, rank, memory location). Returns a string like
   /// "Array2DR8" used as part of the operator variant key.
   static std::string
   getArrayTypeName(ArrayDataType DType, ///< [in] scalar data type
                    I4 Rank,             ///< [in] array rank
                    ArrayMemLoc MemLoc   ///< [in] memory location
   );

}; // end class AnalysisOpFactory

} // end namespace OMEGA

#endif
