# Emulator API Reference

This section provides comprehensive API documentation for
the Emulator Components framework, auto-generated from source
code using Doxygen and mkdoxy.

!!! note
    The API documentation is automatically generated when building the docs.
    If you see broken links, run `mkdocs build` to generate the API reference files.

---

## Browse the API

### Classes

- **[Class List](../emulator/annotated.md)** — All classes with brief descriptions
- **[Class Index](../emulator/classes.md)** — Alphabetical class index
- **[Class Hierarchy](../emulator/hierarchy.md)** — Inheritance relationships

#### Class Members

- [All Members](../emulator/class_members.md)
- [Functions](../emulator/class_member_functions.md)
- [Variables](../emulator/class_member_variables.md)
- [Typedefs](../emulator/class_member_typedefs.md)
- [Enumerations](../emulator/class_member_enums.md)

---

### Namespaces

- **[Namespace List](../emulator/namespaces.md)** — All namespaces

#### Namespace Members

- [All Members](../emulator/namespace_members.md)
- [Functions](../emulator/namespace_member_functions.md)
- [Variables](../emulator/namespace_member_variables.md)
- [Typedefs](../emulator/namespace_member_typedefs.md)
- [Enumerations](../emulator/namespace_member_enums.md)

---

### Other References

- **[Functions](../emulator/functions.md)** — All functions
- **[Variables](../emulator/variables.md)** — All variables
- **[Macros](../emulator/macros.md)** — Preprocessor macros
- **[Files](../emulator/files.md)** — Source file listing
- **[Links](../emulator/links.md)** — External links

---

## Key Classes

Once generated, these are the main entry points:

### Core Framework

| Class | Description |
| ----- | ----------- |
| `emulator::EmulatorComp` | Abstract base class for all emulator components |
| `emulator::EmulatorAtm` | Atmosphere emulator implementation |
| `emulator::EmulatorContext` | Shared state for component lifecycle |

### Inference Backends

| Class | Description |
| ----- | ----------- |
| `emulator::inference::InferenceBackend` | NN inference interface |
| `emulator::inference::LibTorchBackend` | LibTorch/PyTorch C++ backend |
| `emulator::inference::StubBackend` | No-op testing backend |
| `emulator::inference::InferenceConfig` | Inference backend configuration |

### Diagnostics

| Class | Description |
| ----- | ----------- |
| `emulator::DerivedDiagnostic` | Base class for computed diagnostics |
| `emulator::HorizAvgDiagnostic` | Horizontal averaging diagnostic |
| `emulator::VertSliceDiagnostic` | Vertical slicing diagnostic |

---

## Namespace Structure

```cpp
namespace emulator {
    // Core component classes
    class EmulatorComp;     // Base class
    class EmulatorAtm;      // Atmosphere emulator
    
    // Configuration and I/O
    struct EmulatorContext;
    class EmulatorIO;
    
    // Diagnostics
    class DerivedDiagnostic;
    class HorizAvgDiagnostic;
    class VertSliceDiagnostic;
    
    namespace inference {
        // ML backend classes
        class InferenceBackend;
        class LibTorchBackend;
        class StubBackend;
        struct InferenceConfig;
        struct ValidationResult;
    }
}
```

## Documentation Style

We use standard Doxygen comment blocks:

```cpp
/**
 * @brief Short one-line description.
 *
 * @param input Description of input parameter
 * @return Description of return value
 * @see Related classes or functions
 */
```

## Related Documentation

- [Quick Start Guide](../user-guide/quickstart.md) — Getting started with EATM
- [Architecture](../tech-guide/architecture.md) — Design overview
- [Inference Backends](../tech-guide/inference-backends.md) — ML integration
