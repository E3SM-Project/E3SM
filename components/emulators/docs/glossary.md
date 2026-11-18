# Glossary

We use the following terms and definitions to describe E3SM's emulator infrastructure.

- **`CoupledProcess`**: a C++ concept that defines an interface that types can implement to
  implement **processes**.
- **Process**: a type that implements the `CoupledProcess` concept to interact with other processes
  using the emulator framework. Dynamical components, emulators, and intra-component
  parameterizations can all be represented by **processes**.
- **Coupler**: a mechanism responsible for exchanging data between **processes**
- **Imported Field**: a field whose values come from another `process` via the `coupler`
- **Exported Field**: a field sent to the `coupler` to be sent to another `process`
- **Coupler Driver**: object that instantiates and orchestrates `processes` with the `coupler`
- **Coupling Interval**: The number of base timesteps a `process` takes before coupling.
- **Slow Process**: A `process` is slow if its `coupling interval` is larger than another `process`
- **External Process**: `process` not owned by the `CouplerDriver`
- **Harness**: Helper object for binding external `processes` to `Coupler Driver`
- **Prognostic Field**: variable that acts as both an input and output s.t. its computed value at step `t_i` is needed for step `t_{i+1}`
- **Field Registry**: List of fields that are used for any `process` in the simulation
- **Tensor**: a multidimensional array that stores quantities represented in an emulator
- **Layout**: the set of dimensions describing the shape of a **tensor**.
