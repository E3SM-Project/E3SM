# Functions and Methods

### Naming

- Please name functions (methods) to be ***descriptive*** and ***verb***_-ish_ so
  that the name describes what the function *does*.
    - For example, `install_flux_capacitor()` instead of `f_capacitor_method()`
      or `ifc()`.
    - To note, if your function cannot be named verb-ishly, that is probably a
      sign that there is a fundamental issue with your function.
- In general, functions should use `snake_case` and not contain capitalization,
  unless capitalizing that quantity is a standard convention (e.g.,
  `find_Doc_Brown()`).

### Usage

- In situations where multiple class-member variables are passed as function
  arguments, favor passing the object, rather than the individual variables,
  for the sake of shorter lines and function signatures.
- ==FIXME:== Passing as reference vs. value??
