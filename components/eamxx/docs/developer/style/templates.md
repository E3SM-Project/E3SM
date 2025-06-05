# Templating

- Templating and polymorphism are arguably the most powerful features of C++,
  but with great power comes great potential pain.
      - EAMxx makes extensive use of templates, throughout, and this is encouraged
        but should be done judiciously.
      - Adding complexity should only be done in service of improved
        functionality or speed because these typically come at the price of
        clarity and readability.
- Template parameters should mostly obey the same naming conventions of the
  type they correspond to.
      - `lower_snake_case` for variables.
          - E.g., variables of integral types or objects.
      - `UpperCamelCase` for types, classes, or structs.
      - However, terseness that follows EAMxx or Kokkos conventions can improve
        readability, though favoring ***descriptiveness*** should be the default
        case.
          - Take this `update()` function declaration as an example.

              ```c++
              template<int *atom_counter, HostOrDevice HD = Device, typename ST>
              void update(const Field &x, const ST alpha, const ST beta);
              ```

          - The first template parameter is an integer-pointer that self-describes.
          - The second template parameter, `HD` (`enum` type), is
            tersely-named but the type provides sufficient descriptive
            information.
          - The second and third template parameters, including `ST`
            ("scalar type"), are standard conventions used throughout EAMxx, and
            the abbreviation saves an enormous amount of space in aggregate.
