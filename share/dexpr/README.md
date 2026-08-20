# dexpr

A lexer and parser for the diagnostics expression language: the small DSL used
to describe derived diagnostics as expressions over model fields, rather than
as hand-written code per diagnostic.

Given a string like

```text
where(T_mid > 273.15, sum(qc, dims=[lev]))
```

`dexpr` produces an abstract syntax tree. It does not evaluate anything, and it
knows nothing about fields, grids or timesteps -- turning an AST into an actual
diagnostic is the caller's job.

## Layout

| Path | Contents |
| --- | --- |
| `include/dexpr/` | Public headers |
| `src/` | Lexer, parser, AST printing, token and precedence tables |
| `tests/` | Catch2 unit tests |
| `tools/` | `dexpr`, a small command line helper |

The grammar is handled by a hand-written lexer feeding a Pratt (precedence
climbing) parser. AST nodes are a `std::variant`, visited through
`Expression::visit`, so adding a node type is a matter of extending
`ExpressionVariant` and letting the compiler point at what no longer compiles.

The set of callable functions lives in one place, `supported_functions.hpp`.
Nothing else in the library is diagnostics-specific.

## Building

`dexpr` is deliberately standalone. It has no Kokkos, MPI, netCDF or EKAT
dependency, it is not part of the CIME or `csm_share` build, and it requires
only CMake 3.20 and a C++20 compiler:

```shell
./run_tests.sh
```

which is shorthand for

```shell
cmake -S . -B build
cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel
```

Catch2 v3 is used for the tests. An installed copy is used if there is one;
otherwise CMake fetches a pinned version.

### Compiler requirement

C++20, and in practice GCC 11 or newer: numeric literals are parsed and
printed with floating-point `<charconv>`, which is the binding constraint.

### Options

| Option | Default | Effect |
| --- | --- | --- |
| `DEXPR_ENABLE_TESTS` | `ON` | Build the unit tests |
| `DEXPR_ENABLE_TOOL` | `ON` | Build the `dexpr` command line tool |
| `DEXPR_WERROR` | `OFF` | Treat compiler warnings as errors; CI sets this |

## Testing

Every `TEST_CASE` is registered with ctest individually, so a failure names
itself. The same commands run in CI, across gcc and clang in both Debug and
Release, on every pull request touching this directory.

## Provenance

Vendored from [peterdschwartz/e3sm_diags_parser](https://github.com/peterdschwartz/e3sm_diags_parser)
at `d680d50`, MIT licensed; see `LICENSE`. The upstream name, and its `edp`
abbreviation, were dropped in favour of `dexpr` once the code moved here.
