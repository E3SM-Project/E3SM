# Variables

## Naming

- Please name variables as ***descriptive nouns***.
    - That is, for a variable holding data related to *horizontal winds*, choose
      `horizontal_winds` or even `horiz_winds` versus `wh`.
    - In cases for which this may be untenable, add comments explaining
      non-obvious variable names.
- In general, variables should use `snake_case` and be entirely lowercase,
  unless capitalizing is a standard convention (e.g., `T` for temperature).
- We do not currently employ a stylistic distinction between different types of
  variables--e.g., `const`, `static`, member variables, etc.
    - To that end, some developers use a non-standardized convention of
      distinguishing member variables with a leading `m_`
      (`m_var` $\approx$ "my var"), or a trailing/leading underscore
      (`var_` or `_var`).

## Usage

- There is a balance to be struck between introducing extra, or intermediate,
  variables and avoiding excess allocations/de-allocations or minimizing clutter.
    - In line with the hierarchy of correct/fast/understandable, mentioned in
      the [Style Standards Overview](style.md#general-guidelines),
      it can sometimes improve matters to introduce an intermediate variable,
      pointer, subview, etc., rather than triply-indexing into a variable
      that's contained within a class-member struct.
    - This principle applies to calculations as well, as there are situations
      in which readability is improved by breaking complex operations into
      multiple steps.
