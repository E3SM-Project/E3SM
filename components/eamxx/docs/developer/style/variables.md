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
- We do, however, ask that some convention is employed to identify member
  variables of a class.
      - Some commonly used techniques include distinguishing member variables
        with a leading `m_` (`m_var` $\approx$ "my var"), or a trailing/leading
        underscore (`var_` or `_var`).
      - That said, this convention should be ***locally consistent*** within a
        given class, struct, or otherwise.
- As to favoring ***nouns*** for naming, one should avoid using articles or
  verbs in variable names.
      - The key "exception" here would be prefixes like `is_`, `has_`, `can_`
        that are commonly used with boolean variables.

## Usage

### Intermediate Variables

There is a balance to be struck between introducing extra, or intermediate,
variables and avoiding excess allocations/de-allocations or minimizing clutter.

- In line with the hierarchy of correct/fast/understandable, mentioned in
  the [Style Standards Overview](style.md#general-guidelines),
  it often improves matters to introduce an intermediate variable,
  pointer, subview, etc., rather than triply-indexing into a variable
  that's contained within a class-member struct.
- This principle applies to calculations as well, as there are situations
  in which readability is improved by breaking complex operations into
  multiple steps.
- **Generally speaking**, unless the code is highly-optimized and/or
  performance-critical, **usage of intermediate variables is preferred**.
      - Here, readability and bug-safety are the top concerns, so the clarity
        provided by descriptive intermediate variables and resultant shorter
        lines, trumps most other factors.

### Descriptive Power vs. Brevity in Naming

There is obvious a crossover point at which a name becomes **too** descriptive
or simply too long, so as to be cumbersome to use.
Thus, the issue of descriptive naming should be balanced against brevity in
name choice, and we give some examples and rules of thumb.

- Consider synonyms that may shorten a variable name without sacrificing clarity.
- Make use of common word truncations or abbreviations.
      - E.g., `mgr` instead of `manager` or `dev` in place of `device`.
      - However, avoid non-standard abbreviations that may only be common
        parlance within a small group.
            - This includes names that correspond to standard arithmetic variables
              in an equation (`gamma`, `y_hat`), ***unless*** the equation is
              included as a comment and the connection is readily apparent.
- The converse of the previous is avoid making contractions that negatively
  impact clarity.
      - A notorious example of this is removing vowels.
          - E.g., `loop_cnt`, rather than `loop_count` provides little savings,
            at the expense of clarity.[^funfact]

[^funfact]: Much like some [dangerous options-trading strategies](https://www.investopedia.com/ask/answers/050115/what-types-options-positions-create-unlimited-liability.asp),
this practice offers limited upside, yet the potential downside is unlimited.
