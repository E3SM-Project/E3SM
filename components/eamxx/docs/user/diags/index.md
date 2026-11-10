# Online diagnostics

EAMxx has facilities to output optional diagnostics
that are computed during runtime. These diagnostics
are designed generically and composably, and are requestable by users.

## Available diagnostics

- [Field contraction](field_contraction.md)
- [Conditional sampling](conditional_sampling.md)
- [Binary arithmetics](binary_ops.md)
- [Vertical derivative](vert_derivative.md)
- [Previous-timestep field](field_prev.md)
- [Field divided by timestep](field_over_dt.md)
- [Diagnostic expressions](expressions.md)

Derived diagnostics are requested as [expressions](expressions.md), e.g.
`T_500 := T_mid.interp(p_mid=500,units='hPa')`.

More details to follow.
