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

## Diagnostics with several outputs

Some diagnostics compute several fields at once. E.g., `LiqNumberPath`,
`IceNumberPath` and `RainNumberPath` are computed together. Each field is
requested by its own name, and requesting any of them, in one or more output
streams, computes them all together, at most once per step (an instant stream
evaluates its diagnostics on output steps only). Only the requested ones are
written. A request of any of them needs the inputs of all of them (e.g.,
`LiqNumberPath` needs the ice and rain tracers too).
To write one under another name, use `<name> := <field>`, e.g.
`lnp := LiqNumberPath`.

More details to follow.
