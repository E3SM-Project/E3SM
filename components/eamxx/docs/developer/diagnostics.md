# Diagnostics

A diagnostic is a field that is not part of the model state, and cannot
affect it: it is computed from the state, usually only to be written to file.
Diagnostics are created by output streams, from the names users request
(see the [user guide](../user/diags/index.md)), and are computed by the
stream before it writes. All diagnostics derive from `AbstractDiagnostic`,
in `src/share/diagnostics/abstract_diagnostic.hpp`.

## Lifecycle

A diagnostic goes through these steps:

- constructor: it receives a communicator, a parameter list and the grid. All
  diagnostics are created by a factory, so the signature must be exactly that
  of the base class. The constructor fills `m_field_in_names` with the names of
  the input fields.
- `set_input_field`: the customer (e.g., an output stream) sets each input field.
- `initialize`: called once all inputs are set. It calls `initialize_impl`,
  where the diagnostic can do any setup that needs the inputs (e.g., read their
  layout). When it returns, the output fields must be allocated.
- `init_timestep`: optional, called at the start of each step (e.g., `FieldPrev`
  copies its input here).
- `compute`: calls `compute_impl`, unless the inputs and the time stamp are the
  same as at the last call. Once `compute_impl` returns, every output gets the
  time stamp of the evaluation.
- `get`: returns the output field. It is allocated once, and reused at every step.

## Registration and names

A diagnostic class is registered in the `DiagnosticFactory`, in
`src/share/diagnostics/register_diagnostics.hpp`:

```cpp
diag_factory.register_product("WaterPath",&create_diagnostic<WaterPath>);
```

Output streams map a requested name to a diagnostic in `create_diagnostic`
(`src/share/io/eamxx_io_utils.cpp`):

- an expression, such as `T_mid.interp(p_mid=500,units='hPa')`, is lowered
  onto one of the diagnostics above (see `src/share/io/eamxx_dexpr_diags.cpp`);
- a name in the `named_variants` table selects a class, and possibly a parameter
  (e.g., `LiqWaterPath` is `WaterPath` with `water_kind=Liq`);
- any other name is taken as the name of a registered class.

## A diagnostic with one output

Most diagnostics compute one field, which they store in `m_diagnostic_output`,
in the constructor if its layout is known there, or else in `initialize_impl`:

```cpp
WaterPath::
WaterPath (const ekat::Comm& comm, const ekat::ParameterList& params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  // ... set m_qname from params ...
  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back(m_qname);

  FieldIdentifier fid (m_kind + name(), m_grid->get_2d_scalar_layout(), kg/m2, m_grid->name());
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}
```

`compute_impl` then reads the inputs from `m_fields_in` and writes into
`m_diagnostic_output`. Customers get the output with `get()`.

A customer may also set the parameter `output_field_name`, in which case
`initialize` renames the output to that name (e.g., the dexpr lowering names
the output of an expression after the expression). The renamed field is an
alias of the original, so it shares data, time stamp, and mask (whose name does
not change).

## A diagnostic with several outputs

Some diagnostics compute several fields at once, e.g. because they share inputs
and loops. `NumberPath` computes the liquid, ice and rain number paths in one
kernel. Such a diagnostic declares each output with `add_output`, in the
constructor:

```cpp
NumberPath::
NumberPath (const ekat::Comm& comm, const ekat::ParameterList& params,
            const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  m_field_in_names = {"pseudo_density","qc","nc","qi","ni","qr","nr"};

  auto create = [&](const std::string& kind) {
    FieldIdentifier fid(kind + name(), m_grid->get_2d_scalar_layout(), kg/(kg*m2), m_grid->name());
    Field f(fid);
    f.allocate_view();
    add_output(f);
    return f;
  };
  m_liq  = create("Liq");
  m_ice  = create("Ice");
  m_rain = create("Rain");
}
```

The rules for such a diagnostic are:

- all outputs must be declared in the constructor: an output stream needs all
  output names right after construction, to share one instance among requests
  of any of them (see below);
- each output must be allocated, be on the diagnostic grid, and have a name that
  no other output has (case insensitive);
- a diagnostic sets `m_diagnostic_output` or calls `add_output`, not both;
- `get()` errors out: customers use `get(name)` (case insensitive), or
  `get_outputs()`, which lists the outputs in declaration order;
- the outputs cannot be renamed via `output_field_name`. Users who want an
  output under another name write `<name> := <output>`;
- one `compute` evaluates all outputs, and stamps them all.

The inputs are those of all outputs together: a stream that requests only
`LiqNumberPath` still needs `qi`, `ni`, `qr` and `nr`, and still computes the
ice and rain paths.

`get_outputs()` is complete right after construction for a diagnostic with
several outputs. For a diagnostic with one output, it may be empty until
`initialize` is called.

To expose the outputs to users, add each name to the `named_variants` table in
`create_diagnostic`, mapped to the class (with no parameter). `create_diagnostic`
returns the class for any of the names, and the stream picks the output by name.

## How output streams use diagnostics

An output stream keeps the diagnostics it creates in a static repo, keyed by
name, so that two streams requesting the same field share one instance, which
is computed at most once per step (an instant stream evaluates its diagnostics
on output steps only). A diagnostic with several outputs is stored under each
output name, when the first of them is requested. Hence, requesting
`LiqNumberPath` in one stream and `IceNumberPath` in another gives one
`NumberPath`, which computes all three paths in one evaluation.

Each request adds only its own output to the stream. An output that is not
requested is computed, but is not written, and does not affect the stream.

## Masks and avg counts

An output may have a valid mask (see [Fields](field.md)). An averaged stream
then counts, at each entry, the steps where the output is valid, and divides by
that count. The count is shared by diagnostic outputs whose masks share a
name:

- a mask created with `create_valid_mask()` is named `<output>_valid_mask`,
  and its output gets a count of its own;
- a mask with any other name says where a condition holds (e.g.,
  `FieldAtPressureLevel` names it after the pressure level), and all
  diagnostic outputs with a mask of that name share one count.

A masked model field always gets a count of its own, whatever its mask is
named.

So outputs of one diagnostic can share a count by sharing a mask:

```cpp
m_a.create_valid_mask("my_region_mask");
m_z.set_valid_mask(m_a.get_valid_mask());
```

Masks with the same name must hold the same values, since the stream updates
the count once, and the horizontal remapper remaps such masks once. Hence
`initialize` requires that outputs whose masks have the same name share the
mask and have the same layout. An output that may contain fill values must
also call `set_may_be_filled(true)` on its header.

## Testing

- `src/share/diagnostics/tests/multi_output_test.cpp` shows how to test a
  diagnostic directly: create it, set its inputs, initialize, compute, and check
  each output.
- `src/share/io/tests/io_multi_output_diags.cpp` shows how to test a diagnostic
  through output streams: it registers a test diagnostic in the factory, writes
  INSTANT and AVERAGE files and reads them back, and inspects the diagnostics
  of a stream with a class that derives from `AtmosphereOutput`.
- `src/share/io/tests/io_diags.cpp` does the same for a diagnostic with one
  output.

## Converting a diagnostic to several outputs

A diagnostic with one output may call `add_output` once (in the constructor or
in `initialize_impl`) instead of setting `m_diagnostic_output`, with no change
in behavior: replace `m_diagnostic_output = f` with `add_output(f)`, after
allocating `f`, and keep a handle to `f` for `compute_impl`. To merge several
diagnostics that share their inputs into one with several outputs, as was done
for `NumberPath`, move the outputs to the constructor, and add their names to
`named_variants`.
