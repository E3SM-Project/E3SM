# Atmospheric Processes

In EAMxx, `AtmosphereProcess` (AP) is an abstract class representing a portion
of the atmosphere timestep algorithm. In simple terms, an AP is an object that
given certain input fields performs some calculations to compute some output
fields. The concrete AP classes allow to create a buffer layer between
particular packages (e.g., dynamics dycore, physics parametrizations) and the
atmosphere driver (AD), allowing separation of concerns, so that the AD does
not need to know details about the package, and the package does not need to
know about the EAMxx infrastructure.

To enhance this separation of concerns, EAMxx implements two more classes for
handling APs:

- the concrete class `AtmosphereProcessGroup` (APG), which allows to group
  together a set of AP's, which can be seen from outside as a single process;
- the `AtmosphereProcessFactory` class, which allows an APG to create its
  internal processes without any knowledge of what they are.

This infrastructure allows the AD to view the whole atmosphere as a single APG,
and to be completely agnostic to what processes are run, and in which order.
This design allows to have a code base that is cleaner, self-container, and
easy to test via a battery of targeted unit tests.

In EAMxx, we already have a few concrete AP's, interfacing the AD to the
Hommexx non-hydrostatic dycore as well as some physics parametrizations (P3,
SHOC, RRMTPG, etc). In the next section we describe the interfaces of an AP
class, and we show an example of how to write a new concrete AP class.

## Atmosphere process interfaces

An AP has several interfaces, which can be grouped into three categories:

- initialization: these interfaces are used to create the AP, as well as to
  initialize internal data structures;
- run: these interfaces are used to make the AP compute its output fields from
  its input fields;
- finalization: these interfaces are used to perform any clean up operation
  (e.g., release files) before the AP is destroyed.

Among the above, the initialization sequence is the most complex, and consists
of several steps:

- The AD creates the APG corresponding to the whole atmosphere. As mentioned
  above, this phase will make use of a factory, which allows the AD to be
  agnostic to what is actually in the group. All AP's can start performing any
  initialization work that they can, but at this point they are limited to use
  only an MPI communicator as well as a list of runtime parameters (which were
  previously read from an input file).
- The AD passes a `GridsManager` to the AP's, so that they can get information
  about the grids they need. At this point, all AP's have all the information
  they need to establish the layout of the input and output fields they need,
  and can store a list of these "requests"
- After creating all fields (based on AP's requests), the AD passes a copy of
  each input and output field to the AP's. These fields will be divided in
  "required" and "computed", which differ in that the former are only passed
  to the AP's as 'read-only' fields (see the [field](field.md)
  documentation for more details)
- The AP's are queried for how much scratch memory they may need at run time.
  After all AP's communicate their needs, the AD will provide a pointer to
  scratch memory to the AP's. This is memory that can be used to initialize
  temporary views/fields or other internal data structures. All AP's are given
  the same pointer, which means no data persistence should be expected at run
  time between one timestep and the next.
- The AD calls the 'initialize' method on each AP. At this point, all fields
  are set, and AP's can complete any remaining initialization task

While the base AP class provides an (empty) implementation for some methods, in
case derived classes do not need a feature, some methods are purely virtual,
and concrete classes will have to override them. Looking at existing concrete
AP implementations is a good way to have a first idea of what a new AP class
needs to implement. Here, we show go over the possible implementation of these
methods in a hypothetical AP class. The header file may look something like
this

```c++
#include <share/atm_process/atmosphere_process.hpp>

class MyProcess : public AtmosphereProcess
{
public:
  using gm_ptr = std::shared_ptr<const GridsManager>;

  MyProcess(const ekat::Comm& comm, const ekat::ParameterList& pl);

  std::string name () const override { return "my_fancy_process"; }
  void set_grids (const gm_ptr& grids_manager) override;
  size_t requested_buffer_size_in_bytes () const override;
  void init_buffers (const ATMBufferManager& buffer_manager) override;
protected:

  void initialize_impl (const RunType run_type) override;
  void run_impl        (const double dt) override;
  void finalize_impl   () override;

  using view_1d = typename KokkosTypes<DefaultDevice>::view_1d<Real>;
  using view_2d = typename KokkosTypes<DefaultDevice>::view_2d<Real>;

  view_1d m_temp1;
  view_2d m_temp2;

  int m_ncols;
  int m_nlevs;
  bool m_has_blah;
};
```

A few comments:

- we added two views to the class, which are meant to be used to store
  intermediate results during calculations at runtime;
- there are other methods that the class can override (such as additional
  operations when the AD sets a field in the AP), but most AP's only need to
  override only these;
- we strongly encourage to add the keyword `override` when overriding a method;
  in case of small typos (e.g., missing a `&` or a `const`, the compiler will
  be erroring out, since the signature will not match any virtual method in the
  base class;
- `finalize_impl` is often empty; unless the AP is managing external resources,
  everything should be correctly released during destruction;
- the two methods for buffers can be omitted if the AP does not need any
  scratch memory (and the default implementation from the base class will be
  used).

Here is a possible implementation of the methods, with some inline comments to
explain

```c++
MyProcess::MyProcess (const ekat::Comm& comm, const ekat::ParameterList& pl)
 : AtmosphereProcess(comm,pl)
{
  // The base class copies pl into protected member m_params
  m_has_blah = m_params.get<bool>("enable_blah");
}

void MyProcess::set_grids (const std::shared_ptr<GridsManager>& gm)
{
  using namespace ekat::units;
  const auto nondim = Units::nondimensional();

  auto grid = gm->get_grid("Physics");
  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();

  // In these names, 2d refers to "horizontal only", while 3d to "horiz+vert".
  // But the grid stores dofs linearly, so there is only one array dimension
  FieldLayout scalar2d = grid->get_2d_scalar_layout();
  FieldLayout vector3d = grid->get_3d_vector_layout(true,2);

  // Declare fields needed:
  //  - Required: 'input' (read-only)
  //  - Computed: 'output'
  //  - Updated: 'input'+'output'
  // Tell the AD we need 'velocity' to accommodate a Pack scalar type
  add_field<Required>("coeff_2d",scalar2d,nondim,grid->name);
  add_field<Updated>("velocity",vector3d,m/s,grid->name,SCREAM_PACK_SIZE);
}

size_t MyProcess::requested_buffer_size_in_bytes ()
{
  // We use temp2 only if the blah feature is on
  return m_ncols + (m_has_blah ? m_ncols*m_nlev : 0);
}

void MyProcess::init_buffers (const ATMBufferManager& buffer_manager)
{
  auto mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  m_temp1 = view_1d<Real>(mem,m_ncols);
  mem += m_ncols;
  
  if (m_has_blah) {
    m_temp2 = view_2d<Real>(mem,m_ncols,m_nlevs);
    mem += m_ncols*m_nlevs;
  }

  // Make sure we use exactly the mem we said we would
  size_t used_mem = (mem - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
    "Error! Used memory != requested memory for MyProcess."
    "  used memory: " + std::to_string(used_mem) + "\n"
    "  requested: " + std::to_string(requested_buffer_size_in_bytes()) + "\n");
}

void MyProcess::initialize_impl(const RunType run_type)
{
  // Can complete any initialization of the background pkg
  my_process_pkg_init(m_has_blah);
}

void MyProcess:run_impl (const double dt)
{
  using Policy = typename KokkosTypes<DefaultDevice>::TeamPolicy
  using Member = typename KokkosTypes<DefaultDevice>::MemberType
  using PackT  = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  // Create team policy
  Policy policy(m_ncols,m_nlevs,1);

  // Create local copies of class members (or else use KOKKOS_CLASS_LAMBDA)
  auto temp1 = m_temp1;
  auto temp2 = m_temp2;
  auto do_blah = m_has_blah;

  // Get views from fields. We 
  auto coeff2d  = get_field_in("coeff_2d").get_view<const Real*>();
  auto velocity = get_field_out("velocity").get_view<PackT**>();

  // Since we process velocity with a Pack scalar type, find out how many packs
  // we have in each column
  auto nlevs_packs = ekat::PackInfo<SCREAM_PACK_SIZE>::num_packs(m_nlevs);

  // Call some function in the background pkg
  do_some_work (coeff_2d,velocity,temp1,temp2,do_blah);

  // Do some more work here
  auto work = KOKKOS_LAMBDA (const Member& team) {
    int icol = team.league_rank();
    auto col_work = [&](const int ilev) {
      velocity(icol,ilev) *= coeff_2d;
    };
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs_packs),col_work);
  };
  Kokkos::parallel_for(policy,work);
  Kokkos::fence();
}

void MyProcess::finalize_impl ()
{
  // If the background package needs to cleanup something, do it now
  my_process_pkg_cleanup();
}
```
