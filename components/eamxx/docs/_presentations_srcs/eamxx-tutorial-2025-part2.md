---

theme: gaia
_class: lead
paginate: true
header: "<!-- paginate: true -->"
style: |
  section::after {
    content: attr(data-marpit-pagination);
    font-size: 0.6em;
    position: absolute;
    top: 10px;
    right: 10px;
  }
backgroundColor: #fff
backgroundImage: url('./background.png')

---

<!-- markdownlint-disable MD013 MD033 MD041 -->

# 2025 E3SM All-Hands EAMxx Tutorial

## Part 2: Adding a new process to EAMxx

---

## Agenda

Part 1: Running and changing an EAMxx case
Runtime options, controlling output, DPxx, running unit tests

**Part 2: Adding a new process to EAMxx**
**AtmosphereProcess interface, adding a process, running it**

---

## Adding a new process to EAMxx

- EAMxx expects a process interface
  - An instant of AtmosphereProcess
- Must add new process to CMake, MCT, etc.
- Runtime options

---

## POMPEI

POMPEI: Practicing On Manually Adding Processes to EAMxx Interfaces

- Idealized injection of a new volcanic tracer `ash`
- Injection details are defined in `pompei.hpp` and `pompei.cpp` files
- Interfacing with EAMxx via an `AtmosphereProcess` _instance_ (`*_interface.*pp` files)
  - `set_grid`, `initialize_impl`, `run_impl`, `finalize_impl`

---

## Add `pompei` skeleton to EAMxx

In the EAMxx physics directory, create a pompei folder and copy templates

```shell
cd components/eamxx/src/physics
mkdir pompei
cp /global/cfs/cdirs/e3sm/eamxx-tutorial-2025/source/eamxx_interface/eamxx_template_process_interface.*  .
mv eamxx_template_process_interface.cpp eamxx_pompei_process_interface.cpp
mv eamxx_template_process_interface.hpp eamxx_pompei_process_interface.hpp
cp /global/cfs/cdirs/e3sm/eamxx-tutorial-2025/source/pompei/pompei.* ./
```

---

## `pompei.hpp`

```cpp
#ifndef POMPEI_INTERNAL_HPP
#define POMPEI_INTERNAL_HPP

namespace pompei {

double ash_emission_rate(const double days_after_eruption);

}  // namespace pompei

#endif  // POMPEI_INTERNAL_HPP

```

---

## `pompei.cpp`

```cpp
#include "pompei.hpp"
#include <cmath>

namespace pompei {

double ash_emission_rate(const double days_since_eruption) {
  if(days_since_eruption <= 0) return 0;

  double initial_emission = 1e4;
  double decay_rate       = -0.05;
  return initial_emission * exp(days_since_eruption * decay_rate);
}

}  // namespace pompei
```

---

## An `AtmosphereProcess` in EAMxx

- Intantiated inside an interface layer and contains _at a minimum_:
  - `set_grid`
  - `initialize_impl`
  - `run_impl`
  - `finalize_impl`
- `pompei` is an _instance_ of `AtmosphereProcess`

We will walk through the interface layer in detail next

---

## The interface hpp (header) file

- Stores the definitions of all interface methods
- Developer can define internal function and variables

Tasks:

- Replace `AP_TEMPLATE` with `POMPEI`
- Set a desired name for the process
- Add global parameters for the process under `protected` section

See next slide for complete hpp file!

---

```cpp
/*
 * First define this process intereface as a new header file
 */
#ifndef EAMXX_POMPEI_PROCESS_INTERFACE_HPP
#define EAMXX_POMPEI_PROCESS_INTERFACE_HPP
/*
 * Any include statements.  By default you want to include the
 * atmosphere_process.hpp header
 */
#include "share/atm_process/atmosphere_process.hpp"

namespace scream {

/*
 * Class definition for new process interface.  Change "AP_TEMPLATE" to a name
 * that makes sense.
 */
class POMPEI : public AtmosphereProcess {
  // Define public functions
 public:
  POMPEI(const ekat::Comm &comm, const ekat::ParameterList &params);

  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  std::string name() const override { return "my_disastrous_pompei"; }

  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;
// Define the protected functions, usually at least initialize_impl, run_impl
// and finalize_impl, but others could be included.  See
// eamxx_template_process_interface.cpp for definitions of each of these.
#ifndef KOKKOS_ENABLE_CUDA
 protected:
#endif
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  // Could also define here universal variables for the process interface.
  util::TimeStamp m_eruption_start;
  std::shared_ptr<const AbstractGrid> m_grid;
  int m_ncols, m_nlevs;
  Field m_emission_mask;
};

}  // namespace scream

#endif  // EAMXX_POMPEI_PROCESS_INTERFACE_HPP
```

---

## The interface cpp (source) file

Replace `AP_TEMPLATE` with `POMPEI` and then:

- Define a constructor for the process
- Define `set_grids`, which tells EAMxx variables and grid
- Define `initialize_impl`, which initializes the process
- Define `run_impl`, what this process does each time step it is called
- Define `finalize_impl`, and last actions to take at the end of the simulation

---

## Include statements (complete section on next slide)

Top of file, use of #include statements to access any other codes we need

We add the header fo the process and other files we need

---

```cpp
/*
 * Add any #include statements pointing to headers needed by this process.
 * Typical headers would be
 *     eamxx_$PROCESS_process_interface.hpp  (the eamxx interface header file
 * for this process) physics/share/physics_constants.hpp   (a header storing all
 * physical constants in eamxx) $PROCESS.hpp                          (a header
 * file for the process itself)
 */
#include "eamxx_pompei_process_interface.hpp"  // interface header

#include "physics/share/physics_constants.hpp"  // constants in eamxx/src (eamxx/src is in compiler path)
#include "pompei.hpp"                           // emission rate
#include "share/field/field_utils.hpp"  // field utilities

namespace scream {

```

---

## Constructor (complete section on next slide)

First section of code we see in the template

Instantiates an `AtmosphereProcess`

Note that `params.get<T>("name")` only supports
`T=int, bool, double, std::string, std::vector`

---

```cpp
/*-----------------------------------------------------------------------------------------------
 * The Constructor for this interface
 *
 * Inputs:
 *     comm - an EKAT communication group
 *     params - a parameter list of options for the process.
 */

POMPEI::POMPEI(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  // The 'params' variable will hold all runtime options.
  // Note that `params.get<X>("Y") is the syntax to get the parameter labeled
  // "Y" from the list which is of type X.
  //
  // ex:
  // m_Y = params.get<X>("Y");

  m_eruption_start =
      util::str_to_time_stamp(params.get<std::string>("eruption_date"));
}
```

---

## `set_grids`

- EAMxx stores grid info
- FieldManager stores fields
- Fields have:
  - Pointers to actual data
  - Units
  - Data layout
  - Grid
- `add_field`/`add_tracer` calls

---

```cpp
/*-----------------------------------------------------------------------------------------------
 * set_grids(grids_manager)
 *
 * set_grids establishes which grid this process will be run on.
 * It is also where all fields (variables) that the process will need access to.
 *
 * Inputs:
 *     grids_manager - a grids manager object which stores all grids used in
 * this simulation
 */
void POMPEI::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  // Some typical namespaces used in set_grids,
  //   using namespace ShortFieldTagsNames;
  //   using PC = scream::physics::Constants<Real>;

  // Specify which grid this process will act upon, typical options are
  // "Dynamics" or "Physics". auto m_grid = grids_manager->get_grid("Physics");

  // When declaring a field for use we need:
  //
  // to set the units.  These can be found from ekat::units::Units
  // For example,
  //   constexpr auto nondim = ekat::units::Units::nondimensional();
  //
  // to set the layout, which is how the field is configured in space
  // For example,
  //   const auto layout = m_grid->get_3d_scalar_layout(true);

  // We can now add all fields used by this process and designate how they will
  // be used, The format for declaring a field is:
  //   add_field<MODE>(NAME,LAYOUT,UNITS,GRID_NAME);
  // where,
  //   MODE is either "Required", "Computed" or "Updated".
  //     "Required" = input only  (process cannot change this field)
  //     "Computed" = output only (field doesn't have to exist, process will set
  //     values for this field) "Updated"  = input and output (process is able
  //     to update field which should already exist)
  //   NAME is the name of the field as a string
  //   LAYOUT is the physical layout of the field on the grid, e.g. 2D, 3D.  See
  //   layout example above. UNITS are the units of the field, see units example
  //   above GRID_NAME is the name of the grid this process is run on, typically
  //   "m_grid->name()"
  //
  //   *Note, there is a special add_field call for tracers which has a
  //   different signature.  The format is,
  //     add_tracer<MODE>(NAME,GRID,UNITS)
  //   where,
  //     NAME and UNITS are as above
  //     GRID is the actual grid, typically "m_grid"

  m_grid  = grids_manager->get_grid("Physics");
  m_ncols = m_grid->get_num_local_dofs();
  m_nlevs = m_grid->get_num_vertical_levels();

  constexpr auto Pa = ekat::units::Pa;
  constexpr auto kg = ekat::units::kg;
  constexpr auto m  = ekat::units::m;

  const auto layout_3d = m_grid->get_3d_scalar_layout(true);
  const auto layout_2d = m_grid->get_2d_scalar_layout();

  // we need dp to scale layer content --> strictly input (Required)
  add_field<Required>("pseudo_density", layout_3d, Pa, m_grid->name());

  // we are declaring a tracer to be advected (will be init'ed to zero later, so
  // Updated) but note add_tracer
  add_tracer<Updated>("ash", m_grid, kg / kg);

  // we are interested in creating a "diagnostic" (in the context of this
  // process) ash_column (a la LWP) hence, Computed
  add_field<Computed>("ash_column", layout_2d, kg / (m * m), m_grid->name());
}
```

---

## `initialize_impl` (complete section on next slide)

Before simulation begins and happens only once

- Let's create a mask to control where ash is injected into the system
- By far the most "Kokkos" you will do all day!!!
  - Create a mask for the injection region
  - Use Kokkos to parallelize the loop

---

```cpp
/*-----------------------------------------------------------------------------------------------
 * intialize_impl(run_type)
 *
 * called once for each process at initialization of EAMxx.  This impl can be
 * defined with any actions or functions that are needed at initialization.
 *
 * Inputs:
 *     run_type - an enum which describes the run type.  Initial or Restart
 *
 * can also be empty
 */
void POMPEI::initialize_impl(const RunType /* run_type */) {
  // NOTE: run_type tells us if this is an initial or restarted run,
  // NOTE: run_type tells us if this is an initial or restarted run,
  // Some universal constants that we will need
  using PC               = scream::physics::Constants<Real>;
  constexpr Real deg2rad = PC::Pi / 180;
  constexpr Real r_earth = PC::r_earth;

  // We can create our helper "mask" field. It will be 1 where the volcano
  // injection in the atmosphere happens, and 0 elsewhere.
  const std::string mask_name = "emission_mask";
  constexpr auto nondim       = ekat::units::Units::nondimensional();
  const auto layout_3d        = m_grid->get_3d_scalar_layout(true);
  FieldIdentifier rate_fid(mask_name, layout_3d, nondim, m_grid->name());
  m_emission_mask = Field(rate_fid);
  m_emission_mask.allocate_view();
  m_emission_mask.deep_copy(0);  // 0 means "not injecting here"

  // The grid already stores lat/lon, so go ahead and pre-compute the mask field
  auto lat =
      m_grid->get_geometry_data("lat");  // WATCH OUT: it's in degrees, NOT rad
  auto lon =
      m_grid->get_geometry_data("lon");  // WATCH OUT: it's in degrees, NOT rad

  // Target location of volcanic eruption
  auto volcano_lat  = 40.8214 * deg2rad;
  auto volcano_lon  = 14.4260 * deg2rad;
  auto emission_lev = 27;
  auto radius       = m_params.get<double>("plume_radius_in_km");
  EKAT_REQUIRE_MSG(
      radius > 0,
      "Error! Plume radius must be positive. Input value: " << radius << ".\n");

  // Extract Kokkos (device) views from the fields
  // NOTE: views are "just" multi-dimensional arrays, which are accessible on
  // device
  auto emission_view = m_emission_mask.get_view<Real **>();
  auto lat_view      = lat.get_view<const Real *>();
  auto lon_view      = lon.get_view<const Real *>();

  // -------------------------------------------------------------------------------------------
  // // We define a lambda function for computing the mask that we can then use
  // with Kokkos. A lambda is just a function defined "on-the fly".
  // KOKKOS_LAMBDA simply adds some decoration for GPU execution, nothing to
  // worry about for now
  auto compute_mask = KOKKOS_LAMBDA(const int icol) {
    auto lat_rad   = lat_view(icol) * deg2rad;
    auto lon_rad   = lon_view(icol) * deg2rad;
    auto delta_lat = lat_rad - volcano_lat;
    auto delta_lon = lon_rad - volcano_lon;

    auto dist =
        r_earth / 1e3 * sqrt(delta_lat * delta_lat + delta_lon * delta_lon);

    // If the distance between this point and the center of the volcano is
    // within the radius set the MASK value to 1 (true).
    if(dist < radius) {
      emission_view(icol, emission_lev) =
          1;  // this is the actual calculation!!!!
    }
  };

  // A policy tells kokkos how to parallelize the loop. Here, we are doing
  // a single for loop over the range of indices [0,ncols)
  auto policy =
      Kokkos::RangePolicy<Field::device_t::execution_space>(0, m_ncols);

  // Execute the lambda in parallel according to the execution policy
  Kokkos::parallel_for("compute_ash_mask", policy, compute_mask);
}

```

---

## `run_impl` (complete section on next slide)

where the volcano does its thing

- Let's use the emission function we defined in `pompei.cpp`
- Note the use of built-in "field" functions to simplify this code

---

```cpp
/*-----------------------------------------------------------------------------------------------
 * run_impl(dt)
 *
 * The main run call for the process.  This is where most of the interface to
 * the underlying process takes place.  This impl is called every timestep.
 *
 * Inputs:
 *     dt - the timestep for this run step.
 */
void POMPEI::run_impl(const double dt) {
  // Typically here the developer would retrieve data from fields,
  // For example,
  //   const auto X = get_input_field(NAME_X);
  //   const auto Y = get_output_field(NAME_Y);

  // Another typical step is to issue a call to the underlying process code.
  using PC = scream::physics::Constants<Real>;
  auto g   = PC::gravit;

  auto t    = end_of_step_ts();
  auto rate = pompei::ash_emission_rate(t.days_from(m_eruption_start));
  auto mass = dt * rate;

  auto qash     = get_field_out("ash");
  auto qash_col = get_field_out("ash_column");
  auto rho      = get_field_in("pseudo_density");

  qash.scale(rho);
  qash.scale(1 / g);
  qash.update(m_emission_mask, mass, 1.0);
  qash.scale_inv(rho);
  qash.scale(g);

  vert_contraction<Real>(qash_col, qash, rho);
}

```

---

## `finalize_impl` (complete section on next slide)

anything we need to clean up

- Typically empty
- For demonstration purposes, one can print a statement to the e3sm log

---

```cpp
/*-----------------------------------------------------------------------------------------------
 * finalize_impl()
 *
 * Called at the end of the simulation, handles all finalization of the process.
 *
 * In most cases this is left blank, as EAMxx takes care of most finalization
 * steps. But just in case a process has specific needs the option is available.
 */
void POMPEI::finalize_impl() {
  // Usually blank
}
/*-----------------------------------------------------------------------------------------------*/

}  // namespace scream
```

---

## Hooking everything up

- We have the code, now we need to make sure the compiler sees it.
- We will need to create a CMakeLists.txt file for pompei.

- In the pompei folder create a new file: CMakeLists.txt

---

```cmake
# List cpp sources to be compiled
# In this case we have the pompei model code in pompei.cpp
# and the interface code between pompei and eamxx.
set(POMPEI_SRCS
  eamxx_pompei_process_interface.cpp
  pompei.cpp
)

# Create a library of the pompei code.
set(POMPEI_LIBS "pompei")
add_library(pompei ${POMPEI_SRCS})
target_compile_definitions(pompei PUBLIC EAMXX_HAS_POMPEI)
target_link_libraries(pompei physics_share scream_share)

# Add this library to eamxx_physics
target_link_libraries(eamxx_physics INTERFACE pompei)

```

---

## Adding pompei to parent CMakeLists.txt

Need to make CMake aware of our new file: Open ../CMakeLists.txt and add the following:

```cmake
add_subdirectory(pompei)
```

---

## Register pompei with the physics factory

To `components/src/physics/register_physics.hpp`, add:

```cpp
#ifndef EAMXX_HAS_POMPEI
#include "eamxx/pompei/eamxx_pompei_process_interface.hpp"
#endif

// later in the file

#ifndef EAMXX_HAS_POMPEI
  proc_factory.register_product("pompei", &create_atmosphere_process<POMPEI>);
```

---

## Link pompei via MCT coupling layer

in `components/eamxx/src/mct_coupling/CMakeLists.txt` add:

```cmake
set (SCREAM_LIBS
# later
  pompei
)
```

---

## Add pompei to the namelist

In `components/eamxx/cime_config/namelist_defaults_scream.xml`, add:

```xml
    <!-- Pompei -->
    <pompei inherit="atm_proc_base">
      <eruption_date type="string" doc="Start date/time in YYYY-MM-DD-SSSSS format for the eruption">${RUN_REFDATE}-${RUN_REFTOD}</eruption_date>
      <plume_radius_in_km type="real" doc="Radius of eruption plume, in km" constraints="gt 0; le 100000">100</plume_radius_in_km>
    </pompei>
    <!-- later... -->
    <initial_conditions>
    <!-- ... -->
    <!-- Ash initial condition to 0.0 -->
    <ash               >0.0</ash>
```

---

## Ready to launch

- Just need to update our submission to run with pompei turned on.
- Remember the atmchange commands we discussed in part 1.
- Go back to the runscript we used and search for atmchange.  You’ll see a block of  commented out code which will turn pompei on.  Uncomment this and rerun the script.
- Note, you may need to change the case name or delete the old case to get this going.

---

## The runtime parameters in the run script

```bash
# increase SCREAM_NUM_TRACERS from 10 to 11 (note that scream supports both 128 and 72 levels for low-res)
./xmlchange SCREAM_CMAKE_OPTIONS="SCREAM_NP 4 SCREAM_NUM_VERTICAL_LEV 128 SCREAM_NUM_TRACERS 11"

# add pompei to the list of aerosol processes
./atmchange mac_aero_mic::atm_procs_list+=pompei
# if you want to change the eruption date, you can do so with the following command
./atmchange atmosphere_processes::physics::mac_aero_mic::pompei::eruption_date="0001-01-02-00000"
# if you want to change the eruption radius, you can do so with the following command
./atmchange atmosphere_processes::physics::mac_aero_mic::pompei::plume_radius_in_km=1000.0
```
