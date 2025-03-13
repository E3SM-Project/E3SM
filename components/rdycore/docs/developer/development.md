# RDycore Development Process

Here we discuss the practices and tools used to develop RDycore. The RDycore
team is a small interdisciplinary team, and our development process is tailored
to its particular needs. These needs necessarily evolve over the lifecycle of
the project, so this process is likely to change from time to time.

_There is no single process that works best in every situation._ This is as true
for software development as it is for the approximations used in basic science
and mathematical modeling, and for rules of engagement in all human activites.
The working styles, expertise, and needs of individual team members determine
the shape of the process that makes them most productive.

## Git Code Repository and Workflow

Like all E3SM "ecosystem" projects, RDycore stores its source code
in a [GitHub repository](https://www.github.com/RDycore/RDycore). The
[RDycore GitHub Organization](https://www.github.com/RDycore) contains this
and other repositories related to the project.

We develop RDycore using the [Git feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow),
which is among the more popular Git-based development methodologies. In this
workflow:

* the [main branch](https://github.com/RDycore/RDycore/tree/main) contains the
  latest working code that has passed all tests and quality controls
* all work on any new feature or bugfix is performed by one or more team members
  in a "feature branch" created from this main branch
* when work in a feature branch is finished and ready for incorporation into
  the main branch, one of the participating team members creates a pull request
  which, upon successful review, is merged to the `main` branch
    * each pull request to be merged with the `main` branch triggers a set of
      automated tests that build RDycore, run its test suite, and perform some
      basic quality checks (code formatting, test coverage, etc)
    * one or more members of the RDycore team must be assigned to review the
      code changes in the feature branch before it can be merged
    * a code reviewer can ask questions about, request changes to, or approve
      the code in the pull request
    * at least one approval from a reviewer is required to merge the pull
      request (which may or may not require conflicts to be resolved between
      the feature branch and the `main` branch)

### Keep it short and simple!

Because long-lived feature branches are correlated with increasing numbers of
merge conflicts and other inconsistent states within the repository, it is
recommended that feature branches are **short-lived**, **tightly-scoped**, and
merged in a timely fashion. Such branches are also much easier to review for
errors and inconsistencies.

Occasionally, it may be necessary to create a long-lived feature branch. Such a
branch must be managed actively and carefully to accommodate the complexity it
introduces to the development process. Reviewing a long-lived feature branch is
a difficult and time-intensive process, so make sure you know what you're
getting into if you think you need one.

## Third-Party Libraries

RDycore relies on [PETSc](http://petsc.org/release) for much of its
functionality, including most of its third-party libraries such as NetCDF,
HDF5, and Exodus. You can find details for installing the right version of
PETSc with all the necessary third-party libraries in the [installation guide](../common/installation.md).

Additionally, RDycore uses a few other libraries not available from PETSc:

* [libcyaml](https://github.com/tlsa/libcyaml) - A schema-based YAML parser we
  use to parse RDycore's flexible and expressive [input format](../common/input.md)
* [libyaml](https://pyyaml.org/wiki/LibYAML) - the original C YAML parser, which
  is used by `libcyaml` above
* [cmocka](https://cmocka.org/) - a C unit testing framework that we use for
  testing some low-level features

These libraries are imported into the RDycore repo as [Git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules).
Git submodules are like "repos within a repo", which is a half-baked concept
that makes them tricky to use at times, but you won't spend much time thinking
about them unless one of the above libraries is updated. Usually, all you have
to do is type

```
git submodule update --init --recursive
```

from the top level of the RDycore source tree, which ensures that the submodules
in your local Git workspace are consistent with the branch you're working in.

## CMake Build System

Like all E3SM-related projects, RDycore uses [CMake](https://cmake.org/) as a
configuration/build system. CMake is the most common build system for C and C++
projects. It can be very complicated if all of its features are used
indiscriminately, and it's difficult even to understand all of its features.
The RDycore project has tried to confine its CMake usage to the bare essentials,
in the interest of creating a simple and reproducible build and development
environment.

Here we describe the basic structure of our CMake setup. For instructions on
how to build RDycore using this build system, refer to the [installation guide](../common/installation.md).

The build system is composed of a set of `CMakeLists.txt` files that define
various parts of RDycore:

* The `CMakeLists.txt` file in the top-level source directory defines the
  project and configures PETSc and compilers and refers to other directories
  containing their own `CMakeLists.txt` files.
* The `external` directory has its own `CMakeLists.txt` file that builds the
  third-party libraries that aren't available as part of PETSc's distribution.
* The `src` directory contains all the source files, and the `CMakeLists.txt`
  therein builds the RDycore C library. A `tests` subdirectory has a `CMakeLists.txt`
  file that defines unit tests for the library, and the `f90-mod` subdirectory
  has one that builds the Fortran library.
* The `driver` directory contains source files for the C and Fortran standalone
  driver programs and other developer tools, and the `CMakeLists.txt` file
  defines their build configurations. A `tests` subdirectory defines several
  tests for these drivers and tools.
* Various other directories have `CMakeLists.txt` files that perform tasks
  related to their content.

### CMake targets and dependencies

Like Make, CMake defines "targets" to build and expresses dependenceis between
these targets to determine the order in which everything is built. For example,
the RDycore standalone driver (represented by the `rdycore_exe` target defined
in [driver/CMakeLists.txt](https://github.com/RDycore/RDycore/blob/main/driver/CMakeLists.txt))
depends upon the RDycore C library (the `rdycore` target in [src/CMakeLists.txt](https://github.com/RDycore/RDycore/blob/main/src/CMakeLists.txt)),
so the `rydcore` target must be built before the `rdycore_exe` target.

If you're curious about how CMake works in detail, Kitware has a [tutorial](https://cmake.org/cmake/help/latest/guide/tutorial/index.html)
on their CMake website.

### "Modern CMake"

You might hear a CMake enthusiast refer to "Modern CMake". This refers to the
practice of setting build properties on specific targets (with commands like
[target_include_directories](https://cmake.org/cmake/help/latest/command/target_include_directories.html)
and [target_link_libraries](https://cmake.org/cmake/help/latest/command/target_link_libraries.html))
instead of setting these properties globally (with [include_directories](https://cmake.org/cmake/help/latest/command/include_directories.html)
and [link_libraries](https://cmake.org/cmake/help/latest/command/link_libraries.html)).

Using target-specific properties makes it easier to find issues with the build
system. We adhere as much as possible to this practice in the development of
RDycore.

## CTest Automated Testing

The `CMakeLists.txt` files within the test-related directories ([src/tests](https://github.com/RDycore/RDycore/tree/main/src/tests)
and [driver/tests](https://github.com/RDycore/RDycore/tree/main/driver/tests)
define tests that are run when you type `make test`. These tests use [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html),
which is an automated testing system built into CMake. The link in the previous
sentence is a good resource for learning about this testing system, but you can
also just run `make test` after you've successfully built RDycore to see it in
action. The results of this command (testing logs, lists of failed tests, etc)
can be found in the `Testing/Temporary` folder within your build directory.

## GitHub Continuous Integration Environment

Previously, we mentioned that any pull request for a feature branch to be merged
to the `main` branch triggers an automated set of tests and quality checks.
These tests and checks, which we refer to as RDycore's "continuous integration
(CI) environment", are implemented using [GitHub Actions](https://docs.github.com/en/actions).

"Github Actions" are sets of logic that does work on things in repositories.
These Actions can be assembled into "GitHub Actions workflows" that perform tasks such as
building RDycore or running unit tests. GitHub Actions workflows are triggered
by events in a GitHub repository such as pull requests, merges, the creation of
issues, etc. Each workflow can pass or fail, making it convenient to use them
as criteria for accepting new code contributions.

Each GitHub Actions workflow is defined by a YAML file placed in the [.github/workflows](https://github.com/RDycore/RDycore/tree/main/.github/workflows)
directory. The following workflows are defined for RDycore and are triggered by
a pull request to the `main` branch:

* `auto_test`: builds RDycore in a container in which PETSc is installed and
  runs all of its unit tests, generating a code coverage reports. This workflow
  fails if any of these operations cannot be completed successfully.
* `clang-format-check`: checks C source code formatting on selected files,
  failing if the code has not been formatted correctly. For information on
  code formatting, see our [style guide](style.md).
* `gh-pages`: publishes updated documentation generated by [mkdocs](https://www.mkdocs.org/)
  from Markdown files in the `docs` directory. The documentation is published to
  RDycore's [GitHub Page](https://rdycore.github.io/RDycore/), and a preview is
  generated on the pull request page for inspection prior to a merge.

For details on the syntax of the YAML files used to define these workflows, see
the [GitHub Actions documentation](https://docs.github.com/en/actions/using-workflows/workflow-syntax-for-github-actions).
