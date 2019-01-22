.. _adding-components:

===================
Adding components
===================

Here are the steps to add prognostic components to CIME models.  

There are a couple of aspects of a component interface to CIME, the
scripts interface which controls setting up component inputs and
building the component and the run interface which controls connecting
the component to the coupler and through the coupler, the other
components of the CIME based model.

The component should have a subdirectory **cime_config** and this
subdirectory should have two files **buildnml** and **buildlib** The
**buildnml** script is used to build the components instructional,
runtime inputs.  These have traditionally been in the form of fortran
namelists but may also follow other formats.  The **buildnml** may
either be called from the command line or as a python subroutine.  If
buildnml is called from the command line it will be passed the
caseroot directory on the command line.  If it is called as a
subroutine, the subroutine name must be buildnml and it will take
three arguments, a Case object, a caseroot directory and a component
name.  The **buildlib** script will always be called from the command
line, it is called in the case.build step and is expected to build the
the buildlib script will be called with three arguments in order they
are caseroot, libroot (the location of the installed library,
typically EXEROOT/lib) and bldroot, the location of the component
build directory.  Look at the cime internal components such as datm
for an example.

The coupler interface is dependent on which coupler is used, for the mct coupler in cime
the component model must provide NNN_INIT_MCT, NNN_RUN_MCT, NNN_FINAL_MCT where NNN is the
component type of the particular component (eg ATM for an atmosphere, LND for a land model)
these subroutines are expected to be in the component library.
