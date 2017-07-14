.. _adding-components:

===================
Adding components
===================

Here are the steps to add prognostic components to CIME models.  

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
build directory.
