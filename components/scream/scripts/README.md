# The SCREAM toolsuite

This directory contains a python3-based set of tools for SCREAM that compliments the CIME tool
suite and also works with standalone SCREAM. This document will list our tools, roughly in order
of importance, along with a brief description. For more info on a tool, run it with `--help`. We will
always strive to have an informative help dump that contains detailed info on a tool's purpose, args/options,
and some example usages.

## test-all-scream

test-all-scream is our core testing script for standalone SCREAM. When developing on a branch, if this script
runs successfully with the default options, then your branch can be considered to work on the current machine.
When a PR is issued, our continuous integration will run this script (via gather-all-data (see below)) on your
branch on all our core testing machines. If all those tests pass, then your PR will be eligible for merging.

Because this tool is the foundation of our CI testing, you should always expect this tool to work correctly
and have up-to-date documentation.

## gather-all-data

A tool for dispatching jobs to machines, loading the SCREAM env, and doing batch submissions. This tool is
usually used to get test-all-scream jobs running on compute nodes on machines, but it can be used to run
anything you want. Our CI system uses this tool to run test-all-scream. We have used this tool in the past
to relatively easily gather performance data for SCREAM across all the platforms we care about with a single command.

Because this tool is used in our CI testing, you should always expect this tool to work correctly
and have up-to-date documentation.

## scripts-tests

A test suite for this toolsuite. This should be used for any significant developments to the core testing scritpts
(test-all-scream or gather-all-data). This suite is NOT run by our CI system, so it's up to our toolsuite developers
to remember to run this.

Because this testsuite is used in developing our core test scripts, you should always expect this testsuite to pass
and have up-to-date documentation.

## scream-env-cmd

A scream analog for modulecmd, a tool for dumping shell commands that would load the SCREAM approved env for
your machine. Designed for easy integration with bash functions for easy loading of env.

This tool is not used in our core testing tools, but it is simple and actively used by some of our developers
when they work on various machines. It leverages the machines_specs.py library which is essential for the entire
toolsuite and actively maintained, so you should expect this to tool work and load the correct env.

## perf-analysis

The tool we use to gather performance data for SCREAM. It supports scaling experiments, run repitition to reduce
noise, and plot-friendly output among other features. The raw data for many of the plots used in SCREAM presentatons
was generated with this tool.

This tool is not used in our core testing tools and we only occassionally do major perfomance analysis on our code,
so this tool can fall into minor disrepair.

## plot

The tool we use to plot performance data from perf-analysis. It
supports a variety of options for color selection, line formatting
customizations (bold, dashed, etc), font size, etc.

This tool is not used in our core testing tools and we only occassionally do major perfomance analysis on our code,
so this tool can fall into minor disrepair.

## const-maker

A tool for generating CXX and F90 source code for constant definitions where the constant is derived from
a mathematical forumla. This is necessary because cxx constexprs don't always allow for math calls in the
definition of the constant.

This tool is not used in our core testing tools, but it is extremely simple and not coupled to anything else in the repo,
so it should be expected to work.

## gen-boiler

A tool for parsing fortran source files and generating boilerplate for bridging, testing, etc. This is intended to be
used in the situation where we are porting a fortran reference code to CXX like we did for P3 and SHOC.

This tool is not used in our core testing tools and goes through long periods of time where it does not get used much
(IE there are no active ongoing porting efforts). It is likely this tool will not work exactly as expected if it
has not been run in a while or is being used on a package on which it has not been used before.

## cf-xml-to-yaml

Given an XML file containing the CF conventions for standardized field names
(https://cfconventions.org/standard-names.html), this tool generates a YAML
file with the same information.

This tool is not used in our core testing tools, but it is extremely simple and not coupled to anything else in the repo,
so it should be expected to work.
