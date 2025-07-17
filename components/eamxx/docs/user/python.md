# Python support in EAMxx

EAMxx has some limited support for interfacing with external python code.
In particular, we allow calling EAMxx C++ code from a python module, as well
as calling python code from inside an EAMxx atmosphere process. The former
is described in [this page](py2eamxx.md), while the latter is described in
[this page](eamxx2py.md). The two cannot be used at the same time.

NOTE: This feature is currently under development, so details may change in the future.
