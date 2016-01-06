CPAN modules
============

The perl5lib/CPAN directory contains **pure perl** modules from
CPAN. It should **NOT** contain modules that depend on compiled c
code. Editing or modifying these modules runs the risk that changes
will be overwritten if the module is updated from CPAN.

If you add a new module, please update the list with the version,
reason for the addition and url.

Included modules:

 * Test-Class-0.50 - Xunit style test class that provides startup and shutdown methods
 * Test-Exception-0.40 - allow testing for exceptions, e.g. tests that should die
 * MRO-Compat-0.12 - dependency
 * Sub-Uplevel-0.25 - dependency
 * Try-Tiny-0.22 - dependency
 * Module-Runtime-0.014 - dependency

All modules were licensed under the same license as perl (Artistic
License or GNU GPL). As of 2015-06-23 see:

http://cpansearch.perl.org/src/ETHER/Test-Class-0.50/lib/Test/Class.pm
http://cpansearch.perl.org/src/EXODIST/Test-Exception-0.40/lib/Test/Exception.pm
http://cpansearch.perl.org/src/DAGOLDEN/Sub-Uplevel-0.25/lib/Sub/Uplevel.pm
http://cpansearch.perl.org/src/BOBTFISH/MRO-Compat-0.12/lib/MRO/Compat.pm
http://cpansearch.perl.org/src/DOY/Try-Tiny-0.22/lib/Try/Tiny.pm
http://cpansearch.perl.org/src/ZEFRAM/Module-Runtime-0.014/lib/Module/Runtime.pm

add the 3 modules for File::Find::Rule

http://search.cpan.org/~rclamp/Text-Glob-0.09/lib/Text/Glob.pm
http://search.cpan.org/~rclamp/Number-Compare-0.03/lib/Number/Compare.pm
http://search.cpan.org/~rclamp/File-Find-Rule-0.34/lib/File/Find/Rule.pm
