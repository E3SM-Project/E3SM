# wrapper script to build cprnc on bluefire
#!/bin/bash

. /contrib/Modules/${MODULE_VERSION}/init/bash

module load netcdf/4.1.3_seq
gmake
