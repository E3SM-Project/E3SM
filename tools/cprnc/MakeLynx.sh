#wrapper script to build cprnc on lynx
#!/bin/bash

. /opt/modules/${MODULE_VERSION}/init/bash

module load PGI/netcdf4/4.1.3_seq

gmake NETCDF=/contrib/netcdf/4.1.3/netcdf4-pgi
