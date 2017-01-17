#--- setup Macros file
foreach file (Macros.*)
cat >> $file << EOF1

FFLAGS := \$(FFLAGS) -qrealsize=8 -qdpc=e

EOF1
end

