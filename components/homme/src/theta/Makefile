FC = mpif90

../share/physical_constants.o: ../share/physical_constants.F90
	$(FC) -c -Wall ../share/physical_constants.F90 -o ../share/physical_constants.o

../share/kinds.o: ../share/kinds.F90
	$(FC) -c -Wall ../share/kinds.F90 -o ../share/kinds.o

../share/dimensions_mod.o: ../share/dimensions_mod.F90
	$(FC) -c -Wall -DNP=4 -DPLEV=26 ../share/dimensions_mod.F90 -o ../share/dimensions_mod.o

../share/parallel_mod.o: ../share/parallel_mod.F90 ../share/dimensions_mod.o
	$(FC) -c -Wall -D_MPI ../share/parallel_mod.F90 -o ../share/parallel_mod.o

../share/coordinate_systems_mod.o: ../share/coordinate_systems_mod.F90 \
	../share/physical_constants.o ../share/parallel_mod.o
	$(FC) -c -Wall ../share/coordinate_systems_mod.F90 -o ../share/coordinate_systems_mod.o

../share/edgetype_mod.o: ../share/edgetype_mod.F90
	$(FC) -c -Wall ../share/edgetype_mod.F90 -o ../share/edgetype_mod.o

../share/control_mod.o: ../share/control_mod.F90
	$(FC) -c -Wall ../share/control_mod.F90 -o ../share/control_mod.o

../share/params_mod.o: ../share/params_mod.F90
	$(FC) -c -Wall ../share/params_mod.F90 -o ../share/params_mod.o

../share/gridgraph_mod.o: ../share/gridgraph_mod.F90 ../share/control_mod.o \
	../share/params_mod.o
	$(FC) -c -Wall ../share/gridgraph_mod.F90 -o ../share/gridgraph_mod.o

element_state.o: element_state.F90 ../share/kinds.o
	$(FC) -c -Wall element_state.F90 -o element_state.o

../share/element_mod.o: ../share/element_mod.F90 ../share/kinds.o \
	../share/coordinate_systems_mod.o ../share/edgetype_mod.o \
	../share/gridgraph_mod.o element_state.o
	$(FC) -c -Wall ../share/element_mod.F90 -o ../share/element_mod.o

../share/hybvcoord_mod.o: ../share/hybvcoord_mod.F90
	$(FC) -c -Wall ../share/hybvcoord_mod.F90 -o ../share/hybvcoord_mod.o

../share/hybrid_mod.o: ../share/hybrid_mod.F90
	$(FC) -c -Wall ../share/hybrid_mod.F90 -o ../share/hybrid_mod.o

../share/quadrature_mod.o: ../share/quadrature_mod.F90
	$(FC) -c -Wall ../share/quadrature_mod.F90 -o ../share/quadrature_mod.o

../share/derivative_mod_base.o: ../share/derivative_mod_base.F90 ../share/quadrature_mod.o
	$(FC) -c -Wall ../share/derivative_mod_base.F90 -o ../share/derivative_mod_base.o

derivative_mod.o: derivative_mod.F90 ../share/derivative_mod_base.o
	$(FC) -c -Wall -ffree-line-length-none derivative_mod.F90 -o derivative_mod.o

homme_nvector.o: homme_nvector.F90 ../share/element_mod.o ../share/hybvcoord_mod.o \
	../share/hybrid_mod.o derivative_mod.o
	$(FC) -c -Wall homme_nvector.F90 -o homme_nvector.o

homme_nvector_driver: homme_nvector_driver.F90 homme_nvector.o
	$(FC) -c -Wall homme_nvector_driver.F90 -o homme_nvector_driver.o
	$(FC) homme_nvector_driver.o homme_nvector.o -o homme_nvector_driver
