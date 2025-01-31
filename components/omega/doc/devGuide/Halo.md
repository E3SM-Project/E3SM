(omega-dev-halo)=

# Halo Exchanges (Halo)

Omega utilizes domain-based parallelism, where the domain is divided into
partitions that run in parallel as separate tasks distributed across the
resources of the machine. Arrays representing physical quantities that are
defined throughout the domain are divided into overlapping chunks and
distributed among the tasks. As the model evolves, the values stored in these
arrays are updated. Array elements that are needed by a task to evolve the
model forward in time, but are defined on and updated by a neighboring task
are halo elements. It is necessary for each parallel task to regularly send
locally defined array elements that belong to the halos of neighboring tasks
and to receive local halo elements owned by neighboring tasks. In Omega, these
halo exchanges are perfromed using the Halo class which, once constructed,
contains lists of indices to send to and receive from neighboring tasks, and
methods to conduct the exchanges. The halo exchanges are carried out using
non-blocking MPI (Message Passing Interface) communication.

An instance of the Halo class is dependent on a parallel domain decomposition
and the associated parallel machine environment, so those objects need to be
constructed first. The default Halo is created with the call:
```c++
OMEGA::Halo::init();
```
This must be done after the default MachEnv and Decomp have been initialized.
Once initialized, a pointer to the default Halo can be retrieved with:
```c++
OMEGA::Halo *DefHalo = OMEGA::Halo::getDefault();
```
Additional Halo objects can be constructed for defined MachEnv and Decomp
objects by calling the Halo create method with a supplied `std::string` Name:
```c++
OMEGA::Halo::create(Name, NewEnv, NewDecomp);
```
A pointer to each constructed Halo is stored in a `std::map` container, which
can be retrieved by supplying the Name to `OMEGA::Halo::get(Name)`.

The Halo constructor will save the MPI communicator handle and the MPI task ID
of the local task from the MachEnv object, and will extract and save from the
Decomp object indices defining the local and neighboring halo elements for the
cell, edge, and vertex index spaces. This information is organized into objects
of two nested classes of the Halo class: the ExchList class and the Neighbor
class. The ExchList class contains the indices of Halo elements to either send
to or receive from a single neighboring task for a particular index space, and
the Neighbor class contains all the ExchList objects needed to carry out a halo
exchange with a single neighboring task in any index space, as well as the
buffer memory space to communicate with that neighbor. All member variables
and objects and most member methods of the Halo class are declared private
as they are only needed by the Halo class methods to execute an exchange.

The main methods of the Halo class which execute an exchange are
  - startReceives: wrapper to call MPI_Irecv for each neighbor
  - packBuffer: packs halo elements into a buffer to send to a neighbor
  - startSend: wrapper to call MPI_Isend to send a buffer to each neighbor
  - unpackBuffer: unpacks halo elements from a received buffer into an array

The packBuffer and unpackBuffer functions are overloaded to support different
array types.

The exchangeFullArrayHalo function is the main interface to conduct a halo
exchange for any supported array defined in the DataTypes module. The
exchangeFullArrayHalo function is a template which passes the array to the
Halo methods above which perform the exchange. The index space on which
the array is defined must also be passed to this function. This is currently
done using an enum defined in the Halo header file.
```c++
enum MeshElement{OnCell, OnEdge, OnVertex};
```
So, the exchangeFullArrayHalo function is called as follows:
```c++
MyHalo.exchangeFullArrayHalo(SomeCellBasedArray, OMEGA::OnCell);
```
for a constructed Halo named MyHalo and any supported array type in the cell
index space.
