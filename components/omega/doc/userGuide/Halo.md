(omega-user-halo)=

# Halo Exchanges (Halo)

When Omega is run in a parallel machine environment, the computaional domain
is broken into subdomains that are distrubted among the computational nodes
of the machine. It is necessary to regularly exchange data that is located
at or near the interfaces of adjacent subdomains, these regions are called
halos. The Halo class accomplishes these halo exchanges using the MPI
(Message Passing Interface) standard.

The Halo class depends on the MachEnv class and the Decomp class. Currently,
the only configuration parameters that directly impact Halo are those set
to configure the [Decomp](#omega-user-decomp) class.

Once a Halo object is constructed, all the information needed to perform
a halo exchange on any supported array type defined in any index space
is contained in the Halo object. Halo exchanges are executed via the
exchangeFullArrayHalo function.
