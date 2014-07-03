This directory contains the Ross ice shelf intercomparison experiment data in its raw form.
Retrieved from http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html on May 4, 2009

-----------------------
Original README follows
-----------------------
Oct 2nd 1996,


This directory is part of the package described in the
EISMINT paper: "An Ice-shelf model test based on the Ross Ice Shelf"
published in Annals of Glaciology 23.
The following files contain the Ross Ice Shelf dataset for a
finite difference grid.
Note that only a very few of these data have been used in the
intercomparison work described in the paper.

1- readme.txt (this text file)
-------------

2- 111by147Grid.dat
-------------------
The main part of the dataset for a finite difference grid. It contains
the x and y coordinates of the nodes and 10 interpolated fields of data:
the fields are given in the following order:
(In the x direction the grid has 147 grid points and 111 in the y direction,
this corresponds to a 6822m step)

* a mask (mask) which determines whether or not the grid point is in the
Ross Ice Shelf.

* the observed velocity azimuth (azvel) in degrees

* the velocity magnitude (magvel) in m/year
The velocities (uxbar,uybar) in the grid coordinate system are thus given
by the relationship:
uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))

* Ice thickness (H) in meters

* A mask determining where the interpolation is accurate (velobs=1.) and
where it isn't (velobs=0.)

* the seabed depth (seabed) in meters

* a third mask determining where the ice velocity has to be computed
to satisfy the dynamic boundary condition on the ice front but where
the ice doesn't exist in reality. Indeed, one of the major problem
inherent to the finite difference method itself is to follow correctly the
geometry of the ice front.
One solution (the one I use) is to apply the boundary condition on the edge
of the grid and to compute the velocity where the ice doesn't really exist,
using a 1 meter thickness.
Icefront is equal to 1 where I set the thickness equal to 1 meter, and its value
is 0 everywhere else. It is also possible (it depends on the modeling approach)
that you don't need this mask.

* surface accumulation (accumulation) in mm/year.

* The Bbar coefficient related to the ice rheology in Pas^(1/3). If you
write the Glen's flow law (n=3) the following way:

2*viscosity = (At * tau^2)ª-1

Then At is simply Bbar^-3. These values of Bbar are interpolated from
Thomas's data. Their accuracy is not guaranteed.

* the surface temperature in ³C

3- kbc.dat
----------
One of the two files which specifies the kinematic boundary conditions.
This file contains the index of the grid points for which the interpolated
value of velocity (magnitude and azimuth in 111by147Grid.dat) has to be
taken into account.

4- inlets.dat
-------------
also contains the index of the grid points on the "kinematic boundary". However
for these glaciers, a more accurate value than the interpolated velocity has
been found in the literature. The value of the velocity magnitude and
azimuth is given after the x & y index of the corresponding grid point.

5-6) dataset.f & inputfile.f
----------------------------
Fortran routines which read the different files. They probably cannot be
applied straight away; nevertheless they can give a fair idea of the structure
of the input files.

7- petitfront.dat
-----------------
The Ross Ice Shelf front extends on both sides of the Ross Island. One
part of this ice front is about thirty kilometers whereas the second
is more than five hundred kilometers length. This file gives the index
of the grid points on the "30 kilometers length" ice front.

8- RIGGS.dat:
------------
contains the data at the RIGGS stations which have been used for the
interpolation.

9- visco.ross
-------------
contains values of effective viscosity in MPa.a. These values were
obtained by a control method which is described in "Large-scale
rheology of the Ross Ice Shelf computed by a control method" by
Rommelaere and MacAyeal, Annals of Glaciology 24, to be published.
I can send copies of this paper to anyone interested.
To read this file correctly, do the following loop
do i=1,nx
do j=1,ny
read(10,*) eta(i,j)
end do
end do


For more explanations:
---------------------------------------
Vincent Rommelaere
Laboratoire de glaciologie et
geophysique de l'environnement-France
vince@glaciog.grenet.fr
---------------------------------------


