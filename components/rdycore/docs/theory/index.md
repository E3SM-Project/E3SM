# RDycore Theory Guide - Overview

This guide describes the mathematical and numerical approximations used by
RDycore to model compound flooding.

## Notation

### Two types of vectors

We write _multi-component quantities_ in $\mathbf{bold}$, sometimes referring to
them as vectors. For example, we refer this way to the _solution vector_
$\mathbf{U} = [h, hu, hv]^T$. These multi-component quantities are _vectors_ in
the sense that they are used in matrix-vector expressions suitable for solving
systems of equations. By convention, we write them as column vectors.

By contrast, we write _spatial vectors_ (those with $x$ and $y$ components
aligned respectively with corresponding spatial axes) with arrows over their
symbols, like the flow velocity $\vec{u} = (u_x, u_y) = (u, v)$. These are the
vectors familiar to physical scientists. We write these as row vectors.

Sometimes a multi-component quantity itself can be considered a spatial vector,
such as the _flux vector_ $\mathbf{\vec{F}} = (\mathbf{F}_x, \mathbf{F}_y),$
where $\mathbf{F}_x$ and $\mathbf{F}_y$ are themselves multi-component
quantities. We write these "multi-component vector quantities" in bold with
arrows overhead.

Distinguishing between these "types" of vectors allows us to make use of
concepts from vector calculus such as _divergence_ ($\vec{\nabla}\cdot\mathbf{\vec{F}}$)
and _projections_ along vectors normal to surfaces ($\mathbf{\vec{F}}\cdot\vec{n}$).
This distinction is often not made in numerical analysis, which can confuse the
reader who is trying to unpack a complicated expression.

### Geometry

We use some elementary ideas from [set theory](https://en.wikipedia.org/wiki/Set_theory)
to describe the geometry in our model formulation. Specifically:

* We deal exclusively with sets containing points in the plane. Any two sets
$x$ and $y$ are _equal_ $(x = y)$ if they contain exactly the same points.
* $\varnothing$ is the _null set_, which contains no points.
* We represent the 2D _cells_ used by RDycore as [closed sets](https://mathworld.wolfram.com/ClosedSet.html) in the plane. We typically write the cell $i$ as $\Omega_i$, the set of all points contained within that cell, _including those on its boundary_.
* The boundary of a 2D cell consists of the _faces_ (alternatively _edges_) that bound the cell. We write the boundary of the cell $\Omega_i$ as $\partial\Omega_i$, the set of all points in each (piecewise linear) face attached to the cell, _including its endpoints_.
* The vertices of a 2D cell are the isolated points shared by adjacent faces in that cell. A triangular cell has 3 vertices, while a quadrilateral cell has 4.
* $p \in x$ indicates that the point $p$ belongs to the set $x$.
* $x \subset y$ indicates that the set $x$ is a [subset](https://en.wikipedia.org/wiki/Subset) of the set $y$. $x \not\subset y$ indicates that $x$ is _not_ a subset of $y$.
* $x \bigcup y$ indicates the [union](https://en.wikipedia.org/wiki/Union_(set_theory)) of two sets $x$ and $y$.
* $\bigcup_i x_i$ indicates the [union](https://en.wikipedia.org/wiki/Union_(set_theory)) of a number of indexed sets $x_i$.
* $x \bigcap y$ indicates the [intersection](https://en.wikipedia.org/wiki/Intersection_(set_theory)) of two sets $x$ and $y$.
* $\bigcap_i x_i$ indicates the [intersection](https://en.wikipedia.org/wiki/Intersection_(set_theory)) of a number of indexed sets $x_i$.

Cells, faces, and vertices can be easily expressed in terms of these basic ideas:

* The intersection of two cells $i$ and $j$ (written $\Omega_i \bigcap \Omega_j$),
is the face separating those two cells, consisting of all points contained in
both cells.
* We sometimes refer to a face $f \in \partial\Omega_i$, which reads "face $f$,
which belongs to the boundary of cell $i$."
* The intersection of two adjacent faces in a cell is a vertex of that cell,
consisting of the single point common to those faces.

The union of all cells $\Omega_i$ in a grid is the _computational domain_ $\Omega$
over which the grid is defined. The boundary of this domain, which can be written
$\partial\Omega$ but which we also sometimes write as $\Gamma$ for clarity,
is the set of all faces attached to only one cell.
