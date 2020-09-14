# Notes on the maxeler-sombrero interface.

Maxeler provides an implementation of the CG algorithm in terms of
- spinor fields
- gauge fields
- clover term "fields" - which are not currently used in Sombrero.

The interface needs to be able to convert spinor fields and gauge fields
between the two memory layouts. 

## Sombrero Memory layout.

The memory layout for spinor fields and gauge fields in Sombrero can be 
inferred from:
- from the HiRep documentation, the section `Geometry` in `Doc/data_op.tex`,
  where a brief but clear definition of the geometry descriptor structure 
  is available. There is written
  
  > The mapping between the cartesian coordinates of the local lattice 
  > and the array index is given by the macro `ipt(t,x,y,z)=n`.
  
  This seems to be confirmed looking at `LibHR/IO/archive.c` 
  as well, at least for the gauge fields. 
  
  The `ipt` macro, defined in `Include/global.h`, uses the `ipt` 
  lookup table, which is filled in in 
  `LibHR/Geometry/geometry_mpi_eo.c`.
  
-`LibHR/Geometry/geometry_init.c`: some easy-to-understand, basic 
  initialisations.
-`LibHR/Geometry/geometry_mpi_eo.c`: this is where the lookup tables are
  initialised. The logic is cumbersome and the use of unclear variable 
  names,
  In particular it is worrying that the `ipt` lookup table is not used
  directly to set up the `idn` and `iup` lookup tables. Other lookup tables
  are used instead, e.g., `map_overlexi2id` and its inverse, which are 
  filled in using functions that use static variables as counters in a 
  way that requires some thinking to be understood.
-`LibHR/Update/Dphi.c`: To see how the lookup tables, the spinor field and 
  the gauge field really interact.
  From this, we can gather that the indices obtained via the `_FIELD_AT` 
  macro and the `pu_gauge` macros. 
 
### In short

Trusting the documentation, a solution seems to be  to use the `ipt` macro 
in the ways that it is used in the `_FIELD_AT` macro and in the `pu_gauge` macro. 
**It is believed that the links at each site are ordered as TXYZ for 
consistency**.

## MaxCG memory layout. 

Since the mentioning of even/odd is the same everywhere, one can assume that
both the spinor fields and the gauge fields arrays are split in two. 
This is for sure true for the gauge field, for which a function is provided
to interleave even and odd "gauges".
It is also assumed that the ordering is plain lexicographic with a 
XYZT hierarchy.
**It is believed that the links at each site are ordered as XYZT for 
consistency**.


