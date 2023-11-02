## Resizing slabs

Context on the `resize` parameter:

Pymatgen's slab generation algorithm works by replicating the
oriented unit cell(OUC) a number of times to each the minimum slab size.
However, it does not determine the number of layers in the OUC correctly.
This leads to much larger slabs than one asks for, but it ensures that
the terminations of the top and bottom are always complementary, which
can be useful at times.

This tag determines whether to remove layers from the bottom until
the desired thickness is reached. It also modifies the vacuum by
modifying the c parameter so that we have the desired vacuum.
Defaults to 'False'.

Context on the `preserve_terminations` parameter:

Pymatgen sometimes generates more than one unique slab for a given miller
index, meaning there are distinct terminations for each one. This tag
determines whether to preserve the terminations when resizing slabs.
This is done accomplished by removing layers in 'chunks' where each chunk
corresponds to a portion of the oriented unit cell with height equal to
the distance between miller planes for the given conventional bulk structure.
Defaults to 'None'.

Context on the `min_thick_A` parameter:

For certain structures, the layers are very close to each other (depending
on the tol parameter of course), where one could have a slab with many layers
but still thin in the physical space. This could lead to convergence issues
in surface energies. This tag allows the user to set a minimum slab thickness
in units of Angstroms, and the resizing algorithm will not remove any more layers
that will reduce the slab thickness below min_thick_A.
Defaults to 'None'.

For more info about the description of the parameters concerning pymatgen's `SlabGenerator`,
please refer to the documentation of `pymatgen.core.surface.SlabGenerator`.