Using vasprun from the command line
===================================

Here we describe the basic functionality of PyXtal, with descriptions of the most common options a user will need.

Basic view
----------------------

PyXtal allows the user to generate random crystal structures with given symmetry constraints. The main class for this is the `crystal.random_crystal <pyxtal.crystal.html#pyxtal.crystal.random_crystal>`_ class. There are several parameters which can be specified, but only four are necessary: the symmetry group, the types of atoms, the number of each atom in the primitive cell, and the volume factor. Here is a simple example of a 3D carbon crystal:

.. code-block:: Python

  from pyxtal.crystal import random_crystal
  my_crystal = random_crystal(225, ['C'], [3], 1.0)

This would create a crystal structure with space group 225, 3 carbon atoms in the primitive cell, and a volume factor of 1.0. For stoichiometries with more than one type of atom, replace [‘C’] with a list of atomic symbols, and replace [3] with a list of numbers. For example,

.. code-block:: Python

  my_crystal = random_crystal(99, ['Ba','Ti','O'], [1,1,3], 1.0)

would create a random BaTiO3 crystal.

If the generation is successful, the value ``my_crystal.valid`` will be set to True; otherwise, it will be False. The geometric properties of the crystal are stored in ``my_crystal.struct``, which is a `pymatgen.core.structure.Structure <http://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.Structure>`_ object. You can print most of the useful information using the `print_all <pyxtal.crystal.html#pyxtal.crystal.random_crystal.print_all>`_ function:

.. code-block:: Python

  >>> my_crystal.print_all()
  --Random Crystal--
  Dimension: 3
  Group: 99
  Volume factor: 1.0
  Wyckoff sites:
    Ba: [0.         0.         0.87409062] 1a, site symmetry 4 m m
    Ti: [0.         0.         0.36696892] 1a, site symmetry 4 m m
    O: [0.         0.         0.83119703] 1a, site symmetry 4 m m
    O: [0.         0.         0.19427634] 1a, site symmetry 4 m m
    O: [0.5        0.5        0.65656028] 1b, site symmetry 4 m m
  Pymatgen Structure:
  Full Formula (Ba1 Ti1 O3)
  Reduced Formula: BaTiO3
  abc   :   5.867513   5.867513   2.595120
  angles:  90.000000  90.000000  90.000000
  Sites (5)
    #  SP       a     b         c
  ---  ----  ----  ----  --------
    0  Ba    -0    -0    0.874091
    1  Ti     0     0    0.366969
    2  O     -0    -0    0.831197
    3  O     -0    -0    0.194276
    4  O      0.5   0.5  0.65656
  
You can also store the structure to a file for later usage: 
 
.. code-block:: Python

  my_crystal.to_file(filename)

By default, this will create a cif file storing the structure information. Other file types are supported by specifying the value fmt as ‘poscar’, ‘cssr’, or ‘json’.

2D Crystals
~~~~~~~~~~~

PyXtal can also generate subperiodic crystals. To generate a 2d crystal, use the class `crystal.random_crystal_2D <pyxtal.crystal.html#pyxtal.crystal.random_crystal_2D>`_. For example,

.. code-block:: Python

  from pyxtal.crystal import random_crystal_2D
  my_crystal = random_crystal_2D(20, ['C'], [4], 1.0, thickness=2.0)

would generate a 2d crystal with layer group 20, 4 carbon atoms in the primitive cell, a volume factor of 1.0, and a thickness of 2.0 Angstroms. As with the 3d case, for crystals with multiple atom types, you may replace [‘C’] and [4] with lists of the atomic symbols and amounts, respectively. The crystal will be periodic in two directions instead of 3. PyXtal adds 10 Angstroms of vacuum on each side of the 2D lattice, so that optimization may be performed without altering the structure file. However, care should be taken when using the cif file for applications designed for 3D crystals. The axis of non-periodicity can be accessed via my_crystal.PBC; each axis will either be 1 or 0, representing either periodicity or non-periodicity. For example, PBC = [1,1,0] means that the x and y axes are periodic, while the z axis is non-periodic.

Note that the layer group number is different from the international space group number, and ranges between 1 and 80. For a list of the layer groups and their symmetry operations, see `the International Tables of Crystallography, Volume E, part 4 <https://it.iucr.org/Eb/ch4o1v0001/contents/>`_.

By default, PyXtal will automatically generate a value for the thickness of the unit cell, based on the volume. By specifying a value for thickness, you override this behavior. So, if you are testing over a range of volume factors, consider how the shape of the unit cell will be affected, and change the thickness accordingly. Alternatively, you may supply a custom Lattice object, as described below.

