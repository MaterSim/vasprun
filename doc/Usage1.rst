Using vasprun from the command line
===================================
After installation, the executable called ``vasprun`` will be automatically included to the path. This executable provides several convenient shortcut to various functions.

Overview
--------
``$ vasprun -h``

Options:
  -h, --help            show this help message and exit
  -i, --incar           export incar file
  -p poscar_file, --poscar=poscar_file 
                        export poscar file
  -c cif_file, --cif=cif_file
                        export symmetrized cif
  -k, --kpoints         kpoints file
  -d dos_plot, --dosplot=dos_plot
                        export dos plot, options: t, spd, a, a-Si, a-1
  -b band_plot, --bandplot=band_plot
                        export band plot, options: normal or projected
  -v vasprun, --vasprun=vasprun
                        path of vasprun.xml file, default: vasprun.xml
  -f, --showforce       show forces, default: no
  -a, --allparameters   show all parameters
  -e, --eigenvalues     show eigenvalues in valence/conduction band
  -s smearing, --smear=smearing
                        smearing parameter for dos plot, e.g., 0.1 A
  -n figname, --figname=figname
                        dos/band figure name, default: fig.png
  -l lim, --lim=lim     dos/band plot lim, default: -3, 3
  -m max, --max=max     band plot colorbar, default: 0.5


1, A quick run
-----------

``$ vasprun -v vasprun.xml -f yes``

::

  formula :   Sc4C4
  efermi :   5.71888444
  energy :   -65.82108919
  metal :   True
  gap :   -0.4621
  +----+--------------+----------+-----------+
  |    | functional   | labels   |   valence |
  |----+--------------+----------+-----------|
  |  0 | PBE          | Sc_sv    |        11 |
  |  1 | PBE          | C        |         4 |
  +----+--------------+----------+-----------+
  +----+----------------------+-------------------------+
  |    | lattice              | stress                  |
  |----+----------------------+-------------------------|
  |  0 | [3.478641, 0.0, 0.0] | [-4.62998221, 0.0, 0.0] |
  |  1 | [0.0, 4.682565, 0.0] | [0.0, 1.73781963, 0.0]  |
  |  2 | [0.0, 0.0, 6.176701] | [0.0, 0.0, 5.81088683]  |
  +----+----------------------+-------------------------+
  +----+--------------------------+---------------------------------+
  |    | atom                     | force                           |
  |----+--------------------------+---------------------------------|
  |  0 | [0.0, 0.5, 0.606519]     | [0.0, 0.0, 0.02408274]          |
  |  1 | [0.0, 0.0, 0.958081]     | [0.0, 0.0, 0.07723516]          |
  |  2 | [0.5, 0.5, 0.041919]     | [0.0, 0.0, -0.07723516]         |
  |  3 | [0.5, 0.0, 0.393481]     | [0.0, 0.0, -0.02408274]         |
  |  4 | [0.5, 0.83775, 0.744469] | [0.0, 0.0691885, -0.00447732]   |
  |  5 | [0.5, 0.16225, 0.744469] | [-0.0, -0.0691885, -0.00447732] |
  |  6 | [0.0, 0.33775, 0.255531] | [0.0, 0.0691885, 0.00447732]    |
  |  7 | [0.0, 0.66225, 0.255531] | [0.0, -0.0691885, 0.00447732]   |
  +----+--------------------------+---------------------------------+

2, Plot DOS
--------
A number of ways of plotting dos are also supported, some basic options are

- t: total dos
- spd: spd dos

one can just use the -d option to customize the plots. If spin is included in vasprun.xml, the plot will show both up and down spin states separately.

``$ vasprun -v vasprun.xml -f yes``

3, Plot Band structure
-------------------
3.1 band plot with customized energy range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``$ vasprun -v vasprun.xml-band -b normal -l -3,3 -m 0.4 -n band.png``

3.2 colored band based on the occupation of projected DOS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


4, IR intensity analysis
------------------------

5, Look up the eigenvalues by band index
----------------------------------------

