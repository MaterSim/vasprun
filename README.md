# vasprun
This repo used for quick analysis of vasp calculation solely from vasprun.xml. It has the following features:

- band gap calculation
- dos plot (total dos, needs to add partial dos)
- incar/potcar/poscar generation
- force analysys
- eigenvalue analysys (needs to add more)

## Prerequisites
To use it, one must have python 3 installed. It is recommended to install conda if one don't have root access in his computer.
In addition, several packages will be required
- pymatgen
- pandans
- lxml

One could follow the [wikipage](https://github.com/qzhu2017/CMS/wiki/Python-environment-setup) to set up your python environment.
```
$ python vasprun.py -h
```
## Usage
```
Usage: vasprun.py [options]

Options:
  -h, --help            show this help message and exit
  -i, --incar           export incar file
  -p poscar file, --poscar=poscar file
                        export poscar file
  -c cif file, --cif=cif file
                        export symmetrized cif
  -k, --kpoints         kpoints file
  -d dos_plot, --dosplot=dos_plot
                        export dos plot, options: t, spd, a, a-Si, a-1
  -v vasprun, --vasprun=vasprun
                        path of vasprun.xml file, default: vasprun.xml
  -f, --showforce       show forces, default: no
  -a, --allparameters   show all parameters
  -e, --eigenvalues     show eigenvalues in valence/conduction band
  -s smearing, --smear=smearing
                        smearing parameter for dos plot, e.g., 0.1 A
  -n dosfig, --dosfig=dosfig
                        dos figure name, default: dos.png
```
## Force information
```
$ python vasprun.py -v vasprun.xml -f yes
```
```
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

```
## DOS plots

A number of ways of plotting dos are also supported, some basic options are

- t: total dos
- spd: spd dos

one can just use the -d option to customize the plots.
If spin is included in vasprun.xml, the plot will show both up and down spin states separately.
```
$ python vasprun.py -v vasprun.xml -d t+spd -s 0.15 -n dos-spd.png
```
It will geneate a dos-spd.png figure like the following:
![dos.png](https://github.com/qzhu2017/vasprun/blob/master/images/dos.png)
