# vasprun
This code is used for quick analysis of vasp calculation solely from vasprun.xml. It has the following features:

- band gap calculation
- dos plot (total dos / orbital partial dos / atomic partial dos)
- band structure plot (with color map enhancement)
- incar/potcar/poscar generation
- force analysys
- KS orbital eigenvalue analysys
- dynamical matrix (to add)
- elastic constants (to add)
- dielectric constants (to add)

## Prerequisites
To use it, one must have python 3 installed. It is recommended to install conda if one don't have root access in his computer.
In addition, several packages will be required
- pymatgen
- pandans
- lxml

A full documentation has been moved to https://vasprun-xml.readthedocs.io/en/latest/index.html

