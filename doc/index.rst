Introduction
============

``vasprun`` is a python project used for quick analysis of `VASP
<https://www.vasp.at>`_ calculation solely from ``vasprun.xml``. It has the following features:

-  band gap calculation
-  dos plot (total dos / orbital partial dos / atomic partial dos)
-  band structure plot (with color map enhancement)
-  incar/potcar/poscar generation
-  force analysys
-  Kohn-Sham orbital eigenvalue analysys
-  Infrared intensity analysis
-  dielectric constants
-  elastic constants (to add)


Version info
============
The current version is ``1.0`` at `GitHub <https://github.com/qzhu2017/vasprun>`_. 

Expect updates upon request by `Qiang Zhu <http://www.physics.unlv.edu/~qzhu/index.html>`_ at University of Nevada Las Vegas.

Installation and Setup
======================

Dependencies
------------

Versions indicated are those used during development. 

  * `SciPy 1.0.1 <https://www.scipy.org/install.html>`_  
  * `NumPy 1.14.3 <https://www.scipy.org/scipylib/download.html>`_  
  * `Pandas 0.20.3 <https://pandas.pydata.org/getpandas.html>`_  
  * `Pymatgen <http://pymatgen.org/#getting-pymatgen>`_  

Installation
------------

To install it, first install all dependencies, then make a copy of the source code:

``git clone https://github.com/qzhu2017/vasprun``

Then, inside of the downloaded directory, run

``python setup.py install``

This will install the module. The code can be used within Python via

.. code-block:: Python

  import vasprun
  print(vasprun.__version__)

Tutorials
========================

.. toctree::
   Usage1
   Usage2
   :maxdepth: 2
