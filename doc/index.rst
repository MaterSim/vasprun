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
This code is written based on Python 3. Python 2.x won't be supported

Dependencies
------------

Required packages:

- `lxml <https://pypi.org/project/lxml/>`_
- `Matplotlib <https://matplotlib.org>`_
- `NumPy <https://www.scipy.org/scipylib/download.html>`_  
- `Pandas <https://pandas.pydata.org/getpandas.html>`_  

Optional for some features:

- `Pymatgen <http://pymatgen.org/#getting-pymatgen>`_  
- `SciPy <https://www.scipy.org/install.html>`_  

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
