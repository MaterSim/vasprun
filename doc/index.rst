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
The current version is ``1.0.3`` at `GitHub <https://github.com/qzhu2017/vasprun>`_. 

Expect updates upon request by `Qiang Zhu <https://qzhu2017.github.io>`_ at University of Nevada Las Vegas.

Installation and Setup
======================
This code is written based on Python 3. Python 2.x won't be supported

Dependencies
------------

Required packages:

- `lxml>=4.2.5 <https://pypi.org/project/lxml/>`_
- `Matplotlib>=2.0.0 <https://matplotlib.org>`_
- `NumPy>=1.13.3 <https://www.scipy.org/scipylib/download.html>`_  
- `SciPy>1.1.0 <https://www.scipy.org/install.html>`_  
- `Pandas>=0.23.4 <https://pandas.pydata.org/getpandas.html>`_  

Installation
------------

To install it, one can simply type ``pip install vasprun-xml`` or make a copy of the source code, and then install it manually.
::

  git clone https://github.com/qzhu2017/vasprun.git
  cd vasprun
  python setup.py install

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

Modules
========================

.. toctree::
   :maxdepth: 4

   vasprun
