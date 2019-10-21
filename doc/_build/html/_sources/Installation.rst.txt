Installation and Setup
======================

Dependencies
------------

Versions indicated are those used during development. Other versions may be compatible, but have not yet been tested.

  * `SciPy 1.0.1 <https://www.scipy.org/install.html>`_  
  * `NumPy 1.14.3 <https://www.scipy.org/scipylib/download.html>`_  
  * `Pandas 0.20.3 <https://pandas.pydata.org/getpandas.html>`_  
  * `Pymatgen 2017.9.3 <http://pymatgen.org/#getting-pymatgen>`_  
  * `SpgLib for Python 1.9.9.44 <https://atztogo.github.io/spglib/python-spglib.html#installation>`_  

Installation
------------

To install PyXtal, first install all dependencies, then make a copy of the source code:

``git clone https://github.com/qzhu2017/pyxtal``

Then, inside of the downloaded directory, run

``python setup.py install``

This will install the module. The code can be used within Python via

.. code-block:: Python

  import pyxtal
  
More extensive test can be invoked by running

``python $root/pyxtal/test_cases/test_all.py``

Ideally, one should see the completion of all modules in the end.
