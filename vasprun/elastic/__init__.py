"""
   Submodule for calculation of the elastic properties
   from the stres--strain relation based on VASP's ```vasprun.xml``` file.

   Module CMatrix contains a class ```ElasticTensor``` 
   with all 9 different symmetry-based constrains 6x6 elastic tensors.

   Module ```elastic``` contains two classes:
     (i) ```Calculation``` containing the scipy-based engine 
     fitting the elastic tensor C to the stress-strain relation.
     (ii) ```Properites``` with the methods deriving the elastic properties
     from the elastic tensor (within different approximations)

   Module ```readdata``` contains a class ```XMLParser``` extracting
   both stress and strain from the ```vasprun.xml```
"""
