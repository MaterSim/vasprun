"""
   Module ```readdata``` contains a class ```XMLParser``` extracting
   both stress and strain from the ```vasprun.xml```
"""
import numpy as np
import defusedxml.ElementTree as ET

class VasprunError(Exception):
    def __init__(self, data):    
        self.data = data
    def __str__(self):
        return repr(self.data)

class XMLParser:
    def __init__(self,vasprun='vasprun.xml'):
        self.vasprun  = 'vasprun.xml'
        if isinstance(vasprun,str):
            self.vasprun = vasprun
        self.sigmas   = None
        self.epsilons = None

    def parse(self):
        root   = ET.parse(self.vasprun)
        # Checking if it makes sense to even start calculations
        IBRION = 0 # default value
        for line in root.find('incar').iter('i'):
            try:
                if line.attrib['name'] == 'IBRION':
                    IBRION=int(line.text)
            except KeyError:
                continue
        if IBRION < 5:
            raise VasprunError("File "+self.vasprun+" does not contain stain-stress calculations!")
        self.bases  = []
        self.stress = []
        for i,calc in enumerate(root.iter('calculation')):
            basis = []
            strss = []
            for varray in calc.iter('varray'):
                if varray.attrib['name'] == 'stress':
                    for line in varray:
                        strss.append(np.fromstring(line.text,sep=' '))
            self.stress.append(np.array(strss))
            for varray in calc.find('structure').find('crystal').iter('varray'):
                if varray.attrib['name'] == 'basis':
                    for line in varray:
                        basis.append(np.fromstring(line.text,sep=' '))
            self.bases.append(np.array(basis))
        del root

    def get_strain_stress(self):
        if self.sigmas is None or self.epsilons is None:
            self.generate_strain_stress()
        return self.epsilons,self.sigmas


    def generate_strain_stress(self):
        self.sigmas   = []
        self.epsilons = []
        for stress,basis in zip(self.stress[1:],self.bases[1:]):
            strain = np.dot(basis,np.linalg.inv(self.bases[0]))-np.identity(3)
            sigma  = stress - self.stress[0]
            if (np.abs(strain) > 1e-6).any():
                self.sigmas  .append(np.array([ sigma[0,0], sigma[1,1], sigma[2,2], sigma[0,1], sigma[1,2], sigma[2,0]]))
                self.epsilons.append(np.array([strain[0,0],strain[1,1],strain[2,2],strain[0,1],strain[1,2],strain[2,0]]))
