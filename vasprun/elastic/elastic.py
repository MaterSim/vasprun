"""
   Module ```elastic``` contains two classes:
     (i) ```Calculation``` containing the scipy-based engine 
     fiting the elastic tensor C to the stress-strain relation.
     (ii) ```Properites``` with the methods deriving the elastic properties
     from the elastic tensor (within different approximations)
"""
import numpy as np
from vasprun.elastic.CMatrices import ElasticTensor
from vasprun.elastic.readdata import XMLParser
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

class Calculator:
    def __init__(self,options=None,verbosity=0):
        if options is None:
            self.opt   = Options()
            self.opt,_ = self.opt.parse()
        else:
            self.opt = options
        
        data = XMLParser(self.opt.vasprun)
        data.parse()
        self.epsilons,self.sigmas = data.get_strain_stress()
        
        if verbosity > 0:
            print("Using "+self.opt.symmetry+" symmetry")
        self.ET = ElasticTensor(self.opt.symmetry)
        if verbosity > 0:
            print(self.ET.CMatrix[self.opt.symmetry].__doc__)

    def penalty(self,trialState):
        result = 0.0
        for e,s in zip(self.epsilons,self.sigmas):
            result += np.sum(np.power(np.transpose([s])+np.dot(self.ET(*trialState),np.transpose([e])),2))
        return result

    def __call__(self):
        solution        = least_squares(self.penalty,1e3*np.ones(ElasticTensor.DegreesOfFreedom[self.opt.symmetry]))
        self.resultantSystem = Properties(self.ET(*(0.1*solution.x))) #1e-1 -> moving from kBar to GPa
        return self.resultantSystem

    def print(self,approximation=None):
        print("Elastic Tensor C (GPa):")
        np.set_printoptions(suppress=True, precision=0)
        print(self.resultantSystem.Cij)
        print("Compliance Tensor s (1/GPa):")
        np.set_printoptions(suppress=True, precision=5)
        print(self.resultantSystem.sij)
        
        print("Elastic Properties :")
        print("Approximation Bulk_modulus(GPa) Shear_modulus(GPa) Young_modulus(GPa) Poisson_ratio Cauchy_pressure(GPa) Pugh_ratio")
        if approximation is None:
            result = self.get_moduli(self.opt.approximation)
        else:
            result = self.get_moduli(approximation)
        for r in result:
            print("%13s %12.0f %18.0f %18.0f %17.3f %16.0f %14.3f"%(r,*result[r]))

    def get_moduli(self,approximation='all'):
        return self.resultantSystem(approximation)


class Properties:
    maskC11 = np.array([
                       [ 1/3, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 1/3, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 1/3, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       ])
    maskC12 = np.array([
                       [ 0.0, 1/6, 1/6, 0.0, 0.0, 0.0, ],
                       [ 1/6, 0.0, 1/6, 0.0, 0.0, 0.0, ],
                       [ 1/6, 1/6, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       ])
    maskC44 = np.array([
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 1/3, 0.0, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 1/3, 0.0, ],
                       [ 0.0, 0.0, 0.0, 0.0, 0.0, 1/3, ],
                       ])
    def setCij(self,Cij):
        try:
            Cij[5,5] += 0.0
            self.Cij = Cij
        except (TypeError,IndexError):
            try:
                self.Cij = np.reshape(np.array(Cij),(6,6))
            except:
                raise TypeError('elastic matrix must be either numpy.array or an array-like with at least 6x6 fields')
        self.sij = np.linalg.inv(self.Cij)
    
    def __init__(self,Cij):
        self.setCij(Cij)

    def average(self,mask):
        return np.sum(mask*self.Cij)

    def __call__(self,approximation="Hill"):
        self.elastic_properties = {}
        if approximation[0] in "AaVv":
            K,G     = self.voigt()
            E,nu    = Properties.derivative_properties(K,G)
            Cauchy  = self.average(self.maskC12)
            Cauchy -= self.average(self.maskC44)
            Pugh    = G/K
            self.elastic_properties['Voigt'] = (K,G,E,nu,Cauchy,Pugh)
        if approximation[0] in "AaRr":
            K,G     = self.reuss()
            E,nu    = Properties.derivative_properties(K,G)
            Cauchy  = self.average(self.maskC12)
            Cauchy -= self.average(self.maskC44)
            Pugh    = G/K
            self.elastic_properties['Reuss'] = (K,G,E,nu,Cauchy,Pugh)
        if approximation[0] in "AaHh":
            KV,GV     = self.voigt()
            KR,GR     = self.reuss()
            K       = 0.5*(KV+KR)
            G       = 0.5*(GV+GR)
            E,nu    = Properties.derivative_properties(K,G)
            Cauchy  = self.average(self.maskC12)
            Cauchy -= self.average(self.maskC44)
            Pugh    = G/K
            self.elastic_properties['Hill'] = (K,G,E,nu,Cauchy,Pugh)
        return self.elastic_properties

    def voigt(self):
        K = ((self.Cij[0,0] + self.Cij[1,1] + self.Cij[2,2]) 
       + 2.0*(self.Cij[0,1] + self.Cij[1,2] + self.Cij[2,0]))/9.0
        G = ((self.Cij[0,0] + self.Cij[1,1] + self.Cij[2,2]) 
       -     (self.Cij[0,1] + self.Cij[1,2] + self.Cij[2,0])
       + 4.0*(self.Cij[3,3] + self.Cij[4,4] + self.Cij[5,5]))/15.0
        return K,G

    def reuss(self):
        K =      np.reciprocal(
                (self.sij[0,0] + self.sij[1,1] + self.sij[2,2])
          + 2.0*(self.sij[0,1] + self.sij[1,2] + self.sij[2,0]))
        G = 15.0*np.reciprocal(
            4.0*(self.sij[0,0] + self.sij[1,1] + self.sij[2,2]) 
          - 4.0*(self.sij[0,1] + self.sij[1,2] + self.sij[2,0])
          + 3.0*(self.sij[3,3] + self.sij[4,4] + self.sij[5,5]))
        return K,G

    @staticmethod
    def derivative_properties(K,G):
        E = 9.0*K*G/(3*K+G)
        nu= (3*K-2*G)/(6*K+2*G)
        return E,nu
