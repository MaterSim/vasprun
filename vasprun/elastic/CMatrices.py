"""
      Functions allowing for calculation a 6x6 elastic-constant matrix
   for 9 different point-symmetry space groups. Make sure that you use
   the proper ordering of the Cij constants!
"""
import numpy as np

class ElasticTensor:
    DegreesOfFreedom = {
                        'Triclinic'    : 21,
                        'Monoclinic'   : 13,
                        'Orthorhombic' : 9,
                        'Trigonal6'    : 6,
                        'Trigonal7'    : 7,
                        'Hexagonal'    : 5,
                        'Tetragonal6'  : 6,
                        'Tetragonal7'  : 7,
                        'Cubic'        : 3,
                        }
    def __init__(self,symmetry=None):
        self.CMatrix = {
            'Triclinic':    ElasticTensor.CMatrixTriclinic,
            'Monoclinic':   ElasticTensor.CMatrixMonoclinic,
            'Orthorhombic': ElasticTensor.CMatrixOrthorhombic,
            'Trigonal6':    ElasticTensor.CMatrixTrigonal6,
            'Trigonal7':    ElasticTensor.CMatrixTrigonal7,
            'Hexagonal':    ElasticTensor.CMatrixHexagonal,
            'Tetragonal6':  ElasticTensor.CMatrixTetragonal6,
            'Tetragonal7':  ElasticTensor.CMatrixTetragonal7,
            'Cubic':        ElasticTensor.CMatrixCubic
            }

        self.symmetry = None
        if isinstance(symmetry,str):
            self.symmetry = symmetry

    def __call__(self,*args,symmetry=None):
        if symmetry is not None:
            return self.CMatrix[symmetry](*args)
        elif self.symmetry is not None:
            return self.CMatrix[self.symmetry](*args)
        else:
            raise TypeError("You must define symmetry!")

    @staticmethod
    def CMatrixTriclinic(C11,C12,C13,C14,C15,C16,
                             C22,C23,C24,C25,C26,
                                 C33,C34,C35,C36,
                                     C44,C45,C46,
                                         C55,C56,
                                             C66):
        """
        Triclinic cell - 21 constants
        Order:
          C11,C12,C13,C14,C15,C16,
              C22,C23,C24,C25,C26,
                  C33,C34,C35,C36,
                      C44,C45,C46,
                          C55,C56,
                              C66
        returns an array: [
            [ C11, C12, C13, C14, C15, C16 ],
            [ C12, C22, C23, C24, C25, C26 ],
            [ C13, C23, C33, C34, C35, C36 ],
            [ C14, C24, C34, C44, C45, C46 ],
            [ C15, C25, C35, C45, C55, C56 ],
            [ C16, C26, C36, C46, C56, C66 ],
            ]
        """
        return np.array([
            [ C11, C12, C13, C14, C15, C16 ],
            [ C12, C22, C23, C24, C25, C26 ],
            [ C13, C23, C33, C34, C35, C36 ],
            [ C14, C24, C34, C44, C45, C46 ],
            [ C15, C25, C35, C45, C55, C56 ],
            [ C16, C26, C36, C46, C56, C66 ],
            ])
    
    @staticmethod
    def CMatrixMonoclinic(C11,C12,C13,    C15,
                              C22,C23,    C25,
                                  C33,    C35,
                                      C44,    C46,
                                          C55,
                                              C66):
        """
        Monoclinic cell - 13 constants
        Order:
          C11,C12,C13,C15,
              C22,C23,C25,
                  C33,C35,
                      C44,C46,
                          C55,
                              C66
        returns an array: [
            [ C11, C12, C13, 0.0, C15, 0.0 ],
            [ C12, C22, C23, 0.0, C25, 0.0 ],
            [ C13, C23, C33, 0.0, C35, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, C46 ],
            [ C15, C25, C35, 0.0, C55, 0.0 ],
            [ 0.0, 0.0, 0.0, C46, 0.0, C66 ],
            ]
        """
        return np.array([
            [ C11, C12, C13, 0.0, C15, 0.0 ],
            [ C12, C22, C23, 0.0, C25, 0.0 ],
            [ C13, C23, C33, 0.0, C35, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, C46 ],
            [ C15, C25, C35, 0.0, C55, 0.0 ],
            [ 0.0, 0.0, 0.0, C46, 0.0, C66 ],
            ])
    
    @staticmethod
    def CMatrixOrthorhombic(C11,C12,C13,
                                C22,C23,
                            C33,C44,C55,C66):
        """
        Orthorhombic cell - 9 constants
        Order:
          C11,C12,C13,
              C22,C23,
          C33,C44,C55,C66
        returns an array: [
            [ C11, C12, C13, 0.0, 0.0, 0.0 ],
            [ C12, C22, C23, 0.0, 0.0, 0.0 ],
            [ C13, C23, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C55, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C66 ],
            ]
        """
        return np.array([
            [ C11, C12, C13, 0.0, 0.0, 0.0 ],
            [ C12, C22, C23, 0.0, 0.0, 0.0 ],
            [ C13, C23, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C55, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C66 ],
            ])
    
    @staticmethod
    def CMatrixTrigonal7(C11,C12,C13,C14,C15,
                                 C33,C44):
        """
        Trigonal cell - 7 constants
        Order
          C11,C12,C13,C14,C15,
                  C33,C44
        returns an array: [
            [ C11, C12, C13, C14, C15, 0.0 ],
            [ C12, C11, C13,-C14,-C15, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ C14,-C14, 0.0, C44, 0.0,-C15 ],
            [ C15,-C15, 0.0, 0.0, C44, C14 ],
            [ 0.0, 0.0, 0.0,-C15, C14, C66 ],
            ]
        where C66 = 0.5*(C11-C12)
        """
        C66 = 0.5*(C11-C12)
        return np.array([
            [ C11, C12, C13, C14, C15, 0.0 ],
            [ C12, C11, C13,-C14,-C15, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ C14,-C14, 0.0, C44, 0.0,-C15 ],
            [ C15,-C15, 0.0, 0.0, C44, C14 ],
            [ 0.0, 0.0, 0.0,-C15, C14, C66 ],
            ])
    
    @staticmethod
    def CMatrixTrigonal6(C11,C12,C13,C14,
                                 C33,C44):
        """
        Trigonal cell - 6 constants
        Order
          C11,C12,C13,C14,
                  C33,C44
        returns an array: [
            [ C11, C12, C13, C14, 0.0, 0.0 ],
            [ C12, C11, C13,-C14, 0.0, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ C14,-C14, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, C14 ],
            [ 0.0, 0.0, 0.0, 0.0, C14, C66 ],
            ]
        where C66 = 0.5*(C11-C12)
        """
        C66 = 0.5*(C11-C12)
        return np.array([
            [ C11, C12, C13, C14, 0.0, 0.0 ],
            [ C12, C11, C13,-C14, 0.0, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ C14,-C14, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, C14 ],
            [ 0.0, 0.0, 0.0, 0.0, C14, C66 ],
            ])
    
    @staticmethod
    def CMatrixHexagonal(C11,C12,C13,
                                 C33,C44):
        """
        Hexagonal cell - 5 constants
        Order:
          C11,C12,C13,
                  C33,C44
        returns an array: [
            [ C11, C12, C13, 0.0, 0.0, 0.0 ],
            [ C12, C11, C13, 0.0, 0.0, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C66 ],
            ]
        where C66 = 0.5*(C11-C12)
        """
        C66 = 0.5*(C11-C12)
        return np.array([
            [ C11, C12, C13, 0.0, 0.0, 0.0 ],
            [ C12, C11, C13, 0.0, 0.0, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C66 ],
            ])
    
    @staticmethod
    def CMatrixTetragonal7(C11,C12,C13,    C16,
                                   C33,C44,C66):
        """
        Tetragonal cell - 7 constants
        Order:
          C11,C12,C13,    C16,
                  C33,C44,C66
        returns an array: [
            [ C11, C12, C13, 0.0, 0.0, C16 ],
            [ C12, C11, C13, 0.0, 0.0,-C16 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ C16,-C16, 0.0, 0.0, 0.0, C66 ],
            ]
        """
        return np.array([
            [ C11, C12, C13, 0.0, 0.0, C16 ],
            [ C12, C11, C13, 0.0, 0.0,-C16 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ C16,-C16, 0.0, 0.0, 0.0, C66 ],
            ])
    
    @staticmethod
    def CMatrixTetragonal6(C11,C12,C13,
                           C33,C44,C66):
        """
        Tetragonal cell - 6 constants
        Order:
          C11,C12,C13,
          C33,C44,C66
        returns an array: [
            [ C11, C12, C13, 0.0, 0.0, 0.0 ],
            [ C12, C11, C13, 0.0, 0.0, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C66 ],
            ]
        """
        return np.array([
            [ C11, C12, C13, 0.0, 0.0, 0.0 ],
            [ C12, C11, C13, 0.0, 0.0, 0.0 ],
            [ C13, C13, C33, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C66 ],
            ])
    
    @staticmethod
    def CMatrixCubic(C11,C12,C44):
        """
        Cubic cell - 3 constants
        Order:
          C11,C12,C44
        returns an array: [
            [ C11, C12, C12, 0.0, 0.0, 0.0 ],
            [ C12, C11, C12, 0.0, 0.0, 0.0 ],
            [ C12, C12, C11, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C44 ],
            ]
        """
        return np.array([
            [ C11, C12, C12, 0.0, 0.0, 0.0 ],
            [ C12, C11, C12, 0.0, 0.0, 0.0 ],
            [ C12, C12, C11, 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, C44, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, C44, 0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, C44 ],
            ])
