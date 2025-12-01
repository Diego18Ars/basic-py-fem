import numpy as np
from elements.node import Node
from elements.element_1d import Element1D

class Beam2(Element1D):
    '''
    Euler-Bernoulli beam element that incorporates bar behaviour. 
    '''


    def __init__(self, node1: None, node2: Node):
        Element1D.elements.append(self)

        super().__init__(node1, node2)

        self.K_ = None

        # Initialize section and material properties
        self.E = None
        self.A = None
        self.I_zz = None
        self.I_yy = None
        self.t_zz = None
        self.t_yy = None

    def defineSectionProperties(self, A, I_zz, I_yy=None, t_zz=None, t_yy=None):
        self.A = A
        self.I_zz = I_zz
        self.I_yy = I_yy
        self.t_zz = t_zz
        self.t_yy = t_yy

    def defineMaterialProperties(self, E):
        self.E = E

    def resolveLocalStiffness(self):
        K_beam = (2*self.E*self.I_zz/self.L**3) * np.array([
                        [0, 0, 0, 0, 0, 0],
                        [0, 6, -3*self.L, 0, -6, -3*self.L],
                        [0, -3*self.L, 2*self.L**2, 0, 3*self.L, self.L**2],
                        [0, 0, 0, 0, 0, 0],
                        [0, -6, 3*self.L, 0, 6, 3*self.L],
                        [0, -3*self.L, self.L**2, 0, 3*self.L, 2*self.L**2]
                        ])
        
        K_bar = (self.E*self.A/self.L) * np.array([
                        [1, 0, 0, -1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [-1, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0]
                        ])
        
        K_unresolved = K_beam + K_bar
        K_unresolved = K_unresolved @ self.rotate_t
        K_unresolved = self.rotate_t.transpose() @ K_unresolved
        self.K_ = K_unresolved

        return self.K_e

    def getLocalStiffness(self):
        return self.K_e
