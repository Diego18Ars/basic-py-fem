import numpy as np
from elements.node import Node
from elements.element_1d import Element1D

class Bar2(Element1D):
    '''
    Bar element that only resists axial deformation in its nodes.
    '''

    def __init__(self, node1: Node, node2: Node):

        super().__init__(node1, node2)

        # Initialize section and material properties
        self.A = None
        self.E = None

    def defineSectionProperties(self, A):
        self.A = A

    def defineMaterialProperties(self, E):
        self.E = E

    def resolveLocalStiffness(self):
        K_unresolved = (self.E*self.A / self.L) * np.array([
                        [1, 0, 0, -1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [-1, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0]
                        ])
        
        K_unresolved = K_unresolved @ self.rotate_t
        K_unresolved = self.rotate_t.transpose() @ K_unresolved
        self.K_e = K_unresolved

        return self.K_e

    def getLocalStiffness(self):
        return self.K_e
