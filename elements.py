import math as m
import numpy as np

class Node():
    def __init__(self, coords: np.array, DoFs: tuple):
        self.coords = coords
        self.DoFs = DoFs

class Bar():
    # 1st Order Bar Element - Axial Deformation

    def __init__(self, node1: Node, node2: Node):
        # Nodes
        self.node1 = node1
        self.node2 = node2

        self.DoFs = (self.node1.DoFs[0], self.node1.DoFs[1], self.node2.DoFs[0], self.node2.DoFs[1])

        # Bar element lenght
        self.L = m.sqrt((self.node1.coords[0]-self.node2.coords[0])**2 + (self.node1.coords[1]-self.node2.coords[1])**2)

        # Bar angle with +X axis and rotation matrix
        r = self.node2.coords-self.node1.coords
        a = m.acos((r[0]*1 + r[1]*0)/(self.L*1))

        self.rotate_T = np.array([[m.cos(a), m.sin(a), 0, 0],
                                  [-m.sin(a), m.cos(a), 0, 0],
                                  [0, 0, m.cos(a), m.sin(a)],
                                  [0, 0, -m.sin(a), m.cos(a)]])

        # Element shape functions
        self.sf_1 = lambda s: (-1/self.L) * (s) + 1
        self.sf_2 = lambda s: (1/self.L) * (s-self.L) + 1 

        # Section properties
        self.A = None

        # Material Mechanical Properties
        self.E = None

        # Element-wise stiffness matrix for 4 DoFs
        self.K_ = np.array([[1, 0, -1, 0],
                            [0, 0, 0, 0],
                            [-1, 0, 1, 0],
                            [0, 0, 0, 0]])

    def defineSection(self, A):
        self.A = A

    def defineMaterialMechProps(self, E):
        self.E = E

    def resolveLocalStiffness(self):
        K_unresolved = self.K_ @ self.rotate_T
        K_unresolved = self.rotate_T.transpose() @ K_unresolved
        K_unresolved = (self.E*self.A / self.L) * K_unresolved
        self.K_ = K_unresolved

        return self.K_

    def retrieveLocalStiffness(self):
        return self.K_
    
