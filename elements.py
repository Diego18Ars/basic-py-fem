import math as m
import numpy as np

class Node:
    _registry = 0

    def __init__(self, coords: np.array, DoFs: tuple):
        self.coords = coords
        self.DoFs = DoFs
        self._registry += 2 

class Element:
    _registry = []

class Bar(Element):
    # 1st Order Bar Element - Axial Deformation

    def __init__(self, node1: Node, node2: Node):
        # Append element
        self._registry.append(self)

        # Nodes
        self.node1 = node1
        self.node2 = node2

        self.DoFs = (self.node1.DoFs[0], self.node1.DoFs[1], self.node1.DoFs[2], self.node2.DoFs[0], self.node2.DoFs[1], self.node2.DoFs[2])

        # Bar element lenght
        self.L = m.sqrt((self.node1.coords[0]-self.node2.coords[0])**2 + (self.node1.coords[1]-self.node2.coords[1])**2)

        # Bar angle with +X axis and rotation matrix
        r = self.node2.coords-self.node1.coords
        a = m.acos((r[0]*1 + r[1]*0)/(self.L*1))

        self.rotate_T = np.array([[m.cos(a), m.sin(a), 0, 0, 0, 0],
                                  [-m.sin(a), m.cos(a), 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0],
                                  [0, 0, 0, m.cos(a), m.sin(a), 0],
                                  [0, 0, 0, -m.sin(a), m.cos(a), 0],
                                  [0, 0, 0, 0, 0, 1]])

        # Element shape functions
        self.sf_1 = lambda s: (-1/self.L) * (s) + 1
        self.sf_2 = lambda s: (1/self.L) * (s-self.L) + 1 

        # Section properties
        self.A = None

        # Material Mechanical Properties
        self.E = None

        # Element-wise stiffness matrix for 4 DoFs
        self.K_ = None

    def defineSectionProperties(self, A):
        self.A = A

    def defineMaterialProperties(self, E):
        self.E = E

    def resolveLocalStiffness(self):
        K_unresolved = (self.E*self.A/self.L) * np.array([
                        [1, 0, 0, -1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [-1, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0]
                        ])
        
        K_unresolved = K_unresolved @ self.rotate_T
        K_unresolved = self.rotate_T.transpose() @ K_unresolved
        K_unresolved = (self.E*self.A / self.L) * K_unresolved
        self.K_ = K_unresolved

        return self.K_

    def retrieveLocalStiffness(self):
        return self.K_
    
class Beam(Element):
    # Euler-Bernoulli beam with axial deformation as well
    def __init__(self, node1: None, node2: Node):
        self._registry.append(self)

        self.node1 = node1
        self.node2 = node2

        self.DoFs = (self.node1.DoFs[0], self.node1.DoFs[1], self.node1.DoFs[2], self.node2.DoFs[0], self.node2.DoFs[1], self.node2.DoFs[2])

        self.K_ = None

        # Bar element lenght
        self.L = m.sqrt((self.node1.coords[0]-self.node2.coords[0])**2 + (self.node1.coords[1]-self.node2.coords[1])**2)

        # Bar angle with +X axis and rotation matrix
        r = self.node2.coords-self.node1.coords
        a = m.acos((r[0]*1 + r[1]*0)/(self.L*1))

        self.rotate_T = np.array([[m.cos(a), m.sin(a), 0, 0, 0, 0],
                                  [-m.sin(a), m.cos(a), 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0],
                                  [0, 0, 0, m.cos(a), m.sin(a), 0],
                                  [0, 0, 0, -m.sin(a), m.cos(a), 0],
                                  [0, 0, 0, 0, 0, 1]])

        # Initialize geometric and material properties´
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
        K_unresolved = K_unresolved @ self.rotate_T
        K_unresolved = self.rotate_T.transpose() @ K_unresolved
        K_unresolved = K_unresolved
        self.K_ = K_unresolved

        return self.K_

    def retrieveLocalStiffness(self):
        return self.K_
