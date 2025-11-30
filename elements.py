import math as m
import numpy as np

class Node:
    dof_count = 0

    def __init__(self, x, y):
        self.coords = np.array([x, y])
        self.dofs = [Node.dof_count, Node.dof_count+1, Node.dof_count+2]
        Node.dof_count += 3

class Element1D:
    elements = []

    def __init__(self, node1: Node, node2: Node):
        self.coords1 = node1.coords
        self.coords2 = node2.coords
        self.dofs = node1.dofs + node2.dofs

        Element1D.elements.append(self)

        self.L= m.sqrt((self.coords1[0]-self.coords2[0])**2 + (self.coords1[1]-self.coords2[1])**2)

        r = self.coords2-self.coords1
        a = m.acos((r[0]*1 + r[1]*0)/(self.L*1))
        print(a)

        self.rotate_t = np.array([[m.cos(a), m.sin(a), 0, 0, 0, 0],
                                  [-m.sin(a), m.cos(a), 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0],
                                  [0, 0, 0, m.cos(a), m.sin(a), 0],
                                  [0, 0, 0, -m.sin(a), m.cos(a), 0],
                                  [0, 0, 0, 0, 0, 1]])

        # Stiffness matrix intialization
        self.K_ = None


class Bar(Element1D):
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
        self.K_ = K_unresolved

        return self.K_

    def getLocalStiffness(self):
        return self.K_


class Beam(Element1D):
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

        return self.K_

    def getLocalStiffness(self):
        return self.K_
