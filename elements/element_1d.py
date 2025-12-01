import math as m
import numpy as np
from elements.node import Node

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

        self.rotate_t = np.array([[m.cos(a), m.sin(a), 0, 0, 0, 0],
                                  [-m.sin(a), m.cos(a), 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0],
                                  [0, 0, 0, m.cos(a), m.sin(a), 0],
                                  [0, 0, 0, -m.sin(a), m.cos(a), 0],
                                  [0, 0, 0, 0, 0, 1]])

        # Stiffness matrix intialization
        self.K_e = None
