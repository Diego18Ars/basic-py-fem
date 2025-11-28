import math as m
import numpy as np
import matplotlib.pyplot as plt

from elements import Node, Element, Bar, Beam
from boundary_conditions import constraint_dof

n1 = Node(np.array([0, 0]), (1, 2, 3))
n2 = Node(np.array([1, 0]), (4, 5, 6))

A = 1.2*10**-3
E = 210*10**9
I_zz = 60*10**-7

beam1 = Beam(n1, n2)
beam1.defineSectionProperties(A, I_zz)
beam1.defineMaterialProperties(E)

# Assembly of global stiffness matrix
K_glb = np.zeros([6, 6])

for e in Element._registry:
    Ke = e.resolveLocalStiffness()
    e_DoFs = e.DoFs
    for i in range(len(e_DoFs)):
        for j in range(len(e_DoFs)):
            K_glb[e_DoFs[j-1]-1][e_DoFs[i-1]-1] += Ke[j-1][i-1]

F = np.array([0, 0, 0, 0, -1700, 0])

constraint_dof((1, 2, 3), 0, K_glb, F)

u = np.linalg.solve(K_glb, F)

# Basic Postproc, total nodal deformation
print(K_glb)
print(F)
print(u)