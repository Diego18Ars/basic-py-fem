import math as m
import numpy as np

from elements import Node, Bar

n1 = Node(np.array([0, 0]), (1, 2))
n2 = Node(np.array([1, 0]), (3, 4))
n3 = Node(np.array([0, 1]), (5, 6))

A = 1.2*10**-3
E = 210*10**9

bar1 = Bar(n1, n2)
bar1.defineSection(A)
bar1.defineMaterialMechProps(E)

bar2 = Bar(n2, n3)
bar2.defineSection(A)
bar2.defineMaterialMechProps(E)

# Assembly of global stiffness matrix
K1 = bar1.resolveLocalStiffness()
K2 = bar2.resolveLocalStiffness()

K_glb = np.zeros([6, 6])

e1_DoFs = bar1.DoFs

for i in range(len(e1_DoFs)):
    for j in range(len(e1_DoFs)):
        K_glb[e1_DoFs[j-1]-1][e1_DoFs[i-1]-1] += K1[j-1][i-1]

e2_DoFs = bar2.DoFs

for i in range(len(e2_DoFs)):
    for j in range(len(e2_DoFs)):
        K_glb[e2_DoFs[j-1]-1][e2_DoFs[i-1]-1] += K2[j-1][i-1]

K_glb[0][0] += 2.52*10**(8+4)
K_glb[1][1] += 2.52*10**(8+4)
K_glb[4][4] += 2.52*10**(8+4)
K_glb[5][5] += 2.52*10**(8+4)

F = np.array([0, 0, -17000, 0, 0, 0])

u = np.linalg.solve(K_glb, F)

# Debug
print(u)
print(K1)
print(K2)
print(K_glb)