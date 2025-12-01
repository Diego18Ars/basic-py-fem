from elements.node import Node
from elements.element_1d import Element1D
from elements.bar2 import Bar2
from finite_element_method import FiniteElementMethod

n1 = Node(0, 0)
n2 = Node(0, 1)
n3 = Node(1, .5)

A = 1.2*10**-3
E = 210*10**9
I_zz = 60*10**-7

bar1 = Bar2(n1, n2)
bar1.defineSectionProperties(A)
bar1.defineMaterialProperties(E)

bar2 = Bar2(n2, n3)
bar2.defineSectionProperties(A)
bar2.defineMaterialProperties(E)

# Assembly of global stiffness matrix
FEM = FiniteElementMethod(Node.dof_count)
FEM.assembleMatrix(Element1D.elements)
FEM.applyConstraint(n1, 1, 0)
FEM.applyConstraint(n1, 2, 0)
FEM.applyConstraint(n1, 3, 0)
FEM.applyConstraint(n3, 1, 0)
FEM.applyConstraint(n3, 2, 0)
FEM.applyConstraint(n3, 3, 0)
FEM.applyConstraint(n2, 3, 0)

FEM.applyForce(n2, 1, 5000)
FEM.applyForce(n2, 2, -5000)

FEM.solve()

print(FEM.K)
print(FEM.F)
print(FEM.u)
