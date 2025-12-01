import numpy as np

class Node:
    dof_count = 0

    def __init__(self, x, y):
        self.coords = np.array([x, y])
        self.dofs = [Node.dof_count, Node.dof_count+1, Node.dof_count+2]
        Node.dof_count += 3