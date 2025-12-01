import numpy as np

class FiniteElementMethod:
    def __init__(self, dof):
        self.dofs = dof
        self.K = np.zeros([dof, dof])
        self.F = np.zeros([dof, 1])
        self.u = None
        self.penalty = None

    def assembleMatrix(self, elements):
        for elem in elements:
            Ke = elem.resolveLocalStiffness()
            elem_dofs = elem.dofs
            for i in range(len(elem_dofs)):
                for j in range(len(elem_dofs)):
                    self.K[elem_dofs[j-1]][elem_dofs[i-1]] += Ke[j-1][i-1]

        self.penalty = np.max(self.K) * 10**4

        return self.K

    def solve(self):
        self.u = np.linalg.solve(self.K, self.F)

        return self.u
    
    def applyConstraint(self, node, nodal_dof, value):
        '''Contrain Degree of Freedom for node. \n
        node - Node class object that indicates the node to be affected by the constraint \n
        nodal_dof - Index integer indicating which DoF is to be constrained \n
            1 - X translation \n
            2 - Y translation \n
            3 - Rotations \n
        value - Float indicating the displacement value on the input node to be constrained
        '''
        global_dof = node.dofs[nodal_dof-1] # Ready for 0-indexing arrays
        self.K[global_dof][global_dof] += self.penalty
        self.F[global_dof] += value*self.penalty

    def applyForce(self, node, nodal_dof, value):
        '''Apply a point load on node. \n
        node - Node class object that indicates the node to be affected by the constraint \n
        nodal_dof - Index integer indicating the force direction \n
            1 - +X direction \n
            2 - +Y direction \n
            3 - +Anticlockwise moments \n
        value - Force/moment magnitude
        '''
        global_dof = node.dofs[nodal_dof-1] # Ready for 0-indexing arrays
        self.F[global_dof] += value
