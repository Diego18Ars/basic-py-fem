import numpy as np

from elements import Node

def constrain_dof(node: Node, nodal_dof: int, value: float, K: np.array, F: np.array, penalty):
    '''Contrain Degree of Freedom for Node. \n
    node - Node class object that indicates the node to be affected by the constraint \n
    nodal_dof - Index integer indicating which DoF is to be constrained \n
        1 - X translation \n
        2 - Y translation \n
        3 - Rotations \n
    value - Float indicating the displacement value on the input node to be constrained
    '''
    global_dof = node.dofs[nodal_dof-1] # Ready for 0-indexing arrays
    K[global_dof][global_dof] += penalty
    F[global_dof] += value*penalty

def couple_dofs(dofs: tuple, stiff: np.array):
    penalty = np.max(stiff) * 10**4
    for dof1 in dofs:
        for dof2 in dofs:
            if dof1 == dof2:
                stiff[dof1-1][dof2-1] += penalty
            else: 
                stiff[dof1-1][dof2-1] -= penalty