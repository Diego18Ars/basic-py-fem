import math as m
import numpy as np

from elements import Node, Bar

def constraint_dof(dofs: tuple, value: float, stiff: np.array, loads: np.array):
    penalty = np.max(stiff) * 10**4
    for dof in dofs:
        stiff[dof-1][dof-1] += penalty
        loads[dof-1] += value*penalty
