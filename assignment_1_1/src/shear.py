import numpy as np
from geometry import compute_mesh_centroid
# -----------------------------------------------------------------------------
#                               Mesh geometry
# -----------------------------------------------------------------------------

def shear_transformation(V, nu):
    """
    Computes vertices' postion after the shear transformation.

    Input:
    - V : np.array (|V|, 3)
        The array of vertices positions.
        Contains the coordinates of the i-th vertex in i-th row
    - nu : float
        The shear parameter
    
    Output:
    - V1 : np.array (|V|, 3)
        The array of vertices positions after transformation.
    """
    V1 = V.copy()


    # HW 1.3.6
    # enter your code here

    return V1


def shear_equilibrium(x_csl, V, F):
    """
    Shear the input mesh to make it equilibrium.

    Input:
    - x_csl: float
        The x coordinate of the target centroid
    - V : np.array (|V|, 3)
        The array of vertices positions.
        Contains the coordinates of the i-th vertex in i-th row
    - F : np.array (|F|, 3)
        The array of triangle faces.
    Output:
    - V1 : np.array (|V|, 3)
        The array of vertices positions that are equilibrium.
    """
    V1 = V.copy()

    # HW 1.3.7
    # enter your code here

    return V1
