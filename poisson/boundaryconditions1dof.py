import numpy as np
import assemble1dof
from scipy import sparse
import math
from numba import jit


def dirichlet(K, B, mesh, temperature):
    """Apply Dirichlet BC.

    .. note::

        How its done:

        1. Loop over the lines where Dirichlet boundary conditions are
        applied. This is specified as key argument of the dictionary. The
        boundaries are defined on the gmsh nad tagged with a number.

        2. loop over all the nodes at the boundaries.

        3. If those nodes, identified as follows::

            boundary_nodes = [line node1 node2]

            are in the line where dirichlet BC were specified, change the
            stiffness matrix and the B vector based on this node index.

    Args:
        K (2nd order array): Stiffness matrix.
        B (1st order array): Vector with the load and traction.
        temperature (function): Function with 4 components.


    Returns:
        K (2nd order array), B (1st order array): Modified stiffness matrix
        and vector.



    """
    for line in temperature(1,1).keys():
        for n in range(len(mesh.boundary_nodes[:, 0])):
            if line == mesh.boundary_nodes[n, 0]:
                K[mesh.boundary_nodes[n, 1], :] = 0.0
                K[mesh.boundary_nodes[n, 2], :] = 0.0

                K[mesh.boundary_nodes[n, 1], mesh.boundary_nodes[n, 1]] =  1.0
                K[mesh.boundary_nodes[n, 2], mesh.boundary_nodes[n, 2]] = 1.0

                t1 = temperature(mesh.nodes_coord[mesh.boundary_nodes[n,1], 0],
                                mesh.nodes_coord[mesh.boundary_nodes[n, 1], 1],)

                t2 = temperature(mesh.nodes_coord[mesh.boundary_nodes[n, 2], 0],
                                mesh.nodes_coord[mesh.boundary_nodes[n, 2], 1],)

                B[mesh.boundary_nodes[n, 1]] = t1[line]
                B[mesh.boundary_nodes[n, 2]] = t2[line]

    return K, B


@jit
def neumann(mesh, traction):
    """Apply Neumann BC.

    Computes the integral from the weak form with the boundary term.

    .. note::

        How its done:

        1. Define an array with the Gauss points for each path, each path will
        have 2 sets of two gauss points. One of the gp is fixed, which indicates
        that its going over an specific boundary of the element.

        Gauss points are defines as::

            gp = [gp, -1    1st path --> which has 2 possibilities for gp.
                  1,  gp    2nd path
                  gp,  1    3rd
                  -1,  gp]  4th

        2. A loop over the elements on the boundary extracting also the side
        where the boundary is located on this element.

    .. note::

        Edge elements are necessary because we need to extract the nodes
        from the connectivity. Then we can create a T for this element.

    Args:
        traction: Function with the traction and the line where the traction
            is applied.
        mesh: Object with the mesh attributes.

    Returns:
        T: Traction vector with size equals the dof.

    """
    Tele = np.zeros((4, mesh.num_ele))


    gp = np.array([[[-1.0/math.sqrt(3), -1.0],
                    [1.0/math.sqrt(3), -1.0]],
                   [[1.0, -1.0/math.sqrt(3)],
                    [1.0, 1.0/math.sqrt(3)]],
                   [[-1.0/math.sqrt(3), 1.0],
                    [1.0/math.sqrt(3), 1.0]],
                   [[-1.0, -1.0/math.sqrt(3)],
                    [-1.0, 1/math.sqrt(3)]]])


    for line in traction(1,1).keys():
        for ele, side, l in mesh.boundary_elements:
            if l == line:
                for w in range(2):
                    mesh.basisFunction2D(gp[side, w])
                    mesh.eleJacobian(mesh.nodes_coord[mesh.ele_conn[ele, :]])

                    x1_o_e1e2, x2_o_e1e2 = mesh.mapping(ele)
                    t = traction(x1_o_e1e2, x2_o_e1e2)

                    Tele[0, ele] += mesh.phi[0]*t[l]*mesh.ArchLength[side]
                    Tele[1, ele] += mesh.phi[1]*t[l]*mesh.ArchLength[side]
                    Tele[2, ele] += mesh.phi[2]*t[l]*mesh.ArchLength[side]
                    Tele[3, ele] += mesh.phi[3]*t[l]*mesh.ArchLength[side]

    T = assemble1dof.globalVector(Tele, mesh)

    return T


def dirichlet_transient(K, B, mesh, temperature, time):
    """Apply Dirichlet BC.

    .. note::

        How its done:

        1. Loop over the lines where Dirichlet boundary conditions are
        applied. This is specified as key argument of the dictionary. The
        boundaries are defined on the gmsh nad tagged with a number.

        2. loop over all the nodes at the boundaries.

        3. If those nodes, identified as follows::

            boundary_nodes = [line node1 node2]

            are in the line where dirichlet BC were specified, change the
            stiffness matrix and the B vector based on this node index.

    Args:
        K (2nd order array): Stiffness matrix.
        B (1st order array): Vector with the load and traction.
        temperature (function): Function with 4 components.


    Returns:
        K (2nd order array), B (1st order array): Modified stiffness matrix
        and vector.


    """
    for line in temperature(1, 1, 1).keys():
        for n in range(len(mesh.boundary_nodes[:, 0])):
            if line == mesh.boundary_nodes[n, 0]:
                K[mesh.boundary_nodes[n, 1], :] = 0.0
                K[mesh.boundary_nodes[n, 2], :] = 0.0

                K[mesh.boundary_nodes[n, 1], mesh.boundary_nodes[n, 1]] =  1.0
                K[mesh.boundary_nodes[n, 2], mesh.boundary_nodes[n, 2]] = 1.0

                t1 = temperature(mesh.nodes_coord[mesh.boundary_nodes[n,1], 0],
                                mesh.nodes_coord[mesh.boundary_nodes[n, 1],
                                                 1], time)

                t2 = temperature(mesh.nodes_coord[mesh.boundary_nodes[n, 2], 0],
                                mesh.nodes_coord[mesh.boundary_nodes[n, 2],
                                                 1], time)

                B[mesh.boundary_nodes[n, 1]] = t1[line]
                B[mesh.boundary_nodes[n, 2]] = t2[line]

    return K, B



def neumann_transient(mesh, traction, time):
    """Apply Neumann BC.

    Computes the integral from the weak form with the boundary term.

    .. note::

        How its done:

        1. Define an array with the Gauss points for each path, each path will
        have 2 sets of two gauss points. One of the gp is fixed, which indicates
        that its going over an specific boundary of the element.

        Gauss points are defines as::

            gp = [gp, -1    1st path --> which has 2 possibilities for gp.
                  1,  gp    2nd path
                  gp,  1    3rd
                  -1,  gp]  4th

        2. A loop over the elements on the boundary extracting also the side
        where the boundary is located on this element.

    .. note::

        Edge elements are necessary because we need to extract the nodes
        from the connectivity. Then we can create a T for this element.

    Args:
        traction: Function with the traction and the line where the traction
            is applied.
        mesh: Object with the mesh attributes.

    Returns:
        T: Traction vector with size equals the dof.

    """
    Tele = np.zeros((4, mesh.num_ele))


    gp = np.array([[[-1.0/math.sqrt(3), -1.0],
                    [1.0/math.sqrt(3), -1.0]],
                   [[1.0, -1.0/math.sqrt(3)],
                    [1.0, 1.0/math.sqrt(3)]],
                   [[-1.0/math.sqrt(3), 1.0],
                    [1.0/math.sqrt(3), 1.0]],
                   [[-1.0, -1.0/math.sqrt(3)],
                    [-1.0, 1/math.sqrt(3)]]])


    for line in traction(1, 1, 1).keys():
        for ele, side, l in mesh.boundary_elements:
            if l == line:
                for w in range(2):
                    mesh.basisFunction2D(gp[side, w])
                    mesh.eleJacobian(mesh.nodes_coord[mesh.ele_conn[ele, :]])

                    x1_o_e1e2, x2_o_e1e2 = mesh.mapping(ele)
                    t = traction(x1_o_e1e2, x2_o_e1e2, time)

                    Tele[0, ele] += mesh.phi[0]*t[l]*mesh.ArchLength[side]
                    Tele[1, ele] += mesh.phi[1]*t[l]*mesh.ArchLength[side]
                    Tele[2, ele] += mesh.phi[2]*t[l]*mesh.ArchLength[side]
                    Tele[3, ele] += mesh.phi[3]*t[l]*mesh.ArchLength[side]

    T = assemble1dof.globalVector(Tele, mesh)

    return T

