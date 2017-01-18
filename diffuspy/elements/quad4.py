"""Module for quad element with 4 nodes - type 3 in gmsh

"""
from diffuspy.element import Element
import numpy as np


class Quad4(Element):
    """Constructor of a 4-node quadrangle (TYPE 3) element

    """
    def __init__(self, eid, model, material, EPS0):

        super().__init__(eid, model)

        # Nodal coordinates in the natural domain (isoparametric coordinates)
        self.XEZ = np.array([[-1.0, -1.0],
                             [1.0, -1.0],
                             [1.0, 1.0],
                             [-1.0, 1.0]])

        try:
            self.ρ = material.ρ[self.surf]
            self.c = material.c[self.surf]
            self.λ = material.λ[self.surf]
        except AttributeError:
            print('Check if material properties were assigned for all surfaces!')
        except KeyError:
            print('Surface ', self.surf,
                  ' with no material assigned! (Default used)')

        # check if its a boundary element
        if eid in model.bound_ele[:, 0]:

            # index where bound_ele refers to this element
            index = np.where(model.bound_ele[:, 0] == eid)[0]

            # side of the element at the boundary
            self.side_at_boundary = model.bound_ele[index, 1]

            # boundary line where the element side share interface
            self.at_boundary_line = model.bound_ele[index, 2]

        else:
            self.side_at_boundary = []
            self.at_boundary_line = []

    def shape_function(self, xez):
        """Create the basis function and evaluate them at xez coordinates

        Args:
            xez (array): position in the isoparametric coordinate xi, eta, zeta

        Return:
            N (array): shape functions

        """
        # variables in the natural (iso-parametric) domain
        e1 = xez[0]
        e2 = xez[1]

        # Terms of the shape function
        e1_term = 0.5*(1.0 + self.XEZ[:, 0] * e1)
        e2_term = 0.5*(1.0 + self.XEZ[:, 1] * e2)

        # Basis functions
        # N = [ N_1 N_2 N_3 N_4 ]
        N = e1_term*e2_term
        self.N = np.array(N)

        # Derivative of the shape functions
        # dN = [ dN1_e1 dN2_e1 ...
        #         dN1_e2 dN2_e2 ... ]
        self.dN_ei = np.zeros((2, 4))
        self.dN_ei[0, :] = 0.5 * self.XEZ[:, 0] * e2_term
        self.dN_ei[1, :] = 0.5 * self.XEZ[:, 1] * e1_term

        return self.N, self.dN_ei

    def mapping(self, xyz):
        """maps from cartesian to isoparametric.

        """
        x1, x2 = self.N @ xyz
        return x1, x2

    def jacobian(self, xyz, dN_ei):
        """Creates the Jacobian matrix of the mapping between an element

        Args:
        xyz (array of floats): coordinates of element nodes in cartesian
        coordinates
        dN_ei (array of floats): derivative of shape functions

        Return:
        det_jac (float): determinant of the jacobian matrix
        dN_xi (array of floats): derivative of shape function
        with respect to cartesian system
        arch_length (array of floats): arch length for change of variable
        in the line integral
        """
        # Jac = [ x1_e1 x2_e1
        #         x1_e2 x2_e2 ]
        Jac = dN_ei @ xyz

        det_jac = abs((Jac[0, 0]*Jac[1, 1] -
                       Jac[0, 1]*Jac[1, 0]))

        # jac_inv = [ e1_x1 e2_x1
        #            e1_x2 e2_x2 ]
        jac_inv = np.linalg.inv(Jac)

        # Using Chain rule,
        # N_xi = N_eI * eI_xi (2x8 array)
        dN_xi = np.zeros((2, 4))
        dN_xi[0, :] = (dN_ei[0, :]*jac_inv[0, 0] +
                       dN_ei[1, :]*jac_inv[0, 1])

        dN_xi[1, :] = (dN_ei[0, :]*jac_inv[1, 0] +
                       dN_ei[1, :]*jac_inv[1, 1])

        # Length of the transofmation arch
        # Jacobian for line integral-2.
        arch_length = np.array([
            (Jac[0, 0]**2 + Jac[0, 1]**2)**(1/2),
            (Jac[1, 0]**2 + Jac[1, 1]**2)**(1/2),
            (Jac[0, 0]**2 + Jac[0, 1]**2)**(1/2),
            (Jac[1, 0]**2 + Jac[1, 1]**2)**(1/2)
        ])
        return det_jac, dN_xi, arch_length

    def stiffness_matrix(self, t=1):
        """Build the element stiffness matrix

        """
        k = np.zeros((8, 8))

        C = self.c_matrix(t)

        gauss_points = self.XEZ / np.sqrt(3.0)

        for gp in gauss_points:
            _, dN_ei = self.shape_function(xez=gp)
            dJ, dN_xi, _ = self.jacobian(self.xyz, dN_ei)

            B = np.array([
                [dN_xi[0, 0], 0, dN_xi[0, 1], 0, dN_xi[0, 2], 0,
                 dN_xi[0, 3], 0],
                [0, dN_xi[1, 0], 0, dN_xi[1, 1], 0, dN_xi[1, 2], 0,
                 dN_xi[1, 3]],
                [dN_xi[1, 0], dN_xi[0, 0], dN_xi[1, 1], dN_xi[0, 1],
                 dN_xi[1, 2], dN_xi[0, 2], dN_xi[1, 3], dN_xi[0, 3]]])

            k += (B.T @ C @ B)*dJ

        return k

    def mass_matrix(self, t=1):
        """Build element mass matrix

        """
        return None
        

    def c_matrix(self, t=1):
        """Build the element constitutive matrix

        """
        self.C = np.zeros((3, 3))
        self.C[0, 0] = 1.0
        self.C[1, 1] = 1.0
        self.C[1, 0] = self.nu
        self.C[0, 1] = self.nu
        self.C[2, 2] = (1.0 - self.nu)/2.0
        self.C = (self.E/(1.0 - self.nu**2.0))*self.C

        return self.C

    def load_body_vector(self, b_force, t=1):
        """Build the element vector due body forces b_force

        """
        gauss_points = self.XEZ / np.sqrt(3.0)

        pb = np.zeros(8)
        for gp in gauss_points:
            N, dN_ei = self.shape_function(xez=gp)
            dJ, dN_xi, _ = self.jacobian(self.xyz, dN_ei)

            x1, x2 = self.mapping(self.xyz)

            pb[0] += N[0]*b_force(x1, x2, t)[0]*dJ
            pb[1] += N[0]*b_force(x1, x2, t)[1]*dJ
            pb[2] += N[1]*b_force(x1, x2, t)[0]*dJ
            pb[3] += N[1]*b_force(x1, x2, t)[1]*dJ
            pb[4] += N[2]*b_force(x1, x2, t)[0]*dJ
            pb[5] += N[2]*b_force(x1, x2, t)[1]*dJ
            pb[6] += N[3]*b_force(x1, x2, t)[0]*dJ
            pb[7] += N[3]*b_force(x1, x2, t)[1]*dJ

        return pb

    def load_strain_vector(self, t=1):
        """Build the element vector due initial strain

        """
        C = self.c_matrix(t)

        gauss_points = self.XEZ / np.sqrt(3.0)

        pe = np.zeros(8)
        for gp in gauss_points:
            _, dN_ei = self.shape_function(xez=gp)
            dJ, dN_xi, _ = self.jacobian(self.xyz, dN_ei)

            B = np.array([
                [dN_xi[0, 0], 0, dN_xi[0, 1], 0, dN_xi[0, 2], 0,
                 dN_xi[0, 3], 0],
                [0, dN_xi[1, 0], 0, dN_xi[1, 1], 0, dN_xi[1, 2], 0,
                 dN_xi[1, 3]],
                [dN_xi[1, 0], dN_xi[0, 0], dN_xi[1, 1], dN_xi[0, 1],
                 dN_xi[1, 2], dN_xi[0, 2], dN_xi[1, 3], dN_xi[0, 3]]])

            pe += (B.T @ C @ self.eps0)*dJ

        return pe

    def load_traction_vector(self, traction_bc, t=1):
        """Build element load vector due traction_bction boundary condition

        """
        gp = np.array([
            [[-1.0/np.sqrt(3), -1.0],
             [1.0/np.sqrt(3), -1.0]],
            [[1.0, -1.0/np.sqrt(3)],
             [1.0, 1.0/np.sqrt(3)]],
            [[-1.0/np.sqrt(3), 1.0],
             [1.0/np.sqrt(3), 1.0]],
            [[-1.0, -1.0/np.sqrt(3)],
             [-1.0, 1/np.sqrt(3)]]])

        pt = np.zeros(8)

        # loop for specified boundary conditions
        for key in traction_bc(1, 1).keys():
            line = key[1]

            for ele_boundary_line, ele_side in zip(self.at_boundary_line,
                                                   self.side_at_boundary):
                # Check if this element is at the line with traction
                if line == ele_boundary_line:

                    # perform the integral with GQ
                    for w in range(2):
                        N, dN_ei = self.shape_function(xez=gp[ele_side, w])
                        _, _, arch_length = self.jacobian(self.xyz, dN_ei)
                        
                        dL = arch_length[ele_side]
                        x1, x2 = self.mapping(self.xyz)

                        pt[0] += N[0] * traction_bc(x1, x2, t)[key][0] * dL
                        pt[1] += N[0] * traction_bc(x1, x2, t)[key][1] * dL
                        pt[2] += N[1] * traction_bc(x1, x2, t)[key][0] * dL
                        pt[3] += N[1] * traction_bc(x1, x2, t)[key][1] * dL
                        pt[4] += N[2] * traction_bc(x1, x2, t)[key][0] * dL
                        pt[5] += N[2] * traction_bc(x1, x2, t)[key][1] * dL
                        pt[6] += N[3] * traction_bc(x1, x2, t)[key][0] * dL
                        pt[7] += N[3] * traction_bc(x1, x2, t)[key][1] * dL

                else:
                    # Catch element that is not at boundary
                    continue

        return pt
