"""Module for quad element with 4 nodes - type 3 in gmsh

"""
from diffuspy.element import Element
import numpy as np


class Quad4(Element):
    """Constructor of a 4-node quadrangle (TYPE 3) element

    """
    def __init__(self, eid, model, material):

        super().__init__(eid, model)

        # Nodal coordinates in the natural domain (isoparametric coordinates)
        self.XEZ = np.array([[-1.0, -1.0],
                             [1.0, -1.0],
                             [1.0, 1.0],
                             [-1.0, 1.0]])

        # check if conductivity was assigned
        try:
            self.λ = material.λ[self.surf]
        except AttributeError:
            print('Conductivity (λ) not defined!')
        except KeyError:
            print('Surface ', self.surf,
                  ' with no conductivity (λ) assigned!')

        # check if capacitance material properties were assigned
        # if not, just pass because it maybe not a transient analysis
        try:
            self.ρ = material.ρ[self.surf]
            self.c = material.c[self.surf]
        except:
            pass

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

    @staticmethod
    def mapping(N, xyz):
        """maps from cartesian to isoparametric.

        """
        x1, x2 = N @ xyz
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

    def heat_stiffness_matrix(self, t=1):
        """Build the element heat (q) stiffness (k) matrix

        """
        k_q = np.zeros((4, 4))

        gauss_points = self.XEZ / np.sqrt(3.0)

        for gp in gauss_points:
            N, dN_ei = self.shape_function(xez=gp)
            dJ, dN_xi, _ = self.jacobian(self.xyz, dN_ei)

            B = dN_xi

            # Check if condutivity is a function
            if callable(self.λ) is True:
                x1, x2 = self.mapping(N, self.xyz)
                λ = self.λ(x1, x2, t)
            else:
                λ = self.λ

            k_q += λ * (B.T @ B) * dJ

        return k_q

    def heat_capacitance_matrix(self, t=1):
        """Build element matrix (k) due internal thermal energy storage (s)

        """
        k_s = np.zeros((4, 4))

        gauss_points = self.XEZ / np.sqrt(3.0)

        for gp in gauss_points:
            N, dN_ei = self.shape_function(xez=gp)
            dJ, dN_xi, _ = self.jacobian(self.xyz, dN_ei)

            # check if attribute and surface were assigned correctly
            try:
                # Check if specific heat is a function
                if callable(self.c) is True:
                    x1, x2 = self.mapping(N, self.xyz)
                    c = self.c(x1, x2, t)
                else:
                    c = self.c
            except AttributeError:
                print('Specific heat (c) not defined')
            except KeyError:
                print('Surface ', self.surf,
                      ' with no specific heat (c) assigned!')

            try:
                # Check if density is a function
                if callable(self.ρ) is True:
                    x1, x2 = self.mapping(N, self.xyz)
                    ρ = self.ρ(x1, x2, t)
                else:
                    ρ = self.ρ
            except AttributeError:
                print('Density (ρ) not defined')
            except KeyError:
                print('Surface ', self.surf,
                      ' with no density (ρ) assigned!')

            k_s += c*ρ*(np.atleast_2d(N).T @ np.atleast_2d(N))*dJ

        return k_s

    def heat_convection_matrix(self, h, t=1):
        """Build the element matrix (k) due convection boundary (c)

        """
        k_c = np.zeros((4, 4))

        gp = np.array([
            [[-1.0/np.sqrt(3), -1.0],
             [1.0/np.sqrt(3), -1.0]],
            [[1.0, -1.0/np.sqrt(3)],
             [1.0, 1.0/np.sqrt(3)]],
            [[-1.0/np.sqrt(3), 1.0],
             [1.0/np.sqrt(3), 1.0]],
            [[-1.0, -1.0/np.sqrt(3)],
             [-1.0, 1/np.sqrt(3)]]])

        # check if there is convection
        if h is not None:

            # loop for specified boundary conditions
            for key in h(1, 1).keys():
                line = key

                # loop over each boundary line that intersects the element
                # sides
                for ele_boundary_line, ele_side in zip(self.at_boundary_line,
                                                       self.side_at_boundary):
                    # Check if this element is at the line with convection bc
                    if line == ele_boundary_line:

                        # solve the integral with GQ
                        for w in range(2):
                            N, dN_ei = self.shape_function(xez=gp[ele_side, w])
                            _, _, arch_length = self.jacobian(self.xyz, dN_ei)

                            # check if condutance is a function
                            if callable(h) is True:
                                x1, x2 = self.mapping(N, self.xyz)
                                h_v = h(x1, x2, t)[line]
                            else:
                                h_v = h[line]

                            dL = arch_length[ele_side]

                            k_c += h_v * (
                                np.atleast_2d(N).T @ np.atleast_2d(N)) * dL

                    else:
                        # Catch element that is not at boundary
                        continue

        return k_c

    def heat_source_vector(self, σ_q=None, t=1, T_ip=1, dα=1, tol=1e-5):
        """Build the element vector due internal heat (q) source (σ)

        Args:
            T_ip: Nodal temperature form previous iteration i

        """
        gauss_points = self.XEZ / np.sqrt(3.0)

        pq = np.zeros(4)

        # find the average temperature of element
        # check if heat source is nonlinear in T
        if np.size(T_ip) > 1:
            T_avg = np.average(T_ip[self.conn])

        for gp in gauss_points:
            N, dN_ei = self.shape_function(xez=gp)
            dJ, dN_xi, _ = self.jacobian(self.xyz, dN_ei)

            x1, x2 = self.mapping(N, self.xyz)

            if σ_q is not None:
                if 'Reaction Degree' in σ_q.__defaults__:
                    pq[:] += N[:] * σ_q(x1, x2, t=t, dα=dα) * dJ
                elif 'Temperature' in σ_q.__defaults__:
                    pq[:] += N[:] * σ_q(x1, x2, t=t, T=T_avg) * dJ
                else:
                    pq[:] += N[:] * σ_q(x1, x2, t=t) * dJ
        return pq

    def heat_boundary_flux_vector(self, q_bc, t=1):
        """Build element load vector due q_bc boundary condition

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

        p_t = np.zeros(4)

        if q_bc is not None:

            # loop for specified boundary conditions
            for key in q_bc(1, 1).keys():
                line = key

                for ele_boundary_line, ele_side in zip(self.at_boundary_line,
                                                       self.side_at_boundary):
                    # Check if this element is at the line with traction
                    if line == ele_boundary_line:

                        # solve the integral with GQ
                        for w in range(2):
                            N, dN_ei = self.shape_function(xez=gp[ele_side, w])
                            _, _, arch_length = self.jacobian(self.xyz, dN_ei)

                            dL = arch_length[ele_side]
                            x1, x2 = self.mapping(N, self.xyz)

                            p_t[:] += N[:] * q_bc(x1, x2, t)[line] * dL

                    else:
                        # Catch element that is not at boundary
                        continue

        return p_t

    def heat_boundary_convection_vector(self, T_a, h, t=1):
        """Build the element heat vector due convection bc

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

        p_c = np.zeros(4)

        # Try compute the vector due convection
        if h is not None:

            # loop for specified boundary condition line
            for key in h(1, 1).keys():
                line = key

                for ele_boundary_line, ele_side in zip(self.at_boundary_line,
                                                       self.side_at_boundary):
                    # Check if this element is at the line with traction
                    if line == ele_boundary_line:

                        # solve the integral with GQ
                        for w in range(2):
                            N, dN_ei = self.shape_function(xez=gp[ele_side, w])
                            _, _, arch_length = self.jacobian(self.xyz, dN_ei)

                            dL = arch_length[ele_side]

                            # check if condutance is a function
                            if callable(h) is True:
                                x1, x2 = self.mapping(N, self.xyz)
                                h_v = h(x1, x2, t)[line]
                            else:
                                h_v = h[line]

                            # check if the surrounded fluid temperature is
                            # a function
                            if callable(T_a) is True:
                                x1, x2 = self.mapping(N, self.xyz)
                                T_a_v = T_a(x1, x2, t)[line]
                            else:
                                T_a_v = T_a[line]

                            p_c[:] += N[:] * h_v * T_a_v * dL

                    else:
                        # Catch element that is not at boundary
                        continue

        return p_c
