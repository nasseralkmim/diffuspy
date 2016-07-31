import numpy as np
import os
import re


def find_num(string):
    """Find all numbers in a string

    """
    num = re.findall(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)', string)
    return num


class Parse:
    """Parse the .geo and .msh file into dictionaries

    """
    def __init__(self, filename):
        geo_path = os.path.join(filename+'.geo')
        geo_file = open(geo_path, 'r')

        # physical_line_tag: line_tag
        self.physical_line = {}
        # physical_surf_tag: surface_tag
        self.physical_surf = {}
        # surf_tag: line_loop_tag
        self.surf = {}
        # line_loop_tag: [line1_tag line2_tag line3_tag ...]
        self.line_loop = {}
        # line_tag: [node1_tag node2_tag]
        self.line = {}

        for txt_line in geo_file:
            num_list = find_num(txt_line)

            if txt_line.startswith('Physical Line'):
                nl = [int(f) - 1 for f in num_list]
                self.physical_line[nl[0]] = nl[1]

            if txt_line.startswith('Plane Surface'):
                nl = [int(f) - 1 for f in num_list]
                self.surf[nl[0]] = nl[1]

            if txt_line.startswith('Physical Surface'):
                nl = [int(f) - 1 for f in num_list]
                self.physical_surf[nl[0]] = nl[1]

            if txt_line.startswith('Line('):
                nl = [int(f) - 1 for f in num_list]
                self.line[nl[0]] = nl[1:]

            if txt_line.startswith('Line Loop'):
                nl = [abs(int(f)) - 1 for f in num_list]
                self.line_loop[nl[0]] = nl[1:]

        msh_path = os.path.join(filename+'.msh')
        msh_file = open(msh_path, 'r')

        # node_tag: [node1 node2]
        XYZ = {}
        # element_tag: [node1_tag node2_tag node3_tag node4_tag]
        CONN = {}
        # element_tag: physical_surf_tag
        self.surf_of_ele = {}
        # [line_tag node1_tag node2_tag]
        self.nodes_in_bound_line = []

        e_i = 0
        for txt_line in msh_file:
            num_list = find_num(txt_line)

            # nodes coordinates xyz
            if len(num_list) == 4:
                n_tag = int(num_list[0]) - 1
                XYZ[n_tag] = [float(f) for f in num_list[1:3]]

            if len(num_list) == 9:
                conn = [int(f) - 1 for f in num_list[5:]]
                CONN[e_i] = conn
                self.surf_of_ele[e_i] = int(num_list[4]) - 1
                e_i += 1

            if len(num_list) == 7:
                nl = [int(f) - 1 for f in num_list[4:]]
                self.nodes_in_bound_line.append([nl[0],  nl[1], nl[2]])

        self.XYZ = np.array(list(XYZ.values()))
        self.CONN = np.array(list(CONN.values()))
        self.ne = len(CONN)
        self.nn = len(XYZ)

        # Nodal coordinates in the natural domain
        self.chi = np.array([[-1.0, -1.0],
                             [1.0, -1.0],
                             [1.0, 1.0],
                             [-1.0, 1.0]])

        # [ele side_of_ele_at_bound bound_line]
        bound_ele = []
        for e, conn in enumerate(self.CONN):
            for l, n1, n2 in self.nodes_in_bound_line:

                if np.all([n1, n2] == self.CONN[e, 0:2]):
                    bound_ele.append([e, 0, l])

                if np.all([n1, n2] == self.CONN[e, 1:3]):
                    bound_ele.append([e, 1, l])

                if np.all([n1, n2] == self.CONN[e, 2:4]):
                    bound_ele.append([e, 2, l])

                if np.all([n1, n2] == self.CONN[e, ::-3]):
                    bound_ele.append([e, 3, l])

        self.bound_ele = np.array(bound_ele)

    def basis_function(self, natural_coord):
        """Create the basis function and its properties.

        """
        # variables in the natural (iso-parametric) domain
        e1 = natural_coord[0]
        e2 = natural_coord[1]

        # Terms of the shape function
        e1_term = 0.5*(1.0 + self.chi[:, 0] * e1)
        e2_term = 0.5*(1.0 + self.chi[:, 1] * e2)

        # Basis functions
        # phi = [ phi_1 phi_2 phi_3 phi4 ]
        self.phi = e1_term*e2_term
        self.phi = np.array(self.phi)

        # Derivative of the shape functions
        # dphi = [ dphi1_e1 dphi2_e1 ...
        #         dphi1_e2 dphi2_e2 ... ]
        self.dphi_ei = np.zeros((2, 4))
        self.dphi_ei[0, :] = 0.5 * self.chi[:, 0] * e2_term
        self.dphi_ei[1, :] = 0.5 * self.chi[:, 1] * e1_term

    def mapping(self, xyz):
        """maps from cartesian to isoparametric.

        """
        x1, x2 = self.phi @ xyz

        return x1, x2

    def jacobian(self, element_nodes_coord):
        """Creates the Jacobian matrix of the mapping between an element

        """
        # Jac = [ x1_e1 x2_e1
        #         x1_e2 x2_e2]
        self.Jac = np.dot(self.dphi_ei, element_nodes_coord)

        self.detJac = ((self.Jac[0, 0]*self.Jac[1, 1] -
                        self.Jac[0, 1]*self.Jac[1, 0]))

        # JacInv = [ e1_x1 e2_x1
        #            e1_x2 e2_x2 ]
        self.JacInv = ((1.0 / self.detJac) *
                       np.array([[self.Jac[1, 1], -self.Jac[0, 1]],
                                 [-self.Jac[1, 0], self.Jac[0, 0]]]))
        self.JacInv = np.linalg.inv(self.Jac)

        self.dphi_xi2 = np.dot(self.JacInv, self.dphi_ei)

        # Using Chain rule,
        # phi_xi = phi_eI * eI_xi (2x8 array)
        self.dphi_xi = np.zeros((2, 4))
        self.dphi_xi[0, :] = (self.dphi_ei[0, :]*self.JacInv[0, 0] +
                              self.dphi_ei[1, :]*self.JacInv[0, 1])

        self.dphi_xi[1, :] = (self.dphi_ei[0, :]*self.JacInv[1, 0] +
                              self.dphi_ei[1, :]*self.JacInv[1, 1])

        # Length of the transofmation arch
        # Jacobian for line integral-2.
        self.ArchLength = np.array([
            (self.Jac[0, 0]**2. + self.Jac[0, 1]**2.)**(1./2.),
            (self.Jac[1, 0]**2. + self.Jac[1, 1]**2.)**(1./2.),
            (self.Jac[0, 0]**2. + self.Jac[0, 1]**2.)**(1./2.),
            (self.Jac[1, 0]**2. + self.Jac[1, 1]**2.)**(1./2.)
        ])
