import numpy as np
import os


class parse:
    """Produce mesh with gmsh file.

    .. note ::

        parse: analyze into its parts.


    .. note::

        In order to create a mesh for the code follow the steps:

        1. Draw the geometry

            1.1 Draw Points.

            1.2 Connect these points with straight lines.

            1.3 Define a plane surface.

        2. Add physical groups where the BC are going to be applied

            2.1 Add lines in order to get the nodes in that line.

            2.2 Add a Surface in order to get the connectivity.

        3. Change the subdivision algorithm to "all quads" ontools-mesh-general tab.

        4. Click on 2D to create the mesh. It should contain quad elements only.

        5. Save the .msh file with the mesh with the same name as the .geo.


    Args:
        filename (str): Name of the file with the mesh geometry and properties.

    Attributes:
        surfaces (array of int): Index that identifies the surfaces assigned
            in the "gmsh.geo" file.
        boundary_lines (array of int): Index of lines that define the
            boundary globally, i.e, whole geometry.
        nodes_coord (array of flaot): Nodes of the mesh in the format.
            ::
                nodes_coord = [ coordinate x1 of node 0, coordinade x2 of node 0]
        ele_conn (array of int): Nodes that form an quad element in the format
            ::
                ele_conn = [node1, node2, node3, node4]
        ele_surface (array of int): element index and the physical surface
            where it is. The format is::

                ele_surface = [element index, surface index]
        boundary_nodes (array of float): Nodes at the boundary in the format.
            ::

                boundary_nodes = [boundary line, node1, node2]
        boundary_elements (array of int): Elements at the boundary and the
            boundary, which could be 0->first two points on the connectivity;
            1->2nd 2points; 2-> 3rd 2 points; 3-> 4th 2 points. Example::

                ele_conn = [1, 5, 13, 10]
                element boundary 0 = [1, 5]
                element boundary 1 = [5, 13]
                element boundary 2 = [13, 10]
                element boundary 3 = [10, 1]

            the format of this array is::

                boundary_elements = [ele index, element boundary, boundary line]

    """
    def __init__(self, filename):
        path = os.path.join("..\mesh", filename+'.geo')
        geometry_feeder = open(path, 'r')

        surfaces = []
        boundaryLines = []

        self.physicalLine = {}
        self.physicalSurface = {}
        self.lineLoop = {}
        self.line = {}
        #counters for line tag
        i = 0
        for txtLine in geometry_feeder:
            txtLine = txtLine.strip()

            if txtLine.startswith('Physical Line('):
                boundaryLines.append(
                    [int(txtLine[txtLine.find('{')+1:txtLine.find('}')]) - 1, i])
                i += 1

            columns = txtLine.split()

            # physical line means boundary line - [pL tag, line number tag]
            if txtLine.startswith('Physical Line('):
                plTag = int(txtLine[txtLine.find('(')+1:txtLine.find(')')]) - 1

                t = columns[3]
                number = int(t[t.find('{')+1:t.find('}')]) - 1
                self.physicalLine[plTag] = number


            if txtLine.startswith('Physical Surface('):
                surfaces.append(
                    int(txtLine[txtLine.find('(')+1:txtLine.find(')')]) - 1
                )

            # line tag: node 1 node 2
            if txtLine.startswith('Line('):
                lTag = int(txtLine[txtLine.find('(')+1:txtLine.find(')')]) - 1

                tf = columns[2]
                nf = int(tf[tf.find('{')+1:tf.find(',')]) - 1

                tl = columns[3]
                nl = int(tl[:tf.find('}')-1]) - 1

                self.line[lTag] = [nf, nl]

            #line loop tag, line1 line2 ....
            if txtLine.startswith('Line Loop('):
                lpTag = int(txtLine[txtLine.find('(')+1:txtLine.find(')')])
                lpList = []

                # Get the first entry
                tf = columns[3]
                pf = 1
                if tf[1] == '-':
                    pf = 2
                nf= int(tf[tf.find('{')+pf:tf.find(',')]) - 1
                lpList.append(nf)

                for i in range(4,len(columns)-1):
                    t = columns[i]
                    if columns[i][0] == '-':
                        print(t[-1])
                        n = int(t[1:t.find(',')]) - 1
                        lpList.append(n)
                    else:
                        n = int(t[:t.find(',')]) - 1
                        lpList.append(n)

                tl = columns[-1]
                pl = 0
                if tl[0] == '-':
                    pl = 1
                nl = int(tl[pl:tl.find('}')]) - 1
                lpList.append(nl)

                self.lineLoop[lpTag] = lpList

            if txtLine.startswith('Physical Surface('):
                psTag = int(txtLine[txtLine.find('(')+1:txtLine.find(')')]) - 1

                t = columns[3]
                number = int(t[t.find('{')+1:t.find('}')]) - 1
                self.physicalSurface[psTag] = number

        self.boundaryLines = np.asarray(boundaryLines)
        self.surfaces = np.asarray(surfaces)


        path2 = os.path.join("..\mesh", filename+'.msh')
        mesh_feeder = open(path2, 'r')

        nodes_coord = []
        ele_conn = []
        boundary_nodes = []
        ele_surface = []
        i=0
        for line in mesh_feeder:

            line = line.strip()
            columns = line.split()

            if len(columns) == 4:
                nodes_coord.append([float(columns[1]),
                                    float(columns[2])])

            # -1 because python starts lists with 0 index.
            if len(columns) == 9:
                ele_conn.append([int(columns[5]) - 1,
                                 int(columns[6]) - 1,
                                 int(columns[7]) - 1,
                                 int(columns[8]) - 1])

                ele_surface.append([i, int(columns[3]) - 1])
                i += 1

            # boundary nodes = [line, node 1, node 2]
            if len(columns) == 7:
                boundary_nodes.append([int(columns[4]) - 1,
                                       int(columns[5]) - 1,
                                       int(columns[6]) - 1])



        self.boundary_nodes = np.asarray(boundary_nodes)
        self.nodes_coord = np.asarray(nodes_coord)
        self.ele_conn = np.asarray(ele_conn)
        self.ele_surface = np.asarray(ele_surface)

        self.num_ele = len(self.ele_conn[:, 0])
        self.num_nodes = len(self.nodes_coord[:, 0])

        boundary_elements = []

        for ele in range(len(self.ele_conn[:, 0])):
            for no in range(len(self.boundary_nodes[:, 0])):
                if np.all(self.boundary_nodes[no, 1:3] == self.ele_conn[ele,
                                                          0:2]):
                    boundary_elements.append([ele, 0, self.boundary_nodes[no,
                    0]])

                if np.all(self.boundary_nodes[no, 1:3] == self.ele_conn[ele,
                                                          1:3]):
                    boundary_elements.append([ele, 1, self.boundary_nodes[no,
                    0]])

                if np.all(self.boundary_nodes[no, 1:3] == self.ele_conn[ele,
                                                          2:4]):
                    boundary_elements.append([ele, 2, self.boundary_nodes[no,
                    0]])

                if np.all(self.boundary_nodes[no, 1:3] == self.ele_conn[ele,
                                                          ::-3]):
                    boundary_elements.append([ele, 3, self.boundary_nodes[no,
                    0]])

        self.boundary_elements = np.asarray(boundary_elements)

        self.AvgLength = (self.nodes_coord[1, 0] - self.nodes_coord[0, 0])/30.

        # Nodal coordinates in the natural domain
        self.chi = np.array([[-1.0, -1.0],
                             [1.0, -1.0],
                             [1.0, 1.0],
                             [-1.0, 1.0]])

        self.gmsh = 1.0


    def basisFunction2D(self, natural_coord):
        """Create the basis function and its properties.

        The function will be evaluated at the quadrature points that will be
        passed as a natural coordinate.

        Args:
            nodes_cooord (float, float): Nodes coordinates pair in the real
                domain.
            natural_coord (flaot, float): Variables in the natural (
                iso-parametric domain) [e1,e2].

        Attributes:
            phi (float): Function of the natural coordinates with the format.
                ::
                    phi = [phi1, phi2, phi3, phi4]
            dphi (float): Derivative of the shape functions.
                ::
                    dphi = [[dphi1_e1, dphi2_e1, dphi3_e1, dphi4_e1 ],
                        [dphi1_e2, dphi2_e2, dphi3_e2, dphi4_e2 ]]
            Jac (float): Jacobian of the transformation as a function


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
        self.phi = np.asarray(self.phi)

        # Derivative of the shape functions
        # dphi = [ dphi1_e1 dphi2_e1 ...
        #         dphi1_e2 dphi2_e2 ... ]
        self.dphi_ei = np.zeros((2, 4))
        self.dphi_ei[0, :] = 0.5 * self.chi[:, 0] * e2_term
        self.dphi_ei[1, :] = 0.5 * self.chi[:, 1] * e1_term


    def mapping(self, e):
        """maps from cartesian to isoparametric.

        Args:
            e (int): Designated the element in which the calculation is made.

        Returns:
            x1_o_e1e2 , x2_o_e1e2 (float, float) : first and second coordinate
            as a function of [e1, e2].

        """
        x1_o_e1e2 = np.dot(
            self.phi[:], self.nodes_coord[self.ele_conn[e, :], 0])

        x2_o_e1e2 = np.dot(
            self.phi[:], self.nodes_coord[self.ele_conn[e, :], 1])

        return x1_o_e1e2, x2_o_e1e2


    def eleJacobian(self, element_nodes_coord):
        """Creates the Jacobian matrix of the mapping between an element
        and its natural form.

        Args:
            element_nodes_coord ([float, float]: Coordinates of an individual
            element.

        Returns:
            Jac (float): jacobian matrix 2x2
            detJac (float): determinant of the jacobian matrix
            dphi_xi (float): derivative of the shape functions with respect
            the cartesian coordinates

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

