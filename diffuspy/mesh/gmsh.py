import numpy as np
import os
import re


def find_num(string):
    """Find all numbers in a string

    """
    # [+-]?(?:\d+(?:\.\d*)?|\.\d+) this one does't account for scientific notation
    num = re.findall(r'[-+]?\d+[\.]?\d*[eE]?[-+]?\d*', string)
    return num


class Parse(object):
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
        # element TYPE: [e1_type, e2_type ...]
        TYPE = []

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
                TYPE.append(int(num_list[1]))
                e_i += 1

            if len(num_list) == 7:
                nl = [int(f) - 1 for f in num_list[4:]]
                self.nodes_in_bound_line.append([nl[0],  nl[1], nl[2]])

        self.nodes_in_bound_line = np.array(self.nodes_in_bound_line)
        self.XYZ = np.array(list(XYZ.values()))
        self.CONN = np.array(list(CONN.values()))
        self.TYPE = np.array(TYPE)
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
