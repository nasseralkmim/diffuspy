"""Creates an element object with basic attributes

"""
import numpy as np


class Element(object):
    """Build an Element base clase

    """
    def __init__(self, eid, model):
        self.eid = eid
        self.type = model.TYPE[eid]
        self.conn = model.CONN[eid]
        self.xyz = model.XYZ[self.conn]
        self.dof = model.DOF[eid]
        self.surf = model.surf_of_ele[eid]

        self.id_m = np.ix_(self.dof, self.dof)
        self.id_v = self.dof
