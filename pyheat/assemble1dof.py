import numpy as np


def globalMatrix(K_ele, mesh):
    K = np.zeros((mesh.num_nodes, mesh.num_nodes))

    for e in range(mesh.num_ele):
        for i in range(len(mesh.ele_conn[e])):
            for j in range(len(mesh.ele_conn[e])):
                K[mesh.ele_conn[e, i],
                  mesh.ele_conn[e, j]] += K_ele[i, j, e]

    return K


def globalVector(V_ele, mesh):
    V = np.zeros((mesh.num_nodes, 1))

    for e in range(mesh.num_ele):
        for i in range(len(mesh.ele_conn[e])):
            V[mesh.ele_conn[e, i]] += V_ele[i, e]
    return V

def globalVectorAverage(V_ele, mesh):
    V = np.zeros((mesh.num_nodes, 1))

    for e in range(mesh.num_ele):
        for i in range(len(mesh.ele_conn[e])):
            V[mesh.ele_conn[e, i]] += V_ele[i, e]

    return V

def globalVectorAverage2(V_ele, mesh):
    V = np.zeros((mesh.num_nodes, 1))

    for e in range(mesh.num_ele):
        for i in range(len(mesh.ele_conn[e])):
            V[mesh.ele_conn[e, i]] += V_ele[i, e]/((mesh.ele_conn ==
                                                    mesh.ele_conn[e, i]).sum())

    return V
