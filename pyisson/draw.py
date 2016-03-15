import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import networkx as nx


def elements_label(mesh, ax):

    for e in range(len(mesh.ele_conn)):
        x_element = (mesh.nodes_coord[mesh.ele_conn[e, 0], 0] +
                     mesh.nodes_coord[mesh.ele_conn[e, 1], 0] +
                     mesh.nodes_coord[mesh.ele_conn[e, 2], 0] +
                     mesh.nodes_coord[mesh.ele_conn[e, 3], 0])/4.

        y_element = (mesh.nodes_coord[mesh.ele_conn[e, 0], 1] +
                     mesh.nodes_coord[mesh.ele_conn[e, 1], 1] +
                     mesh.nodes_coord[mesh.ele_conn[e, 2], 1] +
                     mesh.nodes_coord[mesh.ele_conn[e, 3], 1])/4.

        ax.annotate(str(e), (x_element, y_element), size=9,
                    color='r')


def surface_label(mesh, ax):

    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]

    G2 = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G2.add_node(i, posxy=(X[i], Y[i]))

    # adjust the node numbering order of an element
    if mesh.gmsh == 1.0:
        temp = np.copy(mesh.ele_conn[:, 2])
        mesh.ele_conn[:, 2] = mesh.ele_conn[:, 3]
        mesh.ele_conn[:, 3] = temp
        mesh.gmsh += 1

    i = 0
    for surface, lpTag in mesh.physicalSurface.items():
        xm = 0.0
        ym = 0.0
        for node in mesh.lineLoop[lpTag]:
            G2.add_edge(mesh.line[node][0],
                        mesh.line[node][1])
            xm += (mesh.nodes_coord[mesh.line[node][0], 0] +
                   mesh.nodes_coord[mesh.line[node][1], 0])
            ym += (mesh.nodes_coord[mesh.line[node][0], 1] +
                   mesh.nodes_coord[mesh.line[node][1], 1])

        xs, ys = xm/(2*len(mesh.lineLoop[lpTag])), ym/(2*len(mesh.lineLoop[
                                                                 lpTag]))
        ax.annotate(str(i), (xs, ys), size=9, color='g')
        i += 1


def edges_label(mesh, ax):
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G.add_node(i, posxy=(X[i], Y[i]))

    bound_middle = {}
    iant = mesh.boundary_nodes[0, 0]

    cont = 0
    for i, e1, e2 in mesh.boundary_nodes:
        if i == iant:
            cont += 1
            bound_middle[i] = cont
        else:
            cont = 1
        iant = i

    edge_labels = {}
    cont = 0
    for i, e1, e2 in mesh.boundary_nodes:
        cont += 1
        if cont == int(bound_middle[i]/2.0):
            edge_labels[e1, e2] = str(i)
        if cont == bound_middle[i]:
            cont = 0

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_edge_labels(G, positions, edge_labels, label_pos=0.5,
                                 font_size=9, font_color='b', ax=ax)


def nodes_label(mesh, ax):
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    label = {}
    for i in range(len(X)):
        label[i] = i
        G.add_node(i, posxy=(X[i], Y[i]))

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_nodes(G, positions, node_color='w', node_size=140,
                           node_shape='s', ax=ax)
    nx.draw_networkx_labels(G, positions, label, font_size=9, ax=ax)


def domain(mesh, ax, color):
    """Draw domain region

    """
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G.add_node(i, posxy=(X[i], Y[i]))

    for plTag, lineTag in mesh.physicalLine.items():
        lineNodes = mesh.line[lineTag]
        G.add_edge(lineNodes[0], lineNodes[1])

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_edges(G, positions, edge_color=color,
                           font_size=0, width=1, origin='lower', ax=ax)


def elements(mesh, ax, color):

    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]

    G2 = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G2.add_node(i, posxy=(X[i], Y[i]))

    if mesh.gmsh == 1.0:
        temp = np.copy(mesh.ele_conn[:, 2])
        mesh.ele_conn[:, 2] = mesh.ele_conn[:, 3]
        mesh.ele_conn[:, 3] = temp
        mesh.gmsh += 1

    for i in range(len(mesh.ele_conn)):
        G2.add_cycle([mesh.ele_conn[i, 0],
                     mesh.ele_conn[i, 1],
                     mesh.ele_conn[i, 3],
                     mesh.ele_conn[i, 2]], )

    edge_line_nodes = {}
    for i in range(len(mesh.boundary_nodes[:, 0])):
        edge_line_nodes[(mesh.boundary_nodes[i, 1], mesh.boundary_nodes[i,
                                                                        2])] \
            = mesh.boundary_nodes[i, 0]

    positions = nx.get_node_attributes(G2, 'posxy')

    nx.draw_networkx(G2, positions, node_size=0, edge_color=color,
                     font_size=0,  width=1)


def tricontourf(mesh, U, ax, cmap, lev):
    """Plot contour with the tricoutour function and the boundary line with
    the boundary node.

    """
    c = mesh.nodes_coord

    bn = mesh.boundary_nodes

    xx, yy, zz = c[:, 0], c[:, 1], U

    ccx = np.append(c[bn[:, 1], 0], c[bn[0, 1], 0])
    ccy = np.append(c[bn[:, 1], 1], c[bn[0, 1], 1])

    triangles = []
    for n1, n2, n3, n4 in mesh.ele_conn:
        triangles.append([n1, n2, n3])
        triangles.append([n1, n3, n4])

    triangles = np.asarray(triangles)

    CS2 = ax.tricontourf(xx, yy, triangles, zz, lev=lev,
                         origin='lower',
                         cmap=cmap, antialiased=True)

    CS3 = ax.tricontour(xx, yy, triangles, zz, lev, colors='black')
    ax.clabel(CS3, fontsize=8, colors='k', fmt='%1.1f')

    ax.patch.set_visible(False)
    ax.axis('off')

    ax.plot(ccx, ccy, '-k')

    plt.colorbar(CS2)
