import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import networkx as nx


def elements_label(model, ax):

    for e in range(len(model.CONN)):
        x_element = (model.XYZ[model.CONN[e, 0], 0] +
                     model.XYZ[model.CONN[e, 1], 0] +
                     model.XYZ[model.CONN[e, 2], 0] +
                     model.XYZ[model.CONN[e, 3], 0])/4.

        y_element = (model.XYZ[model.CONN[e, 0], 1] +
                     model.XYZ[model.CONN[e, 1], 1] +
                     model.XYZ[model.CONN[e, 2], 1] +
                     model.XYZ[model.CONN[e, 3], 1])/4.

        ax.annotate(str(e), (x_element, y_element), size=9,
                    color='r')


def surface_label(model, ax):
    # initalize the surface count
    i = 0
    for phy_surf_tag, surf_tag in model.physical_surf.items():
        # line loop tag
        ll_tag = model.surf[surf_tag]

        xm = 0.0
        ym = 0.0
        for node in model.line_loop[ll_tag]:
            xm += (model.XYZ[model.line[node][0], 0] +
                   model.XYZ[model.line[node][1], 0])
            ym += (model.XYZ[model.line[node][0], 1] +
                   model.XYZ[model.line[node][1], 1])

        xs = xm/(2*len(model.line_loop[ll_tag]))
        ys = ym/(2*len(model.line_loop[ll_tag]))

        ax.annotate(str(i), (xs, ys), size=9, color='g')
        i += 1


def edges_label(model, ax, font_size=6):
    c = model.XYZ

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G.add_node(i, posxy=(X[i], Y[i]))

    bound_middle = {}
    iant = model.nodes_in_bound_line[0][0]

    cont = 0
    for l, n1, n2 in model.nodes_in_bound_line:
        if l == iant:
            cont += 1
            bound_middle[l] = cont
        else:
            cont = 1
        iant = l

    edge_labels = {}
    cont = 0
    for l, n1, n2 in model.nodes_in_bound_line:
        cont += 1
        if cont == int(bound_middle[l]/2.0):
            edge_labels[n1, n2] = str(l)
        if cont == bound_middle[l]:
            cont = 0

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_edge_labels(G, positions, edge_labels, label_pos=0.5, font_color='b',
                                 ax=ax, font_size=font_size)


def nodes_label(model, ax, font_size=6, node=None):
    c = model.XYZ

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    if node is None:
        label = {}
        for i in range(len(X)):
            label[i] = i
            G.add_node(i, posxy=(X[i], Y[i]))
    else:
        label = {}
        label[0] = node
        G.add_node(0, posxy=(X[node], Y[node]))

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_nodes(G, positions, node_color='w',
                           node_shape='s', ax=ax, node_size=70)
    nx.draw_networkx_labels(G, positions, label, ax=ax, font_size=font_size)


def domain(model, ax, color='k'):
    """Draw domain region

    """
    c = model.XYZ

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G.add_node(i, posxy=(X[i], Y[i]))

    for pl_tag, l_tag in model.physical_line.items():
        n = model.line[l_tag]
        G.add_edge(n[0], n[1])

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_edges(G, positions, edge_color=color, node_size=0,
                           font_size=0, width=1, ax=ax)


def elements(model, ax, color='k'):
    c = model.XYZ

    X, Y = c[:, 0], c[:, 1]

    G2 = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G2.add_node(i, posxy=(X[i], Y[i]))

    for n in model.CONN:
        G2.add_cycle(n)

    positions = nx.get_node_attributes(G2, 'posxy')

    nx.draw_networkx_edges(G2, positions, node_size=0, edge_color=color,
                           width=1, ax=ax)


def surface(model, ax, color='k'):
    c = model.XYZ

    X, Y = c[:, 0], c[:, 1]

    G2 = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G2.add_node(i, posxy=(X[i], Y[i]))

    for phy_surf_tag, surf_tag in model.physical_surf.items():
        # line loop tag
        ll_tag = model.surf[surf_tag]

        for l_tag in model.line_loop[ll_tag]:
            n1, n2 = model.line[l_tag]
            G2.add_edge(n1, n2)

    positions = nx.get_node_attributes(G2, 'posxy')
    nx.draw_networkx_edges(G2, positions, node_size=0, edge_color=color,
                           font_size=0, width=1, ax=ax)


def tricontourf(model, U, ax, cmap='plasma', lev=10, cl=True, vmin=None,
                vmax=None, font_size=12):
    """Plot contour with the tricoutour function and the boundary line with
    the boundary node.

    """
    c = model.XYZ

    bn = np.array(model.nodes_in_bound_line)

    xx, yy, zz = c[:, 0], c[:, 1], U

    # ccx = np.append(c[bn[:, 2], 0], c[bn[0, 1], 0])
    # ccy = np.append(c[bn[:, 2], 1], c[bn[0, 1], 1])

    triangles = []
    for n1, n2, n3, n4 in model.CONN:
        triangles.append([n1, n2, n3])
        triangles.append([n1, n3, n4])

    triangles = np.asarray(triangles)

    CS2 = ax.tricontourf(xx, yy, triangles, zz, lev=lev,
                         origin='lower',
                         cmap=cmap, antialiased=True,
                         vmin=vmin, vmax=vmax)

    if cl is True:
        CS3 = ax.tricontour(xx, yy, triangles, zz, lev=lev, colors='k',
                            vmin=vmin, vmax=vmax)
        ax.clabel(CS3, colors='k', fmt='%1.1f')

    # ax.plot(ccx, ccy, '-k')

    return CS2
