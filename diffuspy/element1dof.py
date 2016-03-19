import numpy as np
import math


class Matrices:
    """Build the elemental matrices.

    Creates an object that has as attributes the elemental matrices.

    """
    def __init__(self, objectmesh):
        self.mesh = objectmesh


    def stiffness(self, mat):
        """Build the elemental stiffness matrix.

        Runs over each individual element properties when the object methods
        are called. The object is the mesh and its methods are basisFunction2D
        which defines the phi function and its derivative;  eleJacobian which
        produce a Jacobian matrix and its determinant and also the derivative
        of the basis function with respect the spatial coordinates.

        .. note::

            How it works:

            1. Loop over elements.

            2. loop over the 4 possible combination of 2 GP for each direction.

            3. Call the methods to create the derivative of shape functions  and Jacobian for each GP combination.

            4. Build the stiffness matrix from a matrix multiplication.

        Gauss points from natural nodal coordinates::

             gp = [ [-1, -1]
                    [ 1, -1]
                    [ 1,  1]
                    [-1,  1] ]/ sqrt(3)


        Args:
            mesh: object that includes the mesh atributes and methods for basis
                functions and transformation jacobians.
            k (function) : material properties for the constitutive relation.

        Return:
            K_ele (array of float): 3nd order array 4x4 with the elemental
            stiffness matrix.

        """
        mesh = self.mesh

        self.gp = mesh.chi / math.sqrt(3.0)

        self.K = np.zeros((4, 4, mesh.num_ele))

        for surface, [k, c, rho] in mat.items():
            for e in range(mesh.num_ele):
                if surface == mesh.ele_surface[e, 1]:
                    for w in range(4):
                        mesh.basisFunction2D(self.gp[w])
                        mesh.eleJacobian(mesh.nodes_coord[
                            mesh.ele_conn[e, :]])

                        x1_o_e1e2, x2_o_e1e2 = mesh.mapping(e)

                        B = mesh.dphi_xi

                        self.K[:, :, e] += k*(np.dot(np.transpose(B), B) *
                                            mesh.detJac)


    def load(self, q):
        """Build the load vector thermal diffusivity.

        Args:
            mesh: object that includes the mesh atributes and methods for basis
                functions and transformation jacobians
            q: thermal diffusivity

        """
        mesh = self.mesh
        self.R = np.zeros((4, mesh.num_ele))

        for e in range(mesh.num_ele):
            for w in range(4):
                mesh.basisFunction2D(self.gp[w])
                mesh.eleJacobian(mesh.nodes_coord[
                    mesh.ele_conn[e]])

                x1_o_e1e2, x2_o_e1e2 = mesh.mapping(e)

                load = q(x1_o_e1e2, x2_o_e1e2)

                self.R[0, e] += load*mesh.phi[0]*mesh.detJac
                self.R[1, e] += load*mesh.phi[1]*mesh.detJac
                self.R[2, e] += load*mesh.phi[2]*mesh.detJac
                self.R[3, e] += load*mesh.phi[3]*mesh.detJac



    def mass(self, mat):
        """Build the mass matrix for each element.

        """
        mesh = self.mesh
        self.M = np.zeros((4, 4, mesh.num_ele))

        for surface, [k, c, rho] in mat.items():
            for e in range(mesh.num_ele):
                if surface == mesh.ele_surface[e, 1]:
                    for w in range(4):
                        mesh.basisFunction2D(self.gp[w])
                        mesh.eleJacobian(mesh.nodes_coord[
                            mesh.ele_conn[e]])

                        dA = mesh.detJac

                        self.M[0, :, e] += c*rho*mesh.phi[0]*mesh.phi[:]*dA
                        self.M[1, :, e] += c*rho*mesh.phi[1]*mesh.phi[:]*dA
                        self.M[2, :, e] += c*rho*mesh.phi[2]*mesh.phi[:]*dA
                        self.M[3, :, e] += c*rho*mesh.phi[3]*mesh.phi[:]*dA


    def load_transient(self, q, time):
        """Build the load vector thermal diffusivity.

        Args:
            mesh: object that includes the mesh atributes and methods for basis
                functions and transformation jacobians
            q: thermal diffusivity

        """
        mesh = self.mesh
        self.P0q = np.zeros((4, mesh.num_ele))

        for e in range(mesh.num_ele):
            for w in range(4):
                mesh.basisFunction2D(self.gp[w])
                mesh.eleJacobian(mesh.nodes_coord[
                    mesh.ele_conn[e]])

                x1_o_e1e2, x2_o_e1e2 = mesh.mapping(e)

                load = q(x1_o_e1e2, x2_o_e1e2, time)

                self.P0q[0, e] += load*mesh.phi[0]*mesh.detJac
                self.P0q[1, e] += load*mesh.phi[1]*mesh.detJac
                self.P0q[2, e] += load*mesh.phi[2]*mesh.detJac
                self.P0q[3, e] += load*mesh.phi[3]*mesh.detJac

