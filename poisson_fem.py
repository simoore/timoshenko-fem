from math import floor
import time
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Node(object):
    
    def __init__(self, i, j, dof):
        """
        :i: is the x-index of the node.
        :j: is the y-index of the node.
        :dof: is the dof index associated with the node.
        """
        self.i = i
        self.j = j
        self.dof = dof
        
        
    def __repr__(self):
        
        return 'Node %d: %d, %d' % (self.dof, self.i, self.j)
        
        
class Element(object):
    
    def __init__(self, i, j, nodes_2D):
        """
        :i: is the x-index of the element.
        :j: is the y-index of the element.
        :nodes_2D: list of nodes indexed by (i,j).
        """
        self.i = i
        self.j = j
        self.n1 = nodes_2D[i][j]
        self.n2 = nodes_2D[i + 1][j]
        self.n3 = nodes_2D[i + 1][j + 1]
        self.n4 = nodes_2D[i][j + 1]
        self.dofs = np.array([self.n1.dof, self.n2.dof, self.n3.dof, self.n4.dof])
        
        
class RectangularMesh(object):
    
    def __init__(self, nx, ny):
        """
        :nx: the number of elements in the x-direction.
        :ny: the number of elements in the y-direction.
        """
        # Create the nodes.
        node_index = np.ndindex(nx + 1, ny + 1)
        self.nodes = [Node(i, j, k) for k, (i, j) in enumerate(node_index)]
        self.nodes_2D = [[None for j in range(ny + 1)] for i in range(nx + 1)]
        for i, n in enumerate(self.nodes):
            self.nodes_2D[n.i][n.j] = n
            
        # Create the elements.
        elem_index = np.ndindex(nx, ny)
        self.elements = [Element(i, j, self.nodes_2D) for i, j in elem_index]
        self.elements_2D = [[None for j in range(ny)] for i in range(nx)]
        for e in self.elements:
            self.elements_2D[e.i][e.j] = e
            
        # Number of nodes and elements in the mesh.
        self.n_elem = len(self.elements)
        self.n_node = len(self.nodes)
        
        # Node constraint matrix used by the block diagonal assembler.
        row = np.arange(4 * self.n_elem)
        col = np.empty(4 * self.n_elem)
        for i, e in enumerate(self.elements):
            col[4*i:4*i+4] = e.dofs
        val = np.ones(4 * self.n_elem)
        shape = (4 * self.n_elem, len(self.nodes))
        self.con = sparse.coo_matrix((val, (row, col)), shape=shape).tocsr()
        
        # All, free, and boundary dofs. Boundary nodes are at j == 0.
        self.all_dofs = np.arange(self.n_node)
        self.fixed_dofs = [n.dof for n in self.nodes if n.j == 0]
        self.free_dofs = np.setdiff1d(self.all_dofs, self.fixed_dofs)
        
        
class ElementModel(object):
    
    def __init__(self, a, b):
        """
        :a: half the width of an element in the x-direction.
        :b: half the width of an element in the y-direction.
        """
        points = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / np.sqrt(3)
        jacobian = a*b
        ke = np.zeros((4, 4))
        fe = np.zeros((4, 1))
        for p in points:
            n, dndxi, dndeta = self.shapes(p)
            ke += jacobian * (dndxi.T @ dndxi / a*a + dndeta.T @ dndeta / b*b)
            fe += jacobian * n.T
        self.ke = ke
        self.fe = fe
        
          
    @staticmethod
    def shapes(point):
        """
        :point: the point (xi, eta) to evaluate the shape functions.
        """
        xi, eta = point
        xs = np.array([[-1, 1, 1, -1]])
        es = np.array([[-1, -1, 1, 1]])
        n = 0.25 * (1 + xs * xi) * (1 + es * eta)
        dndxi = xs * 0.25 * (1 + es * eta)
        dndeta = es * 0.25 * (1 + xs * xi) 
        return n, dndxi, dndeta
    
    
class AssemblerBlock(object):
    
    def __init__(self, mesh, model, k, q):  
    
        t0 = time.time()
        k_block = [ke * model.ke for ke in k]
        f_block = [qe * model.fe for qe in q]
        self.sysk1 = mesh.con.T @ sparse.block_diag(k_block, format='csr') @ mesh.con
        self.sysf1 = mesh.con.T @ sparse.coo_matrix(np.vstack(f_block)).tocsr()
        t1 = time.time()
        print('Block Assembley: %g (s)' % (t1-t0))
        
        
    def get_system(self):
        
        return self.sysk1, self.sysf1


class AssemblerIterative(object):
    
    def __init__(self, mesh, model, k, q):

        t0 = time.time()
        ntriplet = 0
        fntriplet = 0
        row = np.zeros(mesh.n_elem * 16)
        col = np.zeros(mesh.n_elem * 16)
        val = np.zeros(mesh.n_elem * 16)
        frow = np.zeros(mesh.n_elem * 4)
        fcol = np.zeros(mesh.n_elem * 4)
        fval = np.zeros(mesh.n_elem * 4)
        for qe, ke, e in zip(q, k, mesh.elements):
            dof = [e.n1.dof, e.n2.dof, e.n3.dof, e.n4.dof]
            for ii, jj in np.ndindex(4, 4):
                row[ntriplet] = dof[ii]
                col[ntriplet] = dof[jj]
                val[ntriplet] = ke * model.ke[ii, jj]
                ntriplet += 1
            for ii in range(4):
                frow[fntriplet] = dof[ii]
                fcol[fntriplet] = 0
                fval[fntriplet] = qe * model.fe[ii]
                fntriplet += 1
                
        shape = (mesh.n_node, mesh.n_node)
        fshape = (mesh.n_node, 1)
        self.sysk2 = sparse.coo_matrix((val, (row, col)), shape=shape).tocsr()
        self.sysf2 = sparse.coo_matrix((fval, (frow, fcol)), shape=fshape).tocsr()
        t1 = time.time()
        print('Iterative Assembley: %g (s)' % (t1-t0))


    def get_system(self):
        
        return self.sysk2, self.sysf2


class AssemblerVector(object):
    """This approach follows the work in the article:
    Effecient topology optimization in MATLAB using 88 lines of code;
    Erik Andreassen, Anders Clausen, Mattias Schevenels, Boyan S. Lazarov,
    Ole Sigmund; Structural and Multidisciplinary Optimization; January 
    2011; Volume 43, Issue 1, pp 1–16;
    """
    def __init__(self, mesh, model, k, q):

        t0 = time.time()
        kr = np.hstack([np.hstack([e.dofs for _ in range(4)]) for e in mesh.elements])
        kc = np.hstack([np.hstack([d*np.ones(4) for d in e.dofs]) for e in mesh.elements])
        fr = np.hstack([e.dofs for e in mesh.elements])
        fc = np.zeros_like(fr)
        kv = np.expand_dims(k, axis=1) @ np.expand_dims(model.ke.ravel(), axis=0)
        kv = np.ravel(kv)
        fv = np.expand_dims(q, axis=1) @ model.fe.T 
        fv = np.ravel(fv)
        
        shape = (mesh.n_node, mesh.n_node)
        fshape = (mesh.n_node, 1)     
        self.sysk3 = sparse.coo_matrix((kv, (kr, kc)), shape=shape).tocsr()
        self.sysf3 = sparse.coo_matrix((fv, (fr, fc)), shape=fshape).tocsr()
        t1 = time.time()
        print('Vector Assembley: %g (s)' % (t1-t0))
        
        
    def get_system(self):
        
        return self.sysk3, self.sysf3
    
    
class FiniteElement(object):
    
    def __init__(self):
        """On initialization, the finite element analysis is executed step
        by step and the solution is saved.
        """
        # Constants defining the domain of the problem.
        self.a = 0.1
        self.b = 0.1
        self.nelx = 20
        self.nely = 20
        
        # Initialize mesh and the element matrices.
        self.mesh = RectangularMesh(self.nelx, self.nely)
        self.model = ElementModel(self.a, self.b)
        
        # Get the parameters of the problem and assemble system matrices.
        k, q = self.define_k_and_q(self.mesh)
        assembler = AssemblerVector(self.mesh, self.model, k, q)
        sysk, sysf = assembler.get_system()
        
        # Apply boundary conditions.
        sysk = sysk[self.mesh.free_dofs, :][:, self.mesh.free_dofs]
        sysf = sysf[self.mesh.free_dofs]
        
        # Solve system of equations.
        self.u = linalg.spsolve(sysk, sysf)
        
        # Add in the boundary dofs.
        self.uall = np.zeros(self.mesh.all_dofs.shape)
        self.uall[self.mesh.free_dofs] = self.u
        
        
    def interpolate(self, x, y):
        """Using the solution stored in 'self.uall' the temperature value is
        interpolated for the point (x, y).
        
        Parameters
        ----------
        :x,y: The coordinate of the point at which the temperature is 
              calculated.
              
        Returns
        -------
        :u: The temperature at the point (x,y).
        """
        i0 = floor(x / (2*self.a))
        j0 = floor(y / (2*self.b))
        x0 = (x - self.a*(2*i0 + 1)) / self.a
        y0 = (y - self.b*(2*j0 + 1)) / self.b
        e = self.mesh.elements_2D[i0][j0]
        ue = self.uall[e.dofs]
        n, _, _ = self.model.shapes((x0, y0))
        u = np.sum(ue * n)
        return u
        
        
    def plot(self):
        """Plots the temperature field calculated using the finite element
        analysis.
        """
        
        # Determine the points over which to evaluate the temperature.
        delta_x = 2 * self.a * 0.1
        delta_y = 2 * self.b * 0.1
        x = np.arange(0, 2*self.nelx*0.1, delta_x)
        y = np.arange(0, 2*self.nely*0.1, delta_y)
        mx, my = np.meshgrid(x, y)
        
        # Determine the temperature at each point for the plot.
        mz = np.zeros(mx.shape)
        nx, ny = my.shape
        for i in range(nx):
            for j in range(ny):
                mz[i, j] = self.interpolate(mx[i, j], my[i, j])

        # Create a contour plot.
        fig, ax = plt.subplots()
        cs = ax.contour(mx*5, my*5, mz, 8, colors='k')
        ax.clabel(cs, inline=1, fontsize=10)
        ax.set_title('Temperature Distribution')
        ax.set_yticks(list(range(20)))
        ax.set_xticks(list(range(20)))
        
        ax.grid(linestyle='-', linewidth='0.5', color='grey')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.tick_params(direction='in')
        
        # Annotate the point of maximum temperature.
        index_array = np.argmax(mz)
        max_x = np.ravel(mx)[index_array]
        max_y = np.ravel(my)[index_array]
        max_z = mz.ravel()[index_array]
        ax.annotate('max: %g' % max_z, 
                    xy=(5*max_x, 5*max_y), xytext=(max_x*5, 5*(max_y-0.4)),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2))
        
        
        # Highlight the area of high and low conductivity.
        k, q = self.define_k_and_q(self.mesh)
        for i, e in enumerate(self.mesh.elements):
            x = e.i * 1
            y = e.j * 1
            if k[i] < 1:
                ax.add_patch(patches.Rectangle((x, y), 1, 1, alpha=0.1))
            else:
                ax.add_patch(patches.Rectangle((x, y), 1, 1, alpha=0.8))
        
        
        
        note = 'Area of low conductivity'
        ax.annotate(note, xy=(2*5, 0.85*5))
        note = 'Area of high conductivity'
        ax.annotate(note, xy=(2*5, 0.60*5))
        note = 'Temperature Boundary (T=0)'
        ax.annotate(note, xy=(2*5, 0.05*5))
        
        fig.savefig('poisson.png', dpi=300)


    @staticmethod
    def define_k_and_q(mesh):
        """
        `k` is the heat conductivity, and `q` is the heat source per unit area. 
        In this code, each are constant for each element.
        
        Parameters
        ----------
        :mesh: The mesh is used to determine where the areas of low 
               conductivity are, and the location of the heat source.
               
        Returns
        -------
        :k: The array of heat conductivity values for each element in an
            order matching that of mesh.elements.
        :q: The heat source on each element.
        """
        k = np.zeros(len(mesh.elements))
        q = np.zeros(len(mesh.elements))
        
        for i, e in enumerate(mesh.elements):
            k[i] = 1e-2 if (3 < e.j < 16) or (e.j > 15 and e.i > 10) else 1
            q[i] = 1 if k[i] == 1 else 0
            
        return k, q
                
    
if __name__ == '__main__':
     
    fem = FiniteElement()
    fem.plot()
            