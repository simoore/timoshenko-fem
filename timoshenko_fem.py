import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg


class Element(object):
    def __init__(self, node_a, node_b):
        self.node_a = node_a
        self.node_b = node_b
        
        
    def get_dof(self):
        return (self.node_a.dof_a, self.node_a.dof_b,
                self.node_b.dof_a, self.node_b.dof_b)
    
    
    def get_boundary(self):
        return (self.node_a.boundary, self.node_a.boundary,
                self.node_b.boundary, self.node_b.boundary)
        
        
class Node(object):
    def __init__(self, coord):
        self.coord = coord
        self.boundary = False
        self.dof_a = 0
        self.dof_b = 0
        
        
def mesh(xa, xb, n_elements):
    n_nodes = n_elements + 1
    xcoords, step = np.linspace(xa, xb, num=n_nodes, retstep=True)
    nodes = [Node(x) for x in xcoords]
    nodes[0].boundary = True
    elements = [None for _ in range(n_elements)]
    for i in range(n_elements):
        elements[i] = Element(nodes[i], nodes[i + 1])
    a = 0.5 * step
    n_dof = 2 * (n_nodes - 1)
    dof_count = 0
    for n in nodes:
        if n.boundary == False:
            n.dof_a = dof_count
            n.dof_b = dof_count + 1
            dof_count += 2
    return elements, a, n_dof
        
    
def shapes(xi, a):
    n = [0.5 * (1 - xi), 0.5 * (1 + xi)]
    dn = [-0.5 / a, 0.5 / a]
    return n, dn


def element_matrices(a, G, E, Iz, rho, A, kappa):
    ke = np.zeros((4, 4))
    me = np.zeros((4, 4))
    
    points = [0.577350269189626, -0.577350269189626]
    for p in points:
        n, dn = shapes(p, a)
        Ba = np.array([[0, dn[0], 0, dn[1]]])
        Bb = np.array([[n[0], 0, n[1], 0]])
        Bc = np.array([[0, n[0], 0, n[1]]])
        ke += E * Iz * Ba.T @ Ba * a
        me += a * rho * A * Bb.T @ Bb + a * rho * Iz * Bc.T @ Bc
    
    point, weight = 0, 2
    n, dn = shapes(point, a)
    Bd = np.array([[dn[0], n[0], dn[1], n[1]]])
    ke += weight * a * kappa * G * A * Bd.T @ Bd
    
    return ke, me


def assemble(elements, ke, me, n_dof):
    ntriplet = 0
    num = 16 * len(elements)
    row, col = np.zeros(num), np.zeros(num)
    kk, mm = np.zeros(num), np.zeros(num)
    for e in elements:
        boundary = e.get_boundary()
        dof = e.get_dof()
        for ii in range(4):
            for jj in range(4):
                if boundary[ii] is False and boundary[jj] is False:
                    row[ntriplet] = dof[ii]
                    col[ntriplet] = dof[jj]
                    mm[ntriplet] = me[ii, jj]
                    kk[ntriplet] = ke[ii, jj]
                    ntriplet += 1
                    
    kks = sparse.coo_matrix((kk, (row, col)), shape=(n_dof, n_dof)).tocsr()
    mms = sparse.coo_matrix((mm, (row, col)), shape=(n_dof, n_dof)).tocsr()
    return kks, mms

                    
def modal_analysis(E, h, poisson, L, wid, rho, n_elements, n_modes):
    Iz = h ** 3 / 12
    kappa = 5 / 6
    rho = 1
    A = wid * h
    G = 0.5 * E / (1 + poisson)

    elements, a, n_dof = mesh(0, L, n_elements)
    ke, me = element_matrices(a, G, E, Iz, rho, A, kappa)
    kks, mms = assemble(elements, ke, me, n_dof)

    w, v = linalg.eigsh(kks, k=n_modes, M=mms, sigma=0, which='LM')
    return w, v


def main():
    E = 30e6
    h = 0.01
    poisson = 0.3
    L = 1
    wid = 1
    rho = 1
    num_elements = 5000
    num_modes = 5
        
    w, v = modal_analysis(E, h, poisson, L, wid, rho, num_elements, num_modes)
    Iz = h ** 3 / 12
    A = wid * h
    normlams = [1.8751, 4.6941, 7.8548, 10.996, 14.137]
    lam = [(l / L) ** 4 * E * Iz / (rho * A) for l in normlams]
    
    print('Modal Frequencies of the Beam (Hz).')
    print('          Analytical      Numerical')
    for i, (m, n) in enumerate(zip(lam, w)):
        print('Mode %d: %12g   %12g' % (i+1, m, n))

    
if __name__ == '__main__':
    main()
    
    