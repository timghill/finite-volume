"""
Module pymesh reads matlab-generated triangular mesh objects into equivalent
python classes.
"""

import scipy.io

class dmesh:
    """
    Class dmesh represents triangular mesh objects generated using the
    matlab `mattri` package (https://bitbucket.org/maurow/mattri/src/master/).
    """
    def __init__(self, meshfile):
        dmesh_mat=scipy.io.loadmat(meshfile)

        self.tri=__tri__(dmesh_mat)

class __tri__:
    """
    Class __tri__ reads the dmesh.tri attributes of mattri-generated objects.
    """
    def __init__(self, mesh):
        self.type=mesh['tri']['type'][0][0][0]

        self.connect=mesh['tri']['connect'][0][0].astype(int)
        self.connect_edge=mesh['tri']['connect_edge'][0][0].astype(int)
        self.connect_el_el=mesh['tri']['connect_el_el'][0][0].astype(int)
        self.connect_el_edge=mesh['tri']['connect_el_edge'][0][0].astype(int)

        self.nodes=mesh['tri']['nodes'][0][0].astype(float)
        self.elements=mesh['tri']['elements'][0][0].astype(float)
        
        self.bmark=mesh['tri']['bmark'][0][0].astype(int)
        self.bmark_edge=mesh['tri']['bmark_edge'][0][0].astype(int)
        self.bmark_el=mesh['tri']['bmark_el'][0][0].astype(int)

        self.n_nodes=int(mesh['tri']['n_nodes'][0][0][0])
        self.n_edges=int(mesh['tri']['n_edges'][0][0][0])
        self.n_elements=int(mesh['tri']['n_elements'][0][0][0])

        self.area=mesh['tri']['area'][0][0].astype(float)
        self.area_nodes=mesh['tri']['area_nodes'][0][0].astype(float)
        
        self.edge_length=mesh['tri']['edge_length'][0][0].astype(float)
        self.edge_midpoints=mesh['tri']['edge_midpoints'][0][0].astype(float)
        
        self.connect_edge_inv=mesh['tri']['connect_edge_inv'][0][0][0]
        
        self.neigh_node=mesh['tri']['neigh_node'][0][0]
        self.neigh_edge_node=mesh['tri']['neigh_edge_node'][0][0]

        self.nx=mesh['tri']['nx'][0][0].astype(float)
        self.ny=mesh['tri']['ny'][0][0].astype(float)

        self.ds=mesh['tri']['ds'][0][0].astype(float)


meshfile='rect_mesh2.mat'
x=dmesh(meshfile)
