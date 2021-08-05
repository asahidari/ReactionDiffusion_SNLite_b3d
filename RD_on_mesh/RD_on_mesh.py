"""
in steps s d=10 n=2
in Du s d=0.87 n=2
in Dv s d=0.23 n=2
in F s d=0.055 n=2
in K s d=0.062 n=2
in dt s d=1.0 n=2
in seed s d=14 n=2
in framenum s d=0 n=2
in verts_in v
in polygons_in s
in selected_verts s
out a_out s
out b_out s
"""

def setup():
    
    import random
    import numpy as np
    import numpy.ctypeslib as npct
    import ctypes as ct
    import os

    # Parameters for Gray-Scott Reaction Diffusion
    class RdParams(ct.Structure):
        _fields_ = [
            ('Du', ct.c_double),
            ('Dv', ct.c_double),
            ('F', ct.c_double),
            ('k', ct.c_double),
            ('dt', ct.c_double),
        ]
        
    # callback function
    def callback_func(step, type, count, data, selfp):
        print("step:", step)

        # get frame data from array data
        v_arr = np.ctypeslib.as_array(ct.POINTER(ct.c_double).from_address(ct.addressof(data)), shape=(count, ))
        arr_stored = v_arr[:count].tolist()
        
        # store the frame data
        instance = ct.cast(selfp, ct.py_object).value
        instance.store_frame(type, step, arr_stored)

    # Reaction Diffusion main class
    class RD_mesh(ct.Structure):

        # Function to get frames
        def get_frame(self, index, number):
            return self.frame_storage_a.get(number) if index == 0 else self.frame_storage_b.get(number) 

        # Function to store frames
        def store_frame(self, index, framestep, data):
            if index == 0:
                self.frame_storage_a[framestep] = data
            else:
                self.frame_storage_b[framestep] = data
            
        # initialize function
        def __init__(self, np, verts=[], polygons=[], steps=1200, seed=14):
            self.frame_storage_a, self.frame_storage_b = {}, {}

            # If any input data is None, do nothing.
            if verts is None or polygons is None:
                print("Not found any verts or polygons.")
                return
            
            num_verts = len(verts)
            np.random.seed(seed)
            
            # Initialize U/V values
            U = np.zeros(num_verts)
            U += 0.05 * np.random.randint(0, 2, num_verts)
            V = np.random.uniform(0, 0.4, num_verts)
            V += 0.05 * np.random.randint(0, 2, num_verts)

            # Use active cells (selected verts) if designated
            if selected_verts is not None:
                for i in range(len(verts)):
                    if i in selected_verts[0]:
                        U[i] = 1.00
            
            # Store neighbor vertex indices
            neighbors = []
            max_neighbor_count = 20
            neighbor_counts = np.zeros(num_verts, dtype=np.int32)
            neighbor_indices = np.zeros((num_verts, max_neighbor_count), dtype=np.int32)
            
            # Store vertex indices of each vertices
            for i in range(num_verts):
                neighbor = []
                for poly in polygons:
                    if i in poly:
                        j = poly.index(i)
                        num_p = len(poly)
                        i_prev = poly[(j-1 + num_p) % num_p]
                        i_next = poly[(j+1) % num_p]
                        if i_prev not in neighbor: neighbor.append(i_prev)
                        if i_next not in neighbor: neighbor.append(i_next)
                neighbors.append(neighbor)
                neighbor_counts[i] = len(neighbor)
                for k in range(len(neighbor)):
                    neighbor_indices[i][k] = neighbor[k]

            # Declare callback function type 
            NpPtr = np.ctypeslib.ndpointer(dtype=ct.c_double, flags='C_CONTIGUOUS')
            cfunc_type = ct.CFUNCTYPE(None, ct.c_int, ct.c_int, ct.c_int, NpPtr, ct.c_void_p)
            
            # Load external C library
            libRDOnMesh = npct.load_library('libRD_on_mesh', os.path.dirname('/Path/to/library/directory/'))

            # Declare process function arguments and result types
            libRDOnMesh.process.argtypes = [
               ct.POINTER(RdParams),
               ct.c_int,
               npct.ndpointer(dtype=ct.c_int, ndim=1, flags='C'),
               npct.ndpointer(dtype=np.uintp, ndim=1, flags='C'),
               npct.ndpointer(dtype=ct.c_double, ndim=1, flags='C'),
               npct.ndpointer(dtype=ct.c_double, ndim=1, flags='C'),
               ct.c_int,
               cfunc_type,
               ct.py_object
               ]
            libRDOnMesh.process.restype = ct.c_int

            # Set Gray-Scott parameters
            rd_params = RdParams()
            rd_params.Du = Du
            rd_params.Dv = Dv
            rd_params.F = F
            rd_params.k = K
            rd_params.dt = dt

            # Convert array types for C function
            p_neighbor_indices = (neighbor_indices.__array_interface__['data'][0] + np.arange(neighbor_indices.shape[0])*neighbor_indices.strides[0]).astype(np.uintp)
            
            # Call C function
            res = libRDOnMesh.process(rd_params, num_verts, neighbor_counts, p_neighbor_indices, U, V, steps, cfunc_type(callback_func), ct.py_object(self))
            
    # Instantiate main class
    RD = RD_mesh(np, verts_in[0], polygons_in[0], steps, seed)

# Return results for each frame
a_out, b_out = [RD.get_frame(0, framenum)], [RD.get_frame(1, framenum)]

