"""
in steps s d=1200 n=2
in seed s d=14 n=2
in framenum s d=0 n=2
in Du s d=0.16 n=2
in Dv s d=0.08 n=2
in F s d=0.035 n=2
in k s d=0.060 n=2
out verts_out v
"""


def setup():
    
    import random
    import numpy as np
    import numpy.ctypeslib as npct
    import os
    import ctypes as ct

    # Parameters for Gray-Scott Reaction Diffusion
    class RdStruct(ct.Structure):
        _fields_ = [
            ('Du', ct.c_double),
            ('Dv', ct.c_double),
            ('F', ct.c_double),
            ('k', ct.c_double),
        ]

    # callback function
    def callback_func(step, length, data, selfp):
        print("step:", step)

        # get frame data from array data
        v_arr = np.ctypeslib.as_array(ct.POINTER(ct.c_double).from_address(ct.addressof(data)), shape=(length, 3))
        arr_stored = v_arr[:length].tolist()
        
        # store the frame data
        instance = ct.cast(selfp, ct.py_object).value
        instance.store_frame(step, arr_stored)

    # Reaction Diffusion main class
    class DiffReact(ct.Structure):

        verts = []
               
        # initialize function
        def __init__(self, np, steps=500, seed=24):
            self.frame_storage = {}
            self.n = n = 256
            random.seed(seed)

            # Initialize U,V array
            U, V = np.zeros((n+2, n+2), dtype=np.double), np.zeros((n+2, n+2), dtype=np.double)
            u, v = U[1:-1, 1:-1], V[1:-1, 1:-1]

            # Initialize rectangle vertices
            r = 20
            u[...] = 1.0
            U[n//2-r:n//2 + r, n//2-r:n//2 + r] = 0.50
            V[n//2-r:n//2 + r, n//2-r:n//2 + r] = 0.25
            u += 0.05 * np.random.random((n, n))
            v += 0.05 * np.random.random((n, n))

            # Declare callback function type
            NpPtr = np.ctypeslib.ndpointer(dtype=ct.c_double, flags='C_CONTIGUOUS')
            cfunc_type = ct.CFUNCTYPE(None, ct.c_int, ct.c_int, NpPtr, ct.c_void_p)
            
            # Load external C library
            libBasicRD2d = npct.load_library('libBasic_RD_2d', os.path.dirname('/Path/to/library/directory/'))

            # Declare process function arguments and result types
            libBasicRD2d.process.argtypes = [
               ct.POINTER(RdStruct),
               npct.ndpointer(dtype=np.uintp, ndim=1, flags='C'),
               npct.ndpointer(dtype=np.uintp, ndim=1, flags='C'),
               ct.c_int,
               ct.c_int,
               cfunc_type,
               ct.py_object
               ]
            libBasicRD2d.process.restype = ct.c_int

            # Set Gray-Scott parameters
            rd_struct = RdStruct()
            rd_struct.Du = Du
            rd_struct.Dv = Dv
            rd_struct.F = F
            rd_struct.k = k

            # Convert array types for C function
            Upp = (U.__array_interface__['data'][0] + np.arange(U.shape[0])*U.strides[0]).astype(np.uintp)
            Vpp = (V.__array_interface__['data'][0] + np.arange(V.shape[0])*V.strides[0]).astype(np.uintp)

            # Call C function
            res = libBasicRD2d.process(rd_struct, Upp, Vpp, n, steps, cfunc_type(callback_func), ct.py_object(self))

        # Function to get frames
        def get_frame(self, number):
            return self.frame_storage.get(number)

        # Function to store frames
        def store_frame(self, framestep, data):
            self.frame_storage[framestep] = data

    # Instantiate main class
    GK = DiffReact(np, steps, seed)

# Return results for each frame
verts = GK.get_frame(framenum)
verts_out.append(verts)


