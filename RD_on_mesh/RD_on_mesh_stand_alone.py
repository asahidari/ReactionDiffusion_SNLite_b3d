"""
in steps s d=1200 n=2
in Du s d=0.84 n=2
in Dv s d=0.41 n=2
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

# Example: Du, Dv, F, K = 0.87, 0.19, 0.065, 0.062

def setup():
    
    import random
    import numpy as np

    class RD_mesh():

        def get_frame(self, index, number):
            return self.frame_storage_a.get(number) if index == 0 else self.frame_storage_b.get(number) 

        def store_frame(self, index, framestep, data):
            if index == 0:
                self.frame_storage_a[framestep] = data
            else:
                self.frame_storage_b[framestep] = data
            
        def __init__(self, np, verts=[], polygons=[], steps=1200, seed=14):
            self.frame_storage_a, self.frame_storage_b = {}, {}

            # If any input data is None, do nothing.
            if verts is None or polygons is None:
                print("Not found any verts or polygons.")
                return
            
            num_verts = len(verts)
            np.random.seed(seed)
            
            # Initialize U/V values
            # U = np.full(len(verts), 1.0)
            U = np.zeros(num_verts)
            U += 0.05 * np.random.randint(0, 2, num_verts)
            
            # V = np.random.randint(0, 2, len(verts)).astype("float")
            # V = np.zeros(len(verts))
            V = np.random.uniform(0, 0.4, num_verts)
            V += 0.05 * np.random.randint(0, 2, num_verts)

            # use active cells if designated
            if selected_verts is not None:
                for i in range(len(verts)):
                    if i in selected_verts[0]:
                        U[i] = 1.00
            # U -= 0.05 * np.random.randint(0, 2, len(U))
            
            # Store neighbor vertex indices
            neighbors = []
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

            # calculate reaction diffusion
            p = 0
            Ubuff, Vbuff = np.zeros(num_verts), np.zeros(num_verts)
            for i in range(steps):
                for j in range(num_verts):
                    du, dv = 0., 0.
                    for neighbor_j in neighbors[j]:
                        du += U[neighbor_j]
                        dv += V[neighbor_j]
                    du = (du / len(neighbors[j])) - U[j]
                    dv = (dv / len(neighbors[j])) - V[j]
                    
                    uvv = U[j] * V[j] * V[j]
                    Ubuff[j] = U[j] + (Du * du - uvv + F * (1. - U[j])) * dt
                    Vbuff[j] = V[j] + (Dv * dv + uvv - (F + K) * V[j]) * dt
                
                for j in range(num_verts):
                    U[j] = Ubuff[j]
                    V[j] = Vbuff[j]
                
                for j in range(num_verts):
                    U[j] = min(max(U[j], 0.0), 1.0)
                    V[j] = min(max(V[j], 0.0), 1.0)

                Usub, Vsub = [], []
                Ud = U.max() - U.min() if U.max() > U.min() else 0.01
                Vd = V.max() - V.min() if V.max() > V.min() else 0.01
                for j in range(num_verts):
                    Usub.append( (U[j] - U.min())/Ud )
                    Vsub.append( (V[j] - V.min())/Vd )
                    
                self.store_frame(0, i, Usub)
                self.store_frame(1, i, Vsub)
                pn = 100 * (i + 1) // steps # percent completed
                if pn != p:
                    p = pn
                    print("%" + str(p).zfill(2))
                    
        
    RD = RD_mesh(np, verts_in[0], polygons_in[0], steps, seed)

a_out, b_out = [RD.get_frame(0, framenum)], [RD.get_frame(1, framenum)]
