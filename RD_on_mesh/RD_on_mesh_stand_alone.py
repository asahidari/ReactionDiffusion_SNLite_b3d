"""
in steps s d=3 n=2
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

        percent = 0
        def get_frame(self, index, number):
            return self.frame_storage_a.get(number) if index == 0 else self.frame_storage_b.get(number)

        def store_frame(self, index, framestep, data):
            if index == 0:
                self.frame_storage_a[framestep] = data
            else:
                self.frame_storage_b[framestep] = data

        def print_progress(self, total, step):
            pn = 100 * (step + 1) // total # percent completed
            if pn != self.percent:
                self.percent = pn
                print("%" + str(self.percent).zfill(2))


        def __init__(self, np, verts=[], polygons=[], selected_verts=[], steps=1200, seed=14):
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

            # use active cells if exists
            np.put(U, selected_verts, 1.0)

            # Store neighbor vertex indices
            neighbors = []
            for i in range(num_verts):
                # Select polygons including the index
                included_polys = np.delete(polygons, np.where([i not in p for p in polygons]), axis=0)

                # Extract neighbor vertex indices
                neighbors_in_poly = np.array([[poly[(np.where(poly == i)[0]-1 + len(poly))%len(poly)], \
                                    poly[(np.where(poly == i)[0]+1)%len(poly)]] \
                                    for poly in included_polys])
                neighbor = np.unique(neighbors_in_poly.reshape(-1))
                neighbors.append(neighbor)

            # calculate reaction diffusion
            p = 0
            Ubuff, Vbuff = np.zeros(num_verts), np.zeros(num_verts)
            for i in range(steps):

                for j in range(num_verts):
                    du, dv = 0., 0.
                    neighbor = neighbors[j]

                    # Calculate laplacian
                    if len(neighbor) > 0:
                        du += np.sum(np.array([U[k] for k in neighbor]))
                        dv += np.sum(np.array([V[k] for k in neighbor]))
                        du = du / len(neighbor) - U[j]
                        dv = dv / len(neighbor) - V[j]

                    uvv = U[j] * V[j] * V[j]
                    Ubuff[j] = U[j] + (Du * du - uvv + F * (1. - U[j])) * dt
                    Vbuff[j] = V[j] + (Dv * dv + uvv - (F + K) * V[j]) * dt

                # Regulate values
                U = np.array([min(ub2, 1.0) for ub2 in [max(ub1, 0.0) for ub1 in Ubuff]])
                V = np.array([min(vb2, 1.0) for vb2 in [max(vb1, 0.0) for vb1 in Vbuff]])

                denominator_u = U.max() - U.min() if U.max() > U.min() else 0.0001
                denominator_v = V.max() - V.min() if V.max() > V.min() else 0.0001
                Usub = (U - U.min())/denominator_u
                Vsub = (V - V.min())/denominator_v

                self.store_frame(0, i, Usub.tolist())
                self.store_frame(1, i, Vsub.tolist())

                self.print_progress(steps, i)

    RD = RD_mesh(np, verts_in[0], polygons_in[0], selected_verts[0], steps, seed)

a_out, b_out = [RD.get_frame(0, framenum)], [RD.get_frame(1, framenum)]

