"""
in steps s d=1200 n=2
in seed s d=14 n=2
in framenum s d=0 n=2
out verts_out v
"""


def setup():
    
    import random
    import numpy as np


    class DiffReact2():

        verts = []
        params = []
        params.append((0.16, 0.08, 0.035, 0.065)) # Bacteria 1
        params.append((0.14, 0.06, 0.035, 0.065)) # Bacteria 2
        params.append((0.16, 0.08, 0.060, 0.062)) # Coral
        params.append((0.19, 0.05, 0.060, 0.062)) # Fingerprint
        params.append((0.10, 0.10, 0.018, 0.050)) # Spirals
        params.append((0.12, 0.08, 0.020, 0.050)) # Spirals Dense
        params.append((0.10, 0.16, 0.020, 0.050)) # Spirals Fast
        params.append((0.16, 0.08, 0.020, 0.055)) # Unstable
        params.append((0.16, 0.08, 0.050, 0.065)) # Worms 1
        params.append((0.16, 0.08, 0.054, 0.063)) # Worms 2
        params.append((0.16, 0.08, 0.035, 0.060)) # Zebrafish

       
        def __init__(self, np, steps=500, seed=24):
            self.frame_storage = {}
            self.n = n = 256//4
            self.imgx = n
            self.imgy = n
            self.imgz = n

            random.seed(seed)
            (Du, Dv, F, k) = random.choice(self.params)

            Z = np.zeros((n+2, n+2, n+2), [('U', np.double), ('V', np.double)])
            U, V = Z['U'], Z['V']
            u, v = U[1:-1, 1:-1, 1:-1], V[1:-1, 1:-1, 1:-1]

            r = 20
            u[...] = 1.0
            U[n//2-r:n//2 + r, n//2-r:n//2 + r, n//2-r:n//2 + r] = 0.50
            V[n//2-r:n//2 + r, n//2-r:n//2 + r, n//2-r:n//2 + r] = 0.25
            u += 0.05 * np.random.random((n, n, n))
            v += 0.05 * np.random.random((n, n, n))

            ######### loop start ############
            p = 0
            for i in range(steps):
                Lu = (        U[1:-1, 1:-1, 0:-2] +
                                        U[0:-2,1:-1, 1:-1] +
                    U[1:-1,0:-2, 1:-1] - 6*U[1:-1,1:-1, 1:-1] + U[1:-1,2:, 1:-1] +
                                        U[2:  ,1:-1, 1:-1] +
                                                 + U[1:-1, 1:-1, 2:] )
                Lv = (        V[1:-1, 1:-1, 0:-2] +
                                        V[0:-2,1:-1, 1:-1] +
                    V[1:-1,0:-2, 1:-1] - 6*V[1:-1,1:-1, 1:-1] + V[1:-1,2:, 1:-1] +
                                        V[2:  ,1:-1, 1:-1] +
                                                 + V[1:-1, 1:-1, 2:] )
                uvv = u*v*v
                u += (Du*Lu - uvv +  F*(1-u))
                v += (Dv*Lv + uvv - (F+k)*v)

                subverts = []
                add_vert = subverts.append

                vMin=V.min(); vMax=V.max()
                for iy in range(self.imgy):
                    for ix in range(self.imgx):
                        for iz in range(self.imgz):
                            w = V[iy, ix, iz]
                            c = int(255 * (w - vMin) / (vMax - vMin))
                            if c > 190:
                                add_vert((ix/40, iy/40, iz/40))
                self.store_frame(i, data=subverts)

                pn = int(100 * (i + 1) / steps) # percent completed
                if pn != p:
                    p = pn
                    print("%" + str(p).zfill(2))

            ######### loop end ############
            
            # label = "Du=" + str(Du) + " Dv=" + str(Dv) + " F=" + str(F) + " k=" + str(k)
            # print(label)

        def get_frame(self, number):
            return self.frame_storage.get(number)

        def store_frame(self, framestep, data):
            self.frame_storage[framestep] = data

    GK = DiffReact2(np, steps, seed)

verts = GK.get_frame(framenum)
verts_out.append(verts)



