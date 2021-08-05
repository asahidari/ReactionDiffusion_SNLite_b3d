#include <stdlib.h>
#include "Basic_RD_3d.h"

int process(struct rd_struct_t *data, double* U, double* V, int n, int steps, _callback_func callback_func, void* pObj) {

    // Allocate buffer memory
    double* u_buf = malloc(sizeof(double) * n * n * n);
    double* v_buf = malloc(sizeof(double) * n * n * n);
    double* res_buf = malloc(sizeof(double) * n * n * n * 3);

    // Loop for each steps
    int np2 = n + 2;
    for (int s = 0; s < steps; s++) {
        double vmin = 1000.0, vmax = -1000.0;
        for (int i = 1; i < n + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                for (int k = 1; k < n + 1; k++) {
                    // Calculate Laplacian U/V
                    double lapU = U[i*np2*np2 + (j-1)*np2 + k] + U[i*np2*np2 + (j+1)*np2 + k]
                        + U[(i-1)*np2*np2 + j*np2 + k] + U[(i+1)*np2*np2 + j*np2 + k]
                        + U[i*np2*np2 + j*np2 + (k-1)] + U[i*np2*np2 + j*np2 + (k+1)] - 6.0 * U[i*np2*np2 + j*np2 + k];
                    double lapV = V[i*np2*np2 + (j-1)*np2 + k] + V[i*np2*np2 + (j+1)*np2 + k]
                        + V[(i-1)*np2*np2 + j*np2 + k] + V[(i+1)*np2*np2 + j*np2 + k]
                        + V[i*np2*np2 + j*np2 + (k-1)] + V[i*np2*np2 + j*np2 + (k+1)] - 6.0 * V[i*np2*np2 + j*np2 + k];

                    // Calculate Reacton Diffusion
                    double uvv = U[i*np2*np2 + j*np2 + k] * V[i*np2*np2 + j*np2 + k] * V[i*np2*np2 + j*np2 + k];
                    u_buf[(i-1)*n*n + (j-1)*n + (k-1)] = U[i*np2*np2 + j*np2 + k] + data->Du * lapU - uvv + data->F * (1.0 - U[i*np2*np2 + j*np2 + k]);
                    v_buf[(i-1)*n*n + (j-1)*n + (k-1)] = V[i*np2*np2 + j*np2 + k] + data->Dv * lapV + uvv - (data->F + data->k) * V[i*np2*np2 + j*np2 + k];

                    // Store minimal value of v_buff
                    if (vmin > v_buf[(i-1)*n*n + (j-1)*n + (k-1)])
                        vmin = v_buf[(i-1)*n*n + (j-1)*n + (k-1)];
                    if (vmax < v_buf[(i-1)*n*n + (j-1)*n + (k-1)])
                        vmax = v_buf[(i-1)*n*n + (j-1)*n + (k-1)];
                }
            }
        }

        int count = 0;
        for (int iy = 0; iy < n; iy++) {
            for (int ix = 0; ix < n; ix++) {
                for (int iz = 0; iz < n; iz++) {
                    // Select V values over a threshold, and normalize them
                    double w = v_buf[iy*n*n + ix*n + iz];
                    int c = (int)(255.0 * (w - vmin) / (vmax - vmin));
                    if (c > 190) {
                        res_buf[count * 3 + 0] = ix/40.0;
                        res_buf[count * 3 + 1] = iy/40.0;
                        res_buf[count * 3 + 2] = iz/40.0;
                        count++;
                    }
                }
            }
        }
        if (count > 0) {
            // Return results for each step
            callback_func(s, count, res_buf, pObj);
        }

        // Pass the calculated data to next step
        for (int i = 1; i < n + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                for (int k = 1; k < n + 1; k++) {
                    U[i*np2*np2 + j*np2 + k] = u_buf[(i-1)*n*n + (j-1)*n + (k-1)];
                    V[i*np2*np2 + j*np2 + k] = v_buf[(i-1)*n*n + (j-1)*n + (k-1)];
                }
            }
        }
    }

    // Release buffer memory
    free(u_buf);
    free(v_buf);
    free(res_buf);

    return 0;
}
