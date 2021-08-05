#include <stdlib.h>
#include "Basic_RD_2d.h"

int process(struct rd_struct_t *data, double** U, double** V, int n, int steps, _callback_func callback_func, void* pObj) {

    // Allocate buffer memory
    double* u_buff = malloc(sizeof(double) * n * n);
    double* v_buff = malloc(sizeof(double) * n * n);
    double* res_buff = malloc(sizeof(double) * n * n * 3);

    // Loop for each steps
    for (int s = 0; s < steps; s++) {

        double vmin = 1000.0, vmax = -1000.0;
        for (int i = 1; i < n + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                // Calculate Laplacian U/V
                double lapU = U[i][j-1] + U[i][j+1] + U[i-1][j] + U[i+1][j] - 4.0 * U[i][j];
                double lapV = V[i][j-1] + V[i][j+1] + V[i-1][j] + V[i+1][j] - 4.0 * V[i][j];
                
                // Calculate Reacton Diffusion
                double uvv = U[i][j] * V[i][j] * V[i][j];
                u_buff[(i-1)*n + (j-1)] = U[i][j] + data->Du * lapU - uvv + data->F * (1.0 - U[i][j]);
                v_buff[(i-1)*n + (j-1)] = V[i][j] + data->Dv * lapV + uvv - (data->F + data->k) * V[i][j];
                
                // Store minimal value of v_buff
                if (vmin > v_buff[(i-1)*n+(j-1)])
                    vmin = v_buff[(i-1)*n+(j-1)];
                if (vmax < v_buff[(i-1)*n+(j-1)])
                    vmax = v_buff[(i-1)*n+(j-1)];
            }
        }

        int count = 0;
        for (int iy = 0; iy < n; iy++) {
            for (int ix = 0; ix < n; ix++) {
                // Select V values over a threshold, and normalize them
                double w = v_buff[iy*n+ix];
                int c = (int)(255.0 * (w - vmin) / (vmax - vmin));
                if (c > 190) {
                    res_buff[count * 3 + 0] = ix/40.0;
                    res_buff[count * 3 + 1] = iy/40.0;
                    res_buff[count * 3 + 2] = 0.0;
                    count++;
                }
            }
        }

        if (count > 0) {
            // Return results for each step
            callback_func(s, count, res_buff, pObj);
        }

        // Pass the calculated data to next step
        for (int i = 1; i < n + 1; i++) {
            for (int j = 1; j < n + 1; j++) {
                U[i][j] = u_buff[(i-1) * n + (j-1)];
                V[i][j] = v_buff[(i-1) * n + (j-1)];
            }
        }
    }

    // Release buffer memory
    free(u_buff);
    free(v_buff);
    free(res_buff);

    return 0;
}
