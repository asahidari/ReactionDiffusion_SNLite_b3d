#include "RD_on_mesh.h"
#include <stdlib.h>

int process(RdParams* params, int num_verts, int* neighbor_counts, int** neighbor_indices, double* U, double* V, int steps, _callback_func cbf, void* pObj) {

    // Allocate buffer memory
    double* u_buff = (double*)malloc(sizeof(double) * num_verts);
    double* v_buff = (double*)malloc(sizeof(double) * num_verts);
    double* a_out = (double*)malloc(sizeof(double) * num_verts);
    double* b_out = (double*)malloc(sizeof(double) * num_verts);

    for (int s = 0; s < steps; s++) {

        double umin = 100.0, umax = -100.0;
        double vmin = 100.0, vmax = -100.0;

        for (int i = 0; i < num_verts; i++) {

            // Calc laplacian from neighbor substances
            double lapU = 0;
            double lapV = 0;
            int nb_count = neighbor_counts[i];
            for (int j = 0; j < nb_count; j++) {
                lapU += U[neighbor_indices[i][j]];
                lapV += V[neighbor_indices[i][j]];
            }
            lapU = (lapU / (double)neighbor_counts[i]) - U[i];
            lapV = (lapV / (double)neighbor_counts[i]) - V[i];

            // Reaction Diffusion
            double uvv = U[i] * V[i] * V[i];
            u_buff[i] = U[i] + (params->Du * lapU - uvv + params->F * (1.0 - U[i])) * params->dt;
            v_buff[i] = V[i] + (params->Dv * lapV + uvv - (params->F + params->k) * V[i]) * params->dt;

            // Store minimal value of u_buff and v_buff
            if (umax < u_buff[i])
                umax = u_buff[i];
            if (umin > u_buff[i])
                umin = u_buff[i];
            if (vmax < v_buff[i])
                vmax = v_buff[i];
            if (vmin > v_buff[i])
                vmin = v_buff[i];
        }

        if (umax > 1.0) umax = 1.0;
        if (umin < 0.0) umin = 0.0;
        if (vmax > 1.0) vmax = 1.0;
        if (vmin < 0.0) vmin = 0.0;

        for (int i = 0; i < num_verts; i++) {

            // Round the values
            if (u_buff[i] > 1.0) u_buff[i] = 1.0;
            if (u_buff[i] < 0.0) u_buff[i] = 0.0;
            if (v_buff[i] > 1.0) v_buff[i] = 1.0;
            if (v_buff[i] < 0.0) v_buff[i] = 0.0;

            // Normalize results
            double u_denominator = (umax > umin) ? (umax - umin) : 0.01;
            double v_denominator = (vmax > vmin) ? (vmax - vmin) : 0.01;
            a_out[i] = (u_buff[i] - umin) / u_denominator;
            b_out[i] = (v_buff[i] - vmin) / v_denominator;

            // Prepare for next step
            U[i] = u_buff[i];
            V[i] = v_buff[i];
        }

        // Send results in each step
        cbf(s, 0, num_verts, a_out, pObj);
        cbf(s, 1, num_verts, b_out, pObj);
    }

    // Release buffer memory
    free(u_buff);
    free(v_buff);
    free(a_out);
    free(b_out);

    return 0;
}
