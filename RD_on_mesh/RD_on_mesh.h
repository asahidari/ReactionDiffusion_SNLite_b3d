
typedef struct _RdParams {
    double Du;
    double Dv;
    double F;
    double k;
    double dt;
} RdParams;

typedef void _callback_func(int step, int type, int count, double* data, void* pObj);

int process(RdParams* params, int num_verts, int* neighbor_counts, int** neighbor_indices, double* U, double* V, int steps, _callback_func cbf, void* pObj);

