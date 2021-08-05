
struct rd_struct_t {
    double Du;
    double Dv;
    double F;
    double k;
};

typedef void _callback_func(int step, int len, double* data, void* pObj);

int process(struct rd_struct_t *data, double** U, double** V, int n, int steps, _callback_func callback_func, void* pObj);
