
#include "global.h"
#include "maxeler/max_cg_API.h"
#include <complex.h>

static inline void _unit_matrix_clover(float complex c[18]){
    // See Eq. 3 of the Maxeler Conjugate Gradient Application on AWS F1 guide
    for(int ic=0; ic<18; ++ic)out[idx].c0[ic] = 0;
    // we set the diagonal terms to 1
    out[idx].c0[ 0] = 1.0f + I*1.0f;
    out[idx].c0[10] = 1.0f + I*1.0f;
    out[idx].c0[16] = 1.0f + I*1.0f;
   
}

void neutral_cloverh(cg_clover* out){
    const int VOLH = T*X*Y*Z/2;
    for(int idx=0; idx<VOLH; ++idx){
        _unit_matrix_clover(out[idx].c0);
        _unit_matrix_clover(out[idx].c1);
    }
}
