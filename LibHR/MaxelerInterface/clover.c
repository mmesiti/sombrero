
#include "global.h"
#include "maxeler/max_cg_API.h"
#include <complex.h>
#include <stdlib.h>

cg_clover* allocate_maxeler_clover_field(){
  const int LOCVOLH = (T*X*Y*Z)/2;
  return malloc( LOCVOLH * sizeof(cg_clover));
}

static inline void _unit_matrix_clover(float complex c[18]){
    // See Eq. 3 of the Maxeler Conjugate Gradient Application on AWS F1 guide
    for(int ic=0; ic<18; ++ic) c[ic] = 0;
    // we set the diagonal terms to 0.5
    c[ 0] = 0.5f + I*0.5f;
    c[10] = 0.5f + I*0.5f;
    c[16] = 0.5f + I*0.5f;
   
}

void neutral_cloverh(cg_clover* out){
    const int VOLH = T*X*Y*Z/2;
    for(int idx=0; idx<VOLH; ++idx){
        _unit_matrix_clover(out[idx].c0);
        _unit_matrix_clover(out[idx].c1);
    }
}
