#include "su3_types.h"
#include "../maxeler/max_cg_API.h"
#include <complex.h>

// sombrero-maxeler element transformations
void sombrero_to_maxeler_gauge(void* _sombrero_gauge,
                               void*  _maxeler_gauge){
    su3*   maxeler_gauge =  _maxeler_gauge;
    suNg* sombrero_gauge = _sombrero_gauge;

#define _stdcmplx(x) ((x).re + I * (x).im)
    maxeler_gauge->c00 = _stdcmplx(sombrero_gauge->c[0]);
    maxeler_gauge->c01 = _stdcmplx(sombrero_gauge->c[1]);

    maxeler_gauge->c10 = _stdcmplx(sombrero_gauge->c[3]);
    maxeler_gauge->c11 = _stdcmplx(sombrero_gauge->c[4]);

    maxeler_gauge->c20 = _stdcmplx(sombrero_gauge->c[6]);
    maxeler_gauge->c21 = _stdcmplx(sombrero_gauge->c[7]);
#undef _stdcmplx

}

void sombrero_to_maxeler_gauge4(void* _sombrero_gauge4,
                                void*  _maxeler_gauge4){
    int mu;
    const int X = 0, Y = 1, Z = 2, T = 3;
    const int SOMBRERO_ORDERING[4] = {T,X,Y,Z};
    const int  MAXELER_ORDERING[4] = {X,Y,Z,T};

    for(mu=0;mu<4;++mu)
        sombrero_to_maxeler_gauge(
            _sombrero_gauge4 + SOMBRERO_ORDERING[mu]*sizeof(suNg),
             _maxeler_gauge4 +  MAXELER_ORDERING[mu]*sizeof( su3)
        );
}
