#include "su3_types.h"
#include "max_cg_API.h"

// sombrero-maxeler element transformations
void sombrero_to_maxeler_gauge(void* _sombrero_gauge,
                               void*  _maxeler_gauge){
    su3*   maxeler_gauge =  _maxeler_gauge;
    suNg* sombrero_gauge = _sombrero_gauge;

    maxeler_gauge.c00 = sombrero_gauge.c[0];
    maxeler_gauge.c01 = sombrero_gauge.c[1];

    maxeler_gauge.c10 = sombrero_gauge.c[3];
    maxeler_gauge.c11 = sombrero_gauge.c[4];

    maxeler_gauge.c20 = sombrero_gauge.c[6];
    maxeler_gauge.c21 = sombrero_gauge.c[7];

}

void sombrero_to_maxeler_gauge4(void* _sombrero_gauge4,
                                void*  _maxeler_gauge4){
    int mu;
    const int X = 0, Y = 1, Z = 2, T = 3;
    const int SOMBRERO_ORDERING[4] = {T,X,Y,Z};
    const int  MAXELER_ORDERING[4] = {X,Y,Z,T};

    for(mu=0;mu<4;++mu)
        sombrero_to_maxeler_gauge(
            _sombrero_gauge4 + SOMBRERO_ORDERING[i]*sizeof(suNg),
             _maxeler_gauge4 +  MAXELER_ORDERING[i]*sizeof( su3)
        );
}
