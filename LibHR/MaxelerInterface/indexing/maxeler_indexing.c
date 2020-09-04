#include "global.h"
#include <assert.h>
// maxeler maps

static int _lexord4d(int   x1, int   x2, int   x3, int   x4,
                     int max1, int max2, int max3, int max4){

    return x1*max2*max3*max4+
                x2*max3*max4+
                     x3*max4+
                          x4;
}

int maxeler_gaugeE_idx(int x, int y, int z, int t){
    return 4*_lexord4d(x, y, z, t,
                       X, Y, Z, T) / 2 ;
}

int maxeler_gaugeO_idx(int x, int y, int z, int t){
    return maxeler_gaugeE_idx(x, y, z, t);
}

int maxeler_spinorE_idx(int x, int y, int z, int t){
    return _lexord4d(x, y, z, t,
                     X, Y, Z, T) / 2;
}

int maxeler_spinorO_idx(int x, int y, int z, int t){
    return maxeler_spinorE_idx(x, y, z, t);
}

int maxeler_parity(int x, int y, int z, int t){
    // TODO: deal with this case
    assert(X%2 == 0 );
    assert(Y%2 == 0 );
    assert(Z%2 == 0 );
    assert(T%2 == 0 );
    return (x+y+z+t) % 2;
}

int maxeler_even(int x, int y, int z, int t){
    return maxeler_parity( x,  y,  z,  t) == 0;
}

int maxeler_odd(int x, int y, int z, int t){
    return maxeler_parity( x,  y,  z,  t) == 1;
}

int maxeler_any(int x, int y, int z, int t){
    return 1;
}

int maxeler_spinor_idx(int x, int y, int z, int t){
    int VOLSITES_HALF = (X*Y*Z*T/2);
    if      (maxeler_even(x,y,z,t)) return maxeler_spinorE_idx(x,y,z,t);
    else if ( maxeler_odd(x,y,z,t)) return maxeler_spinorO_idx(x,y,z,t)
                                           + VOLSITES_HALF;
    else return -1; //
}
