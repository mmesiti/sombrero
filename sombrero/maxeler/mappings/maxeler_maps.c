#include "global.h"
 
// maxeler maps

static int _lexord4d(int   x1, int   x2, int   x3, int   x4,
                     int max1, int max2, int max3, int max4){

    return x1*max2*max3*max4+
                x2*max3*max4+
                     x3*max4+
                          x4;
}

int maxeler_gaugeE_map(int x, int y, int z, int t){
    return 4*_lexord4d(x, y, z, t,
                       X, Y, Z, T) / 2 ;
}

int maxeler_gaugeO_map(int x, int y, int z, int t){
    return maxeler_gaugeE_map(x, y, z, t);
}

int maxeler_spinorE_map(int x, int y, int z, int t){
    return _lexord4d(t, x, y, z,
                     T, X, Y, Z) / 2;
}

int maxeler_spinorO_map(int x, int y, int z, int t){
    return maxeler_spinorE_map(x, y, z, t);
}

int maxeler_parity(int x, int y, int z, int t){
    return (x+y+z+t) % 2;
}

int maxeler_spinor_map(int x, int y, int z, int t){
    int VOLSITES_HALF = (X*Y*Z*T/2);
    if      (maxeler_parity(x,y,z,t) == 0) return maxeler_spinorE_map(x,y,z,t);
    else if (maxeler_parity(x,y,z,t) == 1) return maxeler_spinorO_map(x,y,z,t) + VOLSITES_HALF;
    else return -1; //
}
