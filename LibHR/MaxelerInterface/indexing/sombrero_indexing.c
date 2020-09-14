#include "global.h"

// sombrero maps

int sombrero_gauge_idx(int x, int y, int z, int t){
    return 4*ipt(t,x,y,z);
}

int sombrero_spinor_idx(int x, int y, int z, int t){
    return ipt(t,x,y,z);
}
