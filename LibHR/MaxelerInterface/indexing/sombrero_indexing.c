#include "global.h"

// sombrero maps

int sombrero_gaugequartet_idx(int x, int y, int z, int t){
    return ipt(t,x,y,z);
}

int sombrero_spinor_idx(int x, int y, int z, int t){
    return ipt(t,x,y,z);
}
