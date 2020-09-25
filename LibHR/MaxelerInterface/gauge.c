#include "global.h"
#include "maxeler/max_cg_API.h"
#include "indexing/reindex.h"
#include "indexing/maxeler_indexing.h"
#include "indexing/sombrero_indexing.h"
#include "element_transforms/gauge_transforms.h"
#include <stdlib.h>

static su3* _allocate_maxeler_gauge_field_with_size(int size){
    return malloc(size*sizeof(su3));
}


su3* allocate_maxeler_gauge_field(){
    const int GAUGESIZEH = 4*(T*X*Y*Z)/2;
    return _allocate_maxeler_gauge_field_with_size(GAUGESIZEH);
}

su3* allocate_maxeler_gauge_fieldEO(){
    const int GAUGESIZE = 4*(T*X*Y*Z);
    return _allocate_maxeler_gauge_field_with_size(GAUGESIZE);
}

void sombrero_to_maxeler_gauge_field_E(suNf_field* in, su3* outE){
    suNf* inptr0 = in->ptr+sombrero_gaugequartet_idx(0,0,0,0);
    reindex((void*) inptr0, (void*) outE,
            4*sizeof(suNf),
            4*sizeof(su3),
            0,X,
            0,Y,
            0,Z,
            0,T,
            sombrero_gaugequartet_idx,
            maxeler_gaugequartetE_idx,
            sombrero_to_maxeler_gauge4,
            maxeler_even);
}

void sombrero_to_maxeler_gauge_field_O(suNf_field* in, su3* outO){
    suNf* inptr0 = in->ptr+sombrero_gaugequartet_idx(0,0,0,0);
    reindex((void*) inptr0, (void*) outO,
            4*sizeof(suNf),
            4*sizeof(su3),
            0,X,
            0,Y,
            0,Z,
            0,T,
            sombrero_gaugequartet_idx,
            maxeler_gaugequartetO_idx,
            sombrero_to_maxeler_gauge4,
            maxeler_odd);
}

void maxeler_to_sombrero_gauge_field_E(su3* inE, suNf_field* out){
    suNf* outptr0 = out->ptr+sombrero_gaugequartet_idx(0,0,0,0);
    reindex((void*) inE, (void*) outptr0,
            4*sizeof(su3),
            4*sizeof(suNf),
            0,X,
            0,Y,
            0,Z,
            0,T,
            maxeler_gaugequartetE_idx,
            sombrero_gaugequartet_idx,
            maxeler_to_sombrero_gauge4,
            maxeler_even);
}

void maxeler_to_sombrero_gauge_field_O(su3* inO, suNf_field* out){
    suNf* outptr0 = out->ptr+sombrero_gaugequartet_idx(0,0,0,0);
    reindex((void*) inO, (void*) outptr0,
            4*sizeof(su3),
            4*sizeof(suNf),
            0,X,
            0,Y,
            0,Z,
            0,T,
            maxeler_gaugequartetE_idx,
            sombrero_gaugequartet_idx,
            maxeler_to_sombrero_gauge4,
            maxeler_odd);
}
