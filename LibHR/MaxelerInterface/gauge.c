#include "global.h"
#include "maxeler/max_cg_API.h"
#include "indexing/reindex.h"
#include "indexing/maxeler_indexing.h"
#include "indexing/sombrero_indexing.h"
#include "element_transforms/gauge_transforms.h"

void sombrero_to_maxeler_gauge_field_E(suNg_field* in, su3* outE){
    suNg* inptr0 = in->ptr+sombrero_gauge_idx(0,0,0,0);
    reindex((void*) inptr0, (void*) outE,
            sizeof(suNg),
            sizeof(su3),
            0,X,
            0,Y,
            0,Z,
            0,T,
            sombrero_gauge_idx,
            maxeler_gaugeE_idx,
            sombrero_to_maxeler_gauge,
            maxeler_even);
}

void sombrero_to_maxeler_gauge_field_O(suNg_field* in, su3* outO){
    suNg* inptr0 = in->ptr+sombrero_gauge_idx(0,0,0,0);
    reindex((void*) inptr0, (void*) outO,
            sizeof(suNg),
            sizeof(su3),
            0,X,
            0,Y,
            0,Z,
            0,T,
            sombrero_gauge_idx,
            maxeler_gaugeO_idx,
            sombrero_to_maxeler_gauge,
            maxeler_odd);
}
