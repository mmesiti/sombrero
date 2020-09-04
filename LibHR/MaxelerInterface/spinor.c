#include "global.h"
#include "maxeler/max_cg_API.h"
#include "indexing/reindex.h"
#include "indexing/maxeler_indexing.h"
#include "indexing/sombrero_indexing.h"
#include "element_transforms/spinor_transforms.h"
#include "spinor_field.h"

void sombrero_to_maxeler_spinor_field(const spinor_field* in, cg_spinor* out){

    suNf_spinor* inptr0 = _FIELD_AT(in,sombrero_spinor_idx(0,0,0,0));
    reindex((void*) inptr0, (void*) out,
            sizeof(suNf_spinor),
            sizeof(cg_spinor),
            0,X,
            0,Y,
            0,Z,
            0,T,
            sombrero_spinor_idx,
            maxeler_spinor_idx,
            sombrero_to_maxeler_spinor,
            maxeler_any);

}

void maxeler_to_sombrero_spinor_field(const cg_spinor* in, spinor_field* out){

    suNf_spinor* outptr0 = _FIELD_AT(out,sombrero_spinor_idx(0,0,0,0));
    reindex((void*) in, (void*) outptr0,
            sizeof(cg_spinor),
            sizeof(suNf_spinor),
            0,X,
            0,Y,
            0,Z,
            0,T,
            maxeler_spinor_idx,
            sombrero_spinor_idx,
            maxeler_to_sombrero_spinor,
            maxeler_any);

}
