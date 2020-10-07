#ifdef WITH_CLOVER
#include "global.h"
#include "field_ordering.h"
#include "memory.h"
#define _4FIELD_AT(s,i,mu) (((s)->ptr)+coord_to_index(i-(s)->type->master_shift,mu))
// Fakery for comparison with MaxCG
static void _zero_clover_term_flt(int id)
{
	for(int ij = 0; ij < NF*NF; ij++)
	{
		for(int mu = 0; mu < 4; ++mu ){
			_4FIELD_AT(cl_term_flt,id,mu)->c[ij].re = 0;
			_4FIELD_AT(cl_term_flt,id,mu)->c[ij].im = 0;
		}
	}
}

// Fakery for comparison with MaxCG
void zero_clover_term_flt()
{
	_MASTER_FOR(&glattice,id)
	{
		_zero_clover_term_flt(id);
	}
}

void clover_init_flt_fake(){
	cl_term_flt = alloc_clover_term_flt(&glattice);
}
#endif // #ifdef WITH_CLOVER
