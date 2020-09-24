#ifndef MAXELER_GAUGE_H
#define MAXELER_GAUGE_H

#include "global.h"
#include "maxeler/max_cg_API.h"

su3* allocate_maxeler_gauge_field();
su3* allocate_maxeler_gauge_fieldEO();

void sombrero_to_maxeler_gauge_field_E(suNf_field* in, su3* outE);
void sombrero_to_maxeler_gauge_field_O(suNf_field* in, su3* outO);

void maxeler_to_sombrero_gauge_field_E(su3* inE, suNf_field* out);
void maxeler_to_sombrero_gauge_field_O(su3* inO, suNf_field* out);
#endif
