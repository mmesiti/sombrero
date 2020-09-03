#ifndef MAXELER_GAUGE_H
#define MAXELER_GAUGE_H

#include "global.h"
#include "maxeler/max_cg_API.h"

void sombrero_to_maxeler_gauge_field_E(suNg_field* in, su3* outE);

void sombrero_to_maxeler_gauge_field_O(suNg_field* in, su3* outO);

#endif
