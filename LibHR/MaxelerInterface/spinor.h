#ifndef MAXELER_SPINOR_H
#define MAXELER_SPINOR_H
#include "spinor_field.h"
#include "maxeler/max_cg_API.h"

void sombrero_to_maxeler_spinor_field(const spinor_field* in, cg_spinor* out);

void maxeler_to_sombrero_spinor_field(const cg_spinor* in, spinor_field* out);
#endif
