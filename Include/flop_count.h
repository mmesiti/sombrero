#ifndef FLOP_COUNT_H
#define FLOP_COUNT_H
#include "suN.h"

float site_maxeler_fake_eopre_sq_flt_flops();

float site_g5Cphi_eopre_sq_flops();

float cg_out_of_loop_flops_per_site(float site_operator_flops);

float cg_iteration_flops_per_site(float site_operator_flops);
#endif
