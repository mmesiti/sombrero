#ifndef FLOP_COUNT_H
#define FLOP_COUNT_H
#include "suN.h"

#ifdef REPR_ADJOINT
static const float matrix_mul_flops = 4*NF*NF - 2*NF;
#else
static const float matrix_mul_flops = 8*NF*NF - 2*NF;
#endif


static const float real_sum_flops = 1;
static const float complex_sum_flops = 2;
static const float complex_mul_flops = 6;
static const float complex_mulR_flops = 2; // capital R for legibility (hopefully)
static const float complex_norm_flops = 3; // only re

// single vector
static const float vector_sum_flops = complex_sum_flops*NF;
static const float vector_mul_flops = complex_mulR_flops*NF;
static const float vector_norm_flops = complex_norm_flops*NF+real_sum_flops*NF; // + NF sums (not NF-1)

// single spinor 
static const float spinor_norm_flops = 4*(vector_norm_flops);
static const float spinor_sum_flops = 4*vector_sum_flops;
static const float spinor_mul_flops = 4*vector_mul_flops;
static const float spinor_mul_add_assign_flops = spinor_mul_flops + spinor_sum_flops;

// spinor field related flop count
static const float site_spinor_field_sqnorm_f_flops = 1.0f/2*spinor_norm_flops; // mult + add
static const float site_spinor_field_mul_add_assign_f_flops = 1.0f/2 * spinor_mul_add_assign_flops;
static const float site_spinor_field_mul_f_flops = 1.0f/2 * spinor_mul_flops;
static const float site_spinor_field_sub_f_flops = 1.0f/2*spinor_sum_flops;
static const float site_spinor_field_prod_re_flop = site_spinor_field_sqnorm_f_flops;
static const float site_spinor_field_add_assign_f_flops = 1.0f/2*spinor_sum_flops;


// see Dphi.c
static const float site_Dphi_flops = 1.0f/2*(16*matrix_mul_flops + 45*vector_sum_flops);
static const float site_Cphi_flops = 1.0f/2*( 8*matrix_mul_flops +  4*vector_sum_flops + spinor_mul_add_assign_flops);

// see Dphi_flt.c
static const float site_maxeler_fake_eopre_flt_flops = 2*(site_Dphi_flops + site_Cphi_flops) + site_spinor_field_sub_f_flops;
static const float site_maxeler_fake_eopre_sq_flt_flops = 2* site_maxeler_fake_eopre_flt_flops;

float cg_out_of_loop_flops_per_site(float site_operator_flops);

float cg_iteration_flops_per_site(float site_operator_flops);
#endif
