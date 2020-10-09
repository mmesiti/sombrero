#include "flop_count.h"

static int initialize = 1;

static float matrix_mul_flops;

static float real_op_flop;
static float complex_sum_flops;
static float complex_mul_flops;
static float complex_mulR_flops; // capital R for legibility (hopefully)
static float complex_norm_flops;

// single vector
static float vector_sum_flops;
static float vector_mul_flops;
static float vector_norm_flops;

// single spinor 
static float spinor_norm_flops;
static float spinor_sum_flops;
static float spinor_mul_flops;
static float spinor_mul_add_assign_flops;

// spinor field related flop count
static float site_spinor_field_sqnorm_f_flops;
static float site_spinor_field_mul_add_assign_f_flops;
static float site_spinor_field_mul_f_flops;
static float site_spinor_field_sub_f_flops;
static float site_spinor_field_prod_re_flop ;
static float site_spinor_field_add_assign_f_flops;
static float site_spinor_field_minus_f_flops;
// see Dphi.c
static float site_Dphi_flops;
static float site_Cphi_flops;
static float site_Cphi_assign_flops;
static float site_Cphi_inv_flops;

// see Dphi_flt.c


static void init(){
    if(initialize == 0) return;
#ifdef REPR_ADJOINT
    matrix_mul_flops = 4*NF*NF - 2*NF;
#else
    matrix_mul_flops = 8*NF*NF - 2*NF;
#endif

    real_op_flop = 1;
    complex_sum_flops = 2;
    complex_mul_flops = 6;
    complex_mulR_flops = 2; // capital R for legibility (hopefully)
    complex_norm_flops = 3; // only re

    // single vector
    vector_sum_flops = complex_sum_flops*NF;
    vector_mul_flops = complex_mulR_flops*NF; // mul uses mulR
    vector_norm_flops = complex_norm_flops*NF+real_op_flop*NF; // + NF sums (not NF-1)

    // single spinor 
    spinor_norm_flops = 4*(vector_norm_flops);
    spinor_sum_flops = 4*vector_sum_flops;
    spinor_mul_flops = 4*vector_mul_flops; // uses mulR
    spinor_mul_add_assign_flops = spinor_mul_flops + spinor_sum_flops;

    // spinor field related flop count
    site_spinor_field_sqnorm_f_flops = 1.0f/2 * spinor_norm_flops; // mult + add
    site_spinor_field_mul_add_assign_f_flops = 1.0f/2 * spinor_mul_add_assign_flops;
    site_spinor_field_mul_f_flops = 1.0f/2 * spinor_mul_flops;
    site_spinor_field_sub_f_flops = 1.0f/2 * spinor_sum_flops;
    site_spinor_field_prod_re_flop = site_spinor_field_sqnorm_f_flops;
    site_spinor_field_add_assign_f_flops = 1.0f/2 * spinor_sum_flops;
    site_spinor_field_minus_f_flops = 1.0f/2 * spinor_mul_flops; // as a mul

    // see Dphi.c
    site_Dphi_flops = 1.0f/2*(16*matrix_mul_flops + 45*vector_sum_flops);

    site_Cphi_flops = 1.0f/2*( 8*matrix_mul_flops +
                               4*vector_sum_flops +
                               spinor_mul_add_assign_flops);

    site_Cphi_assign_flops = 1.0f/2*( 8*matrix_mul_flops +
                                      4*vector_sum_flops +
                                      spinor_mul_add_assign_flops +
                                      spinor_sum_assign_flops);

    {
        int N = 2*NF;
        int forward_substitution_loop_flop = N*(N-1)/2 * ((3 * real_op_flop) + // n = i*(i+1)/2+i;
                                                     // _complex_mul_sub_assign();
                                                     (complex_mul_flops + complex_sum_flops) +
                                                     // _complex_mul_sub_assign();
                                                     (complex_mul_flops + complex_sum_flops));

        int backward_substitution_loop_flop = 0
		for(int i = N-1; i >= 0; i--)
		{
            backward_substitution_loop_flop += 3; // n = i*(i+1)/2+i;
			backward_substitution_loop_flop += complex_mulR_flops;// _complex_mulr();
			backward_substitution_loop_flop += complex_mulR_flops;// _complex_mulr();
			for(int k = i+1; k < N; k++)
			{
                backward_substitution_loop_flop += 3; // n = k*(k+1)/2+i;

                // _complex_mul_sub_assign();
				backward_substitution_loop_flop += complex_mul_flops + complex_sum_flops;
                // _complex_mul_sub_assign();
				backward_substitution_loop_flop += complex_mul_flops + complex_sum_flops;
			}
		}

        site_Cphi_inv_flops = 1.0f/2*(forward_substitution_loop_flop +
                                      backward_substitution_loop_flop);

        site_Cphi_inv_assign_flops = 1.0f/2*(forward_substitution_loop_flop +
                                             backward_substitution_loop_flop +
                                             spino_sum_assign_flops);
    }

    // see Dphi_flt.c
    initialize = 0;
}


float site_maxeler_fake_eopre_sq_flt_flops(){
    init();
    float site_maxeler_fake_eopre_flt_flops = ( site_Dphi_flops +
                                                site_Cphi_flops +
                                                site_Dphi_flops +
                                                site_Cphi_flops ) +
                                              site_spinor_field_sub_f_flops;

    return 2*site_maxeler_fake_eopre_flt_flops;
}

float site_g5Cphi_eopre_sq_flops(){
    init();
    float site_g5Cphi_eopre_flops = site_Dphi_flops +
                                    site_Cphi_inf_flops +
                                    site_Dphi_flops +
                                    site_spinor_field_minus_f_flops +
                                    site_Cphi_assign_flops;

    return 2*site_g5Cphi_eopre_flops;

}

float cg_out_of_loop_flops_per_site(float site_operator_flops){
    init();
                                                      // line in cg_test:
    return site_spinor_field_sqnorm_f_flops +         //
           site_operator_flops +                      //
           site_spinor_field_mul_add_assign_f_flops + //
           site_spinor_field_sub_f_flops +            //
           site_spinor_field_sqnorm_f_flops +         //
           site_operator_flops +                      //
           site_spinor_field_mul_add_assign_f_flops+  //
           site_spinor_field_sub_f_flops +            //
           2*site_spinor_field_sqnorm_f_flops;        //
}

float cg_iteration_flops_per_site(float site_operator_flops) {
    init();
                                                      // line in cg_test:
    return site_operator_flops +                      //
           site_spinor_field_prod_re_flop +           //
           site_spinor_field_mul_add_assign_f_flops + //
           site_spinor_field_mul_add_assign_f_flops + //
           site_spinor_field_sqnorm_f_flops +         //
           site_spinor_field_mul_f_flops +            //
           site_spinor_field_add_assign_f_flops +     //
           site_spinor_field_mul_f_flops +            //
           site_spinor_field_add_assign_f_flops;      //
}
