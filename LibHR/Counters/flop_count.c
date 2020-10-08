#include "flop_count.h"

float cg_out_of_loop_flops_per_site(float site_operator_flops){
                                                      // line in cg_test:
    return site_spinor_field_sqnorm_f_flops +         //
           site_operator_flops +                      //
           site_spinor_field_mul_add_assign_f_flops+ //
           site_spinor_field_sub_f_flops +            //
           site_spinor_field_sqnorm_f_flops +         //
           site_operator_flops +                      //
           site_spinor_field_mul_add_assign_f_flops+  //
           site_spinor_field_sub_f_flops +            //
           2*site_spinor_field_sqnorm_f_flops;        //
}

float cg_iteration_flops_per_site(float site_operator_flops) {
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
