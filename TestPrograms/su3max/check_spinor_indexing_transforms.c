#define MAIN_PROGRAM
#include "suN_types.h"
#include "global.h"
#include "LibHR/MaxelerInterface/maxeler/max_cg_API.h"
#include "LibHR/MaxelerInterface/spinor.h"
#include "spinor_field.h"
#include "memory.h"
#include <assert.h>
#include <stdio.h>


static void glb_setup(){
    NP_T=NP_X=NP_Y=NP_Z=1;
    GLB_T=GLB_X=GLB_Y=GLB_Z=8;
    T=X=Y=Z=8;
}

int main(int argc, char* argv[]){

    {
        glb_setup();

        setup_process(&argc,&argv);

        assert(WORLD_SIZE==1);

        /* setup lattice geometry */
        if (geometry_init() == 1) { finalize_process(); return 0; }
        geometry_mpi_eo();

    }

    spinor_field* sombrero_in = alloc_spinor_field_f(1,&glat_even);
    spinor_field* sombrero_out = alloc_spinor_field_f(1,&glat_even);
    cg_spinor* maxeler_intermediate = allocate_maxeler_spinor_field();

    {// round-trip test

        { // sombrero_in initialisation
            const int max_index = T*X*Y*Z // volume
                                  /2      // Even only
                                  *4      // dirac indices
                                  *3      // colors, SU(3)
                                  *2;     // Real/Imaginary
            _MASTER_FOR(sombrero_in->type,ix){
                suNf_spinor* sombrero_in_tmp = _FIELD_AT(sombrero_in,ix);
                for(int d=0;d<4;++d)
                    for(int c=0;c<3;++c){
                        sombrero_in_tmp->c[d].c[c].re = ((float) (ix*4*3*2 +
                                                                      d*3*2 +
                                                                        c*2 +
                                                                          0)
                                                         /max_index);
                        sombrero_in_tmp->c[d].c[c].im = ((float) (ix*4*3*2 +
                                                                      d*3*2 +
                                                                        c*2 +
                                                                          1)
                                                         /max_index);
                    }
            }
        }

        sombrero_to_maxeler_spinor_field(sombrero_in,maxeler_intermediate);

        maxeler_to_sombrero_spinor_field(maxeler_intermediate,sombrero_out);
        {

            _MASTER_FOR(sombrero_in->type,ix){
                suNf_spinor*  sombrero_in_tmp =  _FIELD_AT(sombrero_in,ix);
                suNf_spinor* sombrero_out_tmp = _FIELD_AT(sombrero_out,ix);
                for(int d=0;d<4;++d)
                    for(int c=0;c<3;++c){
                        double  in_re =  sombrero_in_tmp->c[d].c[c].re;
                        double out_re = sombrero_out_tmp->c[d].c[c].re;

                        if(in_re != out_re){
                            printf("Error: ix=%d, d=%d, c=%d, out_re: %.18lf , in_re: %.18lf\n",
                                   ix,d,c,out_re,in_re);
                            assert(0);
                        }

                        double  in_im =  sombrero_in_tmp->c[d].c[c].im;
                        double out_im = sombrero_out_tmp->c[d].c[c].im;

                        if(in_im != out_im){
                            printf("Error: ix=%d, d=%d, c=%d, out_im: %.18lf , in_im: %.18lf\n",
                                   ix,d,c,out_im,in_im);
                            assert(0);
                        }
                    }
            }
        }
    }

    free_spinor_field_f(sombrero_in);
    free_spinor_field_f(sombrero_out);
    free(maxeler_intermediate);

    finalize_process();
    return 0;

}
