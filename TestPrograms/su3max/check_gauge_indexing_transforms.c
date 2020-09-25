#define MAIN_PROGRAM
#include "suN_types.h"
#include <complex.h>
#undef complex // ugly hack, but necessary
#include "hrcomplex.h"
#include "global.h"
#include "LibHR/MaxelerInterface/maxeler/max_cg_API.h"
#include "LibHR/MaxelerInterface/gauge.h"
#include "memory.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "random.h"


static void glb_setup(){
    NP_T=NP_X=NP_Y=NP_Z=1;
    GLB_T=GLB_X=GLB_Y=GLB_Z=8;
    T=X=Y=Z=8;
}

void random_uf(suNf_field* sombrero_in)
{
    _MASTER_FOR(sombrero_in->type,ix){
        for(int mu=0;mu<4;++mu){
            suNf* sombrero_in_tmp = sombrero_in->ptr + coord_to_index(ix,mu);
            random_suNg((suNg*)sombrero_in_tmp);
        }
    }
}

int main(int argc, char* argv[]){

    { // setup
        glb_setup();

        setup_process(&argc,&argv);

        assert(WORLD_SIZE==1);

        /* setup lattice geometry */
        if (geometry_init() == 1) { finalize_process(); return 0; }
        geometry_mpi_eo();

    }

    //su3_generation_checks();

    { // one-way orthonormality check - even

        suNf_field* sombrero_in = alloc_gfield_f(&glattice);
        su3* maxelerE = allocate_maxeler_gauge_field();

        random_uf(sombrero_in);

        printf("Initialisation Ok\n");

        sombrero_to_maxeler_gauge_field_E(sombrero_in,maxelerE);

        printf("Translation Ok\n");

        { // check that all the matrices are ok
            for(int i=0;i<T*X*Y*Z/2*4;++i){
                su3 t  = maxelerE[i];
                float norm0err = crealf(t.c00*conj(t.c00)+ t.c10*conj(t.c10)+ t.c20*conj(t.c20))-1;
                assert(norm0err*norm0err < 1.0e-8);
                float norm1err = crealf(t.c01*conj(t.c01)+ t.c11*conj(t.c11)+ t.c21*conj(t.c21))-1;
                assert(norm1err*norm1err < 1.0e-8);
                float orthoerr = cabsf(t.c01*conj(t.c00)+ t.c11*conj(t.c10)+ t.c21*conj(t.c20));
                assert(orthoerr < 1.0e-7);
            }
        }
        free_gfield_f(sombrero_in);
        free(maxelerE);
        printf("One-way orthonormality check - even: Ok\n");

    }


    { // one-way orthonormality check - odd

        suNf_field* sombrero_in = alloc_gfield_f(&glattice);
        su3* maxelerO = allocate_maxeler_gauge_field();

        random_uf(sombrero_in);

        printf("Initialisation Ok\n");

        sombrero_to_maxeler_gauge_field_O(sombrero_in,maxelerO);

        printf("Translation Ok\n");

        { // check that all the matrices are ok
            for(int i=0;i<T*X*Y*Z/2*4;++i){
                su3 t  = maxelerO[i];
                float norm0err = crealf(t.c00*conj(t.c00)+ t.c10*conj(t.c10)+ t.c20*conj(t.c20))-1;
                assert(norm0err*norm0err < 1.0e-8);
                float norm1err = crealf(t.c01*conj(t.c01)+ t.c11*conj(t.c11)+ t.c21*conj(t.c21))-1;
                assert(norm1err*norm1err < 1.0e-8);
                float orthoerr = cabsf(t.c01*conj(t.c00)+ t.c11*conj(t.c10)+ t.c21*conj(t.c20));
                assert(orthoerr < 1.0e-7);
            }
        }
        free_gfield_f(sombrero_in);
        free(maxelerO);
        printf("One-way orthonormality check - odd: Ok\n");

    }

    {// round-trip test
        suNf_field* sombrero_in = alloc_gfield_f(&glattice);
        suNf_field* sombrero_out = alloc_gfield_f(&glattice);
        su3* maxeler_intermediateE = allocate_maxeler_gauge_field();
        su3* maxeler_intermediateO = allocate_maxeler_gauge_field();


        random_uf(sombrero_in);

        sombrero_to_maxeler_gauge_field_E(sombrero_in,maxeler_intermediateE);
        sombrero_to_maxeler_gauge_field_O(sombrero_in,maxeler_intermediateO);

        maxeler_to_sombrero_gauge_field_E(maxeler_intermediateE,sombrero_out);
        maxeler_to_sombrero_gauge_field_O(maxeler_intermediateO,sombrero_out);

        {
            _MASTER_FOR(sombrero_in->type,ix){
                for(int mu=0;mu<4;++mu){
                    suNf* sombrero_in_tmp = sombrero_in->ptr + coord_to_index(ix,mu);
                    suNf* sombrero_out_tmp = sombrero_out->ptr + coord_to_index(ix,mu);
                    for(int c=0;c<9;++c){
                        {
                            float  in_re =  sombrero_in_tmp->c[c].re;
                            float out_re = sombrero_out_tmp->c[c].re;
                            float  error = in_re - out_re;

                            if(fabs(error)>1.5e-7){
                                printf("Error: ix=%d, c=%d, out_re: %.18lf , in_re: %.18lf diff:%1.5e\n",
                                       ix,c,out_re,in_re,error);
                                assert(0);
                            }
                        }
                        {
                            float  in_im =  sombrero_in_tmp->c[c].im;
                            float out_im = sombrero_out_tmp->c[c].im;
                            float  error = in_im - out_im;

                            if(fabs(error)>1.5e-7){
                                printf("Error: ix=%d, c=%d, out_im: %.18lf , in_im: %.18lf diff:%1.5e\n",
                                       ix,c,out_im,in_im,error);
                                assert(0);
                            }
                        }
                    }
                }
            }
        }

        free_gfield_f(sombrero_in);
        free_gfield_f(sombrero_out);
        free(maxeler_intermediateE);
        free(maxeler_intermediateO);
        printf("Round trip check: Ok\n");
    }

    finalize_process();
    return 0;

}




