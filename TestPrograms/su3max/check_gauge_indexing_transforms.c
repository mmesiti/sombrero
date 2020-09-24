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
#include <stdlib.h>
#include <math.h>

static float random_sign(){
    return (rand()*2 > RAND_MAX)?1:-1;
}


static void random_vector(complex* c){
    for(int i=0;i<3;++i){
        c[i].re = (1.0e-5+(float)rand()/RAND_MAX)*random_sign();
        c[i].im = 0;//(1.0e-5+(float)rand()/RAND_MAX)*random_sign(); //DEBUG
    }
}

static _Complex float hrtoc99(complex c){
    return c.re + I*c.im;
}

static complex c99tohr(_Complex float c){
    complex res = {.re = crealf(c), .im = cimagf(c)};
    return res;
}

static complex scalar_product(const complex* c1, const complex *c2){
    _Complex float res = 0 + I*0;
    for(int i=0;i<3;++i){
        res += conj( hrtoc99(c1[i]) )*
                     hrtoc99(c2[i]);
    }
    return c99tohr(res);
}

static void normalize_vector(complex* c){

    float z = sqrtf((scalar_product(c,c)).re);
    for(int i=0;i<3;++i){
        c[i].re = c[i].re / z;
        c[i].im = c[i].im / z;
    }

    { // check
        float err = scalar_product(c,c).re - 1;
        assert(fabs(err) < 3e-7);
    }
}

static void remove_component(complex* c1, const complex* normalized_c2){

    complex p = scalar_product(normalized_c2,c1);
    for(int i=0;i<3;++i){
        c1[i].re -= (p.re*normalized_c2[i].re - p.im*normalized_c2[i].im);
        c1[i].im -= (p.re*normalized_c2[i].im + p.im*normalized_c2[i].re);
    }

    {// check
        complex pr = scalar_product(c1,normalized_c2);
        assert(pr.im*pr.im < 1.0e-7);
        assert(pr.re*pr.re < 1.0e-7);
    }

}

static void random_u3(complex* c){
    for(int irow=0;irow<3;++irow){
       complex* v = &c[irow*3];
       random_vector(v);
       normalize_vector(v);
       for(int jrow=0;jrow<irow;++jrow){
           complex* w = &c[jrow*3];
           remove_component(v,w);
           normalize_vector(v);
       }
    }
}


static _Complex float det_u3(complex* c){
    _Complex float cc99[3][3];
    for(int row=0;row<3;++row)
        for(int col=0;col<3;++col)
            cc99[row][col] = hrtoc99(c[row*3+col]);

    _Complex float det = 0;

    // plus
    for(int colstart=0;colstart<3;++colstart){
        _Complex float term = 1;
        for(int row=0; row<3; ++row){
            int col = (row+colstart)%3;
            term *= cc99[row][col];

        }
        det += term;
    }

    // minus
    for(int colstart=0;colstart<3;++colstart){
        _Complex float term = 1;
        for(int row=0; row<3; ++row){
            int col_signed = colstart-row;
            int col = col_signed<0?col_signed+3:col_signed;
            term *= cc99[row][col];
        }
        det -= term;
    }
    return det;

}


static void random_su3(complex* c){
    random_u3(c);
    _Complex float det = det_u3(c);
    {// multiply last row so that the determinant is
        for(int i = 6; i<9;++i){
            _Complex float tmp = c[i].re + I* c[i].im;
            tmp /= det;
            c[i].re = creal(tmp);
            c[i].im = cimag(tmp);
        }
    }
}



static void glb_setup(){
    NP_T=NP_X=NP_Y=NP_Z=1;
    GLB_T=GLB_X=GLB_Y=GLB_Z=8;
    T=X=Y=Z=8;
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

    { // pre-checks
        {
            complex c1[]={{.re=0.5,.im=0},{.re=-.5,.im=1.0},{.re=0,.im=-0.5}};
            complex c2[]={{.re=-0.5,.im=0},{.re=-.5,.im=0.5},{.re=1.0,.im=-1.0}};
            complex res = scalar_product(c1,c2);
            float re_err = res.re - 1;
            float im_err = res.im - (0 + .25 + 0.5);
            assert(fabs(im_err) < 1.0e-7);
            assert(fabs(re_err) < 1.0e-7);
        }

        {
            complex c[9];
            random_su3(c);
            assert(cabsf(det_u3(c)-1)<1.0e-6);
        }
        printf("Pre checks Ok\n");
    }

    { // one-way determinant check
        { // even
            suNf_field* sombrero_in = alloc_gfield_f(&glattice);
            su3* maxelerE = allocate_maxeler_gauge_field();

            {// sombrero_in initialisation

                _MASTER_FOR(sombrero_in->type,ix){
                    for(int mu=0;mu<4;++mu){
                        suNf* sombrero_in_tmp = sombrero_in->ptr + coord_to_index(ix,mu);
                        random_su3(sombrero_in_tmp->c);
                    }
                }
            }
            printf("Initialisation Ok\n");

            sombrero_to_maxeler_gauge_field_E(sombrero_in,maxelerE);

            printf("Translation Ok\n");

            { // check that all the matrices are ok
                for(int i=0;i<T*X*Y*Z/2;++i){
                    su3 t  = maxelerE[i];
                    float norm0err = crealf(t.c00*conj(t.c00)+ t.c10*conj(t.c10)+ t.c20*conj(t.c20))-1;
                    assert(norm0err*norm0err < 1.0e-12);
                    float norm1err = crealf(t.c01*conj(t.c01)+ t.c11*conj(t.c11)+ t.c21*conj(t.c21))-1;
                    assert(norm1err*norm1err < 1.0e-12);
                    float orthoerr = cabsf(t.c01*conj(t.c00)+ t.c11*conj(t.c10)+ t.c21*conj(t.c20));
                    assert(orthoerr < 1.0e-12);
                }
            }
        }
        printf("One-way determinant check Ok\n");

    }



    {// round-trip test
        suNf_field* sombrero_in = alloc_gfield_f(&glattice);
        suNf_field* sombrero_out = alloc_gfield_f(&glattice);
        su3* maxeler_intermediateE = allocate_maxeler_gauge_field();
        su3* maxeler_intermediateO = allocate_maxeler_gauge_field();


        {// sombrero_in initialisation

            _MASTER_FOR(sombrero_in->type,ix){
                for(int mu=0;mu<4;++mu){
                    suNf* sombrero_in_tmp = sombrero_in->ptr + coord_to_index(ix,mu);
                    random_su3(sombrero_in_tmp->c);
                }
            }
        }

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
                        double  in_re =  sombrero_in_tmp->c[c].re;
                        double out_re = sombrero_out_tmp->c[c].re;

                        if(in_re != out_re){
                            printf("Error: ix=%d, d=%d, c=%d, out_re: %.18lf , in_re: %.18lf\n",
                                   ix,c,out_re,in_re);
                            assert(0);
                        }

                        double  in_im =  sombrero_in_tmp->c[c].im;
                        double out_im = sombrero_out_tmp->c[c].im;

                        if(in_im != out_im){
                            printf("Error: ix=%d, d=%d, c=%d, out_im: %.18lf , in_im: %.18lf\n",
                                   ix,c,out_im,in_im);
                            assert(0);
                        }
                    }
                }
            }
        }
    }
}




