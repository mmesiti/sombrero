#include "LibHR/MaxelerInterface/maxeler/max_cg_API.h"
#include "LibHR/MaxelerInterface/element_transforms/gauge_transforms.h"
#include "su3_types.h"
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>

int main(){

    {// sombrero -> maxeler

        suNg sombrero;
        su3 maxeler;

        {// sombrero matrix initialisation
            int cr,cc;
            for(cr=0;cr<3;++cr)
                for(cc=0;cc<3;++cc){
                    sombrero.c[cr*3+cc].re = cr*cc+10;
                    sombrero.c[cr*3+cc].im = cr-cc;
                }
        }

        sombrero_to_maxeler_gauge(&sombrero,&maxeler);

        {// checks on maxeler matrix
            assert(maxeler.c00 == (0*0+10) + I*(0-0));
            assert(maxeler.c01 == (0*1+10) + I*(0-1));

            assert(maxeler.c10 == (1*0+10) + I*(1-0));
            assert(maxeler.c11 == (1*1+10) + I*(1-1));

            assert(maxeler.c20 == (2*0+10) + I*(2-0));
            assert(maxeler.c21 == (2*1+10) + I*(2-1));
        }
    }

    { // sombrero -> maxeler -> sombrero
        su3 maxeler;
        suNg sombrero1, sombrero2;

        { // maxeler matrix initialisation
          // columns must be orthonormal
#define _cmplx_assign(c,RE,IM)                  \
            (c).re = RE;                        \
            (c).im = IM

            // V1/2   iV1/3     V1/6
            // V1/2  -iV1/3    -V1/6
            //    0    V1/3  +2iV1/6

            _cmplx_assign(sombrero1.c[0],     sqrtf(0.5),              0);
            _cmplx_assign(sombrero1.c[1],              0, sqrtf(1.0/3.0));
            _cmplx_assign(sombrero1.c[2], sqrtf(1.0/6.0),              0);
            _cmplx_assign(sombrero1.c[3],     sqrtf(0.5),              0);
            _cmplx_assign(sombrero1.c[4],              0,-sqrtf(1.0/3.0));
            _cmplx_assign(sombrero1.c[5],-sqrtf(1.0/6.0),              0);
            _cmplx_assign(sombrero1.c[6],              0,              0);
            _cmplx_assign(sombrero1.c[7], sqrtf(1.0/3.0),              0);
            _cmplx_assign(sombrero1.c[8],              0, sqrtf(2.0/3.0));

#undef _cmplx_assign
        }

        sombrero_to_maxeler_gauge(&sombrero1,&  maxeler);
        maxeler_to_sombrero_gauge(&  maxeler,&sombrero2);

        {// checks

#define _check_diff(i,_GETMACRO)                            \
            {                                               \
                float diff = _GETMACRO(sombrero2.c[i]) -    \
                             _GETMACRO(sombrero1.c[i]);     \
                if(fabs(diff)>1.0e-7){                      \
                    printf("Error: %.18lf vs %.18lf "       \
                           "diff: %.18lf , i: %d\n",        \
                           _GETMACRO(sombrero2.c[i]),       \
                           _GETMACRO(sombrero1.c[i]),       \
                           fabs(diff),                      \
                           i);                              \
                    assert( _GETMACRO(sombrero2.c[i])       \
                            ==                              \
                            _GETMACRO(sombrero1.c[i]));     \
                }                                           \
            }

            { // checking real part
#define _getre(a) (a).re
#define _check_re(i) _check_diff(i,_getre)
                _check_re(0);
                _check_re(1);
                _check_re(2);
                _check_re(3);
                _check_re(5);
                _check_re(6);
                _check_re(7);
                _check_re(8);
#undef _getre
#undef _check_re
            }
            { // checking imaginary part
#define _getim(a) (a).im
#define _check_im(i) _check_diff(i,_getim)
                _check_im(0);
                _check_im(1);
                _check_im(2);
                _check_im(3);
                _check_im(5);
                _check_im(6);
                _check_im(7);
                _check_im(8);
#undef _getim
#undef _check_im
            }
#undef _check_diff
        }

    }

    return 0;

}
