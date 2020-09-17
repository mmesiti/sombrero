#include "LibHR/MaxelerInterface/maxeler/max_cg_API.h"
#include "LibHR/MaxelerInterface/element_transforms/gauge_transforms.h"
#include "su3_types.h"
#include <assert.h>
#include <complex.h>

int main(){

    {

        suNg sombrero;
        su3 maxeler;

        {
            int cr,cc;
            for(cr=0;cr<3;++cr)
                for(cc=0;cc<3;++cc){
                   sombrero.c[cr*3+cc].re = cr*cc+10;
                   sombrero.c[cr*3+cc].im = cr-cc;
                }
        }

        sombrero_to_maxeler_gauge(&sombrero,&maxeler);

        {
            assert(maxeler.c00 == (0*0+10) + I*(0-0));
            assert(maxeler.c01 == (0*1+10) + I*(0-1));

            assert(maxeler.c10 == (1*0+10) + I*(1-0));
            assert(maxeler.c11 == (1*1+10) + I*(1-1));

            assert(maxeler.c20 == (2*0+10) + I*(2-0));
            assert(maxeler.c21 == (2*1+10) + I*(2-1));
        }
    }

    return 0;

}
