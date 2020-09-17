#include "LibHR/MaxelerInterface/maxeler/max_cg_API.h"
#include "LibHR/MaxelerInterface/element_transforms/spinor_transforms.h"
#include "su3_types.h"
#include <assert.h>
#include <complex.h>

int main(){

    {

        suNg_spinor sombrero;
        cg_spinor maxeler;

        {
            int d,c;
            for(d=0;d<4;++d)
                for(c=0;c<3;++c){
                   sombrero.c[d].c[c].re = d*c+10;
                   sombrero.c[d].c[c].im = c-d;
                }
        }

        sombrero_to_maxeler_spinor(&sombrero,&maxeler);

        assert(maxeler.s0c0 == (0*0+10)+ I*(0-0)) ;
        assert(maxeler.s0c1 == (0*1+10)+ I*(1-0)) ;
        assert(maxeler.s0c2 == (0*2+10)+ I*(2-0)) ;

        assert(maxeler.s1c0 == (1*0+10)+ I*(0-1)) ;
        assert(maxeler.s1c1 == (1*1+10)+ I*(1-1)) ;
        assert(maxeler.s1c2 == (1*2+10)+ I*(2-1)) ;

        assert(maxeler.s2c0 == (2*0+10)+ I*(0-2)) ;
        assert(maxeler.s2c1 == (2*1+10)+ I*(1-2)) ;
        assert(maxeler.s2c2 == (2*2+10)+ I*(2-2)) ;
                                  
        assert(maxeler.s3c0 == (3*0+10)+ I*(0-3)) ;
        assert(maxeler.s3c1 == (3*1+10)+ I*(1-3)) ;
        assert(maxeler.s3c2 == (3*2+10)+ I*(2-3)) ;

    }


    {

        suNg_spinor sombrero;
        cg_spinor maxeler;

        maxeler.s0c0 = (0*0+10)+ I*(0-0) ;
        maxeler.s0c1 = (0*1+10)+ I*(1-0) ;
        maxeler.s0c2 = (0*2+10)+ I*(2-0) ;

        maxeler.s1c0 = (1*0+10)+ I*(0-1) ;
        maxeler.s1c1 = (1*1+10)+ I*(1-1) ;
        maxeler.s1c2 = (1*2+10)+ I*(2-1) ;

        maxeler.s2c0 = (2*0+10)+ I*(0-2) ;
        maxeler.s2c1 = (2*1+10)+ I*(1-2) ;
        maxeler.s2c2 = (2*2+10)+ I*(2-2) ;

        maxeler.s3c0 = (3*0+10)+ I*(0-3) ;
        maxeler.s3c1 = (3*1+10)+ I*(1-3) ;
        maxeler.s3c2 = (3*2+10)+ I*(2-3) ;

        maxeler_to_sombrero_spinor(&maxeler,&sombrero);

        {
            int d,c;
            for(d=0;d<4;++d)
                for(c=0;c<3;++c){
                   assert(sombrero.c[d].c[c].re == d*c+10 );
                   assert(sombrero.c[d].c[c].im == c-d );
                }
        }

    }



    return 0;

}
