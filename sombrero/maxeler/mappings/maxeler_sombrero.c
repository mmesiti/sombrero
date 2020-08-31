#include "maxeler_sombrero.h"
#include "su3_types.h"

void _remap(void* in, void* out,
            size_t inelsize,
            size_t outelsize,
            int xmin, int xmax,
            int ymin, int ymax,
            int zmin, int zmax,
            int tmin, int tmax,
            int (*inmap)(int,int,int,int),
            int (*outmap)(int,int,int,int),
            void (*transform)(void*,void*),
            int (*selection)(int,int,int,int)) // e.g., even-odd
{

    void *inel, *outel;
    inel = malloc(inelsize);
    outel = malloc(outelsize);

    {
        int x,y,z,t;
        for(t = tmin; t < tmax; ++t)
            for(z = zmin; z < zmax; ++z)
                for(y = ymin; y < ymax; ++y)
                    for(x = xmin; x < xmax; ++x){
                        if selection(x,y,z,t){
                                void  *inptr =  in+ inelsize* inmap(x,y,z,t);
                                void *outptr = out+outelsize*outmap(x,y,z,t);

                                memcp( inel, inptr, inelsize);
                                transform(outel,inel);
                                memcp(outptr,outel,outelsize);
                            }
                    }
    }

    free(outel);
    free(inel);
}

// sombrero-maxeler element transformations

void sombrero_to_maxeler_spinor(void* _sombrero_spinor,
                                void* _maxeler_spinor){

    suNf_spinor* sombrero_spinor = _sombrero_spinor;
    cg_spinor*    maxeler_spinor =  _maxeler_spinor;

    maxeler_spinor.s0c0 = sombrero_spinor.c[0].c[0];
    maxeler_spinor.s0c1 = sombrero_spinor.c[0].c[1];
    maxeler_spinor.s0c2 = sombrero_spinor.c[0].c[2];

    maxeler_spinor.s1c0 = sombrero_spinor.c[1].c[0];
    maxeler_spinor.s1c1 = sombrero_spinor.c[1].c[1];
    maxeler_spinor.s1c2 = sombrero_spinor.c[1].c[2];

    maxeler_spinor.s2c0 = sombrero_spinor.c[2].c[0];
    maxeler_spinor.s2c1 = sombrero_spinor.c[2].c[1];
    maxeler_spinor.s2c2 = sombrero_spinor.c[2].c[2];

    maxeler_spinor.s3c0 = sombrero_spinor.c[3].c[0];
    maxeler_spinor.s3c1 = sombrero_spinor.c[3].c[1];
    maxeler_spinor.s3c2 = sombrero_spinor.c[3].c[2];

}

void maxeler_to_sombrero_spinor(void* _maxeler_spinor,
                                void* _sombrero_spinor){

    cg_spinor*    maxeler_spinor =  _maxeler_spinor;
    suNf_spinor* sombrero_spinor = _sombrero_spinor;

    sombrero_spinor.c[0].c[0] = maxeler_spinor.s0c0;
    sombrero_spinor.c[0].c[1] = maxeler_spinor.s0c1;
    sombrero_spinor.c[0].c[2] = maxeler_spinor.s0c2;

    sombrero_spinor.c[1].c[0] = maxeler_spinor.s1c0;
    sombrero_spinor.c[1].c[1] = maxeler_spinor.s1c1;
    sombrero_spinor.c[1].c[2] = maxeler_spinor.s1c2;

    sombrero_spinor.c[2].c[0] = maxeler_spinor.s2c0;
    sombrero_spinor.c[2].c[1] = maxeler_spinor.s2c1;
    sombrero_spinor.c[2].c[2] = maxeler_spinor.s2c2;

    sombrero_spinor.c[3].c[0] = maxeler_spinor.s3c0;
    sombrero_spinor.c[3].c[1] = maxeler_spinor.s3c1;
    sombrero_spinor.c[3].c[2] = maxeler_spinor.s3c2;

}

void sombrero_to_maxeler_gauge(void* _sombrero_gauge,
                               void*  _maxeler_gauge){
    su3*   maxeler_gauge =  _maxeler_gauge;
    suNg* sombrero_gauge = _sombrero_gauge;

    maxeler_gauge.c00 = sombrero_gauge.c[0];
    maxeler_gauge.c01 = sombrero_gauge.c[1];

    maxeler_gauge.c10 = sombrero_gauge.c[3];
    maxeler_gauge.c11 = sombrero_gauge.c[4];

    maxeler_gauge.c20 = sombrero_gauge.c[6];
    maxeler_gauge.c21 = sombrero_gauge.c[7];

}

void sombrero_to_maxeler_gauge4(void* _sombrero_gauge4,
                                void*  _maxeler_gauge4){
    int mu;
    const int X = 0, Y = 1, Z = 2, T = 3;
    const int SOMBRERO_ORDERING[4] = {T,X,Y,Z};
    const int  MAXELER_ORDERING[4] = {X,Y,Z,T};

    for(mu=0;mu<4;++mu)
        sombrero_to_maxeler_gauge(
            _sombrero_gauge4 + SOMBRERO_ORDERING[i]*sizeof(suNg),
             _maxeler_gauge4 +  MAXELER_ORDERING[i]*sizeof( su3)
        );
}
