#include "su3_types.h"
#include "../maxeler/max_cg_API.h"


void sombrero_to_maxeler_spinor(void* _sombrero_spinor,
                                void* _maxeler_spinor){

    suNf_spinor* sombrero_spinor = _sombrero_spinor;
    cg_spinor*    maxeler_spinor =  _maxeler_spinor;

#define _stdcmplx(x) ((x).re + I * (x).im)

    maxeler_spinor->s0c0 = _stdcmplx(sombrero_spinor->c[0].c[0]);
    maxeler_spinor->s0c1 = _stdcmplx(sombrero_spinor->c[0].c[1]);
    maxeler_spinor->s0c2 = _stdcmplx(sombrero_spinor->c[0].c[2]);

    maxeler_spinor->s1c0 = _stdcmplx(sombrero_spinor->c[1].c[0]);
    maxeler_spinor->s1c1 = _stdcmplx(sombrero_spinor->c[1].c[1]);
    maxeler_spinor->s1c2 = _stdcmplx(sombrero_spinor->c[1].c[2]);

    maxeler_spinor->s2c0 = _stdcmplx(sombrero_spinor->c[2].c[0]);
    maxeler_spinor->s2c1 = _stdcmplx(sombrero_spinor->c[2].c[1]);
    maxeler_spinor->s2c2 = _stdcmplx(sombrero_spinor->c[2].c[2]);

    maxeler_spinor->s3c0 = _stdcmplx(sombrero_spinor->c[3].c[0]);
    maxeler_spinor->s3c1 = _stdcmplx(sombrero_spinor->c[3].c[1]);
    maxeler_spinor->s3c2 = _stdcmplx(sombrero_spinor->c[3].c[2]);

#undef _stdcmplx

}

void maxeler_to_sombrero_spinor(void* _maxeler_spinor,
                                void* _sombrero_spinor){

    cg_spinor*    maxeler_spinor =  _maxeler_spinor;
    suNf_spinor* sombrero_spinor = _sombrero_spinor;

    sombrero_spinor->c[0].c[0].re = creal(maxeler_spinor->s0c0);
    sombrero_spinor->c[0].c[0].im = cimag(maxeler_spinor->s0c0);
    sombrero_spinor->c[0].c[1].re = creal(maxeler_spinor->s0c1);
    sombrero_spinor->c[0].c[1].im = cimag(maxeler_spinor->s0c1);
    sombrero_spinor->c[0].c[2].re = creal(maxeler_spinor->s0c2);
    sombrero_spinor->c[0].c[2].im = cimag(maxeler_spinor->s0c2);

    sombrero_spinor->c[1].c[0].re = creal(maxeler_spinor->s1c0);
    sombrero_spinor->c[1].c[0].im = cimag(maxeler_spinor->s1c0);
    sombrero_spinor->c[1].c[1].re = creal(maxeler_spinor->s1c1);
    sombrero_spinor->c[1].c[1].im = cimag(maxeler_spinor->s1c1);
    sombrero_spinor->c[1].c[2].re = creal(maxeler_spinor->s1c2);
    sombrero_spinor->c[1].c[2].im = cimag(maxeler_spinor->s1c2);

    sombrero_spinor->c[2].c[0].re = creal(maxeler_spinor->s2c0);
    sombrero_spinor->c[2].c[0].im = cimag(maxeler_spinor->s2c0);
    sombrero_spinor->c[2].c[1].re = creal(maxeler_spinor->s2c1);
    sombrero_spinor->c[2].c[1].im = cimag(maxeler_spinor->s2c1);
    sombrero_spinor->c[2].c[2].re = creal(maxeler_spinor->s2c2);
    sombrero_spinor->c[2].c[2].im = cimag(maxeler_spinor->s2c2);

    sombrero_spinor->c[3].c[0].re = creal(maxeler_spinor->s3c0);
    sombrero_spinor->c[3].c[0].im = cimag(maxeler_spinor->s3c0);
    sombrero_spinor->c[3].c[1].re = creal(maxeler_spinor->s3c1);
    sombrero_spinor->c[3].c[1].im = cimag(maxeler_spinor->s3c1);
    sombrero_spinor->c[3].c[2].re = creal(maxeler_spinor->s3c2);
    sombrero_spinor->c[3].c[2].im = cimag(maxeler_spinor->s3c2);



}
