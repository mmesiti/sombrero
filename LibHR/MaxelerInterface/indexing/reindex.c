#include "su3_types.h"
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

void reindex(void* in, void* out,
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
        for(x = xmin; x < xmax; ++x)
            for(y = ymin; y < ymax; ++y)
                for(z = zmin; z < zmax; ++z)
                    for(t = tmin; t < tmax; ++t)
                    {
                        if(selection(x,y,z,t)){
                                int inidx =  inmap(x,y,z,t);
                                int outidx = outmap(x,y,z,t);
                                void  *inptr =  in+ inelsize*inidx;
                                void *outptr = out+outelsize*outidx;

                                memcpy( inel, inptr, inelsize);
                                transform(inel,outel);
                                memcpy(outptr,outel,outelsize);
                        }
                    }
    }

    free(outel);
    free(inel);
}

