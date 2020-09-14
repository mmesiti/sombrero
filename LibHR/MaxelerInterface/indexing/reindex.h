#ifndef MAXELER_REINDEX_H
#define MAXELER_REINDEX_H

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
             int (*selection)(int,int,int,int)); // e.g., even-odd
#endif
