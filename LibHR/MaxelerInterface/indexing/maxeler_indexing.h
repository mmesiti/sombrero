#ifndef MAXELER_MAXELER_INDEXING_H
#define MAXELER_MAXELER_INDEXING_H
int maxeler_spinor_idx(int x, int y, int z, int t);

int maxeler_gaugeE_idx(int x, int y, int z, int t);
int maxeler_gaugeO_idx(int x, int y, int z, int t);

int maxeler_even(int x, int y, int z, int t);
int maxeler_odd(int x, int y, int z, int t);
int maxeler_any(int x, int y, int z, int t);

#endif
