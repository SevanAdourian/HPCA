#include <stdio.h>
#include <stdlib.h>

float *allocate_mem(float*** arr, int n, int m) {
  *arr = malloc(n * sizeof **arr);
  float *arr_data = malloc(n * m * sizeof *arr_data);
  for (size_t i = 0; i < n; ++i) {
    (*arr)[i] = arr_data + i * m ;
  }
  return arr_data;
}
void deallocate_mem(float*** arr, float* arr_data) {
  free(arr_data);
  free(*arr);
}

struct edks_type {
  char prefix[128], date[32], version[32], comment[128];
  int nlayer, ndepth, ndista;
  float depthmin, depthmax, distamin, distamax;
  float *rho, *alpha, *beta, *thickness, *depths, *distas;
  float ***zrtdsx;
};

extern void read_edks_(char *, struct edks_type *);
extern void read_receivers_(char *, unsigned int *, float *, float *);
extern void read_patch_(char *, unsigned int *, float *, float *, float *, float *, float *, float *, float *, float *, float *);
extern void patch2pts_(float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned int *, unsigned int *, float *, float *, float *, float *);
extern void mom2disp_(float **, float **, float **, float *, float *, float *, float *, float *, struct edks_type *, float *);
extern void write_disp_(char *, unsigned int *, unsigned int *, float **, float **, float **);

int main(int argc, char const *argv[]) {
  float **uxs, **uys, **uzs;
  float **ux, **uy, **uz;
  float *xr, *yr;
  float *xsp, *ysp, *zsp;
  float *x, *y, *z;
  float *strike, *dip, *strike_r, *dip_r;
  float *rake, *W, *L, *slip;
  float dw, dy, M;
  float wg, yg, zg;
  unsigned int npw, npy, nrec, np, nspp, irec;
  char elem_filename[256], prefix[256];
  struct edks_type edks;

  if (argc != 7) {
    printf("Usage: sum_layered edks_name geom_prefix #receivers #patches NPW NPY\n");
    exit(0);
  }

  sscanf(argv[1], "%s", elem_filename);
  sscanf(argv[2], "%s", prefix);
  sscanf(argv[3], "%d", &nrec);
  sscanf(argv[4], "%d", &np);
  sscanf(argv[5], "%d", &npw);
  sscanf(argv[6], "%d", &npy);

  nspp = npw * npy; // number of  sources per patch

  strike = malloc(np * sizeof *strike);
  dip = malloc(np * sizeof *dip);
  strike_r = malloc(np * sizeof *strike_r);
  dip_r = malloc(np * sizeof *dip_r);
  rake = malloc(np * sizeof *rake);
  W = malloc(np * sizeof *W);
  L = malloc(np * sizeof *L);
  slip = malloc(np * sizeof *slip);
  x = malloc(np * sizeof *x);
  y = malloc(np * sizeof *y);
  z = malloc(np * sizeof *z);
  xsp = malloc(nspp * sizeof *xsp);
  ysp = malloc(nspp * sizeof *ysp);
  zsp = malloc(nspp * sizeof *zsp);
  float *ux_m = allocate_mem(&ux, nrec, np);
  float *uy_m = allocate_mem(&uy, nrec, np);
  float *uz_m = allocate_mem(&uz, nrec, np);
  float *uxs_m = allocate_mem(&uxs, nrec, nspp);
  float *uys_m = allocate_mem(&uys, nrec, nspp);
  float *uzs_m = allocate_mem(&uzs, nrec, nspp);
  xr = malloc(nrec * sizeof *xr);
  yr = malloc(nrec * sizeof *yr);

  read_edks_(elem_filename, &edks);

  /*read_receivers_(prefix, &nrec, xr, yr);

  read_patch_(prefix, &np, x, y, z, strike, dip, rake, W, L, slip);

  for (size_t ip = 0; ip < np; ip++) {
  patch2pts_(strike + ip, dip + ip, rake + ip, slip + ip, W + ip, L + ip, x + ip, y + ip, z + ip, &npw, &npy, xsp, ysp, zsp, &M);
  mom2disp_(uxs, uys, uzs, xr, yr, xsp, ysp, zsp, &edks, &M);

  // need to sum pts back to patch contributions where arrays are
  // indexed as ux(patch,receiver), uxs(pt_source,receiver)
  for (size_t ii = 0; ii < nrec; ii++) {
  for (size_t jj = 0; jj < nspp; jj++) {
  ux[ii][ip] = ux[ii][ip] + uxs[ii][jj];
  uy[ii][ip] = uy[ii][ip] + uys[ii][jj];
  uz[ii][ip] = uz[ii][ip] + uys[ii][jj];
}
}
}

write_disp_(prefix, &np, &nrec, ux, uy, uz);*/

deallocate_mem(&ux, ux_m);
deallocate_mem(&uy, uy_m);
deallocate_mem(&uz, uz_m);
deallocate_mem(&uxs, uxs_m);
deallocate_mem(&uys, uys_m);
deallocate_mem(&uzs, uzs_m);

return 0;
}
