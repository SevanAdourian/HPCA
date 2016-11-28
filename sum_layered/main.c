#include <stdio.h>
#include <stdlib.h>

struct edks_type {
  char prefix[128];
  char date[32];
  char version[32];
  char comment[128];
  unsigned int nlayer;
  int ndepth;
  int ndista;
  float *rho;
  float *alpha;
  float *beta;
  float *thickness;
  float depthmin;
  float depthmax;
  float distamin;
  float distamax;
  float *depths;
  float *distas;
  float ***zrtdsx;
};

extern void read_edks_(char *, struct edks_type *);
/* extern void read_receivers_(char *, unsigned int *, float *, float *); */
/* extern void read_patch_(char *, unsigned int *, float *, float *, float *, float *, float *, float *, float *, float *, float *); */
/* extern void patch2pts_(float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned int *, unsigned int *, float *, float *, float *, float *); */
/* extern void mom2disp_(float [][], float [][], float [][], float *, float *, float *, float *, float *, struct edks_type *, float *); */
/* extern void write_disp_(char *, unsigned int *, unsigned int *, float [][], float [][], float [][]); */

int main(int argc, char *argv[]) {
  /*float *uxs, *uys, *uzs;
    float *ux, *uy, *uz;*/
  float *xr, *yr;
  float *xsp, *ysp, *zsp;
  float *x, *y, *z;
  float *strike, *dip, *strike_r, *dip_r;
  float *rake, *W, *L, *slip;
  float dw, dy, M;
  float wg, yg, zg;
  unsigned int npw, npy, nrec, np, nspp, irec;
  char *elem_filename, prefix[256];
  struct edks_type edks;

  if (argc != 7) {
    printf("Usage: sum_layered edks_name geom_prefix #receivers #patches NPW NPY\n");
    exit(0);
  }

  elem_filename = malloc(20 * sizeof(char));
  elem_filename = argv[1];
  sscanf(argv[2], "%s", prefix);
  sscanf(argv[3], "%d", &nrec);
  sscanf(argv[4], "%d", &np);
  sscanf(argv[5], "%d", &npw);
  sscanf(argv[6], "%d", &npy);

  nspp = npw * npy; // number of  sources per patch

  float ux[nrec][np] ;
  float uy[nrec][np] ;
  float uz[nrec][np] ;
  float uxs[nrec][nspp] ;
  float uys[nrec][nspp] ;
  float uzs[nrec][nspp] ;

  strike = (float *) malloc(np * sizeof(float));
  dip = (float *) malloc(np * sizeof(float));
  strike_r = (float *) malloc(np * sizeof(float));
  dip_r = (float *) malloc(np * sizeof(float));
  rake = (float *) malloc(np * sizeof(float));
  W = (float *) malloc(np * sizeof(float));
  L = (float *) malloc(np * sizeof(float));
  slip = (float *) malloc(np * sizeof(float));
  x = (float *) malloc(np * sizeof(float));
  y = (float *) malloc(np * sizeof(float));
  z = (float *) malloc(np * sizeof(float));
  xsp = (float *) malloc(nspp * sizeof(float));
  ysp = (float *) malloc(nspp * sizeof(float));
  zsp = (float *) malloc(nspp * sizeof(float));
  /*ux = (float *) malloc(np * nrec * sizeof(float));
  uy = (float *) malloc(np * nrec * sizeof(float));
  uz = (float *) malloc(np * nrec * sizeof(float));
  uxs = (float *) malloc(nspp * nrec * sizeof(float));
  uys = (float *) malloc(nspp * nrec * sizeof(float));
  uzs = (float *) malloc(nspp * nrec * sizeof(float));*/
  xr = (float *) malloc(nrec * sizeof(float));
  yr = (float *) malloc(nrec * sizeof(float));

  read_edks_(elem_filename, &edks);

  /* read_receivers_(prefix, &nrec, xr, yr); */

  /* read_patch_(prefix, &np, x, y, z, strike, dip, rake, W, L, slip); */

  /* for (size_t ip = 0; ip < np; ip++) { */
  /*   patch2pts_(&strike[ip], &dip[ip], &rake[ip], &slip[ip], &W[ip], &L[ip], &x[ip], &y[ip], &z[ip], &npw, &npy, xsp, ysp, zsp, &M); */
  /*   mom2disp_(uxs, uys, uzs, xr, yr, xsp, ysp, zsp, &edks, &M); */

  /*   // need to sum pts back to patch contributions where arrays are */
  /*   // indexed as ux(patch,receiver), uxs(pt_source,receiver) */
  /*   for (size_t ii = 0; ii < nrec; ii++) { */
  /*     for (size_t jj = 0; jj < nspp; jj++) { */
  /* 	ux[ii][ip] = ux[ii][ip] + uxs[ii][jj]; */
  /* 	uy[ii][ip] = uy[ii][ip] + uys[ii][jj]; */
  /* 	uz[ii][ip] = uz[ii][ip] + uys[ii][jj]; */
  /*     } */
  /*   } */
  /* } */

  /* write_disp_(prefix, &np, &nrec, ux, uy, uz); */

  return 0;
}
