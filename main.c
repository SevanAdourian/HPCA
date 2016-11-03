#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
  float *uxs, *uys, *uzs;
  float *ux, *uy, *uz;
  float *xr, *yr;
  float *xsp, *ysp, *zsp;
  float *x, *y, *z;
  float *strike, *dip, *strike_r, *dip_r;
  float *rake, *W, *L, *slip;
  float dw, dy, M;
  float wg, yg, zg;
  int npw, npy, nrec, np, nspp, irec, i, ip;
  char elem_filename[256], prefix[256];
  char argum[32];

  if (argc != 7)
  {
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
  ux = (float *) malloc(np * nrec * sizeof(float));
  uy = (float *) malloc(np * nrec * sizeof(float));
  uz = (float *) malloc(np * nrec * sizeof(float));
  uxs = (float *) malloc(nspp * nrec * sizeof(float));
  uys = (float *) malloc(nspp * nrec * sizeof(float));
  uzs = (float *) malloc(nspp * nrec * sizeof(float));
  xr = (float *) malloc(nrec * sizeof(float));
  yr = (float *) malloc(nrec * sizeof(float));

  /*! read in the edks (the greens functions)
  call read_edks(elem_filename, edks)

  ! read in the receiver locations
  call read_receivers(prefix, nrec, xr, yr)

  ! read in the patch information
  call read_patch(prefix, np, x, y, z, strike, dip, rake, W, L, slip)

  ! need to loop over all patches, get local coordinates and fill single arrays
  do ip = 1, np

     call patch2pts(strike(ip), dip(ip), rake(ip), slip(ip), W(ip), L(ip), &
                    x(ip), y(ip), z(ip), npw, npy, xsp, ysp, zsp, M)

     call mom2disp(uxs, uys, uzs, xr, yr, xsp, ysp, zsp, edks, M)

     ! need to sum pts back to patch contributions where arrays are
     ! indexed as ux(patch,receiver), uxs(pt_source,receiver)
     ux(ip, :) = sum(uxs, dim=1)
     uy(ip, :) = sum(uys, dim=1)
     uz(ip, :) = sum(uzs, dim=1)

  end do

  ! write out displacement arrays

  call write_disp(prefix, np, nrec, ux, uy, uz)

  end program sum_layered*/

  return 0;
}
