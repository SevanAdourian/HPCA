#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "model.h"

int N, Np, M; /* A mettre dans le main ? */

float *
acotab(float depth, float htol, struct model_struct struct_model) {
  //int N, Np, M;
  float h, dh;
  float *vp, *vs, *rh, *th, *lambda, *mu, *delta, *gamma;

  h = depth;
  N = sizeof(struct_model.vp0) / sizeof(float);

  if (N < 1){
    printf("[acotab] Error: Nb of layers <= 0 - %s\n", strerror(errno));
    exit(1);
  }
  
  if (h < 0){
    printf("[acotab] Error: Source depth < 0 - %s\n", strerror(errno));
    exit(1);
  }

  if (htol < 0){
    printf("[acotab] Warning: h-tolerance < 0\n using 100m instead - %s\n",
	   strerror(errno));
    htol = 100;
  }

  int flag = 0;
  if (h < htol){
    M = 0;
    flag = -1;
  }
  else if (N == 1){
    M = 1;
    flag = 1;
  }
  else {
    for(M = 1; M < N-1; M++){
      dh = h - struct_model.th0[M];
      if (dh > htol){
	h = dh;
      }
      else if (dh > -htol){
	flag = -1; /* The source is at one pre-existing interface */
	break;
      }
      else{
	flag = 1;
	break;
      }
    }
    if (flag == 0){
      /* the source is in the half space */
      M = N;
      flag = 1;
    }
  }

  if (flag = 1){
    vp = malloc(N+1 * sizeof(float));
    vs = malloc(N+1 * sizeof(float));
    rh = malloc(N+1 * sizeof(float));
    th = malloc(N * sizeof(float));
    memcpy(vp, struct_model.vp0, N*sizeof(float));
    memcpy(vs, struct_model.vs0, N*sizeof(float));
    memcpy(rh, struct_model.rh0, N*sizeof(float));
    if (N > 1){
      memcpy(th, struct_model.th0, (N-1)*sizeof(float));
    }
    if (M < N){
      if (M < N-1){
	/* Code de maniere naive */
	for (int ii = M + 1; ii < N - 1; ii++){
	  th[ii] = th[ii -1];
	}
      }
      th[M] = th[M-1] - h;
    }
    th[M - 1] = h;
    
    for (int jj = M; jj < N; jj++){
      vp[jj] = vp [jj-1];
      vs[jj] = vs [jj-1];
      rh[jj] = rh [jj-1];
      N = N + 1;
    }
  }
  else {
    vp = malloc(N * sizeof(float));
    vs = malloc(N * sizeof(float));
    rh = malloc(N * sizeof(float));
    th = malloc((N - 1) * sizeof(float));
    memcpy(vp, struct_model.vp0, N*sizeof(float));
    memcpy(vs, struct_model.vs0, N*sizeof(float));
    memcpy(rh, struct_model.rh0, N*sizeof(float));
    memcpy(th, struct_model.th0, (N-1)*sizeof(float));
  }

  lambda = malloc(N * sizeof(float));
  mu = malloc(N * sizeof(float));
  delta = malloc(N * sizeof(float));
  gamma = malloc(N * sizeof(float));

  for (int ij = 0; ij < N - 1; ij++){
    lambda[ij] = rh[ij] * (powf(vp[ij], 2) - 2*powf(vs[ij], 2));
    mu[ij] = rh[ij] * powf(vs[ij], 2);
    delta[ij] = 1 - (2 * powf(vs[ij], 2)) /
      (powf(vp[ij], 2) + powf(vs[ij], 2));
    gamma[ij] = (lambda[ij] + mu[ij]) / (lambda[ij] + 2 * mu[ij]);
  }
  Np = N - 1;
  
  return *vp, *vs, *rh, *th, *lambda, *mu, *delta, *gamma;
}
  
float *
discotabE(int nsource, float Pp, float Pm, float SVp, float SVm,
	  float SHp, float SHm, float Ep, float Em){

  float dd, eps, epd, tmp;

  dd = delta[nsource - 1];
  tmp = (-1 + 4 * dd)/2;

  eps = 1;
  epd = eps * dd;

  Pp[0] = tmp / (1+dd); Pp[1] = -epd / (1+dd); Pp[2] = -0.5 / (1+dd);
  SVp[0] = -1.5 / (1+dd); SVp[1] = eps / (1+dd); SVp[2] = -0.5 / (1+dd);
  SHp[0] = 0; SHp[1] = eps; SHp[2] = -1;
  Ep = 1;

  eps = -1;
  epd = eps * dd;
  Pm[0] = tmp / (1+dd); Pm[1] = -epd / (1+dd); Pm[2] = 0.5 / (1+dd);
  SVm[0] = -1.5 / (1+dd); SVm[1] = eps / (1+dd); SVm[2] = -0.5 / (1+dd);
  SHm[0] = 0; SHm[1] = eps / (1+dd); SHm[2] = -1.0 / (1+dd);
  Ep = 1;

}    

float *
calcqwvE(int nk, float *Pp, float *Pm, float *SVp, float *SVm, float *SHp,
	 float *SHm, float Ep, float Em, float *q, float *w, float *v,
	 float *qE, float *wE) {
  float *RUFS, *TURS, *RDRS, *TDRS;
  float RUFSsh, TURSsh, RDRSsh, TDRSsh;
  float *RUSL, *TUSL, *RDSL, *TDSL;
  float RUSLsh, TUSLsh, RDSLsh, TDSLsh;
  float *REV, *RTIL;
  float REVsh,  RTILsh, MOUTsh;
  float *MOUT, *TEMP1, *TEMP2;
  int jm, jk;
  float  *VD, *qw;
  float *VDE, *qwE;
  float k, VDsh;

  float *TUFS, *RDFS, *TDFS, *RURS;
  float *TUFSsh, *RDFSsh, *TDFSsh, *RURSsh;


  /* Malloc à gogo. Creer un h file de toute urgence... */
  RUFS  = (float *)malloc(2*2 * sizeof(float));
  TURS  = (float *)malloc(2*2 * sizeof(float));
  RDRS  = (float *)malloc(2*2 * sizeof(float));
  TDRS  = (float *)malloc(2*2 * sizeof(float));

  RUSL  = (float *)malloc(2*2 * sizeof(float));
  TUSL  = (float *)malloc(2*2 * sizeof(float));
  RDSL  = (float *)malloc(2*2 * sizeof(float));
  TDSL  = (float *)malloc(2*2 * sizeof(float));

  REV   = (float *)malloc(2*2 * sizeof(float));
  RTIL  = (float *)malloc(2*2 * sizeof(float));

  MOUT  = (float *)malloc(2*2 * sizeof(float));
  TEMP1 = (float *)malloc(2*2 * sizeof(float));
  TEMP2 = (float *)malloc(2*2 * sizeof(float));

  VD    = (float *)malloc(2*1 * sizeof(float));
  qw    = (float *)malloc(2*1 * sizeof(float));
  VDE   = (float *)malloc(2*1 * sizeof(float));
  qwE   = (float *)malloc(2*1 * sizeof(float));

  TUFS  = (float *)malloc(2*2 * sizeof(float));
  RDFS  = (float *)malloc(2*2 * sizeof(float));
  TDFS  = (float *)malloc(2*2 * sizeof(float));
  RURS  = (float *)malloc(2*2 * sizeof(float));
	  
  int M1  = M+1;
  int NP1 = NP+1;
  
  for (int jk = 0; jk < nk-1; jk++){
    k = kv[jk];
    
    kenpsv_(0, &M1, k, RUFS, TUFS, RDFS, TDFS);
    kensh_(0, &M1, k, RUFSsh, TUFSsh, RDFSsh, TDFSsh);

    kenpsv_(1, &M1, k, RURS, TURS, RDRS, TDRS);
    kensh_(1, &M1, k, RURSsh, TURSsh, RDRSsh, TDRSsh);

    kenpsv_(&M1, &Np1, k, RUSL, TUSL, RDSL, TDSL);
    kensh_(&M1, &Np1, k, RUSLsh, TUSLsh, RDSLsh, TDSLsh);

    RERPSV_(1, REV, RTIL);
    REVsh = 2;
    RTILsh = 1;

    TEMP1 = matmul22(RDRS, RTIL);
    TEMP2 = matmul22(RDSL, RUFS);

    MOUT = matmul22(matmul22(matmul22(REV, inv22(Id2 - TEMP1)), TURS),
		    inv22(Id2 - TEMP)); /* definir Id2, Ze2 */

    /* et hop, deux de moins */
    free(TEMP1);
    free(TEMP2);
    
    MOUTsh = REVsh / (1 - RDRSsh * RTILsh) * TURSsh / (1 - RDSLsh * RUFSsh);

    for (int jm = 0; jm < 2; jm++) {
      float TEMP1[2];
      TEMP1[0] = Pp[jm] + Pm[jm]; TEMP1[1] = SVp[jm] + SVm[jm];
     
      VD = matmul21(RDSL, TEMP1);
      VDsh = RDSLsh*SHp[jm] + SHm[jm];

      qw = matmul22(MOUT, VD);
      q[jm][jk] = qw[0];
      w[jm][jk] = qw[1];
      v[jm][jk] = MOUTsh * VDsh;
    }

    float Ep0[2] = {0};
    float Em0[2] = {0};
    Ep0[0] = Ep; Ep[1] = 0;
    Em0[0] = Em; Em[1] = 0;
    
    VDE = matmul21(RDSL, Ep0);
    qwE = matmul22(MOUT, VDE);
    qE[jk] = qwE[0];
    wE[jk] = qwE[0];

    free(RUFS); free(TURS); free(RDRS); free(TDRS);
    free(RUSL); free(TUSL); free(RDSL); free(TDSL);
    free(REV); free(RTIL); free(MOUT);
    free(VD); free(qw); free(VDE); free(qwE);
    free(TUFS);free(RDFS); free(TDFS); free(RURS);
    
    return *q, *w, *v, *qE, *wE;
    }
  }
}

float *
inv22(float a) {

  /* Computes the invert of a 2x2 matrix. See if something smarter can be done. */
  float det;
  float ia[2][2];
  
  det = a[0][0] * a[1][1] - a[2][1] * a[1][2];
  if (det < DET_MIN) {
    printf("Error: inv22 called with a singular matrix\n");
    exit(EXIT_ERROR);
  }
  
  ia[0][0] = a[1][1] / det;
  ia[0][1] = -a[1][0] / det;
  ia[1][0] = -a[0][1] / det;
  ia[1][1] = a[0][0] / det;

  return ia;
}
      
float *
matmul22(a, b) {
  /* Computes the product of 2 2x2 matrix. See if something smarter can be done. */
  /* Peut être à déclarer comme pointeur */
  float *ab[2][2];

  /* For safety purpose, should check if they are 2x2. But heck it, user's responsability. */
  ab[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
  ab[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
  ab[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
  ab[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];

  return *ab;
}

float *
matmul21(a, b) {
  /* Computes the product of a 2x2 matrix by a 2x1 vector. See if something smarter can be done. */
  /* Peut être à déclarer comme pointeur */
  float *ab[2];

  /* For safety purpose, should check if they are 2x2 and 2x1. But heck it, user's responsability. */
  ab[0] = a[0][0] * b[0] + a[0][1] * b[1];
  ab[1] = a[1][0] * b[0] + a[1][1] * b[1];
  
  return *ab;
}