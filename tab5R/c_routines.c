#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "model.h"

#define DET_MIN 1e-9

int N, Np, M; /* A mettre dans le main ? */
int Id2[2][2];
/* Id2[0][0] = 1; */
/* Id2[0][1] = 1; */
/* Id2[1][0] = 1; */
/* Id2[1][1] = 1; */

int Ze2[2][2];
/* Ze2[0][0] = 0; */
/* Ze2[0][1] = 0; */
/* Ze2[1][0] = 0; */
/* Ze2[1][1] = 0; */

typedef struct {
  float mat[2][2];
} matrix2by2;

typedef struct {
  float mat[2];
} matrix2by1;

typedef struct {
  float *vp;
  float *vs;
  float *rh;
  float *th;
  float *lambda;
  float *mu;
  float *delta;
  float *gamma;
} in_calcqwE;

typedef struct {
  matrix2by2 q;
  matrix2by2 w;
  matrix2by2 v;
  matrix2by1 qE;
  matrix2by1 wE;
} out_calcqwE;

matrix2by2
inv22(matrix2by2 a) {

  /* Computes the invert of a 2x2 matrix. See if something smarter can be done. */
  float det;
  matrix2by2 ia;

  //ia = malloc(2*2 * sizeof(float));
  
  det = a.mat[0][0] * a.mat[1][1] - a.mat[2][1] * a.mat[1][2];
  if (det < DET_MIN) {
    printf("Error: inv22 called with a singular matrix\n");
    exit(EXIT_FAILURE);
  }
  
  ia.mat[0][0] = a.mat[1][1] / det;
  ia.mat[0][1] = -a.mat[1][0] / det;
  ia.mat[1][0] = -a.mat[0][1] / det;
  ia.mat[1][1] = a.mat[0][0] / det;

  return ia;
}
      
matrix2by2
matmul22(matrix2by2 a, matrix2by2 b) {
  /* Computes the product of 2 2x2 matrix. See if something smarter can be done. */
  /* Peut être à déclarer comme pointeur */
  matrix2by2 ab;

  //ab = malloc(2*2 * sizeof(float));

  /* For safety purpose, should check if they are 2x2. But heck it, user's responsability. */
  ab.mat[0][0] = a.mat[0][0] * b.mat[0][0] + a.mat[0][1] * b.mat[1][0];
  ab.mat[1][0] = a.mat[1][0] * b.mat[0][0] + a.mat[1][1] * b.mat[1][0];
  ab.mat[0][1] = a.mat[0][0] * b.mat[0][1] + a.mat[0][1] * b.mat[1][1];
  ab.mat[1][1] = a.mat[1][0] * b.mat[0][1] + a.mat[1][1] * b.mat[1][1];

  return ab;
}

matrix2by1
matmul21(matrix2by2 a, matrix2by1 b) {
  /* Computes the product of a 2x2 matrix by a 2x1 vector. See if something smarter can be done. */
  /* Peut être à déclarer comme pointeur */
  matrix2by1 ab;

  //ab = malloc(2*1 * sizeof(float));

  /* For safety purpose, should check if they are 2x2 and 2x1. But heck it, user's responsability. */
  ab.mat[0] = a.mat[0][0] * b.mat[0] + a.mat[0][1] * b.mat[1];
  ab.mat[1] = a.mat[1][0] * b.mat[0] + a.mat[1][1] * b.mat[1];
  
  return ab;
}

matrix2by2
id_mat2by2(matrix2by2 a) {
	matrix2by2 b;
	b.mat[0][0] = 1 - a.mat[0][0];
	b.mat[0][1] = - a.mat[0][1];
	b.mat[1][0] = - a.mat[1][0];
	b.mat[1][1] = 1 - a.mat[1][1];
	return b;
}

in_calcqwE
acotab(float depth, float htol, struct model_struct struct_model) {
  //int N, Np, M;
  float h, dh;
  in_calcqwE ret_val;

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

  if (flag == 1){
    ret_val.vp = malloc(N+1 * sizeof(float));
    ret_val.vs = malloc(N+1 * sizeof(float));
    ret_val.rh = malloc(N+1 * sizeof(float));
    ret_val.th = malloc(N * sizeof(float));
    memcpy(ret_val.vp, struct_model.vp0, N*sizeof(float));
    memcpy(ret_val.vs, struct_model.vs0, N*sizeof(float));
    memcpy(ret_val.rh, struct_model.rh0, N*sizeof(float));
    if (N > 1){
      memcpy(ret_val.th, struct_model.th0, (N-1)*sizeof(float));
    }
    if (M < N){
      if (M < N-1){
	/* Code de maniere naive */
	for (int ii = M + 1; ii < N - 1; ii++){
	  ret_val.th[ii] = ret_val.th[ii -1];
	}
      }
      ret_val.th[M] = ret_val.th[M-1] - h;
    }
    ret_val.th[M - 1] = h;
    
    for (int jj = M; jj < N; jj++){
      ret_val.vp[jj] = ret_val.vp [jj-1];
      ret_val.vs[jj] = ret_val.vs [jj-1];
      ret_val.rh[jj] = ret_val.rh [jj-1];
      N = N + 1;
    }
  }
  else {
    ret_val.vp = malloc(N * sizeof(float));
    ret_val.vs = malloc(N * sizeof(float));
    ret_val.rh = malloc(N * sizeof(float));
    ret_val.th = malloc((N - 1) * sizeof(float));
    memcpy(ret_val.vp, struct_model.vp0, N*sizeof(float));
    memcpy(ret_val.vs, struct_model.vs0, N*sizeof(float));
    memcpy(ret_val.rh, struct_model.rh0, N*sizeof(float));
    memcpy(ret_val.th, struct_model.th0, (N-1)*sizeof(float));
  }

  ret_val.lambda = malloc(N * sizeof(float));
  ret_val.mu = malloc(N * sizeof(float));
  ret_val.delta = malloc(N * sizeof(float));
  ret_val.gamma = malloc(N * sizeof(float));

  for (int ij = 0; ij < N - 1; ij++){
    ret_val.lambda[ij] = ret_val.rh[ij] * (powf(ret_val.vp[ij], 2)
					   - 2*powf(ret_val.vs[ij], 2));
    ret_val.mu[ij] = ret_val.rh[ij] * powf(ret_val.vs[ij], 2);
    ret_val.delta[ij] = 1 - (2 * powf(ret_val.vs[ij], 2)) /
      (powf(ret_val.vp[ij], 2) + powf(ret_val.vs[ij], 2));
    ret_val.gamma[ij] = (ret_val.lambda[ij] + ret_val.mu[ij]) /
      (ret_val.lambda[ij] + 2 * ret_val.mu[ij]);
  }
  Np = N - 1;
  
  return ret_val;
}

void
discotabE(int nsource, float *Pp, float *Pm, float *SVp, float *SVm,
	  float *SHp, float *SHm, float Ep, float Em, float *delta){

  float dd, eps, epd, tmp;

  dd = delta[nsource - 1];
  tmp = (-1 + 4 * dd)/2;

  eps = 1.;
  epd = eps * dd;

  Pp[0] = tmp / (1+dd); Pp[1] = -epd / (1+dd); Pp[2] = -0.5 / (1+dd);
  SVp[0] = -1.5 / (1+dd); SVp[1] = eps / (1+dd); SVp[2] = -0.5 / (1+dd);
  SHp[0] = 0; SHp[1] = eps; SHp[2] = -1;
  Ep = 1.;

  eps = -1.;
  epd = eps * dd;
  Pm[0] = tmp / (1+dd); Pm[1] = -epd / (1+dd); Pm[2] = 0.5 / (1+dd);
  SVm[0] = -1.5 / (1+dd); SVm[1] = eps / (1+dd); SVm[2] = -0.5 / (1+dd);
  SHm[0] = 0; SHm[1] = eps / (1+dd); SHm[2] = -1.0 / (1+dd);
  Ep = 1.;

}    


/* void */
/* calcqwvE(int nk, float *Pp, float *Pm, float *SVp, float *SVm, float *SHp, */
/* 	 float *SHm, float Ep, float Em, float *q, float *w, float *v, */
/* 	 float *qE, float *wE) { */
out_calcqwE
calcqwvE(int nk, float *Pp, float *Pm, float *SVp, float *SVm, float *SHp,
	 float *SHm, float Ep, float Em) {
  //  float *RUFS, *TURS, *RDRS, *TDRS;
  matrix2by2 RUFS, TURS, RDRS, TDRS;
  float RUFSsh, TURSsh, RDRSsh, TDRSsh;
  //  float *RUSL, *TUSL, *RDSL, *TDSL;
  matrix2by2 RUSL, TUSL, RDSL, TDSL;
  float RUSLsh, TUSLsh, RDSLsh, TDSLsh;
  //  float *REV, *RTIL;
  matrix2by2 REV, RTIL;
  float REVsh,  RTILsh, MOUTsh;
  matrix2by2 MOUT, TEMP1, TEMP2;
  //  float *MOUT, *TEMP1, *TEMP2;

  int jm, jk;
  matrix2by1 VD, qw, VDE, qwE;
  //  float  *VD, *qw;
  //  float *VDE, *qwE;
  float k, VDsh;

  matrix2by2 TUFS, RDFS, TDFS, RURS;
  //  float *TUFS, *RDFS, *TDFS, *RURS;
  float TUFSsh, RDFSsh, TDFSsh, RURSsh;

  out_calcqwE qwvE;
  /* Malloc à gogo. Creer un h file de toute urgence... */
  /* RUFS  = (float *)malloc(2*2 * sizeof(float)); */
  /* TURS  = (float *)malloc(2*2 * sizeof(float)); */
  /* RDRS  = (float *)malloc(2*2 * sizeof(float)); */
  /* TDRS  = (float *)malloc(2*2 * sizeof(float)); */

  /* RUSL  = (float *)malloc(2*2 * sizeof(float)); */
  /* TUSL  = (float *)malloc(2*2 * sizeof(float)); */
  /* RDSL  = (float *)malloc(2*2 * sizeof(float)); */
  /* TDSL  = (float *)malloc(2*2 * sizeof(float)); */

  /* REV   = (float *)malloc(2*2 * sizeof(float)); */
  /* RTIL  = (float *)malloc(2*2 * sizeof(float)); */

  /* MOUT  = (float *)malloc(2*2 * sizeof(float)); */
  /* TEMP1 = (float *)malloc(2*2 * sizeof(float)); */
  /* TEMP2 = (float *)malloc(2*2 * sizeof(float)); */

  /* VD    = (float *)malloc(2*1 * sizeof(float)); */
  /* qw    = (float *)malloc(2*1 * sizeof(float)); */
  /* VDE   = (float *)malloc(2*1 * sizeof(float)); */
  /* qwE   = (float *)malloc(2*1 * sizeof(float)); */

  /* TUFS  = (float *)malloc(2*2 * sizeof(float)); */
  /* RDFS  = (float *)malloc(2*2 * sizeof(float)); */
  /* TDFS  = (float *)malloc(2*2 * sizeof(float)); */
  /* RURS  = (float *)malloc(2*2 * sizeof(float)); */
	  
  int M1  = M+1;
  int Np1 = Np+1;
  
  for (jk = 0; jk < nk-1; jk++){
    //k = kv[jk];
    
    kenpsv_(0, &M1, &k, &RUFS, &TUFS, &RDFS, &TDFS);
    kensh_(0, &M1, &k, &RUFSsh, &TUFSsh, &RDFSsh, &TDFSsh);

    kenpsv_(1, &M1, &k, &RURS, &TURS, &RDRS, &TDRS);
    kensh_(1, &M1, &k, &RURSsh, &TURSsh, &RDRSsh, &TDRSsh);

    kenpsv_(&M1, &Np1, &k, &RUSL, &TUSL, &RDSL, &TDSL);
    kensh_(&M1, &Np1, &k, &RUSLsh, &TUSLsh, &RDSLsh, &TDSLsh);

    RERPSV_(1, &REV, &RTIL);
    REVsh = 2;
    RTILsh = 1;

    TEMP1 = matmul22(RDRS, RTIL);
    TEMP2 = matmul22(RDSL, RUFS);

    MOUT = matmul22(matmul22(matmul22(REV, inv22(id_mat2by2(TEMP1))), TURS),
		    inv22(id_mat2by2(TEMP2))); /* definir Id2, Ze2 */

    /* et hop, deux de moins */
    /* free(TEMP1); */
    /* free(TEMP2); */
    
    MOUTsh = REVsh / (1 - RDRSsh * RTILsh) * TURSsh / (1 - RDSLsh * RUFSsh);

    for (jm = 0; jm < 2; jm++) {
      matrix2by1 TEMP3;
      TEMP3.mat[0] = Pp[jm] + Pm[jm]; TEMP1.mat[1] = SVp[jm] + SVm[jm];
     
      VD = matmul21(RDSL, TEMP3);
      VDsh = RDSLsh*SHp[jm] + SHm[jm];

      qw = matmul21(MOUT, VD); /* VERIFIER DIMENSIONNNNNS */
      qwvE.q.mat[jm][jk] = qw.mat[0];
      qwvE.w.mat[jm][jk] = qw.mat[1];
      qwvE.v.mat[jm][jk] = MOUTsh * VDsh;
    }

    matrix2by1 Ep0;
    matrix2by1 Em0;
    Ep0.mat[0] = Ep; Ep0.mat[1] = 0;
    Em0.mat[0] = Em; Em0.mat[1] = 0; /* Sert a rien ?*/ 
    
    VDE = matmul21(RDSL, Ep0);
    qwE = matmul21(MOUT, VDE);
    qwvE.qE.mat[jk] = qwE.mat[0];
    qwvE.wE.mat[jk] = qwE.mat[1]; /* A verifier */

    /* free(RUFS); free(TURS); free(RDRS); free(TDRS); */
    /* free(RUSL); free(TUSL); free(RDSL); free(TDSL); */
    /* free(REV); free(RTIL); free(MOUT); */
    /* free(VD); free(qw); free(VDE); free(qwE); */
    /* free(TUFS);free(RDFS); free(TDFS); free(RURS); */

    /* free(vp); free(vs); free(rh); free(th); */
    /* free(lambda); free(mu); free(delta); free(gamma); */
    
    return qwvE;
  }
}

