#define c_4pi 12.5663706144



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
  float *q;
  float *w;
  float *v;
  float *qE;
  float *wE;
} out_calcqwE;

matrix2by2 inv22(matrix2by2);      

matrix2by2 matmul22(matrix2by2, matrix2by2);

matrix2by1 matmul21(matrix2by2, matrix2by1);

matrix2by2 id_mat2by2(matrix2by2);

in_calcqwE acotab(float, float, struct model_struct);

void discotabE(int, float *, float *, float *, float *, float *, float *, float, float, float *);

out_calcqwE calcqwvE(int, float *, float *, float *, float *, float *, float *, float, float);

void calc_BCD0E_(float *, float *, float *, float *, float *, float *, float *, float *, float, float, float *);

void calc_coeffsE_(float *, float *, float *, float *, float *, float, float *, float *, float *, float, float, float, float *);

