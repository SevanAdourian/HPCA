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

matrix2by2 inv22(matrix2by2);      
matrix2by2 matmul22(matrix2by2, matrix2by2);
matrix2by1 matmul21(matrix2by2, matrix2by1);
matrix2by2 id_mat2by2(matrix2by2);
in_calcqwE acotab(float, float, struct model_struct);
void discotabE(int, float *, float *, float *, float *, float *, float *, float, float, float *);
out_calcqwE calcqwvE(int, float *, float *, float *, float *, float *, float *, float, float);
