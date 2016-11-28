#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#define MODE_MAX 100

struct model_struct{
  float rh0[MODE_MAX];
  float vp0[MODE_MAX];
  float vs0[MODE_MAX];
  float th0[MODE_MAX];
};

struct param_struct{
  float rh0[MODE_MAX];
  float vp0[MODE_MAX];
  float vs0[MODE_MAX];
  float th0[MODE_MAX];
}; // le flou absolu
