#ifndef __FFT_FUNCTION_SET__
#define __FFT_FUNCTION_SET__

#include "Parameters.h"
#include "Parameters_input.h"
#include "LLG_kernel.cu"
#include <cmath>
//#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <curand.h>
#include <curand_kernel.h>


using namespace std;
#define   dx         delta_x
#define   dy         delta_y
#define   dz         delta_z


static bool G_matrix(int, int, int);

extern int G_tensor(int, int, int);

extern int Hms(int,        int,             int, 
               dim3, dim3);

#endif
