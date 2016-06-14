// LLG_CUDA.cpp : Defines the entry point for the console application.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <cassert>
#include <ctime>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <vector>
#include <algorithm>

using namespace std;
#include "Parameters.h"
#include "Parameters_input.h"
#include "LLG_kernel.cu"
#include "FFT_Function_Set.h"
#include "Moving_Head.h"
#include "Ha_field.h"
#include "GrainsIndicator.h"
#include "MiscFunctions.h"
#include "Outputs.h"
#include "utilities.h"

double  *Mx = NULL, *My = NULL, *Mz = NULL; 
double  *Mx_bar1 = NULL, *My_bar1 = NULL, *Mz_bar1 = NULL, 
        *Mx_bar2 = NULL, *My_bar2 = NULL, *Mz_bar2 = NULL,
		*Mx_bar3 = NULL, *My_bar3 = NULL, *Mz_bar3 = NULL,
		*Mx_bar4 = NULL, *My_bar4 = NULL, *Mz_bar4 = NULL,
		*Mx_bar5 = NULL, *My_bar5 = NULL, *Mz_bar5 = NULL,
		*Mx_bar6 = NULL, *My_bar6 = NULL, *Mz_bar6 = NULL,
		*Mx_bar = NULL,  *My_bar = NULL,  *Mz_bar = NULL,
		*Mz_CGC = NULL;
int     *indicator1 = NULL, *indicator2 = NULL, *indicator3 = NULL, *indicator7 = NULL;

// Variables defined for recording simulation
double *memTemp0 = NULL, *memTemp1 = NULL, *memTemp2 = NULL, *memTemp3 = NULL,
       *Happl_x_temp = NULL, *Happl_y_temp = NULL, *Happl_z_temp = NULL,
	   *Happl_x = NULL, *Happl_y = NULL, *Happl_z = NULL,
	   *FP_theta = NULL, *FP_mag = NULL;
double *theta0 = NULL;
double *mGrid[3];
double *fGrid[3];
double *FP_inp[3], *FP_trail[3], *FP[3];
double *Mag0[3];
double *fSeq = NULL;
int    *idx_f = NULL, *idx_t = NULL;
int    *Seq_f = NULL;

// Variables defined for cluster size calculation
double *rrt_cor = NULL;
double **Mz_tt;


// FFT
double *Hd_x_1d = NULL, *Hd_y_1d = NULL, *Hd_z_1d = NULL;

double *Hd_x_1d_shift = NULL, *Hd_y_1d_shift = NULL, *Hd_z_1d_shift = NULL;

complex<float> *Hd_x_1d_cmplx = NULL, *Hd_y_1d_cmplx = NULL, *Hd_z_1d_cmplx = NULL;

complex<float> *Mx_1d_cmplx = NULL,   *My_1d_cmplx = NULL,   *Mz_1d_cmplx = NULL;

complex<float> *Gxx_1d_cmplx = NULL,  *Gxy_1d_cmplx = NULL,  *Gxz_1d_cmplx = NULL, 
               *Gyx_1d_cmplx = NULL,  *Gyy_1d_cmplx = NULL,  *Gyz_1d_cmplx = NULL, 
               *Gzx_1d_cmplx = NULL,  *Gzy_1d_cmplx = NULL,  *Gzz_1d_cmplx = NULL;

cufftComplex *dev_Gxx_cufft = NULL,  *dev_Gxy_cufft = NULL,  *dev_Gxz_cufft = NULL,
             *dev_Gyx_cufft = NULL,  *dev_Gyy_cufft = NULL,  *dev_Gyz_cufft = NULL,
             *dev_Gzx_cufft = NULL,  *dev_Gzy_cufft = NULL,  *dev_Gzz_cufft = NULL;

cufftComplex *dev_Mx_cufft = NULL,   *dev_My_cufft = NULL,   *dev_Mz_cufft = NULL;

cufftComplex *dev_Hd_x_cufft = NULL, *dev_Hd_y_cufft = NULL, *dev_Hd_z_cufft = NULL;







bool AllocateDeviceMemory(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad, int mNx, int mNy, int mNz, int DEG_FREEDOM)
{
	/* Allocate memory on device */
	cudaMalloc((void **)&dev_GasArray, DEG_FREEDOM * sizeof(double));             
	if (dev_GasArray == NULL) { printf("cudaMalloc() fail.\n");  return false; }          
	cudaMalloc((void **)&dev_theta, mNx * mNy * mNz * sizeof(double));			
	if (dev_theta == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_phi, mNx * mNy * mNz * sizeof(double));				
	if (dev_phi == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_a_theta, mNx * mNy * mNz * sizeof(double));			
	if (dev_a_theta == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_a_phi, mNx * mNy * mNz * sizeof(double));					
	if (dev_a_phi == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_b_theta, mNx * mNy * mNz * sizeof(double));			
	if (dev_b_theta == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_b_phi, mNx * mNy * mNz * sizeof(double));				
	if (dev_b_phi == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_c_theta, mNx * mNy * mNz * sizeof(double));				
	if (dev_c_theta == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_c_phi, mNx * mNy * mNz * sizeof(double));					
	if (dev_c_phi == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_d_theta, mNx * mNy * mNz * sizeof(double));			
	if (dev_d_theta == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_d_phi, mNx * mNy * mNz * sizeof(double));				
	if (dev_d_phi == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hth_x, mNx * mNy * mNz * sizeof(double));				
	if (dev_Hth_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hth_y, mNx * mNy * mNz * sizeof(double));			
	if (dev_Hth_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hth_z, mNx * mNy * mNz * sizeof(double));				
	if (dev_Hth_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }	
	cudaMalloc((void **)&dev_Hk_x, mNx * mNy * mNz * sizeof(double));				
	if (dev_Hk_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hk_y, mNx * mNy * mNz * sizeof(double));				
	if (dev_Hk_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hk_z, mNx * mNy * mNz * sizeof(double));				
	if (dev_Hk_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }


	cudaMalloc((void **)&dev_Aex_XP, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Aex_XP == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Aex_XM, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Aex_XM == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Aex_YP, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Aex_YP == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Aex_YM, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Aex_YM == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Aex_ZP, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Aex_ZP == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Aex_ZM, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Aex_ZM == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_XP, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Ms_XP == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_XM, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Ms_XM == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_YP, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Ms_YP == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_YM, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Ms_YM == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_ZP, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Ms_ZP == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_ZM, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_Ms_ZM == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	
	cudaMalloc((void **)&dev_H1_x, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H1_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H1_y, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H1_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H1_z, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H1_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H2_x, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H2_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H2_y, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H2_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H2_z, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H2_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H3_x, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H3_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H3_y, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H3_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H3_z, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H3_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H4_x, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H4_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H4_y, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H4_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_H4_z, mNx * mNy * mNz * sizeof(double));						  // 18 zyliu
	if (dev_H4_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }






	cudaMalloc((void **)&dev_Ha_x, mNx * mNy * mNz * sizeof(double));						  // 18
	if (dev_Ha_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ha_y, mNx * mNy * mNz * sizeof(double));						  // 19
	if (dev_Ha_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ha_z, mNx * mNy * mNz * sizeof(double));						  // 20	
	if (dev_Ha_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hal_x, mNx * mNy * mNz * sizeof(double));						  // 18
	if (dev_Hal_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hal_y, mNx * mNy * mNz * sizeof(double));						  // 19
	if (dev_Hal_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hal_z, mNx * mNy * mNz * sizeof(double));						  // 20	
	if (dev_Hal_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hd_x, mNx * mNy * mNz * sizeof(double));						  // 21
	if (dev_Hd_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hd_y, mNx * mNy * mNz * sizeof(double));						  // 22
	if (dev_Hd_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hd_z, mNx * mNy * mNz * sizeof(double));						  // 23
	if (dev_Hd_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hint_x, mNx * mNy * mNz * sizeof(double));						  // 21
	if (dev_Hint_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hint_y, mNx * mNy * mNz * sizeof(double));						  // 22
	if (dev_Hint_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hint_z, mNx * mNy * mNz * sizeof(double));						  // 23
	if (dev_Hint_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }


	// Head Field and Temperature Profile.//
	cudaMalloc((void **)&dev_Happl_x, mNx * mNy * mNz * sizeof(double));					  // 24
	if (dev_Happl_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Happl_y, mNx * mNy * mNz * sizeof(double));					  // 25
	if (dev_Happl_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Happl_z, mNx * mNy * mNz * sizeof(double));					  // 26
	if (dev_Happl_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_D, mNx * mNy * mNz * sizeof(double));							  // 27
	if (dev_D == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_T, mNx * mNy * mNz * sizeof(double));							  // 27
	if (dev_T == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	//-----------------------------------//


	cudaMalloc((void **)&dev_d_theta_d_t, mNx * mNy * mNz * sizeof(double));				  // 28
	if (dev_d_theta_d_t == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_d_phi_d_t, mNx * mNy * mNz * sizeof(double));					  // 29
	if (dev_d_phi_d_t == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Mx, mNx * mNy * mNz * sizeof(double));						  // 30
	if (dev_Mx == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_My, mNx * mNy * mNz * sizeof(double));						  // 31
	if (dev_My == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Mz, mNx * mNy * mNz * sizeof(double));						  // 32
	if (dev_Mz == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	
	// Indicator
	cudaMalloc((void **)&dev_indicator1, mNx * mNy * mNz * sizeof(int));					  // 33
	if (dev_indicator1 == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_indicator1_temp, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(int));	  // 34
	if (dev_indicator1_temp == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_indicator2, mNx * mNy * mNz * sizeof(int));					  // 35
	if (dev_indicator2 == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_indicator2_temp, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(int));	  // 36
	if (dev_indicator2_temp == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_indicator3, mNx * mNy * mNz * sizeof(int));					  // 33
	if (dev_indicator3 == NULL) { printf("cudaMalloc() fail.\n");  return false; }

	// M_temp
	cudaMalloc((void **)&dev_M_temp_x, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));		  // 37
	if (dev_M_temp_x == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_M_temp_y, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));		  // 38
	if (dev_M_temp_y == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_M_temp_z, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));		  // 39
	if (dev_M_temp_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	
	// Watch
	cudaMalloc((void **)&dev_watch2, (mNx) * (mNy) * (mNz) * sizeof(double));				  // 40
	if (dev_watch2 == NULL) { printf("cudaMalloc() fail.\n");  return false; }

	// Input magnetic parameters
	cudaMalloc((void **)&dev_Aex, (mNx) * (mNy) * (mNz) * sizeof(double));				      // 41
	if (dev_Aex == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Aex_temp, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));			  // 42
	if (dev_Aex_temp == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms, (mNx) * (mNy) * (mNz) * sizeof(double));				      // 43
	if (dev_Ms == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ms_temp, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));				  // 44
	if (dev_Ms_temp == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ku, (mNx) * (mNy) * (mNz) * sizeof(double));					  // 45
	if (dev_Ku == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Ku_temp, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));				  // 46
	if (dev_Ku_temp == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_alpha, (mNx) * (mNy) * (mNz) * sizeof(double));				  // 47
	if (dev_alpha == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_alpha_temp, (mNx+2) * (mNy+2) * (mNz+2) * sizeof(double));			  // 48
	if (dev_alpha_temp == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_gamma, (mNx) * (mNy) * (mNz) * sizeof(double));			  // 48
	if (dev_gamma == NULL) { printf("cudaMalloc() fail.\n");  return false; }

	// CUFFT
	cudaMalloc((void **)&dev_Gxx_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gxx_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gxy_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gxy_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gxz_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gxz_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gyx_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gyx_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gyy_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gyy_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gyz_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gyz_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gzx_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gzx_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gzy_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gzy_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Gzz_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Gzz_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Mx_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Mx_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_My_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_My_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Mz_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Mz_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hd_x_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Hd_x_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hd_y_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Hd_y_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }
	cudaMalloc((void **)&dev_Hd_z_cufft, (lx_zero_pad) * (ly_zero_pad) * (lz_zero_pad) * sizeof(cufftComplex));			  // 48
	if (dev_Hd_z_cufft == NULL) { printf("cudaMalloc() fail.\n");  return false; }

	//recording simulation
	cudaMalloc((void **)&dev_Mag0_z, (mNx) * (mNy) * (mNz) * sizeof(double));			  // 48
	if (dev_Mag0_z == NULL) { printf("cudaMalloc() fail.\n");  return false; }

	return true;
}

//---------- Allocate Memory (including real and imaginary parts) ----------//
bool AllocateHostMemory(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad, int mNx, int mNy, int mNz, int DEG_FREEDOM, int TOTAL_TIME)
{

	Hd_x_1d = (double*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad, sizeof(double));
	Hd_y_1d = (double*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad, sizeof(double));
	Hd_z_1d = (double*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad, sizeof(double));

	Hd_x_1d_cmplx = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Hd_y_1d_cmplx = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Hd_z_1d_cmplx = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Mx_1d_cmplx   = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	My_1d_cmplx   = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Mz_1d_cmplx   = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gxx_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gxy_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gxz_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gyx_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gyy_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gyz_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gzx_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gzy_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	Gzz_1d_cmplx  = (complex<float>*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad,  sizeof(complex<float>));
	
	Hd_x_1d_shift =  (double*) calloc(mNx * mNy * mNz,  sizeof(double));
	Hd_y_1d_shift =  (double*) calloc(mNx * mNy * mNz,  sizeof(double));
	Hd_z_1d_shift =  (double*) calloc(mNx * mNy * mNz,  sizeof(double));

	host_Aex   = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	host_Ku    =  (double*) calloc(mNx * mNy * mNz, sizeof(double));
	host_Ms    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	host_alpha = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	host_gamma = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	indicator1 = (int*) calloc(mNx * mNy * mNz, sizeof(int));
	indicator2 = (int*) calloc(mNx * mNy * mNz, sizeof(int));
	indicator3 = (int*) calloc(mNx * mNy * mNz, sizeof(int));
	indicator4 = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	indicator5 = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	indicator6 = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	indicator7 = (int*) calloc(mNx * mNy * mNz, sizeof(int));
	T          = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	D          = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Happl_x    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Happl_y    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Happl_z    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Hint_x    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Hint_y    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Hint_z    = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Mx         = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	My         = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Mz         = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	watch100   = (double*) calloc((mNx+2) * (mNy+2) * (mNz+2), sizeof(double));
	watch2     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	watch3     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	watch4     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	watch5     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	watch100_int = (int*) calloc((mNx+2) * (mNy+2) * (mNz+2), sizeof(int));
	watch2_int = (int*) calloc(mNx * mNy * mNz, sizeof(int));
	watch3_int = (int*) calloc(mNx * mNy * mNz, sizeof(int));
	Mx_float   = (float*) calloc(mNx * mNy * mNz, sizeof(float));
	My_float   = (float*) calloc(mNx * mNy * mNz, sizeof(float));
	Mz_float   = (float*) calloc(mNx * mNy * mNz, sizeof(float));
	Mx_t_bar   = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	My_t_bar   = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	Mz_t_bar   = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	M_bar      = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar     = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar     = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar     = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar1    = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar1    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar1    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar2    = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar2    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar2    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar3    = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar3    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar3    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar4    = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar4    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar4    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar5    = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar5    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar5    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mx_bar6    = (double*) calloc(TOTAL_TIME, sizeof(double));
	My_bar6    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_bar6    = (double*) calloc(TOTAL_TIME, sizeof(double));
	Mz_CGC     = (double*) calloc(TOTAL_TIME, sizeof(double));
	
	Mz_LayerLayer = (double*) calloc(mNz, sizeof(double));
	
	Hint12_x_SingleGrain = (double*) calloc(mNx*mNy, sizeof(double));
	Hint12_y_SingleGrain = (double*) calloc(mNx*mNy, sizeof(double));
	Hint12_z_SingleGrain = (double*) calloc(mNx*mNy, sizeof(double));
	std_Ku     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	std_Aex     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	std_Tc     = (double*) calloc(mNx * mNy * mNz, sizeof(double));
	gasarray   = (double*) calloc(DEG_FREEDOM, sizeof(double));

	rrt_cor = (double*) calloc(mNy, sizeof(double));
	Mz_tt = (double**)calloc(mNx*mNy, sizeof(double));
	for (int ij=0; ij<mNx*mNy; ij++) { Mz_tt[ij] = (double*)calloc(TOTAL_TIME, sizeof(double)); }


	// Allocate Host memory for recording simulation
	theta0 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	mGrid[0] = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	mGrid[1] = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	mGrid[2] = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	fGrid[0] = (double*) calloc(fNx*fNy*fNz, sizeof(double));
	fGrid[1] = (double*) calloc(fNx*fNy*fNz, sizeof(double));
	fGrid[2] = (double*) calloc(fNx*fNy*fNz, sizeof(double));
	Mag0[0] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //Pre-pattern Mx
	Mag0[1] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //Pre-pattern My
	Mag0[2] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //Pre-pattern Mz
	FP_inp[0] = (double*) calloc(fNx*fNy*fNz, sizeof(double)); //field profile directly from Input file
	FP_inp[1] = (double*) calloc(fNx*fNy*fNz, sizeof(double)); //field profile directly from Input file
	FP_inp[2] = (double*) calloc(fNx*fNy*fNz, sizeof(double)); //field profile directly from Input file
	FP_theta   = (double*) calloc(fNx*fNy*fNz, sizeof(double));
	FP_mag     = (double*) calloc(fNx*fNy*fNz, sizeof(double));
	FP[0] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //field profile Hx
	FP[1] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //field profile Hy
	FP[2] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //field profile Hz
	FP_trail[0] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //field profile Hx
	FP_trail[1] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //field profile Hy
	FP_trail[2] = (double*) calloc(mNx*mNy*mNz, sizeof(double)); //field profile Hz
	fSeq = (double*) calloc(sfNc, sizeof(double));
	idx_f = (int*) calloc(1, sizeof(int));
	idx_t = (int*) calloc(1, sizeof(int));
	Seq_f = (int*) calloc(1, sizeof(int));
	Happl_x_temp = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	Happl_y_temp = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	Happl_z_temp = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	memTemp0 = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	memTemp1 = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	memTemp2 = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	memTemp3 = (double*)calloc(mNx*mNy*mNz, sizeof(double));
	printf("Host temp memory allocation succeeded! \n");
	


	return true;
}
	//--------------------------------------------------------------------------//

void ReleaseMemory(void)
{
	// Release memory on device
	if (dev_GasArray)   cudaFree(dev_GasArray);
	if (dev_theta)      cudaFree(dev_theta);				 
	if (dev_phi)        cudaFree(dev_phi);					 
	if (dev_a_theta)    cudaFree(dev_a_theta);				
	if (dev_a_phi)      cudaFree(dev_a_phi);				
	if (dev_b_theta)    cudaFree(dev_a_theta);				
	if (dev_b_phi)      cudaFree(dev_a_phi);				
	if (dev_c_theta)    cudaFree(dev_a_theta);				 
	if (dev_c_phi)      cudaFree(dev_a_phi);				
	if (dev_d_theta)    cudaFree(dev_a_theta);				
	if (dev_d_phi)      cudaFree(dev_a_phi);				
	if (dev_Hth_x)      cudaFree(dev_Hth_x);			    
	if (dev_Hth_y)      cudaFree(dev_Hth_y);				
	if (dev_Hth_z)      cudaFree(dev_Hth_z);	
	if (dev_Hk_x)		cudaFree(dev_Hk_x);                 
	if (dev_Hk_y)		cudaFree(dev_Hk_y);					
	if (dev_Hk_z)       cudaFree(dev_Hk_z);				       
	
	if (dev_Aex_XP)       cudaFree(dev_Aex_XP);   // zyliu
	if (dev_Aex_XM)       cudaFree(dev_Aex_XM);   // zyliu
	if (dev_Aex_YP)       cudaFree(dev_Aex_YP);   // zyliu
	if (dev_Aex_YM)       cudaFree(dev_Aex_YM);   // zyliu
	if (dev_Aex_ZP)       cudaFree(dev_Aex_ZP);   // zyliu
	if (dev_Aex_ZM)       cudaFree(dev_Aex_ZM);   // zyliu
	
	if (dev_Ms_XP)       cudaFree(dev_Ms_XP);     // zyliu
	if (dev_Ms_XM)       cudaFree(dev_Ms_XM);     // zyliu
	if (dev_Ms_YP)       cudaFree(dev_Ms_YP);     // zyliu
	if (dev_Ms_YM)       cudaFree(dev_Ms_YM);     // zyliu
	if (dev_Ms_ZP)       cudaFree(dev_Ms_ZP);     // zyliu
	if (dev_Ms_ZM)       cudaFree(dev_Ms_ZM);     // zyliu

	if (dev_H1_x)       cudaFree(dev_H1_x);     // zyliu
	if (dev_H1_y)       cudaFree(dev_H1_y);     // zyliu
	if (dev_H1_z)       cudaFree(dev_H1_z);     // zyliu
	if (dev_H2_x)       cudaFree(dev_H2_x);     // zyliu
	if (dev_H2_y)       cudaFree(dev_H2_y);     // zyliu
	if (dev_H2_z)       cudaFree(dev_H2_z);     // zyliu
	if (dev_H3_x)       cudaFree(dev_H3_x);     // zyliu
	if (dev_H3_y)       cudaFree(dev_H3_y);     // zyliu
	if (dev_H3_z)       cudaFree(dev_H3_z);     // zyliu
	if (dev_H4_x)       cudaFree(dev_H4_x);     // zyliu
	if (dev_H4_y)       cudaFree(dev_H4_y);     // zyliu
	if (dev_H4_z)       cudaFree(dev_H4_z);     // zyliu
	
	if (dev_Ha_x)       cudaFree(dev_Ha_x);				
	if (dev_Ha_y)       cudaFree(dev_Ha_y);					 
	if (dev_Ha_z)       cudaFree(dev_Ha_z);
	if (dev_Hal_x)       cudaFree(dev_Hal_x);				
	if (dev_Hal_y)       cudaFree(dev_Hal_y);					 
	if (dev_Hal_z)       cudaFree(dev_Hal_z);	
	if (dev_Hd_x)       cudaFree(dev_Hd_x);					
	if (dev_Hd_y)       cudaFree(dev_Hd_y);					
	if (dev_Hd_z)       cudaFree(dev_Hd_z);	
	if (dev_Hint_x)     cudaFree(dev_Hint_x);					
	if (dev_Hint_y)     cudaFree(dev_Hint_y);					
	if (dev_Hint_z)     cudaFree(dev_Hint_z);	
	if (dev_Happl_x)    cudaFree(dev_Happl_x);				
	if (dev_Happl_y)    cudaFree(dev_Happl_y);				
	if (dev_Happl_z)    cudaFree(dev_Happl_z);			
	if (dev_D)          cudaFree(dev_D);
	if (dev_T)          cudaFree(dev_T);
	if (dev_d_theta_d_t)cudaFree(dev_d_theta_d_t);			
	if (dev_d_phi_d_t)  cudaFree(dev_d_phi_d_t);	
	if (dev_Mx)         cudaFree(dev_Mx);				
	if (dev_My)         cudaFree(dev_My);			
	if (dev_Mz)         cudaFree(dev_Mz);	

	if (dev_indicator1)		  cudaFree(dev_indicator1);		 // 33
	if (dev_indicator1_temp)  cudaFree(dev_indicator1_temp); // 34
	if (dev_indicator2)       cudaFree(dev_indicator2);		 // 35
	if (dev_indicator2_temp)  cudaFree(dev_indicator2_temp); // 36
	if (dev_indicator3)       cudaFree(dev_indicator3);

	if (dev_M_temp_x)   cudaFree(dev_M_temp_x);				 // 37
	if (dev_M_temp_y)   cudaFree(dev_M_temp_y);				 // 38
	if (dev_M_temp_z)   cudaFree(dev_M_temp_z);				 // 39

	if (dev_watch2)     cudaFree(dev_watch2);				 // 40

	if (dev_Aex)         cudaFree(dev_Aex);					 // 41
	if (dev_Aex_temp)    cudaFree(dev_Aex_temp);		     // 42
	if (dev_Ms)          cudaFree(dev_Ms);				     // 43
	if (dev_Ms_temp)     cudaFree(dev_Ms_temp);				 // 44
	if (dev_Ku)          cudaFree(dev_Ku);				     // 45
	if (dev_Ku_temp)     cudaFree(dev_Ku_temp);				 // 46
	if (dev_alpha)       cudaFree(dev_alpha);		     	 // 47
	if (dev_alpha_temp)  cudaFree(dev_alpha_temp);		     // 48
	if (dev_gamma)       cudaFree(dev_gamma);
	if (dev_Gxx_cufft)   cudaFree(dev_Gxx_cufft);
	if (dev_Gxy_cufft)   cudaFree(dev_Gxy_cufft);
	if (dev_Gxz_cufft)   cudaFree(dev_Gxz_cufft);
	if (dev_Gyx_cufft)   cudaFree(dev_Gyx_cufft);
	if (dev_Gyy_cufft)   cudaFree(dev_Gyy_cufft);
	if (dev_Gyz_cufft)   cudaFree(dev_Gyz_cufft);
	if (dev_Gzx_cufft)   cudaFree(dev_Gzx_cufft);
	if (dev_Gzy_cufft)   cudaFree(dev_Gzy_cufft);
	if (dev_Gzz_cufft)   cudaFree(dev_Gzz_cufft);
	if (dev_Mx_cufft)    cudaFree(dev_Mx_cufft);
	if (dev_My_cufft)    cudaFree(dev_My_cufft);
	if (dev_Mz_cufft)    cudaFree(dev_Mz_cufft);
	if (dev_Hd_x_cufft)  cudaFree(dev_Hd_x_cufft);
	if (dev_Hd_y_cufft)  cudaFree(dev_Hd_y_cufft);
	if (dev_Hd_z_cufft)  cudaFree(dev_Hd_z_cufft);
	if (dev_Mag0_z)      cudaFree(dev_Mag0_z);

	//---- free memory -----//
	if (Hd_x_1d)        free (Hd_x_1d);
	if (Hd_y_1d)        free(Hd_y_1d);
	if (Hd_z_1d)        free(Hd_z_1d);
	if (Hd_x_1d_cmplx)  free(Hd_x_1d_cmplx);
	if (Hd_y_1d_cmplx)  free(Hd_y_1d_cmplx);
	if (Hd_z_1d_cmplx)  free(Hd_z_1d_cmplx);
	if (Mx_1d_cmplx)    free(Mx_1d_cmplx);
	if (My_1d_cmplx)    free(My_1d_cmplx);
	if (Mz_1d_cmplx)    free(Mz_1d_cmplx);
	if (Gxx_1d_cmplx)   free(Gxx_1d_cmplx);
	if (Gxy_1d_cmplx)   free(Gxy_1d_cmplx);
	if (Gxz_1d_cmplx)   free(Gxz_1d_cmplx);
	if (Gyx_1d_cmplx)   free(Gyx_1d_cmplx);
	if (Gyy_1d_cmplx)   free(Gyy_1d_cmplx);
	if (Gyz_1d_cmplx)   free(Gyz_1d_cmplx);
	if (Gzx_1d_cmplx)   free(Gzx_1d_cmplx);
	if (Gzy_1d_cmplx)   free(Gzy_1d_cmplx);
	if (Gzz_1d_cmplx)   free(Gzz_1d_cmplx);
	if (Hd_x_1d_shift)  free(Hd_x_1d_shift);
	if (Hd_y_1d_shift)  free(Hd_y_1d_shift);
	if (Hd_z_1d_shift)  free(Hd_z_1d_shift);
	if (host_Aex)       free(host_Aex);
	if (host_Ku)        free(host_Ku);
	if (host_Ms)        free(host_Ms);
	if (host_alpha)     free(host_alpha);
	if (host_gamma)     free(host_gamma);
	if (indicator1)     free(indicator1);
	if (indicator2)     free(indicator2);
	if (indicator3)     free(indicator3);
	if (indicator4)     free(indicator4);
	if (indicator5)     free(indicator5);
	if (indicator6)     free(indicator6);
	if (T)              free(T);
	if (D)              free(D);
	if (Happl_x)        free(Happl_x);
	if (Happl_y)        free(Happl_y);
	if (Happl_z)        free(Happl_z);
	if (Mx)             free(Mx);
	if (My)             free(My);
	if (Mz)             free(Mz);
	if (watch100)         free(watch100);
	if (watch2)         free(watch2);
	if (watch3)         free(watch3);
	if (watch4)         free(watch4);
	if (watch5)         free(watch5);
	if (watch100_int)     free(watch100_int);
	if (watch2_int)     free(watch2_int);
	if (watch3_int)     free(watch3_int);
	if (Mx_float)       free(Mx_float);
	if (My_float)       free(My_float);
	if (Mz_float)       free(Mz_float);
	if (Mx_t_bar)       free(Mx_t_bar);
	if (My_t_bar)       free(My_t_bar);
	if (Mz_t_bar)       free(Mz_t_bar);
	if (M_bar)          free(M_bar);
	if (Mx_bar)         free(Mx_bar);
	if (My_bar)         free(My_bar);
	if (Mz_bar)         free(Mz_bar);
	if (Mx_bar1)        free(Mx_bar1);
	if (My_bar1)        free(My_bar1);
	if (Mz_bar1)        free(Mz_bar1);
	if (Mx_bar2)        free(Mx_bar2);
	if (My_bar2)        free(My_bar2);
	if (Mz_bar2)        free(Mz_bar2);
	if (Mx_bar3)        free(Mx_bar3);
	if (My_bar3)        free(My_bar3);
	if (Mz_bar3)        free(Mz_bar3);
	if (Mx_bar4)        free(Mx_bar4);
	if (My_bar4)        free(My_bar4);
	if (Mz_bar4)        free(Mz_bar4);
	if (Mx_bar5)        free(Mx_bar5);
	if (My_bar5)        free(My_bar5);
	if (Mz_bar5)        free(Mz_bar5);
	if (Mx_bar6)        free(Mx_bar6);
	if (My_bar6)        free(My_bar6);
	if (Mz_bar6)        free(Mz_bar6);
	if (Mz_CGC)         free(Mz_CGC);
	if (Mz_LayerLayer)  free(Mz_LayerLayer);
	if (M1_z_SingleGrain)      free(M1_z_SingleGrain);
	if (Hint12_x_SingleGrain)  free(Hint12_x_SingleGrain);
	if (Hint12_y_SingleGrain)  free(Hint12_y_SingleGrain);
	if (Hint12_z_SingleGrain)  free(Hint12_z_SingleGrain);
	if (std_Ku)         free(std_Ku);
	if (std_Aex)         free(std_Aex);
	if (std_Tc)         free(std_Tc);
	if (gasarray)       free(gasarray);

	//cluster
	if (rrt_cor)        free(rrt_cor);
	if (Mz_tt)          free(Mz_tt);

	if (FP_theta)       free(FP_theta);
	if (FP_mag)       free(FP_mag);
	
}

const std::string currentDateTime(){
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	return buf;
}

int main(int argc, char* argv[])
{
	int itt = -1;
	float mm;
	int m,  lx_zero_pad, ly_zero_pad, lz_zero_pad;
	curandGenerator_t gen;
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	

	// Initial input
	if (!Input()) printf("Input() failed!\n");

	// Preliminary Settings for FFT
	lx_zero_pad = 2*Rnd_upto_pow2(mNx);
	ly_zero_pad = 2*Rnd_upto_pow2(mNy);
	lz_zero_pad = 2*Rnd_upto_pow2(mNz);

	
	// Allocate memory
	cudaSetDevice(deviceID);
	
	if (!AllocateDeviceMemory(lx_zero_pad, ly_zero_pad, lz_zero_pad, mNx, mNy, mNz, DEG_FREEDOM)) { printf("AllocateDeviceMemory() failed!\n"); }
	
	if (!AllocateHostMemory(lx_zero_pad, ly_zero_pad, lz_zero_pad, mNx, mNy, mNz, DEG_FREEDOM, TOTAL_TIME)) { printf("AllocateHostMemory() failed!\n"); }
	
	
	// Calculate media cell indicator
	if (!GrainsIndicator()) { printf("GrainsIndicator() failed!\n"); }

	if (!Input_Indicator_1D_Format(mNx, mNy, mNz, indicator1, indicator2, indicator3, indicator4, indicator5, indicator6, indicator7, "indicator.inp"))
	{ printf("Input_Indicator_1D_Format() failed!\n"); }
	

	// Define arrays for anaylyzing Hc for all of the grains
	#ifndef __POST_ANALYSIS_
	double *Happl1_sweep = (double*)calloc(ceil(TOTAL_TIME/FieldSweepTimeStep)+1, sizeof(double));
	vector< vector<double> > M1_z_SingleGrain_field(0);
	int NumOfGrains;
	
	if (!VORO_GRAIN){
		NumOfGrains=mNx*mNy;
		M1_z_SingleGrain_field.resize(NumOfGrains);
		for (int i = 0; i < NumOfGrains; i++) M1_z_SingleGrain_field[i].resize(ceil(TOTAL_TIME/FieldSweepTimeStep)+1, 0.0);
		M1_z_SingleGrain = (double*) calloc(NumOfGrains, sizeof(double));
	}
	
	// Analyze how many voronoi grains
	vector<vector<coor2d_t> > grain_coor(0);
	if (VORO_GRAIN){
		vector<int> sorted_indicator2(mNx*mNy, 0);
		for (int idx=0; idx<mNx*mNy; idx++) sorted_indicator2[idx] = indicator2[idx];
		sort(sorted_indicator2.begin(), sorted_indicator2.end());
		NumOfGrains = *(sorted_indicator2.end()-1);
		M1_z_SingleGrain_field.resize(NumOfGrains);
		for (int i = 0; i < NumOfGrains; i++) M1_z_SingleGrain_field[i].resize(ceil(TOTAL_TIME/FieldSweepTimeStep)+1, 0.0);
		M1_z_SingleGrain = (double*) calloc(NumOfGrains, sizeof(double));

		// Assign each grain grid coordinates
		grain_coor.resize(NumOfGrains);
		for (int i_grain=1; i_grain<=NumOfGrains; i_grain++){
			for (int j=0; j<mNy; j++){
				for (int i=0; i<mNx; i++){
					int idx = i + j*mNx;
					if (indicator2[idx] == i_grain){
						coor2d_t tmp;
						tmp.x = i;
						tmp.y = j;
						grain_coor[i_grain-1].push_back(tmp);
					}
				}
			}
		}
	}
	#endif


	
	
	
	
	




	
	


	




	// Write to files
	FILE *wfile3, *wfile4, *wfile5, *wfile6, *wfile7, *wfile8, *wfile10, *wfile11, *wfile12, *wfile15;
	wfile3  = fopen("Mz_LayerLayer_Happl.out", "w");
	wfile4  = fopen("Mz_LayerAve_Happl.out", "w");
	wfile5  = fopen("Mz1_SingleGrain_Happl.out", "w");
	wfile6  = fopen("Hint12_x_SingleGrain_Happl.out", "w");
	wfile7  = fopen("Hint12_y_SingleGrain_Happl.out", "w");
	wfile8  = fopen("Hint12_z_SingleGrain_Happl.out", "w");
	wfile10 = fopen("Output_Spatial_MxMyMz1.out", "w");
	wfile11 = fopen("Output_Spatial_MxMyMz2.out", "w");
	wfile12 = fopen("Output_Spatial_MxMyMz3.out", "w");
	wfile15 = fopen("Runlog.out", "w");
	//-----------------------------------------------------//
	

	



	//--------I/O------------------------------------------//
	if (!Input_FieldProfile(EXT_FP_PROFILE)) { printf("Input_FieldProfile() failed!\n"); }
	if (EXT_FP_PROFILE) {if (!FP_inp_theta_n_Magnitude()){ printf("FP_inp_theta_n_Magnitude() failed!"); } }
	if (EXT_FP_PROFILE) {if (!Output_Float_3D_Format(fNx, fNy, fNz, FP_mag, "FP_inp_mag.out")) { printf("Output_Float_3D_Format() failed!\n"); } }
	if (EXT_FP_PROFILE) {if (!Output_Float_3D_Format(fNx, fNy, fNz, FP_theta, "FP_inp_theta.out")) { printf("Output_Float_3D_Format() failed!\n"); } }
	if (!Input_FieldSeq(EXT_FP_PROFILE))     { printf("Input_FieldSeq() failed!\n"); }
	if (!Input_Mag0(PRE_DEFINED_PATTERN))    { printf("Input_Mag0() failed!\n"); }
	if (!MH_LOOP){ if (!Output_Float_3D_Format(fNx, fNy, fNz, FP_inp[2], "Field_z_Inp.out")) { printf("OutputFloat_3D_format() failed!\n"); } }
	if (!MH_LOOP){ if (!Output_Float_3D_Format(fNx, fNy, fNz, FP_inp[1], "Field_y_Inp.out")) { printf("OutputFloat_3D_format() failed!\n"); } }
	if (!MH_LOOP){ if (!Output_Float_3D_Format(fNx, fNy, fNz, FP_inp[0], "Field_x_Inp.out")) { printf("OutputFloat_3D_format() failed!\n"); } }
	if (!Set_delta_K_Aex_Tc(std_Ku, std_Aex, std_Tc, indicator4, indicator5, indicator6)) { printf("Set_delta_K_Aex_Tc() failed!\n"); }  // Set delta_K, delta_Aex, delta_Tc
	if (!Output_Int_3D_Format(mNx, mNy, mNz, indicator1, "LayersFlag.out")) { printf("Output_Int_3D_Format() failed!\n"); }
	if (!Output_Int_3D_Format(mNx, mNy, mNz, indicator2, "Output_VoroCell_Map.out")) { printf("Output_Int_3D_Format() failed!\n"); }
	if (!Output_Float_3D_Format(mNx, mNy, mNz, indicator4, "Output_Std_Map.out")) { printf("Output_Float_3D_Format() failed!\n"); }
	if (!Output_Int_3D_Format(mNx, mNy, mNz, indicator3, "Output_Ini_State.out")) { printf("Output_Int_3D_Format() failed!\n"); }
	if (!Output_Int_3D_Format(mNx, mNy, mNz, indicator7, "Output_Defect_Map.out")) { printf("Output_Int_3D_Format() failed!\n"); }
	if (DT_Rec_Analysis) { if (!Set_Happl_DT()) { printf("Set_Happl() failed!\n"); } }
	if (CT_Rec_Analysis) { if (!Set_Happl_CT()) { printf("Set_Happl() failed!\n"); } }
	if (!Output_Float_3D_Format(mNx, mNy, mNz, Happl_z, "Happl_z.out")) { printf("OutputFloat_3D_format() failed!\n"); }
	if (!Output_Float_3D_Format(mNx, mNy, mNz, Happl_y, "Happl_y.out")) { printf("OutputFloat_3D_format() failed!\n"); }
	if (!Output_Float_3D_Format(mNx, mNy, mNz, Happl_x, "Happl_x.out")) { printf("OutputFloat_3D_format() failed!\n"); }
	//-----------------------------------------------------//
	
	





	int count;
	if (!VORO_GRAIN){
		count = 0;
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				if (indicator7[i+j*mNx+(mNz_6-1)*mNx*mNy] == 0) count++;
			}
		}
		printf("count=%d\n", count);
	}

	if (VORO_GRAIN){
		/*count = 0;
		for (int i=0; i<grain_coor.size(); i++){
			printf("%d\n",i);
			coor2d_t tmp=grain_coor[i][0];
			printf("%d\n",i);
			int idx = tmp.x + tmp.y*mNx + (mNz_6-1)*mNx*mNy;
			printf("%d\n",idx);
			if (indicator7[idx] == 0) count++;
		}
		printf("count=%d\n", count);*/
		count = grain_coor.size();
	}
	
	
	clock_t start, stop;
	double time=0.0;
	assert((start = clock())!=-1);

	


	
	//-------- Setup kernel configuration-------------------//
	dim3 blocks1(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z); // for real medium dimensions
	int  grid1_x = (mNx % BLOCK_SIZE_X) ? (mNx/BLOCK_SIZE_X + 1) : (mNx/BLOCK_SIZE_X),
		 grid1_y = (mNy % BLOCK_SIZE_Y) ? (mNy/BLOCK_SIZE_Y + 1) : (mNy/BLOCK_SIZE_Y),
		 grid1_z = (mNz % BLOCK_SIZE_Z) ? (mNz/BLOCK_SIZE_Z + 1) : (mNz/BLOCK_SIZE_Z);
	dim3 grids1(grid1_x, grid1_y, grid1_z);
	
	dim3 blocks2(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLK_SZ_Z);  // for CUFFT dimensions
	int  grid2_x = ((2*mNx) % BLOCK_SIZE_X) ? ((2*mNx)/BLOCK_SIZE_X + 1) : ((2*mNx)/BLOCK_SIZE_X),
		 grid2_y = ((2*mNy) % BLOCK_SIZE_Y) ? ((2*mNy)/BLOCK_SIZE_Y + 1) : ((2*mNy)/BLOCK_SIZE_Y),
		 grid2_z = (lz_zero_pad % BLK_SZ_Z) ? (lz_zero_pad/BLK_SZ_Z + 1) : (lz_zero_pad/BLK_SZ_Z);
	dim3 grids2(grid2_x, grid2_y, grid2_z);
	
	dim3 blocks3(2, 2, mNz+2); // medium-with-apron dimensions
	int  grid3_x = ((mNx+2) % 2) ? ((mNx+2)/2 + 1) : ((mNx+2)/2),
		 grid3_y = ((mNy+2) % 2) ? ((mNy+2)/2 + 1) : ((mNy+2)/2),
		 grid3_z = ((mNz+2) % (mNz+2)) ? ((mNz+2)/(mNz+2) + 1) : ((mNz+2)/(mNz+2));
	dim3 grids3(grid3_x, grid3_y, grid3_z);
	//-----------------------------------------------------//


	//-------------- Initial condition---------------------------------//
	for (int k = 0; k < mNz; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				int idx = i + j*mNx + k*mNx*mNy;
				Hd_x_1d_shift[idx] = 0.0;
				Hd_y_1d_shift[idx] = 0.0;
				Hd_z_1d_shift[idx] = 0.0;
			}
		}
	}
	for (int k = 0; k < lz_zero_pad; k++){
		for (int j = 0; j < ly_zero_pad; j++){
			for (int i = 0; i < lx_zero_pad; i++){
				int idx = i + j*lx_zero_pad + k*lx_zero_pad*ly_zero_pad;
				Gxx_1d_cmplx[idx] = 0.0;
				Gxy_1d_cmplx[idx] = 0.0;
				Gxz_1d_cmplx[idx] = 0.0;
				Gyx_1d_cmplx[idx] = 0.0;
				Gyy_1d_cmplx[idx] = 0.0;
				Gyz_1d_cmplx[idx] = 0.0;
				Gzx_1d_cmplx[idx] = 0.0;
				Gzy_1d_cmplx[idx] = 0.0;
				Gzz_1d_cmplx[idx] = 0.0;
			}
		}
	}



	Happl = field1*IniScalFact;
	for (unsigned long tt = 0; tt < TOTAL_TIME; tt++)
	{
		Mx_bar1[tt] = 0.0;
		My_bar1[tt] = 0.0;
		Mz_bar1[tt] = 0.0;
		Mx_bar2[tt] = 0.0;
		My_bar2[tt] = 0.0;
		Mz_bar2[tt] = 0.0;
		Mx_bar3[tt] = 0.0;
		My_bar3[tt] = 0.0;
		Mz_bar3[tt] = 0.0;
		Mx_bar4[tt] = 0.0;
		My_bar4[tt] = 0.0;
		Mz_bar4[tt] = 0.0;
		Mx_bar5[tt] = 0.0;
		My_bar5[tt] = 0.0;
		Mz_bar5[tt] = 0.0;
		Mx_bar6[tt] = 0.0;
		My_bar6[tt] = 0.0;
		Mz_bar6[tt] = 0.0;
		Mx_bar[tt] = 0.0;
		My_bar[tt] = 0.0;
		Mz_bar[tt] = 0.0;
		M_bar[tt] = 0.0;
	}
	M_bar_t = 0;
	cudaMemcpy(dev_indicator1, indicator1, (mNx)*(mNy)*(mNz)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_indicator2, indicator2, (mNx)*(mNy)*(mNz)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_indicator3, indicator3, (mNx)*(mNy)*(mNz)*sizeof(int), cudaMemcpyHostToDevice);
	
	Kernel_Initialization<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
		                                       dev_theta, dev_phi,
											   dev_a_theta, dev_b_theta, dev_c_theta, dev_d_theta,
											   dev_a_phi, dev_b_phi, dev_c_phi, dev_d_phi,
											   dev_Ha_x, dev_Ha_y, dev_Ha_z,
											   dev_Hal_x, dev_Hal_y, dev_Hal_z,
											   dev_Hth_x, dev_Hth_y, dev_Hth_z,
											   dev_Hk_x, dev_Hk_y, dev_Hk_z,
											   dev_Hd_x, dev_Hd_y, dev_Hd_z,
											   dev_d_theta_d_t, dev_d_phi_d_t,
											   dev_indicator1, dev_indicator3);

	if (AFC==1) Kernel_Initialization_AFC<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
															   AF_layer_label,
															   dev_theta, dev_indicator1, dev_indicator3);

	if (PRE_DEFINED_PATTERN){ cudaMemcpy(dev_theta, theta0, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice); }
	Kernel_Initialization_Cplx<<<grids2, blocks2>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
		                                            dev_Mx_cufft,   dev_My_cufft,   dev_Mz_cufft,
                                                    dev_Gxx_cufft,  dev_Gxy_cufft,  dev_Gxz_cufft,
													dev_Gyx_cufft,  dev_Gyy_cufft,  dev_Gyz_cufft,
													dev_Gzx_cufft,  dev_Gzy_cufft,  dev_Gzz_cufft,
												    dev_Hd_x_cufft, dev_Hd_y_cufft, dev_Hd_z_cufft);
	Kernel_Ini_indicator_temp<<<grids3, blocks3>>>(2, 2, 2, mNx, mNy,
		                                           dev_indicator1_temp, dev_indicator2_temp, 
												   dev_Ms_temp, dev_M_temp_x, dev_M_temp_y, dev_M_temp_z, dev_Aex_temp);
	if (FULL_REC) { if (!External_HeadField(EXT_FP_PROFILE)) { printf("!External_HeadField() failed!\n"); } } // Set Initial Field Profile to the Medium
	//---------------------------------------------------------------//

	// Demag Tensor
	if (DEMAG){	if (!G_tensor(lx_zero_pad, ly_zero_pad, lz_zero_pad)) { printf("!G_tensor() failed!\n"); }	}
	


	Kernel_dev_indicator_with_apron<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy, mNz,
		                                                 dev_indicator1, dev_indicator1_temp, 
												         dev_indicator2, dev_indicator2_temp);
	

    #ifndef __WATCH__
	// Watch
	/*std::ofstream pfile1;
	pfile1.open("watch.out");
	for (int idx=0; idx<lx_zero_pad*ly_zero_pad*lz_zero_pad; ++idx){
		if (idx%lx_zero_pad==0 && idx!=0) pfile1 << "" <<std::endl;
		if (idx%(lx_zero_pad*ly_zero_pad)==0 && idx!=0) pfile1 << "" <<std::endl << "" <<std::endl;
		pfile1 << std::setw(15) << Gxx_1d_cmplx[idx].real();
	}
	pfile1.close();*/
	#endif
	
	// Create pseudo-random number generator
	if (curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT) != CURAND_STATUS_SUCCESS){ printf("CURAND_CREATE() failed!\n"); };
	// Set random seed
	if (curandSetPseudoRandomGeneratorSeed(gen, iseed) != CURAND_STATUS_SUCCESS){ printf("CURAND_SEED() failed!\n"); };



		
	// Start time evolution
	*idx_f = -1;
	*idx_t = -1;
	double SweepFieldStep1;
	SweepFieldStep1 = SweepFieldStep;
	if (FULL_REC) { if (!Moving_Head(0)) { printf("Moving_Head() failed!\n"); } }
	if (MH_LOOP) { if (!NonMoving_Head(0, SweepFieldStep1)) { printf("NonMoving_Head() failed!\n"); } }
	if (!Set_devMxyz(host_Ms, dev_theta, dev_phi, dev_Mx, dev_My, dev_Mz)) printf("Set_devMxyz() failed!\n");
	
	#ifndef __WATCH__
	//---------Watch---------//
	/*double* watch1 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	cudaMemcpy(watch1, dev_theta, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
	FILE *pfile1 = fopen("watch1.out", "w");
	for (int k = 0; k < (mNz); k++){
		for (int j = 0; j < (mNy); j++){
			for (int i = 0; i < (mNx); i++){
				fprintf(pfile1, "%10.3f", watch1[i+j*(mNx)+k*(mNx)*(mNy)]);
			}
			fprintf(pfile1, "\n");
		}
		fprintf(pfile1, "\n\n");
	}
	fclose(pfile1);
	printf("watch1\n");*/
	////
	//double* watch2 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	//cudaMemcpy(watch2, dev_Mx, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
	//FILE *pfile2 = fopen("watch2.out", "w");
	//for (int k = 0; k < (mNz); k++){
	//	for (int j = 0; j < (mNy); j++){
	//		for (int i = 0; i < (mNx); i++){
	//			fprintf(pfile2, "%15.3e", watch2[i+j*(mNx)+k*(mNx)*(mNy)]);
	//		}
	//		fprintf(pfile2, "\n");
	//	}
	//	fprintf(pfile2, "\n\n");
	//}
	//fclose(pfile2);
	//printf("watch2\n");
	////---------Watch---------//
	#endif

	if (DT_Rec_Analysis) { if (!DT_Analysis_Head(0)) { printf("2D_Analysis_Head() failed!\n"); } }
	if (CT_Rec_Analysis) { if (!CT_Analysis_Head(0)) { printf("2D_Analysis_Head() failed!\n"); } }
	if (FULL_REC){ if (!Output_Float_3D_Format(mNx, mNy, mNz, FP[2], "FPz.out")) { printf("OutputFloat_3D_format() failed!\n"); } } //write to output file
	if (FULL_REC){ if (!Output_Float_3D_Format(mNx, mNy, mNz, FP_trail[2], "FPz_trail.out")) { printf("OutputFloat_3D_format() failed!\n"); } } //write to output file
	if (FULL_REC){ if (!Output_Float_3D_Format(mNx, mNy, mNz, Happl_z, "Output_Spatial_Hz_temp.out")) { printf("OutputFloat_3D_format() failed!\n"); } } //write to output file

    // Start iteration
	for (unsigned long tt = 0; tt < TOTAL_TIME; tt++){

		// Moving the head
		if (tt > 0 && tt < EQUI_START_TIME) { 
			/*if (MH_LOOP){
				if (tt > FieldSweepTimeStep && (abs(Mz_CGC[tt-1]-Mz_CGC[(tt-1)-FieldSweepTimeStep]) >= abs(Mz_CGC[(tt-1)-FieldSweepTimeStep])*0.02) && (tt%FieldSweepTimeStep == 0)){
					SweepFieldStep1 = SweepFieldStep/ScalFactFieldStep;
				}
			}*/
			if (FULL_REC) { if (!Moving_Head(tt)) { printf("Moving_Head() failed!\n"); }; }
			if (MH_LOOP) { if (!NonMoving_Head(tt, SweepFieldStep1)) { printf("NonMoving_Head() failed!\n"); }; }
		}
		else FieldReset();
	



		// Demag. Field Calculated by FFT (slow)
		
		#ifndef __WATCH__
		//---------Watch---------//
		/*double* watch1 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		cudaMemcpy(watch1, dev_Mz, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		FILE *pfile1 = fopen("watch1.out", "w");
		for (int k = 0; k < (mNz); k++){
			for (int j = 0; j < (mNy); j++){
				for (int i = 0; i < (mNx); i++){
					fprintf(pfile1, "%15.3e", watch1[i+j*(mNx)+k*(mNx)*(mNy)]);
				}
				fprintf(pfile1, "\n");
			}
			fprintf(pfile1, "\n\n");
		}
		fclose(pfile1);
		printf("watch1\n");*/
		////
		//double* watch2 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch2, dev_Mx, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile2 = fopen("watch2.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile2, "%15.3e", watch2[i+j*(mNx)+k*(mNx)*(mNy)]);
		//		}
		//		fprintf(pfile2, "\n");
		//	}
		//	fprintf(pfile2, "\n\n");
		//}
		//fclose(pfile2);
		//printf("watch2\n");
		////---------Watch---------//

		//---------Watch---------//
		//double* watch1 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		////cudaMemcpy(watch1, Mz, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile1 = fopen("watch1.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile1, "%15.3e", Mz[i+j*(mNx)+k*(mNx)*(mNy)]);
		//		}
		//		fprintf(pfile1, "\n");
		//	}
		//	fprintf(pfile1, "\n\n");
		//}
		//fclose(pfile1);
		//printf("watch1\n");
		////
		//double* watch2 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		////cudaMemcpy(watch2, dev_phi, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile2 = fopen("watch2.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile2, "%15.3e", Mx[i+j*(mNx)+k*(mNx)*(mNy)]);
		//		}
		//		fprintf(pfile2, "\n");
		//	}
		//	fprintf(pfile2, "\n\n");
		//}
		//fclose(pfile2);
		//printf("watch2\n");
		//---------Watch---------//
		#endif
		
		#ifndef __DEMAG__
		// Demag. Field Calculated by FFT (slow)
		if (DEMAG && tt%DmagINT==0){
			for (int k = 0; k < mNz; k++){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){
						int idx = i + j*mNx + k*mNx*mNy;
						int idxx = i + j * lx_zero_pad + k * lx_zero_pad*ly_zero_pad;
						// Exclude demag field contribution from break layers
						if (indicator1[idx]==12 || indicator1[idx]==23 || indicator1[idx]==34 || indicator1[idx]==45 || indicator1[idx]==56){
							Mx_1d_cmplx[idxx].real((float)pow(3,0.5));
							My_1d_cmplx[idxx].real((float)pow(3,0.5));
							Mz_1d_cmplx[idxx].real((float)pow(3,0.5));
							Mx_1d_cmplx[idxx].imag(0.);
							My_1d_cmplx[idxx].imag(0.);
							Mz_1d_cmplx[idxx].imag(0.);
						}
						else {
							Mx_1d_cmplx[idxx].real((float)Mx[idx]);
							My_1d_cmplx[idxx].real((float)My[idx]);
							Mz_1d_cmplx[idxx].real((float)Mz[idx]);
							Mx_1d_cmplx[idxx].imag(0.);
							My_1d_cmplx[idxx].imag(0.);
							Mz_1d_cmplx[idxx].imag(0.);
						}
					}
				}
			}

			////---------Watch---------//
			//double* watch1 = (double*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad, sizeof(double));
			////cudaMemcpy(watch1, Mz, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
			//FILE *pfile1 = fopen("watch1.out", "w");
			//for (int k = 0; k < (lz_zero_pad); k++){
			//	for (int j = 0; j < (ly_zero_pad); j++){
			//		for (int i = 0; i < (lx_zero_pad); i++){
			//			fprintf(pfile1, "%15.3e", Mz_1d_cmplx[i+j*(lx_zero_pad)+k*(lx_zero_pad)*(ly_zero_pad)].real());
			//		}
			//		fprintf(pfile1, "\n");
			//	}
			//	fprintf(pfile1, "\n\n");
			//}
			//fclose(pfile1);
			//printf("watch1\n");
			////
			//double* watch2 = (double*) calloc(lx_zero_pad*ly_zero_pad*lz_zero_pad, sizeof(double));
			////cudaMemcpy(watch2, dev_phi, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
			//FILE *pfile2 = fopen("watch2.out", "w");
			//for (int k = 0; k < (lz_zero_pad); k++){
			//	for (int j = 0; j < (ly_zero_pad); j++){
			//		for (int i = 0; i < (lx_zero_pad); i++){
			//			fprintf(pfile2, "%15.3e", Mz_1d_cmplx[i+j*(lx_zero_pad)+k*(lx_zero_pad)*(ly_zero_pad)].imag());
			//		}
			//		fprintf(pfile2, "\n");
			//	}
			//	fprintf(pfile2, "\n\n");
			//}
			//fclose(pfile2);
			//printf("watch2\n");
			////---------Watch---------//

			#ifndef __WATCH__
			// Watch
			/*std::ofstream pfile1;
			pfile1.open("watch.out");
			for (int idx=0; idx<lx_zero_pad*ly_zero_pad*lz_zero_pad; ++idx){
				if (idx%lx_zero_pad==0 && idx!=0) pfile1 << "" <<std::endl;
				if (idx%(lx_zero_pad*ly_zero_pad)==0 && idx!=0) pfile1 << "" <<std::endl << "" <<std::endl;
				pfile1 << std::setw(15) << Mz_1d_cmplx[idx].real();
			}
			pfile1.close();*/
			#endif
			
			cudaMemcpy(dev_Mx_cufft, Mx_1d_cmplx, (lx_zero_pad*ly_zero_pad*lz_zero_pad)*sizeof(cufftComplex), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_My_cufft, My_1d_cmplx, (lx_zero_pad*ly_zero_pad*lz_zero_pad)*sizeof(cufftComplex), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Mz_cufft, Mz_1d_cmplx, (lx_zero_pad*ly_zero_pad*lz_zero_pad)*sizeof(cufftComplex), cudaMemcpyHostToDevice);
			Hms(lx_zero_pad,   ly_zero_pad,   lz_zero_pad,
		        grids2, blocks2);
			cudaMemcpy(dev_Hd_x, Hd_x_1d_shift, mNx*mNy*mNz*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Hd_y, Hd_y_1d_shift, mNx*mNy*mNz*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Hd_z, Hd_z_1d_shift, mNx*mNy*mNz*sizeof(double), cudaMemcpyHostToDevice);
		}
		#endif
		
		#ifndef __WATCH__
		// Watch
		/*std::ofstream pfile1;
		pfile1.open("watch.out");
		for (int idx=0; idx<mNx*mNy*mNz; ++idx){
			if (idx%mNx==0 && idx!=0) pfile1 << "" <<std::endl;
			if (idx%(mNx*mNy)==0 && idx!=0) pfile1 << "" <<std::endl << "" <<std::endl;
			pfile1 << std::setw(15) << Hd_z_1d_shift[idx];
		}
		pfile1.close();*/
		#endif
	
		#ifndef __THERMAL__
		// Generate random floats on device
		if (THERMAL){
			if (curandGenerateNormalDouble(gen, dev_GasArray, DEG_FREEDOM, 0, 1) != CURAND_STATUS_SUCCESS){ printf("CURAND_GENERATE() failed!\n"); };
			Kernel_Hth_field<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                      dev_GasArray, dev_Hth_x, dev_Hth_y, dev_Hth_z, dev_D);
		}
		#endif

        #ifndef __WATCH__
		////---------Watch---------//
		//// Watch1
		//double* watch1 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch1, dev_Ms, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile1 = fopen("watch1.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile1, "%15.3e", watch1[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile1, "\n");
		//	}
		//	fprintf(pfile1, "\n\n");
		//}
		//fclose(pfile1);
		//printf("watch1\n");
		//// Watch2
		//double* watch2 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch2, dev_Ku, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile2 = fopen("watch2.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile2, "%15.3e", watch2[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile2, "\n");
		//	}
		//	fprintf(pfile2, "\n\n");
		//}
		//fclose(pfile2);
		//printf("watch2\n");
		//// Watch3
		//double* watch3 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch3, dev_Aex, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile3 = fopen("watch3.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile3, "%15.3e", watch3[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile3, "\n");
		//	}
		//	fprintf(pfile3, "\n\n");
		//}
		//fclose(pfile3);
		//printf("watch3\n");
		//// Watch4
		//double* watch4 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch4, dev_alpha, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile4 = fopen("watch4.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile4, "%15.3e", watch4[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile4, "\n");
		//	}
		//	fprintf(pfile4, "\n\n");
		//}
		//fclose(pfile4);
		//printf("watch4\n");
		//// Watch5
		//double* watch5 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch5, dev_D, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile5 = fopen("watch5.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile5, "%15.3e", watch5[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile5, "\n");
		//	}
		//	fprintf(pfile5, "\n\n");
		//}
		//fclose(pfile5);
		//printf("watch5\n");
		//// Watch6
		//double* watch6 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch6, dev_gamma, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile6 = fopen("watch6.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile6, "%15.3e", watch6[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile6, "\n");
		//	}
		//	fprintf(pfile6, "\n\n");
		//}
		//fclose(pfile6);
		//printf("watch6\n");
		//// Watch7
		//double* watch7 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch7, dev_T, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile7 = fopen("watch7.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile7, "%15.3e", watch7[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile7, "\n");
		//	}
		//	fprintf(pfile7, "\n\n");
		//}
		//fclose(pfile7);
		//printf("watch7\n");
		////---------Watch---------//
		#endif

	
		#ifndef __RUNGE_KUTTA__
        // Here a_theta, a_phi make no effect 
		if (!Ha_field(ZERO, dev_a_theta, dev_a_phi, grids1, blocks1)) { printf("Ha_field(ZERO) failed!\n"); }
	
		Kernel_d_theta_phi_d_t<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                        dev_d_theta_d_t, dev_d_phi_d_t,
											        dev_theta,  dev_phi, 
												    dev_a_theta, dev_a_phi,
												    dev_Ha_x, dev_Ha_y, dev_Ha_z,
												    dev_Hth_x, dev_Hth_y, dev_Hth_z,
												    dev_Hd_x, dev_Hd_y, dev_Hd_z,
												    dev_Happl_x, dev_Happl_y, dev_Happl_z,
												    dev_Ku, dev_Ms, dev_alpha, dev_gamma);

		//// Watch1
		//double* watch1 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch1, dev_theta, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile1 = fopen("watch1.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile1, "%16.6e", watch1[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile1, "\n");
		//	}
		//	fprintf(pfile1, "\n\n");
		//}
		//fclose(pfile1);
		//printf("watch1\n");


		Kernel_a_theta_phi<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                    dev_a_theta, dev_a_phi,
											    dev_d_theta_d_t, dev_d_phi_d_t);

		if (!Ha_field(delta_t/2, dev_a_theta, dev_a_phi, grids1, blocks1)) { printf("Ha_field_a(delta_t/2) failed!\n"); }
		Kernel_b_theta_phi<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                    dev_b_theta, dev_b_phi,                                   
											    dev_theta, dev_phi,                                      
											    dev_a_theta, dev_a_phi,
											    dev_Ha_x, dev_Ha_y, dev_Ha_z,
											    dev_Hth_x, dev_Hth_y, dev_Hth_z,
											    dev_Hd_x, dev_Hd_y, dev_Hd_z,
											    dev_Happl_x, dev_Happl_y, dev_Happl_z,
											    dev_Ku, dev_Ms, dev_alpha, dev_gamma, delta_t/2);
		
			
		if (!Ha_field(delta_t/2, dev_b_theta, dev_b_phi, grids1, blocks1)) { printf("Ha_field_a(delta_t/2) failed!\n"); }
		Kernel_c_theta_phi<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                    dev_c_theta, dev_c_phi,
											    dev_theta, dev_phi, 
	 										    dev_b_theta, dev_b_phi,
											    dev_Ha_x, dev_Ha_y, dev_Ha_z,
											    dev_Hth_x, dev_Hth_y, dev_Hth_z,
											    dev_Hd_x, dev_Hd_y, dev_Hd_z,
											    dev_Happl_x, dev_Happl_y, dev_Happl_z,
											    dev_Ku, dev_Ms, dev_alpha, dev_gamma, delta_t/2);
		
		if (!Ha_field(delta_t,   dev_c_theta, dev_c_phi, grids1, blocks1)) { printf("Ha_field_a(delta_t) failed!\n"); }		
		Kernel_d_theta_phi<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                    dev_d_theta, dev_d_phi,
										    	dev_theta, dev_phi, 
											    dev_c_theta, dev_c_phi,
											    dev_Ha_x, dev_Ha_y, dev_Ha_z,
											    dev_Hth_x, dev_Hth_y, dev_Hth_z,
											    dev_Hd_x, dev_Hd_y, dev_Hd_z,
											    dev_Happl_x, dev_Happl_y, dev_Happl_z,
											    dev_Ku, dev_Ms, dev_alpha, dev_gamma, delta_t);
		



		
		// Time evolution
		Kernel_time_increment<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                       dev_theta, dev_phi,
		     									   dev_a_theta, dev_b_theta, dev_c_theta, dev_d_theta,   
												   dev_a_phi,  dev_b_phi, dev_c_phi, dev_d_phi,   
												   dev_Mx, dev_My, dev_Mz,
												   dev_indicator1, dev_indicator3, CGC_DEF,
											       delta_t, dev_Ms);
		#endif

		Kernel_Hk_field<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
			                                 dev_theta, dev_phi, 
											 dev_Hk_x, dev_Hk_y, dev_Hk_z, 
											 dev_Ku, dev_Ms);

		cudaMemcpy(Mx, dev_Mx, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(My, dev_My, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Mz, dev_Mz, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);

		#ifndef __WATCH__
		////---------Watch---------//
		////// Watch1
		//double* watch1 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch1, dev_Mx, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile1 = fopen("watch1.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile1, "%15.3e", watch1[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile1, "\n");
		//	}
		//	fprintf(pfile1, "\n\n");
		//}
		//fclose(pfile1);
		//printf("watch1\n");
		////// Watch2
		//double* watch2 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch2, dev_My, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile2 = fopen("watch2.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile2, "%15.3e", watch2[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile2, "\n");
		//	}
		//	fprintf(pfile2, "\n\n");
		//}
		//fclose(pfile2);
		//printf("watch2\n");
		////// Watch3
		//double* watch3 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch3, dev_Mz, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile3 = fopen("watch3.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile3, "%15.3e", watch3[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile3, "\n");
		//	}
		//	fprintf(pfile3, "\n\n");
		//}
		//fclose(pfile3);
		//printf("watch3\n");
		////// Watch4
		//double* watch4 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch4, dev_theta, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile4 = fopen("watch4.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile4, "%15.3e", watch4[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile4, "\n");
		//	}
		//	fprintf(pfile4, "\n\n");
		//}
		//fclose(pfile4);
		//printf("watch4\n");
		////// Watch5
		//double* watch5 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch5, dev_phi, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile5 = fopen("watch5.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile5, "%15.3e", watch5[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile5, "\n");
		//	}
		//	fprintf(pfile5, "\n\n");
		//}
		//fclose(pfile5);
		//printf("watch5\n");
		//// Watch6
		//double* watch6 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch6, dev_gamma, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile6 = fopen("watch6.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile6, "%15.3e", watch6[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile6, "\n");
		//	}
		//	fprintf(pfile6, "\n\n");
		//}
		//fclose(pfile6);
		//printf("watch6\n");
		//// Watch7
		//double* watch7 = (double*) calloc(mNx*mNy*mNz, sizeof(double));
		//cudaMemcpy(watch7, dev_T, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
		//FILE *pfile7 = fopen("watch7.out", "w");
		//for (int k = 0; k < (mNz); k++){
		//	for (int j = 0; j < (mNy); j++){
		//		for (int i = 0; i < (mNx); i++){
		//			fprintf(pfile7, "%15.3e", watch7[i+j*(mNx)+k*(mNx*mNy)]);
		//		}
		//		fprintf(pfile7, "\n");
		//	}
		//	fprintf(pfile7, "\n\n");
		//}
		//fclose(pfile7);
		//printf("watch7\n");
		////---------Watch---------//
		#endif

		//-------Calculate spacial avg of magnetization------------------------//
		if (tt >= 0)  if (!CalcMagLayerAve(tt, count, Mz_CGC, grain_coor, NumOfGrains)) printf("CalcMagLayerAverage() failed!\n");
		M_bar[tt] = pow((pow(Mx_bar1[tt], 2.0) + pow(My_bar1[tt], 2.0) + pow(Mz_bar1[tt],2.0)), 0.5);
		////---------------------------------------------------------------------//

		//------ Calculate time avg of magnetization -------------------------//
		if (tt >= EQUI_START_TIME) {
			M_bar_t = M_bar_t + M_bar[tt]/(TOTAL_TIME-EQUI_START_TIME);
			for(int idx = 0; idx < mNx*mNy*mNz; idx++){
				Mx_t_bar[idx] = Mx_t_bar[idx] + Mx[idx]/(TOTAL_TIME-EQUI_START_TIME);
				My_t_bar[idx] = My_t_bar[idx] + My[idx]/(TOTAL_TIME-EQUI_START_TIME);
				Mz_t_bar[idx] = Mz_t_bar[idx] + Mz[idx]/(TOTAL_TIME-EQUI_START_TIME);
			}
		}
		//---------------------------------------------------------------------//
		theta = acos(Mz_bar1[tt]/M_bar[tt])*360/2/PI;

		//-----------------I/O------------------------------------------------//
		if (tt>=0 && (tt%(FieldSweepTimeStep*1))==0){
			if (MH_LOOP){ 
				printf("t=%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f%8.2f\n", 
				       tt, Mz_bar1[tt], Mz_bar2[tt], Mz_bar3[tt], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl, SweepFieldStep1);
			    fprintf(wfile15, "t=%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f%8.2f\n", 
				        tt, Mz_bar1[tt], Mz_bar2[tt], Mz_bar3[tt], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl, SweepFieldStep1);
				
				// For checking ZZ's hypothesis
				/*printf("t=%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f%8.2f\n", 
				       tt, Mx[14+14*(mNx)+11*(mNx*mNy)], My[14+14*(mNx)+11*(mNx*mNy)], Mz[14+14*(mNx)+11*(mNx*mNy)], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl, SweepFieldStep1);
			    fprintf(wfile15, "t=%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f%8.2f\n", 
				        tt, Mx[14+14*(mNx)+11*(mNx*mNy)], My[14+14*(mNx)+11*(mNx*mNy)], Mz[14+14*(mNx)+11*(mNx*mNy)], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl, SweepFieldStep1);*/
			}
			if (!MH_LOOP){ 
				printf("t=%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f%5d%3d\n", 
				       tt, Mz_bar1[tt], Mz_bar2[tt], Mz_bar3[tt], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl, f_x0, *idx_f);
				fprintf(wfile15, "t=%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f%5d%3d\n", 
				        tt, Mz_bar1[tt], Mz_bar2[tt], Mz_bar3[tt], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl, f_x0, *idx_f);
			}
			if (!Output_Float_3D_Format(mNx, mNy, mNz, Happl_z, "Output_Spatial_Hz_temp.out")) { printf("OutputFloat_3D_format() failed!\n"); } //write to output file
			if (!Output_Float_3D_Format(mNx, mNy, mNz, Hd_z_1d_shift, "Output_Spatial_Hdz_temp.out")) { printf("OutputFloat_3D_format() failed!\n"); } //write to output file
			if (!Output_Float_3D_Format(mNx, mNy, mNz, Mz, "Output_Spatial_Mz_temp.out")) { printf("OutputFloat_3D_format() failed!\n"); }
			if (THERMAL) { if (!Output_Float_3D_Format(mNx, mNy, mNz, T, "Output_Spatial_T_temp.out")) { printf("OutputFloat_3D_format() failed!\n"); } }
		}

		if (MH_LOOP){ 

			#ifndef __UNIFORM_GRAIN_MEDIA__
			if (!VORO_GRAIN){
				if ((tt >= 0) && (tt%FieldSweepTimeStep == 0)){
					fprintf(wfile4, "%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%10.2f%12.2f\n", Mz_bar1[tt], Mz_bar2[tt], Mz_bar3[tt], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl);
					for (int k = 0; k < mNz; k++){
						Mz_LayerLayer[k] = 0.;
						for (int j = 0; j < mNy; j++){
							for (int i = 0; i < mNx; i++){
								int idx = i + j*mNx + k*mNx*mNy;
								Mz_LayerLayer[k] = Mz_LayerLayer[k] + Mz[idx]/(mNx*mNy);
							}
						}
					}
					itt++;
					Happl1_sweep[itt] = Happl;
					for (int i = 0; i < mNx; i++){
						for (int j = 0; j < mNy; j++){
							int idxx = i + j*mNx;
							M1_z_SingleGrain[idxx] = 0.;
								for (int k = 0; k < mNz_1; k++){  //important
								int idx = i + j*mNx + k*mNx*mNy;
								M1_z_SingleGrain[idxx] = M1_z_SingleGrain[idxx] + Mz[idx]/(mNz_1);  //important 
							}
							M1_z_SingleGrain_field[idxx][itt] = M1_z_SingleGrain[idxx];
							fprintf(wfile5, "%8.2f", M1_z_SingleGrain[idxx]);   // "Mz1_SingleGrain_Happl.out"
						}
					}
					fprintf(wfile5, "%12.2f\n", Happl);
				}

				//cluster
				for (int i=0; i<mNx; i++){
					for (int j=0; j<mNy; j++){
						Mz_tt[i+j*mNx][tt] = Mz[i+j*mNx];
					}
				}
			}
			#endif

			#ifndef __VORO_GRAIN_MEDIA__
			if (VORO_GRAIN){
				
				if ((tt >= 0) && (tt%FieldSweepTimeStep == 0)){
					
					fprintf(wfile4, "%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%12.2f\n", Mz_bar1[tt], Mz_bar2[tt], Mz_bar3[tt], Mz_bar4[tt], Mz_bar5[tt], Mz_bar6[tt], Mz_bar[tt], Happl);
					for (int k = 0; k < mNz; k++){
						Mz_LayerLayer[k] = 0.0;
						int cnt=0;;
						for (int j = 0; j < mNy; j++){
							for (int i = 0; i < mNx; i++){
								int idx = i + j*mNx + k*mNx*mNy;
								if (indicator1[idx]!=0){
									cnt++;
									Mz_LayerLayer[k] = Mz_LayerLayer[k] + Mz[idx];
								}
							}
						}
						Mz_LayerLayer[k] = Mz_LayerLayer[k]/cnt;
					}
					itt++;
					Happl1_sweep[itt] = Happl;

					// Average Mz (bottom-most) over each voronoi grian
					for (int i=0; i<grain_coor.size(); ++i){
						M1_z_SingleGrain[i] = 0.0;
						for (int ii=0; ii<grain_coor[i].size(); ++ii){
							int idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx;
							M1_z_SingleGrain[i] = M1_z_SingleGrain[i] + Mz[idx]/grain_coor[i].size();
						}
						M1_z_SingleGrain_field[i][itt] = M1_z_SingleGrain[i];
						fprintf(wfile5, "%8.2f", M1_z_SingleGrain[i]); //"Mz1_SingleGrain_Happl.out"
					}
					fprintf(wfile5, "%12.2f\n", Happl); //"Mz1_SingleGrain_Happl.out"

				}
			}
			#endif
			//--------------------------------------------------------------------//
		}
	} // Finish iteration
	curandDestroyGenerator(gen);
    
	/////// Watch /////
	if (DT_Rec_Analysis){ if (!Jitter_Calc()) { printf("Jitter_Calc() failed!\n"); }; }
	if (CT_Rec_Analysis){ if (!WPE_Calc()) { printf("Jitter_Calc() failed!\n"); }; }
	if (MH_LOOP){ if (!sigHc_Calc(M1_z_SingleGrain_field, Happl1_sweep, count, NumOfGrains, grain_coor)) { printf("sigHc_Calc() failed!\n"); }; }
		
		
		
	//-----------------I/O------------------------------------------------//
	/*for (int i = 0; i < mNx*mNy*mNz; i++){
		if((i)%(mNx) == 0 && i != 0){ 
			fprintf(wfile10, "\n"); 
			fprintf(wfile11, "\n"); 
			fprintf(wfile12, "\n"); 				
		}
		if((i)%(mNx*mNy) == 0 && i != 0){ 
			fprintf(wfile10, "\n\n"); 
			fprintf(wfile11, "\n\n");
			fprintf(wfile12, "\n\n");
		}
		fprintf(wfile10, "%12.3f", Mx_t_bar[i]);
		fprintf(wfile11, "%12.3f", My_t_bar[i]);
		fprintf(wfile12, "%12.3f", Mz_t_bar[i]);
	}*/
	//-----------------I/O------------------------------------------------//


	#ifndef __WATCH__
	/////// Watch /////
	/*printf("M_bar_t = %10.4f \n\n", M_bar_t);
	printf("T = %5.0f \n\n", T[0]);
	fprintf(wfile10, "M_bar_t = %10.4f \n\n", M_bar_t);*/
	/////// Watch /////
	#endif
	//--------------------------------------------------------------------//


	//------------Postprocessing-------------------//
	/*if (MH_LOOP){
		double *M1_Hc_gbyg = (double*) calloc(mNx*mNy, sizeof(double));
		double M1_Hc_mean = 0, M1_Hc_var = 0, M1_Hc_std = 0;
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				idxx = i + j*mNx;
				itt = 1;
				while (M1_z_SingleGrain_field[idxx][itt]*M1_z_SingleGrain_field[idxx][itt-1] > 0 && itt < ceil(TOTAL_TIME/FieldSweepTimeStep)+1){
					M1_Hc_gbyg[idxx] = Happl1_sweep[itt-1];
					itt++;
				}
			}
		}
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				if (indicator3[i+j*mNx+(mNz-1)*mNx*mNy] == 1){
					M1_Hc_mean = M1_Hc_mean + abs(M1_Hc_gbyg[i+j*mNx])/(count);
				}
			}
		}
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				if (indicator3[i+j*mNx+(mNz-1)*mNx*mNy] != 1){
					M1_Hc_var = M1_Hc_var + pow(abs(M1_Hc_gbyg[i+j*mNx])-M1_Hc_mean, 2);
				}
			}
		}
		M1_Hc_std = pow(M1_Hc_var/(count-1), 0.5);
		printf("\nM1_SFD%=%7.4f% \n", M1_Hc_std/M1_Hc_mean);
		printf("M1_Hc_mean=%9.2f \n", M1_Hc_mean);
		printf("M1_Hc_std=%8.2f \n", M1_Hc_std);
		printf("BL1=%9.5f\n\n", BL12);
		printf("BL2=%9.5f\n\n", BL23);
		printf("BL3=%9.5f\n\n", BL34);
		printf("BL4=%9.5f\n\n", BL45);
		cout << currentDateTime() << std::endl;
	
		FILE *output1;
		timeinfo = localtime(&rawtime);
		fopen_s(&output1, "output.out", "w");
		fprintf(output1,"%s", asctime(timeinfo));
		fprintf(output1,"------------------------\n");
		fprintf(output1,"M1_SFD%=%7.4f% \n", M1_Hc_std/M1_Hc_mean);
		fprintf(output1,"M1_Hc_mean=%9.2f \n", M1_Hc_mean);
		fprintf(output1,"M1_Hc_std=%8.2f \n", M1_Hc_std);
		fprintf(output1,"BL1=%9.5f\n", BL12);
		fprintf(output1,"BL2=%9.5f\n", BL23);
		fprintf(output1,"BL3=%9.5f\n", BL34);
		fprintf(output1,"BL4=%9.5f\n", BL45);
		fclose(output1);
	}*/







	stop = clock();
	time = (double) (stop-start)/CLOCKS_PER_SEC;
	printf("Run time: %f\n\n", time);
	fprintf(wfile15, "Run time: %f\n\n", time);
	printf("Done!\n\n");
	
	
	fclose(wfile3);
	fclose(wfile4);
	fclose(wfile5);
	fclose(wfile6);
	fclose(wfile7);
	fclose(wfile8);
	fclose(wfile10);
	fclose(wfile11);
	fclose(wfile12);
	fclose(wfile15);

	
	

	return 0;
}

