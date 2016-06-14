#include "Parameters.h"
#include "Parameters_input.h"
#include <cufft.h>


#ifndef __LLG_KERNEL_H__
#define __LLG_KERNEL_H__


__device__ __host__ static double GetElement(double* pData, int i, int j, int k, int pitch_x, int pitch_y)
{
	return pData[k*pitch_x*pitch_y + j*pitch_x + i];
}
__device__ __host__ static int GetElement_Int(int* pData, int i, int j, int k, int pitch_x, int pitch_y)
{
	return pData[k*pitch_x*pitch_y + j*pitch_x + i];
}
__device__ __host__ static cufftComplex GetElement_Complex(cufftComplex* pData, int i, int j, int k, int pitch_x, int pitch_y)
{
	return pData[k*pitch_x*pitch_y + j*pitch_x + i];
}
__device__ __host__ static void SetElement(double* pData, int i, int j, int k, int pitch_x, int pitch_y, double value)
{
	pData[k*pitch_x*pitch_y + j*pitch_x + i] = value;
}
__device__ __host__ static void SetElement_Int(int* pData, int i, int j, int k, int pitch_x, int pitch_y, int value)
{
	pData[k*pitch_x*pitch_y + j*pitch_x + i] = value;
}
__device__ __host__ static void SetElement_Complex(cufftComplex* pData, int i, int j, int k, int pitch_x, int pitch_y, cufftComplex value)
{
	pData[k*pitch_x*pitch_y + j*pitch_x + i] = value;
}
// Complex multiplication
__device__ __host__ static cufftComplex ComplexMul(cufftComplex* a, cufftComplex* b, int i, int j, int k, int pitch_x, int pitch_y)
{
        cufftComplex c;
        c.x = a[k*pitch_x*pitch_y + j*pitch_x + i].x * b[k*pitch_x*pitch_y + j*pitch_x + i].x - a[k*pitch_x*pitch_y + j*pitch_x + i].y * b[k*pitch_x*pitch_y + j*pitch_x + i].y;
        c.y = a[k*pitch_x*pitch_y + j*pitch_x + i].x * b[k*pitch_x*pitch_y + j*pitch_x + i].y + a[k*pitch_x*pitch_y + j*pitch_x + i].y * b[k*pitch_x*pitch_y + j*pitch_x + i].x;
        return c;
}
//__device__ __host__ static cufftComplex ComplexAdd(cufftComplex* a, cufftComplex* b, int i, int j, int k, int pitch_x, int pitch_y)
//{
//    cufftComplex c;
//    c.x = a[k*pitch_x*pitch_y + j*pitch_x + i].x + b[k*pitch_x*pitch_y + j*pitch_x + i].x;
//    c.y = a[k*pitch_x*pitch_y + j*pitch_x + i].y + b[k*pitch_x*pitch_y + j*pitch_x + i].y;
//    return c;
//}

__global__ static void Kernel_Initialization(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                         double* dev_theta, double* dev_phi,
											 double* dev_a_theta, double* dev_b_theta, double* dev_c_theta, double* dev_d_theta,
											 double* dev_a_phi, double* dev_b_phi, double* dev_c_phi, double* dev_d_phi,
											 double* dev_Ha_x, double* dev_Ha_y, double* dev_Ha_z,
											 double* dev_Hal_x, double* dev_Hal_y, double* dev_Hal_z,
											 double* dev_Hth_x, double* dev_Hth_y, double* dev_Hth_z,
											 double* dev_Hk_x, double* dev_Hk_y, double* dev_Hk_z,
											 double* dev_Hd_x, double* dev_Hd_y, double* dev_Hd_z,
											 double* dev_d_theta_d_t, double* dev_d_phi_d_t,
											 int* dev_indicator1, int* dev_indicator3/*, int AFC*/)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;


	/*if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == 3 && GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 1)
	{SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Up);}
	if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == 3 && GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 2)
	{SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);}
	if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == 3 && GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 0)
	{SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);}*/
	if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 1) SetElement(dev_theta,          x, y, z, mNx, mNy, Ini_THETA_Up);
	if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 2) SetElement(dev_theta,          x, y, z, mNx, mNy, Ini_THETA_Down);
	if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == 0) SetElement(dev_theta,          x, y, z, mNx, mNy, PI/2);
	if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == 7) SetElement(dev_theta,          x, y, z, mNx, mNy, PI/2);


	// For checking ZZ's hypothesis
	/*if (x%2 == 0) SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Up);
	else          SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);*/

	
	//if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 1)
	//{SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Up);} //Defect
	//if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 2)
	//{SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);}   //Defect


	SetElement(dev_phi,         x, y, z, mNx, mNy, 0.);
	SetElement(dev_a_theta,     x, y, z, mNx, mNy, 0.);
	SetElement(dev_b_theta,     x, y, z, mNx, mNy, 0.);
	SetElement(dev_c_theta,     x, y, z, mNx, mNy, 0.);
	SetElement(dev_d_theta,     x, y, z, mNx, mNy, 0.);
	SetElement(dev_a_phi,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_b_phi,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_c_phi,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_d_phi,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_Ha_x,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Ha_y,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Ha_z,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hal_x,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hal_y,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hal_z,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hth_x,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hth_y,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hth_z,       x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hk_x,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hk_y,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hk_z,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hd_x,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hd_y,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hd_z,        x, y, z, mNx, mNy, 0.);
	SetElement(dev_d_theta_d_t, x, y, z, mNx, mNy, 0.);
	SetElement(dev_d_phi_d_t,   x, y, z, mNx, mNy, 0.);

	__syncthreads();
}

__global__ static void Kernel_Initialization_AFC(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
												 int AF_layer_label,
	                                             double* dev_theta, int* dev_indicator1, int* dev_indicator3)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 1){
		if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == AF_layer_label){
			SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);
		}
		if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy)==GetElement_Int(dev_indicator1, x, y, z-1, mNx, mNy)  &&  GetElement_Int(dev_indicator1, x, y, z+1, mNx, mNy)==AF_layer_label){
			SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);
		}
		if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy)==GetElement_Int(dev_indicator1, x, y, z+1, mNx, mNy)  &&  GetElement_Int(dev_indicator1, x, y, z-1, mNx, mNy)==AF_layer_label){
			SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Down);
		}
	}
	if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 2){
		if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == AF_layer_label){
			SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Up);
		}
		if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy)==GetElement_Int(dev_indicator1, x, y, z-1, mNx, mNy)  &&  GetElement_Int(dev_indicator1, x, y, z+1, mNx, mNy)==AF_layer_label){
			SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Up);
		}
		if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy)==GetElement_Int(dev_indicator1, x, y, z+1, mNx, mNy)  &&  GetElement_Int(dev_indicator1, x, y, z-1, mNx, mNy)==AF_layer_label){
			SetElement(dev_theta,       x, y, z, mNx, mNy, Ini_THETA_Up);
		}
	}
	
	__syncthreads();

}
__global__ static void Kernel_Ini_indicator_temp(int BLK_SIZE_X, int BLK_SIZE_Y, int BLK_SIZE_Z, int mNx, int mNy,
	                                             int* dev_indicator1_temp, int* dev_indicator2_temp, 
												 double* dev_Ms_temp, double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z, double* dev_Aex_temp)
{
	int x = blockIdx.x * BLK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLK_SIZE_Z + threadIdx.z;
	SetElement_Int(dev_indicator1_temp, x, y, z, mNx+2, mNy+2, 0);
	SetElement_Int(dev_indicator2_temp, x, y, z, mNx+2, mNy+2, 0);
    SetElement(dev_Ms_temp, x, y, z, mNx+2, mNy+2, 0.);
	SetElement(dev_M_temp_x, x, y, z, mNx+2, mNy+2, 0.);
	SetElement(dev_M_temp_y, x, y, z, mNx+2, mNy+2, 0.);
	SetElement(dev_M_temp_z, x, y, z, mNx+2, mNy+2, 0.);
	SetElement(dev_Aex_temp, x, y, z, mNx+2, mNy+2, 0.);
	__syncthreads();
}
__global__ static void Kernel_Initialization_Cplx(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                              cufftComplex* dev_Mx_cufft,   cufftComplex* dev_My_cufft,   cufftComplex* dev_Mz_cufft,
												  cufftComplex* dev_Gxx_cufft,  cufftComplex* dev_Gxy_cufft,  cufftComplex* dev_Gxz_cufft,
												  cufftComplex* dev_Gyx_cufft,  cufftComplex* dev_Gyy_cufft,  cufftComplex* dev_Gyz_cufft,
												  cufftComplex* dev_Gzx_cufft,  cufftComplex* dev_Gzy_cufft,  cufftComplex* dev_Gzz_cufft,
												  cufftComplex* dev_Hd_x_cufft, cufftComplex* dev_Hd_y_cufft, cufftComplex* dev_Hd_z_cufft)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	cufftComplex tmp;
	tmp.x = 0;
	tmp.y = 0;

	SetElement_Complex(dev_Mx_cufft,       x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_My_cufft,       x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Mz_cufft,       x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gxx_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gxy_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gxz_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gyx_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gyy_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gyz_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gzx_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gzy_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Gzz_cufft,      x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Hd_x_cufft,     x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Hd_y_cufft,     x, y, z, 2*mNx, 2*mNy, tmp);
	SetElement_Complex(dev_Hd_z_cufft,     x, y, z, 2*mNx, 2*mNy, tmp);

	__syncthreads();
}

__global__ static void Kernel_Hth_field(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                    double* dev_GasArray, double* dev_Hth_x, double* dev_Hth_y, double* dev_Hth_z, double* dev_D)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int index = z * mNx * mNy + y * mNx + x;
	double temp_x, temp_y, temp_z;

	temp_x = dev_GasArray[index*3  ] * pow(GetElement(dev_D, x, y, z, mNx, mNy), 0.5);
	SetElement(dev_Hth_x, x, y, z, mNx, mNy, temp_x);

	temp_y = dev_GasArray[index*3+1] * pow(GetElement(dev_D, x, y, z, mNx, mNy), 0.5);
	SetElement(dev_Hth_y, x, y, z, mNx, mNy, temp_y);

	temp_z = dev_GasArray[index*3+2] * pow(GetElement(dev_D, x, y, z, mNx, mNy), 0.5);
	SetElement(dev_Hth_z, x, y, z, mNx, mNy, temp_z);

	__syncthreads();
}

__global__ static void Kernel_Hk_field(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                   double* dev_theta, double* dev_phi,
									   double* dev_Hk_x,  double* dev_Hk_y, double* dev_Hk_z,
									   double* dev_Ku, double* dev_Ms)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	double temp_z;
	double Ku_temp = GetElement(dev_Ku, x, y, z, mNx, mNy), Ms_temp = GetElement(dev_Ms, x, y, z, mNx, mNy);

	temp_z = 2 * Ku_temp / Ms_temp * cos(GetElement(dev_theta, x, y, z, mNx, mNy));
	SetElement(dev_Hk_x, x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hk_y, x, y, z, mNx, mNy, 0.);
	SetElement(dev_Hk_z, x, y, z, mNx, mNy, temp_z);

	__syncthreads();
}

__global__ static void Kernel_initiate_dev_M_temp(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy, int mNz,
	                                              double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z, 
												  int* dev_indicator1_temp, int* dev_indicator2_temp)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
		
		SetElement(dev_M_temp_x, x, y, z, mNx+2, mNy+2, 0.);
		SetElement(dev_M_temp_y, x, y, z, mNx+2, mNy+2, 0.);
		SetElement(dev_M_temp_z, x, y, z, mNx+2, mNy+2, 0.);

		SetElement_Int(dev_indicator1_temp, x, y, z, mNx+2, mNy+2, 0);
		SetElement_Int(dev_indicator2_temp, x, y, z, mNx+2, mNy+2, 0);
		__syncthreads();
		
		if (x == mNx-1) {
			SetElement(dev_M_temp_x, x+1, y, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+2, y, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+1, y, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+2, y, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+1, y, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+2, y, z, mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x+2, y, z, mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x+2, y, z, mNx+2, mNy+2, 0);
		}
		if (y == mNy-1) {
			SetElement(dev_M_temp_x, x, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y+2, z, mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x, y+1, z, mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x, y+2, z, mNx+2, mNy+2, 0);
		}
		if (z == mNz-1) {
			SetElement(dev_M_temp_x, x, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y, z+2, mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x, y, z+1, mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x, y, z+2, mNx+2, mNy+2, 0);
		}
		if (x == mNx-1 && y == mNy-1) {
			SetElement(dev_M_temp_x, x+1, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+2, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+1, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+2, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+1, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+2, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+1, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+2, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+1, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+2, y+1, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+1, y+2, z, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+2, y+2, z, mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x+1, y+2, z, mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x+2, y+2, z, mNx+2, mNy+2, 0);
		}
		if (x == mNx-1 && z == mNz-1) {
			SetElement(dev_M_temp_x, x+1, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+2, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+1, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x+2, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+1, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+2, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+1, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+2, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+1, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+2, y, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+1, y, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+2, y, z+2, mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x+1, y, z+2, mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x+2, y, z+2, mNx+2, mNy+2, 0);
		}
		if (y == mNy-1 && z == mNz-1) {
			SetElement(dev_M_temp_x, x, y+1, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x, y+2, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x, y+1, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_x, x, y+2, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y+1, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y+2, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y+1, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x, y+2, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y+1, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y+2, z+1, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y+1, z+2, mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x, y+2, z+2, mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x, y+1, z+2, mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x, y+2, z+2, mNx+2, mNy+2, 0);
		}
		__syncthreads();
		if (x < 2 && y < 2 && z < 2) {
			SetElement(dev_M_temp_x, x+mNx,   y+mNy,   z+mNz,   mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_y, x+mNx,   y+mNy,   z+mNz,   mNx+2, mNy+2, 0.);
			SetElement(dev_M_temp_z, x+mNx,   y+mNy,   z+mNz,   mNx+2, mNy+2, 0.);

			SetElement_Int(dev_indicator1_temp, x+mNx,   y+mNy,   z+mNz,   mNx+2, mNy+2, 0);
			SetElement_Int(dev_indicator2_temp, x+mNx,   y+mNy,   z+mNz,   mNx+2, mNy+2, 0);
		}
		__syncthreads();
	



	/*int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
		
		SetElement(dev_M_temp_x, x, y, z, mNx+2, mNy+2, 0.);
		SetElement(dev_M_temp_y, x, y, z, mNx+2, mNy+2, 0.);
		SetElement(dev_M_temp_z, x, y, z, mNx+2, mNy+2, 0.);
		SetElement_Int(dev_indicator1_temp, x, y, z, mNx+2, mNy+2, 0);
		SetElement_Int(dev_indicator2_temp, x, y, z, mNx+2, mNy+2, 0);

		__syncthreads();*/

}

__global__ static void Kernel_dev_indicator_with_apron(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy, int mNz,
	                                                   int* dev_indicator1, int* dev_indicator1_temp, 
												       int* dev_indicator2, int* dev_indicator2_temp)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int temp1, temp2;
	
	temp1 = GetElement_Int(dev_indicator1, x, y, z, mNx, mNy);
	temp2 = GetElement_Int(dev_indicator2, x, y, z, mNx, mNy);
	__syncthreads();

	// Inner cube
	SetElement_Int(dev_indicator1_temp, x+1, y+1, z+1, mNx+2, mNy+2, temp1);
	SetElement_Int(dev_indicator2_temp, x+1, y+1, z+1, mNx+2, mNy+2, temp2);
	__syncthreads();

	// up & down faces
	if (z == mNz-1) {
		SetElement_Int(dev_indicator1_temp, x+1, y+1, 0, mNx+2, mNy+2, 0);
		SetElement_Int(dev_indicator2_temp, x+1, y+1, 0, mNx+2, mNy+2, 0);
	}
	if (z == 0) {
		SetElement_Int(dev_indicator1_temp, x+1, y+1, mNz+1, mNx+2, mNy+2, 0);
		SetElement_Int(dev_indicator2_temp, x+1, y+1, mNz+1, mNx+2, mNy+2, 0);
	}
	__syncthreads();

	// left and right faces
	if (x == mNx-1) {
		SetElement_Int(dev_indicator1_temp, 0, y+1, z+1, mNx+2, mNy+2, temp1);
		SetElement_Int(dev_indicator2_temp, 0, y+1, z+1, mNx+2, mNy+2, temp2);
	}
	if (x == 0) {
		SetElement_Int(dev_indicator1_temp, mNx+1, y+1, z+1, mNx+2, mNy+2, temp1);
		SetElement_Int(dev_indicator2_temp, mNx+1, y+1, z+1, mNx+2, mNy+2, temp2);
	}
	__syncthreads();

	// front and back faces
	if (y == mNy-1) {
		SetElement_Int(dev_indicator1_temp, x+1, 0, z+1, mNx+2, mNy+2, temp1);
		SetElement_Int(dev_indicator2_temp, x+1, 0, z+1, mNx+2, mNy+2, temp2);
	}
	if (y == 0) {
		SetElement_Int(dev_indicator1_temp, x+1, mNy+1, z+1, mNx+2, mNy+2, temp1);
		SetElement_Int(dev_indicator2_temp, x+1, mNy+1, z+1, mNx+2, mNy+2, temp2);
	}
	__syncthreads();

}
__global__ static void Kernel_M_A_with_apron(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy, int mNz,
	                                         double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z, 
										     double* dev_Ms_temp, double* dev_Aex_temp,
										     double* dev_theta,    double* dev_x_theta,  double* dev_phi, double* dev_x_phi,
										     double* dev_Ms, double* dev_Aex, double hh, int CGC_DEF, int* dev_indicator3)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	double temp_x, temp_y, temp_z;
	
	double Ms_temp = GetElement(dev_Ms, x, y, z, mNx, mNy), Aex_temp = fabs(GetElement(dev_Aex, x, y, z, mNx, mNy));


	temp_x = Ms_temp * sin(GetElement(dev_theta, x, y, z, mNx, mNy) + hh * GetElement(dev_x_theta, x, y, z, mNx, mNy)) * 
		               cos(GetElement(dev_phi,   x, y, z, mNx, mNy) + hh * GetElement(dev_x_phi,   x, y, z, mNx, mNy));
	temp_y = Ms_temp * sin(GetElement(dev_theta, x, y, z, mNx, mNy) + hh * GetElement(dev_x_theta, x, y, z, mNx, mNy)) * 
		               sin(GetElement(dev_phi,   x, y, z, mNx, mNy) + hh * GetElement(dev_x_phi,   x, y, z, mNx, mNy));
	temp_z = Ms_temp * cos(GetElement(dev_theta, x, y, z, mNx, mNy) + hh * GetElement(dev_x_theta, x, y, z, mNx, mNy));

	__syncthreads();
	

	// inner cube of M_temp is given values
	SetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2, temp_x);
	SetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2, temp_y);
	SetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2, temp_z);
	SetElement(dev_Ms_temp,  x+1, y+1, z+1, mNx+2, mNy+2, Ms_temp);
	SetElement(dev_Aex_temp, x+1, y+1, z+1, mNx+2, mNy+2, Aex_temp);
	
	__syncthreads();


	// up & down faces are given values
	if (z == mNz-1) {
	//	temp_x = GetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_x, x+1, y+1, 0, mNx+2, mNy+2, 1/pow(3.,0.5));
	//	temp_y = GetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_y, x+1, y+1, 0, mNx+2, mNy+2, 1/pow(3.,0.5));
	//	temp_z = GetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_z, x+1, y+1, 0, mNx+2, mNy+2, 1/pow(3.,0.5));
		SetElement(dev_Ms_temp, x+1, y+1, 0, mNx+2, mNy+2, 1.);
		SetElement(dev_Aex_temp, x+1, y+1, 0, mNx+2, mNy+2, 0.);
	}
	if (z == 0) {
	//	temp_x = GetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_x, x+1, y+1, mNz+1, mNx+2, mNy+2, 1/pow(3.,0.5));
	//	temp_y = GetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_y, x+1, y+1, mNz+1, mNx+2, mNy+2, 1/pow(3.,0.5));
	//	temp_z = GetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_z, x+1, y+1, mNz+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		SetElement(dev_Ms_temp, x+1, y+1, mNz+1, mNx+2, mNy+2, 1.);
		SetElement(dev_Aex_temp, x+1, y+1, mNz+1, mNx+2, mNy+2, 0.);
	}
	__syncthreads();

	// left and right faces are given values
	if (x == mNx-1) {
		/*temp_x = GetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_x, 0, y+1, z+1, mNx+2, mNy+2, temp_x);*/
		SetElement(dev_M_temp_x, 0, y+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_y = GetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_y, 0, y+1, z+1, mNx+2, mNy+2, temp_y);*/
		SetElement(dev_M_temp_y, 0, y+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_z = GetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_z, 0, y+1, z+1, mNx+2, mNy+2, temp_z);*/
		SetElement(dev_M_temp_z, 0, y+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		SetElement(dev_Ms_temp, 0, y+1, z+1, mNx+2, mNy+2, /*GetElement(dev_Ms_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/1.);
		SetElement(dev_Aex_temp, 0, y+1, z+1, mNx+2, mNy+2, /*GetElement(dev_Aex_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/0.);
	}
	if (x == 0) {
		/*temp_x = GetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_x, mNx+1, y+1, z+1, mNx+2, mNy+2, temp_x);*/
		SetElement(dev_M_temp_x, mNx+1, y+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_y = GetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_y, mNx+1, y+1, z+1, mNx+2, mNy+2, temp_y);*/
		SetElement(dev_M_temp_y, mNx+1, y+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_z = GetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_z, mNx+1, y+1, z+1, mNx+2, mNy+2, temp_z);*/
		SetElement(dev_M_temp_z, mNx+1, y+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		SetElement(dev_Ms_temp, mNx+1, y+1, z+1, mNx+2, mNy+2, /*GetElement(dev_Ms_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/1.);
		SetElement(dev_Aex_temp, mNx+1, y+1, z+1, mNx+2, mNy+2, /*GetElement(dev_Aex_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/0.);
	}
	__syncthreads();

	// front and back faces are given values
	if (y == mNy-1) {
		/*temp_x = GetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_x, x+1, 0, z+1, mNx+2, mNy+2, temp_x);*/
		SetElement(dev_M_temp_x, x+1, 0, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_y = GetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_y, x+1, 0, z+1, mNx+2, mNy+2, temp_y);*/
		SetElement(dev_M_temp_y, x+1, 0, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_z = GetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_z, x+1, 0, z+1, mNx+2, mNy+2, temp_z);*/
		SetElement(dev_M_temp_z, x+1, 0, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		SetElement(dev_Ms_temp, x+1, 0, z+1, mNx+2, mNy+2, /*GetElement(dev_Ms_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/1.);
		SetElement(dev_Aex_temp, x+1, 0, z+1, mNx+2, mNy+2, /*GetElement(dev_Aex_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/0.);
	}
	if (y == 0) {
		/*temp_x = GetElement(dev_M_temp_x, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_x, x+1, mNy+1, z+1, mNx+2, mNy+2, temp_x);*/
		SetElement(dev_M_temp_x, x+1, mNy+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_y = GetElement(dev_M_temp_y, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_y, x+1, mNy+1, z+1, mNx+2, mNy+2, temp_y);*/
		SetElement(dev_M_temp_y, x+1, mNy+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		/*temp_z = GetElement(dev_M_temp_z, x+1, y+1, z+1, mNx+2, mNy+2);
		SetElement(dev_M_temp_z, x+1, mNy+1, z+1, mNx+2, mNy+2, temp_z);*/
		SetElement(dev_M_temp_z, x+1, mNy+1, z+1, mNx+2, mNy+2, 1/pow(3.,0.5));
		SetElement(dev_Ms_temp, x+1, mNy+1, z+1, mNx+2, mNy+2, /*GetElement(dev_Ms_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/1.);
		SetElement(dev_Aex_temp, x+1, mNy+1, z+1, mNx+2, mNy+2, /*GetElement(dev_Aex_temp, x+1, y+1, z+1, mNx+2, mNy+2)*/0.);
	}
	__syncthreads();

}

// For obtaining exchange field Ha
__global__ static void Kernel_Ha_with_apron(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                        double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z, 
											double* dev_Ms_temp, double* dev_Aex_temp,
											double* dev_Aex_XP, double* dev_Aex_XM, double* dev_Aex_YP, double* dev_Aex_YM, double* dev_Aex_ZP, double* dev_Aex_ZM, 
											double* dev_Ms_XP, double* dev_Ms_XM, double* dev_Ms_YP, double* dev_Ms_YM, double* dev_Ms_ZP, double* dev_Ms_ZM,
											int* dev_indicator1_temp, int* dev_indicator2_temp, 
						                    double* dev_Ms, double* dev_Aex,
											float L1_Hex_l, float L2_Hex_l, float L3_Hex_l, float L4_Hex_l, float L5_Hex_l, float L6_Hex_l,
											float BL12_Hex_l, float BL23_Hex_l, float BL34_Hex_l, float BL45_Hex_l, float BL56_Hex_l, int AFC, int AF_layer_label)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int i_tmp = x + 1, 
		j_tmp = y + 1, 
		k_tmp = z + 1;
	double Ms = GetElement(dev_Ms, x, y, z, mNx, mNy), 
		   Aex = fabs(GetElement(dev_Aex, x, y, z, mNx, mNy));
	double Ms_xp = GetElement(dev_Ms_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2), Aex_xp = GetElement(dev_Aex_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2),
		   Ms_xm = GetElement(dev_Ms_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2), Aex_xm = GetElement(dev_Aex_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2),
		   Ms_yp = GetElement(dev_Ms_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2), Aex_yp = GetElement(dev_Aex_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2),
		   Ms_ym = GetElement(dev_Ms_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2), Aex_ym = GetElement(dev_Aex_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2),
		   Ms_zp = GetElement(dev_Ms_temp, i_tmp, j_tmp, k_tmp+1, mNx+2, mNy+2), Aex_zp = GetElement(dev_Aex_temp, i_tmp, j_tmp, k_tmp+1, mNx+2, mNy+2),
		   Ms_zm = GetElement(dev_Ms_temp, i_tmp, j_tmp, k_tmp-1, mNx+2, mNy+2), Aex_zm = GetElement(dev_Aex_temp, i_tmp, j_tmp, k_tmp-1, mNx+2, mNy+2);

	/* Write to Ha field. */
	double Ms_XP, Ms_XM, 
		   Ms_YP, Ms_YM,
		   Ms_ZP, Ms_ZM;
	double Aex_XP, Aex_XM,   // Grain boundary exchange (x/y/z's plus/minus)
		   Aex_YP, Aex_YM,   // Grain boundary exchange (x/y/z's plus/minus)
		   Aex_ZP, Aex_ZM;   // Grain boundary exchange (x/y/z's plus/minus)

	#ifdef __6_LAYER__
	
	// Break the exchange when confronting grain boundaries
	// Geometic means of exchange stiffness constant and saturation magnetization are obtained 
	if (GetElement_Int(dev_indicator1_temp, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2) != 0){
		switch (GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2)){
			case 1:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) !=
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)L1_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)L1_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)L1_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)L1_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 12:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)BL12_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)BL12_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)BL12_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)BL12_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 2:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)L2_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)L2_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)L2_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)L2_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 23:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)BL23_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)BL23_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)BL23_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)BL23_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 3:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)L3_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)L3_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)L3_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)L3_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 34:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)BL34_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)BL34_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)BL34_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)BL34_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 4:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)L4_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)L4_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)L4_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)L4_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 45:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)BL45_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)BL45_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)BL45_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)BL45_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 5:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)L5_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)L5_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)L5_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)L5_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 56:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)BL56_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)BL56_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)BL56_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)BL56_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
			case 6:
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XP = ((double)L6_Hex_l)*pow(Aex*Aex_xp,0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				else {
					Aex_XP = pow(Aex*Aex_xp, 0.5);
					Ms_XP  = pow(Ms*Ms_xp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)){
					Aex_XM = ((double)L6_Hex_l)*pow(Aex*Aex_xm,0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				else {
					Aex_XM = pow(Aex*Aex_xm, 0.5);
					Ms_XM  = pow(Ms*Ms_xm, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)){
					Aex_YP = ((double)L6_Hex_l)*pow(Aex*Aex_yp,0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				else {
					Aex_YP = pow(Aex*Aex_yp, 0.5);
					Ms_YP  = pow(Ms*Ms_yp, 0.5);
				}
				if (GetElement_Int(dev_indicator2_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) != 
					GetElement_Int(dev_indicator2_temp, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)){
					Aex_YM = ((double)L6_Hex_l)*pow(Aex*Aex_ym,0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				else {
					Aex_YM = pow(Aex*Aex_ym, 0.5);
					Ms_YM  = pow(Ms*Ms_ym, 0.5);
				}
				break;
		}

		if (AFC){
			if ((GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2)==int(AF_layer_label*10+(AF_layer_label+1)))        && 
				(GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp+1, mNx+2, mNy+2)==int(AF_layer_label*10+(AF_layer_label+1)))){
				Aex_ZP = -pow(Aex*Aex_zp, 0.5);
				Aex_ZM =  pow(Aex*Aex_zm, 0.5);
			}
			else if ((GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2)==int(AF_layer_label*10+(AF_layer_label+1)))     && 
				     (GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp-1, mNx+2, mNy+2)==int(AF_layer_label*10+(AF_layer_label+1)))){
				Aex_ZM = -pow(Aex*Aex_zm, 0.5);
				Aex_ZP = pow(Aex*Aex_zp, 0.5);
			}
			else if ((GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2) == int((AF_layer_label-1)*10+AF_layer_label))    && 
				     (GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp+1, mNx+2, mNy+2) == int((AF_layer_label-1)*10+AF_layer_label))){
				Aex_ZP = -pow(Aex*Aex_zp, 0.5);
				Aex_ZM = pow(Aex*Aex_zm, 0.5);
			}
			else if ((GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp, mNx+2, mNy+2)==int((AF_layer_label-1)*10+AF_layer_label))    && 
				     (GetElement_Int(dev_indicator1_temp, i_tmp,   j_tmp, k_tmp-1, mNx+2, mNy+2)==int((AF_layer_label-1)*10+AF_layer_label))){
				Aex_ZM = -pow(Aex*Aex_zm, 0.5);
				Aex_ZP = pow(Aex*Aex_zp, 0.5);
			}
			else {
				Aex_ZP = pow(Aex*Aex_zp, 0.5);
				Aex_ZM = pow(Aex*Aex_zm, 0.5);
			}
		}
		else {
			Aex_ZP = pow(Aex*Aex_zp, 0.5);
			Aex_ZM = pow(Aex*Aex_zm, 0.5);
		}
		Ms_ZP =  pow(Ms*Ms_zp, 0.5);
		Ms_ZM =  pow(Ms*Ms_zm, 0.5);
	}
	// If the cell happens to be grain boundary
	else
	{
		double Aex_bnd = 0.0;
		double Ms_bnd = 1.0;
		Aex_XP = Aex_bnd;
		Aex_XM = Aex_bnd;
		Aex_YP = Aex_bnd;
		Aex_YM = Aex_bnd;
		Aex_ZP = Aex_bnd;
		Aex_ZM = Aex_bnd;
		Ms_XP = Ms_bnd;
		Ms_XM = Ms_bnd;
		Ms_YP = Ms_bnd;
		Ms_YM = Ms_bnd;
		Ms_ZP = Ms_bnd;
		Ms_ZM = Ms_bnd;
	}
	#endif

	__syncthreads();
	
	// Set calculated Aex_XP to dev_Aex_XP
	SetElement(dev_Aex_XP, x, y, z, mNx, mNy, Aex_XP);
	SetElement(dev_Aex_XM, x, y, z, mNx, mNy, Aex_XM);
	SetElement(dev_Aex_YP, x, y, z, mNx, mNy, Aex_YP);
	SetElement(dev_Aex_YM, x, y, z, mNx, mNy, Aex_YM);
	SetElement(dev_Aex_ZP, x, y, z, mNx, mNy, Aex_ZP);
	SetElement(dev_Aex_ZM, x, y, z, mNx, mNy, Aex_ZM);
	
	// Set calculated Ms_XP to dev_Ms_XP 
	SetElement(dev_Ms_XP, x, y, z, mNx, mNy, Ms_XP);
	SetElement(dev_Ms_XM, x, y, z, mNx, mNy, Ms_XM);
	SetElement(dev_Ms_YP, x, y, z, mNx, mNy, Ms_YP);
	SetElement(dev_Ms_YM, x, y, z, mNx, mNy, Ms_YM);
	SetElement(dev_Ms_ZP, x, y, z, mNx, mNy, Ms_ZP);
	SetElement(dev_Ms_ZM, x, y, z, mNx, mNy, Ms_ZM);
	
	__syncthreads();
	
}


__global__ static void Kernel_Ha_with_Left(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
										   double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z,
										   double* dev_H1_x,     double* dev_H1_y,     double* dev_H1_z,
										   double delta_x,       double* dev_Aex_XM,   double* dev_Ms_XM, int* dev_indicator1)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int i_tmp = x + 1, 
		j_tmp = y + 1, 
		k_tmp = z + 1;
	
	double phitt, Ms_subject, Ms_adjacent, Hn, Hn1, Hn2, Hnx, Hny, Hnz;   //zyliu
	double H1_x, H1_y, H1_z;
	double sstt1, sstt2, sstt3, aatt1, aatt2, aatt3;            //zyliu

	double Ms_XM = GetElement(dev_Ms_XM, x, y, z, mNx, mNy);
	double Aex_XM = GetElement(dev_Aex_XM, x, y, z, mNx, mNy);

	Ms_subject=pow((pow(GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
			        pow(GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
		            pow(GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject == 0.0){
		sstt1 = 0.0;
		sstt2 = 0.0;
		sstt3 = 0.0;
	}
	else {
		sstt1 = GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt2 = GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt3 = GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
	}
	////////////////////////////////********i-1********//////////////////////////////////////////////////////////////////////////////
    Ms_adjacent = pow((pow(GetElement(dev_M_temp_x, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
					   pow(GetElement(dev_M_temp_y, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
					   pow(GetElement(dev_M_temp_z, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject*Ms_adjacent == 0.0){
		H1_x = 0.0;
		H1_y = 0.0;
		H1_z = 0.0;
	}
	else {
		phitt = (GetElement(dev_M_temp_x, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
			     GetElement(dev_M_temp_y, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
			     GetElement(dev_M_temp_z, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2))/(Ms_subject*Ms_adjacent);
	
		if(phitt > 0.999999999999) phitt = 0.999999999999;
		if(phitt < -0.999999999999) phitt = -0.999999999999;
		phitt = acos(phitt);

		if (fabs(phitt) > 2.00*PI) phitt = fmod(phitt,2*PI);
	
		if (phitt < -1.00*PI) phitt = 2*PI+phitt;
		if (phitt > PI) phitt = 2*PI-phitt;

		if (phitt == 0.0){
			H1_x = 2*Aex_XM* GetElement(dev_M_temp_x, i_tmp-1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XM,2.0)/pow(delta_x, 2.0);
			H1_y = 2*Aex_XM* GetElement(dev_M_temp_y, i_tmp-1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XM,2.0)/pow(delta_x, 2.0);
			H1_z = 2*Aex_XM* GetElement(dev_M_temp_z, i_tmp-1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XM,2.0)/pow(delta_x, 2.0);
		}
		else {
			Hnx = 2*Aex_XM* GetElement(dev_M_temp_x, i_tmp-1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XM,2.0)/pow(delta_x, 2.0);
			Hny = 2*Aex_XM* GetElement(dev_M_temp_y, i_tmp-1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XM,2.0)/pow(delta_x, 2.0);
			Hnz = 2*Aex_XM* GetElement(dev_M_temp_z, i_tmp-1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XM,2.0)/pow(delta_x, 2.0);
			Hn = pow((pow(Hnx, 2.0) + pow(Hny, 2.0) + pow(Hnz, 2.0)),0.5)*pow((1+pow(phitt, 2.0)), 0.5);
			
			if ( (phitt-atan(phitt)) < PI/2.0 ){
				Hn1 = Hn*tan(phitt-atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(phitt-atan(phitt));
			}
			else if ( (phitt-atan(phitt)) > PI/2.0 ){
				Hn2 = Hn*tan(atan(phitt))/(tan(PI-phitt+atan(phitt))-tan(atan(phitt)))/cos(PI-phitt+atan(phitt));
				Hn1 = Hn2*sin(PI-phitt+atan(phitt))/sin(atan(phitt));
			}
			else {
				Hn1 = Hn/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt));
			}
			if (Ms_adjacent == 0.0){
				aatt1 = 0.0;
				aatt2 = 0.0;
				aatt3 = 0.0;
			}
			else {
				aatt1 = GetElement(dev_M_temp_x, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt2 = GetElement(dev_M_temp_y, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt3 = GetElement(dev_M_temp_z, i_tmp-1, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
			}
			H1_x = Hn1*sstt1+Hn2*aatt1;
			H1_y = Hn1*sstt2+Hn2*aatt2;
			H1_z = Hn1*sstt3+Hn2*aatt3;
		}
	}

	SetElement(dev_H1_x, x, y, z, mNx, mNy, H1_x);
	SetElement(dev_H1_y, x, y, z, mNx, mNy, H1_y);
	SetElement(dev_H1_z, x, y, z, mNx, mNy, H1_z);
	
	__syncthreads();
}

__global__ static void Kernel_Ha_with_right(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
											double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z,
											double* dev_H2_x,     double* dev_H2_y,     double* dev_H2_z,
											double delta_x,       double* dev_Aex_XP,   double* dev_Ms_XP, int* dev_indicator1)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int i_tmp = x + 1, 
		j_tmp = y + 1, 
		k_tmp = z + 1;
	
	double phitt, Ms_subject, Ms_adjacent, Hn, Hn1, Hn2, Hnx, Hny, Hnz;   //zyliu
	double H2_x, H2_y, H2_z;
	double sstt1, sstt2, sstt3, aatt1, aatt2, aatt3;            //zyliu

	double Ms_XP = GetElement(dev_Ms_XP, x, y, z, mNx, mNy);
	double Aex_XP = GetElement(dev_Aex_XP, x, y, z, mNx, mNy);

	Ms_subject = pow((pow(GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
			          pow(GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
		              pow(GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject == 0.0){
		sstt1 = 0.0;
		sstt2 = 0.0;
		sstt3 = 0.0;
	}
	else {
		sstt1 = GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt2 = GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt3 = GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
	}
	////////////////////////////////********i+1********//////////////////////////////////////////////////////////////////////////////
    Ms_adjacent = pow((pow(GetElement(dev_M_temp_x, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
					 pow(GetElement(dev_M_temp_y, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
					 pow(GetElement(dev_M_temp_z, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject*Ms_adjacent*Ms_XP == 0.0){
		H2_x = 0.0;
		H2_y = 0.0;
		H2_z = 0.0;
	}
	else {
		phitt = (GetElement(dev_M_temp_x, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
			     GetElement(dev_M_temp_y, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
			     GetElement(dev_M_temp_z, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2))/(Ms_subject*Ms_adjacent);
	
		if (phitt > 0.999999999999) phitt = 0.999999999999;
		if (phitt < -0.999999999999) phitt = -0.999999999999;
		phitt = acos(phitt);

		if (fabs(phitt) > 2.0*PI) phitt = fmod(phitt,2*PI);
	
		if (phitt < -1.0*PI) phitt = 2*PI+phitt;
		if (phitt > PI) phitt = 2*PI-phitt;

		if(phitt == 0){
			H2_x = 2*Aex_XP* GetElement(dev_M_temp_x, i_tmp+1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XP,2.0)/pow(delta_x, 2.0);
			H2_y = 2*Aex_XP* GetElement(dev_M_temp_y, i_tmp+1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XP,2.0)/pow(delta_x, 2.0);
			H2_z = 2*Aex_XP* GetElement(dev_M_temp_z, i_tmp+1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XP,2.0)/pow(delta_x, 2.0);
		}
		else{
			Hnx = 2*Aex_XP* GetElement(dev_M_temp_x, i_tmp+1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XP,2.0)/pow(delta_x, 2.0);
			Hny = 2*Aex_XP* GetElement(dev_M_temp_y, i_tmp+1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XP,2.0)/pow(delta_x, 2.0);
			Hnz = 2*Aex_XP* GetElement(dev_M_temp_z, i_tmp+1, j_tmp,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_XP,2.0)/pow(delta_x, 2.0);
			Hn = pow((pow(Hnx, 2.0)+pow(Hny, 2.0)+pow(Hnz, 2.0)),0.5)*pow((1+pow(phitt, 2.0)), 0.5);
			
			if ( (phitt-atan(phitt)) < PI/2.0 ){
				Hn1 = Hn*tan(phitt-atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(phitt-atan(phitt));
			}
			else if ( (phitt-atan(phitt)) > PI/2.0 ){
				Hn2 = Hn*tan(atan(phitt))/(tan(PI-phitt+atan(phitt))-tan(atan(phitt)))/cos(PI-phitt+atan(phitt));
				Hn1 = Hn2*sin(PI-phitt+atan(phitt))/sin(atan(phitt));
			}
			else {
				Hn1 = Hn/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt));
			}
			if (Ms_adjacent == 0.0){
				aatt1 = 0.0;
				aatt2 = 0.0;
				aatt3 = 0.0;
			}
			else {
				aatt1 = GetElement(dev_M_temp_x, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt2 = GetElement(dev_M_temp_y, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt3 = GetElement(dev_M_temp_z, i_tmp+1, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
			}
			H2_x = Hn1*sstt1+Hn2*aatt1;
			H2_y = Hn1*sstt2+Hn2*aatt2;
			H2_z = Hn1*sstt3+Hn2*aatt3;
		}
	}

	SetElement(dev_H2_x, x, y, z, mNx, mNy, H2_x);
	SetElement(dev_H2_y, x, y, z, mNx, mNy, H2_y);
	SetElement(dev_H2_z, x, y, z, mNx, mNy, H2_z);
	
	__syncthreads();
}
__global__ static void Kernel_Ha_with_up(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
									     double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z,
									     double* dev_H3_x,     double* dev_H3_y,     double* dev_H3_z,
									     double delta_y,       double* dev_Aex_YM,   double* dev_Ms_YM, int* dev_indicator1)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int i_tmp = x + 1, 
		j_tmp = y + 1, 
		k_tmp = z + 1;
	
	double phitt, Ms_subject, Ms_adjacent, Hn, Hn1, Hn2, Hnx, Hny, Hnz;   //zyliu
	double H3_x, H3_y, H3_z;
	double sstt1, sstt2, sstt3, aatt1, aatt2, aatt3;            //zyliu

	double Ms_YM = GetElement(dev_Ms_YM, x, y, z, mNx, mNy);
	double Aex_YM = GetElement(dev_Aex_YM, x, y, z, mNx, mNy);

	Ms_subject = pow((pow(GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
			          pow(GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
		              pow(GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject == 0.0){
		sstt1 = 0.0;
		sstt2 = 0.0;
		sstt3 = 0.0;
	}
	else {
		sstt1 = GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt2 = GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt3 = GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
	}
	////////////////////////////////********j-1********//////////////////////////////////////////////////////////////////////////////
    Ms_adjacent = pow((pow(GetElement(dev_M_temp_x, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2), 2.0)+
					   pow(GetElement(dev_M_temp_y, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2), 2.0)+
					   pow(GetElement(dev_M_temp_z, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject*Ms_adjacent*Ms_YM == 0.0){
		H3_x = 0.0;
		H3_y = 0.0;
		H3_z = 0.0;
	}
	else {
		phitt = (GetElement(dev_M_temp_x, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
				 GetElement(dev_M_temp_y, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
				 GetElement(dev_M_temp_z, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2))/(Ms_subject*Ms_adjacent);
	
		if (phitt > 0.999999999999) phitt = 0.999999999999;
		if (phitt < -0.999999999999) phitt = -0.999999999999;
		phitt = acos(phitt);

		if (fabs(phitt) > 2.0*PI) phitt = fmod(phitt,2*PI);
	
		if (phitt < -1.0*PI) phitt = 2*PI+phitt;
		if (phitt > PI) phitt = 2*PI-phitt;

		if (phitt == 0){
			H3_x = 2*Aex_YM* GetElement(dev_M_temp_x, i_tmp, j_tmp-1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YM,2.0)/pow(delta_y, 2.0);
			H3_y = 2*Aex_YM* GetElement(dev_M_temp_y, i_tmp, j_tmp-1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YM,2.0)/pow(delta_y, 2.0);
			H3_z = 2*Aex_YM* GetElement(dev_M_temp_z, i_tmp, j_tmp-1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YM,2.0)/pow(delta_y, 2.0);
		}
		else {
			Hnx = 2*Aex_YM* GetElement(dev_M_temp_x, i_tmp, j_tmp-1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YM,2.0)/pow(delta_y, 2.0);
			Hny = 2*Aex_YM* GetElement(dev_M_temp_y, i_tmp, j_tmp-1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YM,2.0)/pow(delta_y, 2.0);
			Hnz = 2*Aex_YM* GetElement(dev_M_temp_z, i_tmp, j_tmp-1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YM,2.0)/pow(delta_y, 2.0);
			Hn = pow((pow(Hnx, 2.0)+pow(Hny, 2.0)+pow(Hnz, 2.0)),0.5)*pow((1+pow(phitt, 2.0)), 0.5);
			
			if ( (phitt-atan(phitt)) < PI/2.0 ){
				Hn1 = Hn*tan(phitt-atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(phitt-atan(phitt));
			}
			else if ( (phitt-atan(phitt)) > PI/2.0 ){
				Hn2 = Hn*tan(atan(phitt))/(tan(PI-phitt+atan(phitt))-tan(atan(phitt)))/cos(PI-phitt+atan(phitt));
				Hn1 = Hn2*sin(PI-phitt+atan(phitt))/sin(atan(phitt));
			}
			else {
				Hn1 = Hn/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt));
			}
		
			if (Ms_adjacent == 0.0){
				aatt1 = 0.0;
				aatt2 = 0.0;
				aatt3 = 0.0;
			}
			else {
				aatt1 = GetElement(dev_M_temp_x, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt2 = GetElement(dev_M_temp_y, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt3 = GetElement(dev_M_temp_z, i_tmp, j_tmp-1, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
			}
			H3_x = Hn1*sstt1+Hn2*aatt1;
			H3_y = Hn1*sstt2+Hn2*aatt2;
			H3_z = Hn1*sstt3+Hn2*aatt3;
		}
	}

	SetElement(dev_H3_x, x, y, z, mNx, mNy, H3_x);
	SetElement(dev_H3_y, x, y, z, mNx, mNy, H3_y);
	SetElement(dev_H3_z, x, y, z, mNx, mNy, H3_z);
	
	__syncthreads();
}
__global__ static void Kernel_Ha_with_down(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	 									   double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z,
										   double* dev_H4_x,     double* dev_H4_y,     double* dev_H4_z,
										   double delta_y,       double* dev_Aex_YP,   double* dev_Ms_YP, int* dev_indicator1)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int i_tmp = x + 1, 
		j_tmp = y + 1, 
		k_tmp = z + 1;
	
	double phitt, Ms_subject, Ms_adjacent, Hn, Hn1, Hn2, Hnx, Hny, Hnz;   //zyliu
	double H4_x, H4_y, H4_z;
	double sstt1, sstt2, sstt3, aatt1, aatt2, aatt3;            //zyliu

	double Ms_YP = GetElement(dev_Ms_YP, x, y, z, mNx, mNy);
	double Aex_YP = GetElement(dev_Aex_YP, x, y, z, mNx, mNy);

	Ms_subject = pow((pow(GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
			          pow(GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)+
		              pow(GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2), 2.0)),0.5);

	if (Ms_subject == 0.0){
		sstt1 = 0.0;
		sstt2 = 0.0;
		sstt3 = 0.0;
	}
	else {
		sstt1 = GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt2 = GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
		sstt3 = GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)/Ms_subject;
	}
	////////////////////////////////********j+1********//////////////////////////////////////////////////////////////////////////////
    Ms_adjacent = pow((pow(GetElement(dev_M_temp_x, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2), 2.0)+
					   pow(GetElement(dev_M_temp_y, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2), 2.0)+
					   pow(GetElement(dev_M_temp_z, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2), 2.0)),0.5);
	
	if (Ms_subject*Ms_adjacent*Ms_YP == 0.0){
		H4_x = 0.0;
		H4_y = 0.0;
		H4_z = 0.0;
	}
	else {
		phitt = (GetElement(dev_M_temp_x, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
				 GetElement(dev_M_temp_y, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2)+
				 GetElement(dev_M_temp_z, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)*GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp, mNx+2, mNy+2))/(Ms_subject*Ms_adjacent);
	
		if (phitt > 0.999999999999) phitt = 0.999999999999;
		if (phitt < -0.999999999999) phitt = -0.999999999999;
		phitt = acos(phitt);

		if (fabs(phitt) > 2.0*PI) phitt = fmod(phitt,2*PI);
	
		if (phitt < -1.0*PI) phitt = 2*PI+phitt;
		if (phitt > PI) phitt = 2*PI-phitt;

		if (phitt == 0){
			H4_x = 2*Aex_YP* GetElement(dev_M_temp_x, i_tmp, j_tmp+1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YP,2.0)/pow(delta_y, 2.0);
			H4_y = 2*Aex_YP* GetElement(dev_M_temp_y, i_tmp, j_tmp+1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YP,2.0)/pow(delta_y, 2.0);
			H4_z = 2*Aex_YP* GetElement(dev_M_temp_z, i_tmp, j_tmp+1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YP,2.0)/pow(delta_y, 2.0);
		}
		else {
			Hnx = 2*Aex_YP* GetElement(dev_M_temp_x, i_tmp, j_tmp+1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YP,2.0)/pow(delta_y, 2.0);
			Hny = 2*Aex_YP* GetElement(dev_M_temp_y, i_tmp, j_tmp+1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YP,2.0)/pow(delta_y, 2.0);
			Hnz = 2*Aex_YP* GetElement(dev_M_temp_z, i_tmp, j_tmp+1,   k_tmp,   mNx+2, mNy+2)/ pow(Ms_YP,2.0)/pow(delta_y, 2.0);
			Hn = pow((pow(Hnx, 2.0)+pow(Hny, 2.0)+pow(Hnz, 2.0)),0.5)*pow((1+pow(phitt, 2.0)), 0.5);
			
			if ( (phitt-atan(phitt)) < PI/2.0 ){
				Hn1 = Hn*tan(phitt-atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt))/(tan(atan(phitt))+tan(phitt-atan(phitt)))/cos(phitt-atan(phitt));
			}
			else if( (phitt-atan(phitt)) > PI/2.0 ){
				Hn2 = Hn*tan(atan(phitt))/(tan(PI-phitt+atan(phitt))-tan(atan(phitt)))/cos(PI-phitt+atan(phitt));
				Hn1 = Hn2*sin(PI-phitt+atan(phitt))/sin(atan(phitt));
			}
			else {
				Hn1 = Hn/cos(atan(phitt));
				Hn2 = Hn*tan(atan(phitt));
			}
			if (Ms_adjacent == 0.0){
				aatt1 = 0.0;
				aatt2 = 0.0;
				aatt3 = 0.0;
			}
			else {
				aatt1 = GetElement(dev_M_temp_x, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt2 = GetElement(dev_M_temp_y, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
				aatt3 = GetElement(dev_M_temp_z, i_tmp, j_tmp+1, k_tmp, mNx+2, mNy+2)/Ms_adjacent;
			}
		
			H4_x = Hn1*sstt1+Hn2*aatt1;
			H4_y = Hn1*sstt2+Hn2*aatt2;
			H4_z = Hn1*sstt3+Hn2*aatt3;
		}
	}
	
	SetElement(dev_H4_x, x, y, z, mNx, mNy, H4_x);
	SetElement(dev_H4_y, x, y, z, mNx, mNy, H4_y);
	SetElement(dev_H4_z, x, y, z, mNx, mNy, H4_z);
	__syncthreads();
	
}
__global__ static void Kernel_Ha_with_att(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
										  double* dev_M_temp_x, double* dev_M_temp_y, double* dev_M_temp_z,
										  double* dev_H1_x, double* dev_H1_y, double* dev_H1_z,
										  double* dev_H2_x, double* dev_H2_y, double* dev_H2_z,
										  double* dev_H3_x, double* dev_H3_y, double* dev_H3_z,
										  double* dev_H4_x, double* dev_H4_y, double* dev_H4_z, 
										  double* dev_Ha_x,     double* dev_Ha_y,     double* dev_Ha_z,
										  double* dev_Hal_x,    double* dev_Hal_y,     double* dev_Hal_z,
										  double delta_z, double* dev_Aex_ZP, double* dev_Aex_ZM, double* dev_Ms_ZP, double* dev_Ms_ZM)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	int i_tmp = x + 1, 
		j_tmp = y + 1, 
		k_tmp = z + 1;

	double temp_x,  temp_y,  temp_z;
	double temp_xx, temp_yy, temp_zz;
	double Ms_ZP_sqr, Ms_ZM_sqr;

	double Ms_ZP  = GetElement(dev_Ms_ZP, x, y, z, mNx, mNy);
	double Ms_ZM  = GetElement(dev_Ms_ZM, x, y, z, mNx, mNy);
	double Aex_ZP = GetElement(dev_Aex_ZP, x, y, z, mNx, mNy);
	double Aex_ZM = GetElement(dev_Aex_ZM, x, y, z, mNx, mNy);

	Ms_ZP_sqr = pow(Ms_ZP, 2.0);
	Ms_ZM_sqr = pow(Ms_ZM, 2.0);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	temp_x = GetElement(dev_H1_x, x, y, z, mNx, mNy)+ 
		     GetElement(dev_H2_x, x, y, z, mNx, mNy)+ 
			 GetElement(dev_H3_x, x, y, z, mNx, mNy)+
		     GetElement(dev_H4_x, x, y, z, mNx, mNy)+ 
			 2*((Aex_ZP* GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp+1, mNx+2, mNy+2)/ Ms_ZP_sqr + 
			     Aex_ZM* GetElement(dev_M_temp_x, i_tmp, j_tmp, k_tmp-1, mNx+2, mNy+2)/ Ms_ZM_sqr)/pow(delta_z, 2.0));
	temp_xx = GetElement(dev_H1_x, x, y, z, mNx, mNy)+
			  GetElement(dev_H2_x, x, y, z, mNx, mNy)+
			  GetElement(dev_H3_x, x, y, z, mNx, mNy)+
			  GetElement(dev_H4_x, x, y, z, mNx, mNy);

	SetElement(dev_Ha_x, x, y, z, mNx, mNy, temp_x);
	SetElement(dev_Hal_x, x, y, z, mNx, mNy, temp_xx);

	temp_y = GetElement(dev_H1_y, x, y, z, mNx, mNy)+
			 GetElement(dev_H2_y, x, y, z, mNx, mNy)+
			 GetElement(dev_H3_y, x, y, z, mNx, mNy)+
			 GetElement(dev_H4_y, x, y, z, mNx, mNy)+ 
			 2*((Aex_ZP* GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp+1, mNx+2, mNy+2)/ Ms_ZP_sqr + 
			     Aex_ZM* GetElement(dev_M_temp_y, i_tmp, j_tmp, k_tmp-1, mNx+2, mNy+2)/ Ms_ZM_sqr)/pow(delta_z, 2.0));
	temp_yy = GetElement(dev_H1_y, x, y, z, mNx, mNy)+
			  GetElement(dev_H2_y, x, y, z, mNx, mNy)+
			  GetElement(dev_H3_y, x, y, z, mNx, mNy)+
			  GetElement(dev_H4_y, x, y, z, mNx, mNy);

	SetElement(dev_Ha_y, x, y, z, mNx, mNy, temp_y);
	SetElement(dev_Hal_y, x, y, z, mNx, mNy, temp_yy);

	temp_z =  GetElement(dev_H1_z, x, y, z, mNx, mNy)+
			  GetElement(dev_H2_z, x, y, z, mNx, mNy)+
			  GetElement(dev_H3_z, x, y, z, mNx, mNy)+
			  GetElement(dev_H4_z, x, y, z, mNx, mNy)+ 
			  2*((Aex_ZP* GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp+1, mNx+2, mNy+2)/ Ms_ZP_sqr + 
				  Aex_ZM* GetElement(dev_M_temp_z, i_tmp, j_tmp, k_tmp-1, mNx+2, mNy+2)/ Ms_ZM_sqr)/pow(delta_z, 2.0));
	temp_zz = GetElement(dev_H1_z, x, y, z, mNx, mNy)+
			  GetElement(dev_H2_z, x, y, z, mNx, mNy)+
			  GetElement(dev_H3_z, x, y, z, mNx, mNy)+
		      GetElement(dev_H4_z, x, y, z, mNx, mNy);

	SetElement(dev_Ha_z, x, y, z, mNx, mNy, temp_z);
	SetElement(dev_Hal_z, x, y, z, mNx, mNy, temp_zz);
	
	__syncthreads();
}

__device__ static double dtheta_dt(int mNx, int mNy,
	                               int i, int j, int k, 
								   double* dev_theta,	double* dev_phi, 
								   double* dev_x_theta, double* dev_x_phi,
								   double* dev_Ha_x,	double* dev_Ha_y, double* dev_Ha_z,
								   double* dev_Hth_x,	double* dev_Hth_y, double* dev_Hth_z,
								   double* dev_Hd_x,	double* dev_Hd_y, double* dev_Hd_z,
								   double* dev_Happl_x, double* dev_Happl_y, double* dev_Happl_z,	
								   double hh, double* dev_Ku, double* dev_Ms, double* dev_alpha, double* dev_gamma)
{
	double theta_tmp = GetElement(dev_theta, i, j, k, mNx, mNy);
	double phi_tmp   = GetElement(dev_phi,   i, j, k, mNx, mNy);
	double Ms = GetElement(dev_Ms, i, j, k, mNx, mNy), 
		   Ku = GetElement(dev_Ku, i, j, k, mNx, mNy), 
		   alpha = GetElement(dev_alpha, i, j, k, mNx, mNy),
		   gamma = GetElement(dev_gamma, i, j, k, mNx, mNy);

	if (fabs(theta_tmp) > 2.0 * PI)
		theta_tmp = fmod(fabs(theta_tmp), 2.0 * PI);
	else if (fabs(theta_tmp) > PI && fabs(theta_tmp) <= 2.0 * PI)
		theta_tmp = 2.0 * PI - fabs(theta_tmp);
	else
		theta_tmp = fabs(theta_tmp);

	phi_tmp = phi_tmp - (int)((phi_tmp/2.0/PI) + 0.5) * 2.0 * PI;


	double theta_add_dtheta = theta_tmp + hh * GetElement(dev_x_theta, i, j, k, mNx, mNy);
	double phi_add_dphi     = phi_tmp   + hh * GetElement(dev_x_phi  , i, j, k, mNx, mNy);

	return  gamma * (alpha * cos(theta_add_dtheta) * cos(phi_add_dphi) - sin(phi_add_dphi))*
		    (GetElement(dev_Ha_x, i, j, k, mNx, mNy) + GetElement(dev_Hth_x, i, j, k, mNx, mNy) + GetElement(dev_Hd_x, i, j, k, mNx, mNy) + GetElement(dev_Happl_x, i, j, k, mNx, mNy)) + 
		    gamma * (alpha * cos(theta_add_dtheta) * sin(phi_add_dphi) + cos(phi_add_dphi))*
		    (GetElement(dev_Ha_y, i, j, k, mNx, mNy) + GetElement(dev_Hth_y, i, j, k, mNx, mNy) + GetElement(dev_Hd_y, i, j, k, mNx, mNy) + GetElement(dev_Happl_y, i, j, k, mNx, mNy)) -
	   	    gamma * alpha * sin(theta_add_dtheta) * (2*Ku/Ms*cos(theta_add_dtheta) + GetElement(dev_Ha_z, i, j, k, mNx, mNy) + GetElement(dev_Hth_z, i, j, k, mNx, mNy) + GetElement(dev_Hd_z, i, j, k, mNx, mNy) + GetElement(dev_Happl_z, i, j, k, mNx, mNy));   

	__syncthreads();
}

__device__ static double dphi_dt(int mNx, int mNy,
	                             int i, int j, int k, 
								 double* dev_theta, double* dev_phi, 
								 double* dev_x_theta, double* dev_x_phi,
								 double* dev_Ha_x, double* dev_Ha_y, double* dev_Ha_z,
								 double* dev_Hth_x, double* dev_Hth_y, double* dev_Hth_z,
								 double* dev_Hd_x, double* dev_Hd_y, double* dev_Hd_z,
								 double* dev_Happl_x, double* dev_Happl_y, double* dev_Happl_z,
								 double hh, double* dev_Ku, double* dev_Ms, double* dev_alpha, double* dev_gamma)
{
	double theta_tmp = GetElement(dev_theta, i, j, k, mNx, mNy);
	double phi_tmp   = GetElement(dev_phi,   i, j, k, mNx, mNy);
	double Ms = GetElement(dev_Ms, i, j, k, mNx, mNy), 
		   Ku = GetElement(dev_Ku, i, j, k, mNx, mNy), 
		   alpha = GetElement(dev_alpha, i, j, k, mNx, mNy),
		   gamma = GetElement(dev_gamma, i, j, k, mNx, mNy);

	if (fabs(theta_tmp) > 2.0 * PI)
		theta_tmp = fmod(fabs(theta_tmp), 2.0 * PI);
	else if (fabs(theta_tmp) > PI && fabs(theta_tmp) <= 2.0 * PI)
		theta_tmp = 2.0 * PI - fabs(theta_tmp);
	else
		theta_tmp = fabs(theta_tmp);

	phi_tmp = phi_tmp - (int)((phi_tmp/2.0/PI) + 0.5) * 2.0 * PI;


	double theta_add_dtheta = theta_tmp + hh * GetElement(dev_x_theta, i, j, k, mNx, mNy);
	double phi_add_dphi     = phi_tmp   + hh * GetElement(dev_x_phi  , i, j, k, mNx, mNy);

	if (theta_add_dtheta < 0.000001) { theta_add_dtheta = 0.000001; }
	if (theta_add_dtheta > PI-0.000001) { theta_add_dtheta = PI-0.000001; }
	return  gamma * (-alpha*sin(phi_add_dphi)/sin(theta_add_dtheta) - cos(theta_add_dtheta)/sin(theta_add_dtheta)*cos(phi_add_dphi)) * 
		            (GetElement(dev_Ha_x, i, j, k, mNx, mNy) + GetElement(dev_Hth_x, i, j, k, mNx, mNy) + GetElement(dev_Hd_x, i, j, k, mNx, mNy) + GetElement(dev_Happl_x, i, j, k, mNx, mNy)) + 
  		    gamma * ( alpha*cos(phi_add_dphi)/sin(theta_add_dtheta) - cos(theta_add_dtheta)/sin(theta_add_dtheta)*sin(phi_add_dphi)) *
	   	            (GetElement(dev_Ha_y, i, j, k, mNx, mNy) + GetElement(dev_Hth_y, i, j, k, mNx, mNy) + GetElement(dev_Hd_y, i, j, k, mNx, mNy) + GetElement(dev_Happl_y, i, j, k, mNx, mNy)) +
 		    gamma * (2*Ku/Ms*cos(theta_add_dtheta) + GetElement(dev_Ha_z, i, j, k, mNx, mNy) + GetElement(dev_Hth_z, i, j, k, mNx, mNy) + GetElement(dev_Hd_z, i, j, k, mNx, mNy) + GetElement(dev_Happl_z, i, j, k, mNx, mNy));  

	__syncthreads();
}


__global__ static void Kernel_d_theta_phi_d_t(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                          double* dev_d_theta_d_t, double* dev_d_phi_d_t,
											  double* dev_theta,       double* dev_phi, 
											  double* dev_a_theta,     double* dev_a_phi,
											  double* dev_Ha_x,        double* dev_Ha_y,    double* dev_Ha_z,
									          double* dev_Hth_x,       double* dev_Hth_y,   double* dev_Hth_z,
											  double* dev_Hd_x,        double* dev_Hd_y,    double* dev_Hd_z,
											  double* dev_Happl_x,     double* dev_Happl_y, double* dev_Happl_z,
									          double* dev_Ku,          double* dev_Ms,      double* dev_alpha,  	double* dev_gamma)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	SetElement(dev_d_theta_d_t, x, y, z, mNx, mNy, dtheta_dt(mNx, mNy,
		                                                     x, y, z, 
		                                                     dev_theta, dev_phi, 
													 	     dev_a_theta, dev_a_phi,
														     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                         dev_Hth_x, dev_Hth_y, dev_Hth_z,
														     dev_Hd_x, dev_Hd_y, dev_Hd_z,
														     dev_Happl_x, dev_Happl_y, dev_Happl_z,
														     ZERO, dev_Ku, dev_Ms, dev_alpha, dev_gamma));
	SetElement(dev_d_phi_d_t,   x, y, z, mNx, mNy, dphi_dt  (mNx, mNy,
		                                                     x, y, z, 
		                                                     dev_theta, dev_phi, 
														     dev_a_theta, dev_a_phi,
														     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                         dev_Hth_x, dev_Hth_y, dev_Hth_z,
														     dev_Hd_x, dev_Hd_y, dev_Hd_z,
														     dev_Happl_x, dev_Happl_y, dev_Happl_z,
														     ZERO, dev_Ku, dev_Ms, dev_alpha, dev_gamma));

	__syncthreads();
}


__global__ static void Kernel_a_theta_phi(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                      double* dev_a_theta, double* dev_a_phi,
								          double* dev_d_theta_d_t, double* dev_d_phi_d_t)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	SetElement(dev_a_theta, x, y, z, mNx, mNy, GetElement(dev_d_theta_d_t, x, y, z, mNx, mNy));
	SetElement(dev_a_phi,   x, y, z, mNx, mNy, GetElement(dev_d_phi_d_t,   x, y, z, mNx, mNy));

	__syncthreads();
}


__global__ static void Kernel_b_theta_phi(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                      double* dev_b_theta, double* dev_b_phi,
	                                      double* dev_theta, double* dev_phi, 
										  double* dev_a_theta, double* dev_a_phi,
										  double* dev_Ha_x, double* dev_Ha_y, double* dev_Ha_z,
									      double* dev_Hth_x, double* dev_Hth_y, double* dev_Hth_z,
										  double* dev_Hd_x, double* dev_Hd_y, double* dev_Hd_z,
										  double* dev_Happl_x, double* dev_Happl_y, double* dev_Happl_z,
									      double* dev_Ku, double* dev_Ms, double* dev_alpha, double* dev_gamma, double h)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	SetElement(dev_b_theta, x, y, z, mNx, mNy, dtheta_dt(mNx, mNy,
		                                                 x, y, z, 
		                                                 dev_theta, dev_phi, 
													     dev_a_theta, dev_a_phi,
													     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                     dev_Hth_x, dev_Hth_y, dev_Hth_z,
													     dev_Hd_x, dev_Hd_y, dev_Hd_z,
													     dev_Happl_x, dev_Happl_y, dev_Happl_z,
													     h/2, dev_Ku, dev_Ms, dev_alpha, dev_gamma));
	SetElement(dev_b_phi,   x, y, z, mNx, mNy, dphi_dt  (mNx, mNy,
		                                                 x, y, z, 
		                                                 dev_theta, dev_phi, 
													     dev_a_theta, dev_a_phi,
													     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                     dev_Hth_x, dev_Hth_y, dev_Hth_z,
													     dev_Hd_x, dev_Hd_y, dev_Hd_z,
													     dev_Happl_x, dev_Happl_y, dev_Happl_z,
													     h/2, dev_Ku, dev_Ms, dev_alpha, dev_gamma));

	__syncthreads();
}

__global__ static void Kernel_c_theta_phi(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                      double* dev_c_theta, double* dev_c_phi,
	                                      double* dev_theta, double* dev_phi, 
										  double* dev_b_theta, double* dev_b_phi,
										  double* dev_Ha_x, double* dev_Ha_y, double* dev_Ha_z,
									      double* dev_Hth_x, double* dev_Hth_y, double* dev_Hth_z,
										  double* dev_Hd_x, double* dev_Hd_y, double* dev_Hd_z,
										  double* dev_Happl_x, double* dev_Happl_y, double* dev_Happl_z,
									      double* dev_Ku, double* dev_Ms, double* dev_alpha, double* dev_gamma, double h)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	SetElement(dev_c_theta, x, y, z, mNx, mNy, dtheta_dt(mNx, mNy,
		                                                 x, y, z, 
		                                                 dev_theta, dev_phi, 
													     dev_b_theta, dev_b_phi,
													     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                     dev_Hth_x, dev_Hth_y, dev_Hth_z,
													     dev_Hd_x, dev_Hd_y, dev_Hd_z,
													     dev_Happl_x, dev_Happl_y, dev_Happl_z,
													     h/2, dev_Ku, dev_Ms, dev_alpha, dev_gamma));
	SetElement(dev_c_phi,   x, y, z, mNx, mNy, dphi_dt  (mNx, mNy,
		                                                 x, y, z, 
		                                                 dev_theta, dev_phi, 
													     dev_b_theta, dev_b_phi,
													     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                     dev_Hth_x, dev_Hth_y, dev_Hth_z,
													     dev_Hd_x, dev_Hd_y, dev_Hd_z,
													     dev_Happl_x, dev_Happl_y, dev_Happl_z,
													     h/2, dev_Ku, dev_Ms, dev_alpha, dev_gamma));

	__syncthreads();
}

__global__ static void Kernel_d_theta_phi(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                      double* dev_d_theta, double* dev_d_phi,
   									      double* dev_theta, double* dev_phi, 
										  double* dev_c_theta, double* dev_c_phi,
										  double* dev_Ha_x, double* dev_Ha_y, double* dev_Ha_z,
									      double* dev_Hth_x, double* dev_Hth_y, double* dev_Hth_z,
										  double* dev_Hd_x, double* dev_Hd_y, double* dev_Hd_z,
										  double* dev_Happl_x, double* dev_Happl_y, double* dev_Happl_z,
									      double* dev_Ku, double* dev_Ms, double* dev_alpha, double* dev_gamma, double h)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	SetElement(dev_d_theta, x, y, z, mNx, mNy, dtheta_dt(mNx, mNy,
		                                                 x, y, z, 
		                                                 dev_theta, dev_phi, 
													     dev_c_theta, dev_c_phi,
													     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                     dev_Hth_x, dev_Hth_y, dev_Hth_z,
												         dev_Hd_x, dev_Hd_y, dev_Hd_z,
													     dev_Happl_x, dev_Happl_y, dev_Happl_z,
													     h, dev_Ku, dev_Ms, dev_alpha, dev_gamma));
	SetElement(dev_d_phi,   x, y, z, mNx, mNy, dphi_dt  (mNx, mNy,
		                                                 x, y, z, 
		                                                 dev_theta, dev_phi, 
													     dev_c_theta, dev_c_phi,
													     dev_Ha_x, dev_Ha_y, dev_Ha_z,
									                     dev_Hth_x, dev_Hth_y, dev_Hth_z,
													     dev_Hd_x, dev_Hd_y, dev_Hd_z,
													     dev_Happl_x, dev_Happl_y, dev_Happl_z,
													     h, dev_Ku, dev_Ms, dev_alpha, dev_gamma));

	__syncthreads();
}


__global__ static void Kernel_time_increment(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy,
	                                         double* dev_theta, double* dev_phi,
											 double* dev_a_theta, double* dev_b_theta, double* dev_c_theta, double* dev_d_theta,   
											 double* dev_a_phi, double* dev_b_phi, double* dev_c_phi, double* dev_d_phi,   
											 double* dev_Mx, double* dev_My, double* dev_Mz,
											 int* dev_indicator1, int* dev_indicator3, int CGC_DEF,
											 double h, double* dev_Ms)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;
	
	double Ms_temp = GetElement(dev_Ms, x, y, z, mNx, mNy);

	// Rescaling theta and then carry out time increment for updating dev_theta
	double temp;
	temp = GetElement(dev_theta, x, y, z, mNx, mNy);
	if (fabs(temp) > 2.0*PI) SetElement(dev_theta, x, y, z, mNx, mNy, fmod(fabs(temp), 2.0*PI));
	else if (fabs(temp) > PI && fabs(temp) <= 2.0*PI) SetElement(dev_theta, x, y, z, mNx, mNy, 2.0*PI - fabs(temp));
	else SetElement(dev_theta, x, y, z, mNx, mNy, fabs(temp));
	
	temp = GetElement(dev_theta, x, y, z, mNx, mNy) + h/6.0*(GetElement(dev_a_theta, x, y, z, mNx, mNy)   + 
													         GetElement(dev_b_theta, x, y, z, mNx, mNy)*2.0 + 
													         GetElement(dev_c_theta, x, y, z, mNx, mNy)*2.0 +
													         GetElement(dev_d_theta, x, y, z, mNx, mNy));
	SetElement(dev_theta, x, y, z, mNx, mNy, temp);

	// Carry out time increment for updating dev_phi
	temp = GetElement(dev_phi, x, y, z, mNx, mNy);
	SetElement(dev_phi, x, y, z, mNx, mNy, temp - ((int)(temp/2.0/PI+0.5))*2.0*PI);
	
	temp = GetElement(dev_phi, x, y, z, mNx, mNy) + h/6.0*(GetElement(dev_a_phi, x, y, z, mNx, mNy)   + 
		                                                   GetElement(dev_b_phi, x, y, z, mNx, mNy)*2.0 + 
					                                       GetElement(dev_c_phi, x, y, z, mNx, mNy)*2.0 +
					                                       GetElement(dev_d_phi, x, y, z, mNx, mNy));
	SetElement(dev_phi, x, y, z, mNx, mNy, temp);
    
	// Rescaling theta and phi
	temp = GetElement(dev_theta, x, y, z, mNx, mNy);
	if (fabs(temp) > 2.0*PI)				
		SetElement(dev_theta, x, y, z, mNx, mNy, fmod(fabs(temp), 2.0*PI));
	else if (fabs(temp) > PI && fabs(temp) <= 2.0*PI)
		SetElement(dev_theta, x, y, z, mNx, mNy, 2.0*PI - fabs(temp));
	else
		SetElement(dev_theta, x, y, z, mNx, mNy, fabs(temp));

	temp = GetElement(dev_phi, x, y, z, mNx, mNy);
	SetElement(dev_phi, x, y, z, mNx, mNy, temp - ((int)(temp/2.0/PI+0.5))*2.0*PI);
	
	// Set magnetization on grain boundary to zero
	if (GetElement_Int(dev_indicator1, x, y, z, mNx, mNy) == 0){
		SetElement(dev_Mx, x, y, z, mNx, mNy, 0.0);
		SetElement(dev_My, x, y, z, mNx, mNy, 0.0);
		SetElement(dev_Mz, x, y, z, mNx, mNy, 0.0);
	}
	else {
		SetElement(dev_Mx, x, y, z, mNx, mNy, Ms_temp*sin(GetElement(dev_theta, x, y, z, mNx, mNy))*cos(GetElement(dev_phi, x, y, z, mNx, mNy)));
		SetElement(dev_My, x, y, z, mNx, mNy, Ms_temp*sin(GetElement(dev_theta, x, y, z, mNx, mNy))*sin(GetElement(dev_phi, x, y, z, mNx, mNy)));
		SetElement(dev_Mz, x, y, z, mNx, mNy, Ms_temp*cos(GetElement(dev_theta, x, y, z, mNx, mNy)));
	}

	// Defects
	/*if (GetElement_Int(dev_indicator3, x, y, z, mNx, mNy) == 3 && CGC_DEF == 2)
	{SetElement(dev_theta, x, y, z, mNx, mNy, Ini_THETA_Down);}*/

	__syncthreads();
}

__global__ static void Kernel_CUFFT_M_times_G(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int BLK_SZ_Z, int lx_zero_pad, int ly_zero_pad,
	                                          cufftComplex* dev_Gxx_cufft, cufftComplex* dev_Gxy_cufft,   cufftComplex* dev_Gxz_cufft,
											  cufftComplex* dev_Gyx_cufft,  cufftComplex* dev_Gyy_cufft,  cufftComplex* dev_Gyz_cufft,
											  cufftComplex* dev_Gzx_cufft,  cufftComplex* dev_Gzy_cufft,  cufftComplex* dev_Gzz_cufft,
											  cufftComplex* dev_Mx_cufft,   cufftComplex* dev_My_cufft,   cufftComplex* dev_Mz_cufft,
											  cufftComplex* dev_Hd_x_cufft, cufftComplex* dev_Hd_y_cufft, cufftComplex* dev_Hd_z_cufft)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLK_SZ_Z     + threadIdx.z;
	cufftComplex adder1, adder2, adder3,
                 Hd_x, Hd_y, Hd_z;

	adder1 = ComplexMul(dev_Gxx_cufft, dev_Mx_cufft, x, y, z, lx_zero_pad, ly_zero_pad); 
	adder2 = ComplexMul(dev_Gxy_cufft, dev_My_cufft, x, y, z, lx_zero_pad, ly_zero_pad);
	adder3 = ComplexMul(dev_Gxz_cufft, dev_Mz_cufft, x, y, z, lx_zero_pad, ly_zero_pad);
	__syncthreads();
	Hd_x.x = adder1.x + adder2.x + adder3.x;
	Hd_x.y = adder1.y + adder2.y + adder3.y;
	
	adder1 = ComplexMul(dev_Gyx_cufft, dev_Mx_cufft, x, y, z, lx_zero_pad, ly_zero_pad); 
	adder2 = ComplexMul(dev_Gyy_cufft, dev_My_cufft, x, y, z, lx_zero_pad, ly_zero_pad);
	adder3 = ComplexMul(dev_Gyz_cufft, dev_Mz_cufft, x, y, z, lx_zero_pad, ly_zero_pad);
	__syncthreads();
	Hd_y.x = adder1.x + adder2.x + adder3.x;
	Hd_y.y = adder1.y + adder2.y + adder3.y;
	
	adder1 = ComplexMul(dev_Gzx_cufft, dev_Mx_cufft, x, y, z, lx_zero_pad, ly_zero_pad); 
	adder2 = ComplexMul(dev_Gzy_cufft, dev_My_cufft, x, y, z, lx_zero_pad, ly_zero_pad);
	adder3 = ComplexMul(dev_Gzz_cufft, dev_Mz_cufft, x, y, z, lx_zero_pad, ly_zero_pad);
	__syncthreads();
	Hd_z.x = adder1.x + adder2.x + adder3.x;
	Hd_z.y = adder1.y + adder2.y + adder3.y;
	
	SetElement_Complex(dev_Hd_x_cufft, x, y, z, lx_zero_pad, ly_zero_pad, Hd_x);
	SetElement_Complex(dev_Hd_y_cufft, x, y, z, lx_zero_pad, ly_zero_pad, Hd_y);
	SetElement_Complex(dev_Hd_z_cufft, x, y, z, lx_zero_pad, ly_zero_pad, Hd_z);

	__syncthreads();
}

__global__ static void Kernel_Hint_field(int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int BLOCK_SIZE_Z, int mNx, int mNy, 
										 double* dev_Hint_x, double* dev_Hint_y, double* dev_Hint_z, 
										 double* dev_Hal_x, double* dev_Hal_y, double* dev_Hal_z,
										 double* dev_Hd_x, double* dev_Hd_y, double* dev_Hd_z)
{
	int x = blockIdx.x * BLOCK_SIZE_X + threadIdx.x,
	    y = blockIdx.y * BLOCK_SIZE_Y + threadIdx.y,
		z = blockIdx.z * BLOCK_SIZE_Z + threadIdx.z;

	SetElement(dev_Hint_x, x, y, z, mNx, mNy, GetElement(dev_Hal_x, x, y, z, mNx, mNy) + GetElement(dev_Hd_x, x, y, z, mNx, mNy));
	SetElement(dev_Hint_y, x, y, z, mNx, mNy, GetElement(dev_Hal_y, x, y, z, mNx, mNy) + GetElement(dev_Hd_y, x, y, z, mNx, mNy));
	SetElement(dev_Hint_z, x, y, z, mNx, mNy, GetElement(dev_Hal_z, x, y, z, mNx, mNy) + GetElement(dev_Hd_z, x, y, z, mNx, mNy));
	
	__syncthreads();
}
#endif