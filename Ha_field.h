//#include "cutil.h"
#include "Parameters.h"
#include "Parameters_input.h"
#include "LLG_kernel.cu"


bool Ha_field(double hh, double* dev_x_theta, double* dev_x_phi, dim3 grids1, dim3 blocks1)
{
	
	Kernel_M_A_with_apron<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy, mNz,
		                                       dev_M_temp_x, dev_M_temp_y, dev_M_temp_z, 
										       dev_Ms_temp, dev_Aex_temp,
									   	       dev_theta, dev_x_theta, dev_phi, dev_x_phi,
									  	       dev_Ms, dev_Aex, hh, CGC_DEF, dev_indicator3);

	Kernel_Ha_with_apron<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
		                                      dev_M_temp_x, dev_M_temp_y, dev_M_temp_z,
											  dev_Ms_temp, dev_Aex_temp,
											  dev_Aex_XP, dev_Aex_XM, dev_Aex_YP, dev_Aex_YM, dev_Aex_ZP, dev_Aex_ZM, 
											  dev_Ms_XP, dev_Ms_XM, dev_Ms_YP, dev_Ms_YM, dev_Ms_ZP, dev_Ms_ZM,
											  dev_indicator1_temp, dev_indicator2_temp,
											  dev_Ms, dev_Aex,
											  L1_Hex_l, L2_Hex_l, L3_Hex_l, L4_Hex_l, L5_Hex_l, L6_Hex_l,
											  BL12_Hex_l, BL23_Hex_l, BL34_Hex_l, BL45_Hex_l, BL56_Hex_l, AFC, AF_layer_label);

	/////// Watch /////
	/*double* watch1 = (double*) calloc((mNx+2)*(mNy+2)*(mNz+2), sizeof(double));
	cudaMemcpy(watch1, dev_Ms_temp, (mNx+2)*(mNy+2)*(mNz+2)*sizeof(double), cudaMemcpyDeviceToHost);
	FILE *pfile1 = fopen("watch1.out", "w");
	for (int k = 0; k < (mNz+2); k++){
		for (int j = 0; j < (mNy+2); j++){
			for (int i = 0; i < (mNx+2); i++){
				fprintf(pfile1, "%15.6e", watch1[i+j*(mNx+2)+k*(mNx+2)*(mNy+2)]);
			}
			fprintf(pfile1, "\n");
		}
		fprintf(pfile1, "\n\n");
	}
	fclose(pfile1);
	printf("watch1\n");*/
	
	/*double* watch2 = (double*) calloc((mNx)*(mNy)*(mNz), sizeof(double));
	cudaMemcpy(watch2, dev_theta, mNx*mNy*mNz*sizeof(double), cudaMemcpyDeviceToHost);
	FILE *pfile2 = fopen("watch2.out", "w");
	for (int k = 0; k < (mNz); k++){
		for (int j = 0; j < (mNy); j++){
			for (int i = 0; i < (mNx); i++){
				fprintf(pfile2, "%15.6e", watch2[i+j*(mNx)+k*(mNx)*(mNy)]);
			}
			fprintf(pfile2, "\n");
		}
		fprintf(pfile2, "\n\n");
	}
	fclose(pfile2);
	printf("watch2\n");*/
	/////// Watch /////

	Kernel_Ha_with_Left<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
											 dev_M_temp_x, dev_M_temp_y, dev_M_temp_z,
											 dev_H1_x,     dev_H1_y,     dev_H1_z,
											 delta_x,      dev_Aex_XM,   dev_Ms_XM, dev_indicator1);
	
	Kernel_Ha_with_right<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
											  dev_M_temp_x, dev_M_temp_y, dev_M_temp_z,
											  dev_H2_x,     dev_H2_y,     dev_H2_z,
											  delta_x,      dev_Aex_XP,   dev_Ms_XP, dev_indicator1);

	Kernel_Ha_with_up<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
										   dev_M_temp_x, dev_M_temp_y, dev_M_temp_z,
										   dev_H3_x,     dev_H3_y,     dev_H3_z,
										   delta_y,      dev_Aex_YM,   dev_Ms_YM, dev_indicator1);

	Kernel_Ha_with_down<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
											 dev_M_temp_x, dev_M_temp_y, dev_M_temp_z,
											 dev_H4_x,     dev_H4_y,     dev_H4_z,
											 delta_y,      dev_Aex_YP,   dev_Ms_YP, dev_indicator1);
	Kernel_Ha_with_att<<<grids1, blocks1>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, mNx, mNy,
											dev_M_temp_x, dev_M_temp_y, dev_M_temp_z,
											dev_H1_x, dev_H1_y, dev_H1_z,
											dev_H2_x, dev_H2_y, dev_H2_z,
											dev_H3_x, dev_H3_y, dev_H3_z,
											dev_H4_x, dev_H4_y, dev_H4_z, 
											dev_Ha_x, dev_Ha_y, dev_Ha_z,
											dev_Hal_x,dev_Hal_y,dev_Hal_z,
											delta_z,  dev_Aex_ZP,dev_Aex_ZM, dev_Ms_ZP, dev_Ms_ZM);
	
	return true;
}


//bool Hk_field(void)
//{
//	// Setup kernel configuration
//	dim3 blocks1(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);  // dimension of block
//	int  grid_x = (mNx % BLOCK_SIZE_X) ? (mNx/BLOCK_SIZE_X + 1) : (mNx/BLOCK_SIZE_X),
//		 grid_y = (mNy % BLOCK_SIZE_Y) ? (mNy/BLOCK_SIZE_Y + 1) : (mNy/BLOCK_SIZE_Y),
//		 grid_z = (mNz % BLOCK_SIZE_Z) ? (mNz/BLOCK_SIZE_Z + 1) : (mNz/BLOCK_SIZE_Z);
//	dim3 grids1(grid_x, grid_y, grid_z);
//
//	Kernel_Hk_field<<<grids1, blocks1>>>(dev_theta, dev_phi, dev_Hk_x,  dev_Hk_y, dev_Hk_z, Ku, Ms);
//
//	return true;
//}