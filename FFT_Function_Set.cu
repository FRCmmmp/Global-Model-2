#include "FFT_Function_Set.h"


bool G_matrix(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad); 

int G_tensor(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad);

int Hms(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad,
		dim3 grids2, dim3 blocks2);




static bool G_matrix(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad)
{
	int idx;
	double x, y, z;
	double xmd, xpd, ymd, ypd, zmd, zpd; // xmd: x minus dx/2; xpd: x plus dx/2

	for (int k = 0; k < lz_zero_pad-1; k++){
		for (int j = 0; j < ly_zero_pad-1; j++){
			for (int i = 0; i < lx_zero_pad-1; i++){
				x = ((i + 1) - lx_zero_pad/2) * dx;
				y = ((j + 1) - ly_zero_pad/2) * dy;
				z = ((k + 1) - lz_zero_pad/2) * dz;

				xmd = pow(x-dx/2, 2);
				xpd = pow(x+dx/2, 2);
				ymd = pow(y-dy/2, 2);
				ypd = pow(y+dy/2, 2);
				zmd = pow(z-dz/2, 2);
				zpd = pow(z+dz/2, 2);
				
				idx = i + j*(lx_zero_pad) + k*(lx_zero_pad)*(ly_zero_pad);				
				Gxx_1d_cmplx[idx] =  ((atan( (y-dy/2) * (z-dz/2)/((x-dx/2) * pow(xmd + ymd + zmd, 0.5)) ) + 
									   atan( (y+dy/2) * (z+dz/2)/((x-dx/2) * pow(xmd + ypd + zpd, 0.5)) )) - 
					                  (atan( (y-dy/2) * (z+dz/2)/((x-dx/2) * pow(xmd + ymd + zpd, 0.5)) ) + 
									   atan( (y+dy/2) * (z-dz/2)/((x-dx/2) * pow(xmd + ypd + zmd, 0.5)) ))) - 
							         ((atan( (y-dy/2) * (z-dz/2)/((x+dx/2) * pow(xpd + ymd + zmd, 0.5)) ) + 
					                   atan( (y+dy/2) * (z+dz/2)/((x+dx/2) * pow(xpd + ypd + zpd, 0.5)) )) - 
					                  (atan( (y-dy/2) * (z+dz/2)/((x+dx/2) * pow(xpd + ymd + zpd, 0.5)) ) + 
					                   atan( (y+dy/2) * (z-dz/2)/((x+dx/2) * pow(xpd + ypd + zmd, 0.5)) )));
				Gxy_1d_cmplx[idx] =  ( log( 4 * (-(z-dz/2) + pow(xmd + ymd + zmd, 0.5) ) * (-(z+dz/2) + pow(xpd + ymd + zpd, 0.5) )) - 
					                   log( 4 * (-(z-dz/2) + pow(xpd + ymd + zmd, 0.5) ) * (-(z+dz/2) + pow(xmd + ymd + zpd, 0.5) ))) - 
					                 ( log( 4 * (-(z-dz/2) + pow(xmd + ypd + zmd, 0.5) ) * (-(z+dz/2) + pow(xpd + ypd + zpd, 0.5) )) - 
					                   log( 4 * (-(z-dz/2) + pow(xpd + ypd + zmd, 0.5) ) * (-(z+dz/2) + pow(xmd + ypd + zpd, 0.5) )));
				
				Gxz_1d_cmplx[idx] = ( log( 4 * (-(y-dy/2) + pow(xmd + ymd + zmd, 0.5) ) * (-(y+dy/2) + pow(xpd + ypd + zmd, 0.5) )) - 
					                  log( 4 * (-(y-dy/2) + pow(xpd + ymd + zmd, 0.5) ) * (-(y+dy/2) + pow(xmd + ypd + zmd, 0.5) ))) - 
					                ( log( 4 * (-(y-dy/2) + pow(xmd + ymd + zpd, 0.5) ) * (-(y+dy/2) + pow(xpd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(y-dy/2) + pow(xpd + ymd + zpd, 0.5) ) * (-(y+dy/2) + pow(xmd + ypd + zpd, 0.5) )));
				Gyx_1d_cmplx[idx] = ( log( 4 * (-(z-dz/2) + pow(xmd + ymd + zmd, 0.5) ) * (-(z+dz/2) + pow(xmd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(z-dz/2) + pow(xmd + ypd + zmd, 0.5) ) * (-(z+dz/2) + pow(xmd + ymd + zpd, 0.5) ))) - 
					                ( log( 4 * (-(z-dz/2) + pow(xpd + ymd + zmd, 0.5) ) * (-(z+dz/2) + pow(xpd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(z-dz/2) + pow(xpd + ypd + zmd, 0.5) ) * (-(z+dz/2) + pow(xpd + ymd + zpd, 0.5) )));
				Gyy_1d_cmplx[idx] = ((atan( (x-dx/2) * (z-dz/2)/((y-dy/2) * pow(xmd + ymd + zmd, 0.5)) ) + 
					                  atan( (x+dx/2) * (z+dz/2)/((y-dy/2) * pow(xpd + ymd + zpd, 0.5)) )) - 
					                 (atan( (x-dx/2) * (z+dz/2)/((y-dy/2) * pow(xmd + ymd + zpd, 0.5)) ) +
					                  atan( (x+dx/2) * (z-dz/2)/((y-dy/2) * pow(xpd + ymd + zmd, 0.5)) ))) - 
					                ((atan( (x-dx/2) * (z-dz/2)/((y+dy/2) * pow(xmd + ypd + zmd, 0.5)) ) + 
					                  atan( (x+dx/2) * (z+dz/2)/((y+dy/2) * pow(xpd + ypd + zpd, 0.5)) )) - 
					                 (atan( (x-dx/2) * (z+dz/2)/((y+dy/2) * pow(xmd + ypd + zpd, 0.5)) ) + 
					                  atan( (x+dx/2) * (z-dz/2)/((y+dy/2) * pow(xpd + ypd + zmd, 0.5)) )));
				Gyz_1d_cmplx[idx] = ( log( 4 * (-(x-dx/2) + pow(xmd + ymd + zmd, 0.5) ) * (-(x+dx/2) + pow(xpd + ypd + zmd, 0.5) )) - 
					                  log( 4 * (-(x-dx/2) + pow(xmd + ypd + zmd, 0.5) ) * (-(x+dx/2) + pow(xpd + ymd + zmd, 0.5) ))) - 
					                ( log( 4 * (-(x-dx/2) + pow(xmd + ymd + zpd, 0.5) ) * (-(x+dx/2) + pow(xpd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(x-dx/2) + pow(xmd + ypd + zpd, 0.5) ) * (-(x+dx/2) + pow(xpd + ymd + zpd, 0.5) )));
				Gzx_1d_cmplx[idx] = ( log( 4 * (-(y-dy/2) + pow(xmd + ymd + zmd, 0.5) ) * (-(y+dy/2) + pow(xmd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(y-dy/2) + pow(xmd + ymd + zpd, 0.5) ) * (-(y+dy/2) + pow(xmd + ypd + zmd, 0.5) ))) - 
					                ( log( 4 * (-(y-dy/2) + pow(xpd + ymd + zmd, 0.5) ) * (-(y+dy/2) + pow(xpd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(y-dy/2) + pow(xpd + ymd + zpd, 0.5) ) * (-(y+dy/2) + pow(xpd + ypd + zmd, 0.5) )));
				Gzy_1d_cmplx[idx] = ( log( 4 * (-(x-dx/2) + pow(xmd + ymd + zmd, 0.5) ) * (-(x+dx/2) + pow(xpd + ymd + zpd, 0.5) )) - 
					                  log( 4 * (-(x-dx/2) + pow(xmd + ymd + zpd, 0.5) ) * (-(x+dx/2) + pow(xpd + ymd + zmd, 0.5) ))) - 
					                ( log( 4 * (-(x-dx/2) + pow(xmd + ypd + zmd, 0.5) ) * (-(x+dx/2) + pow(xpd + ypd + zpd, 0.5) )) - 
					                  log( 4 * (-(x-dx/2) + pow(xmd + ypd + zpd, 0.5) ) * (-(x+dx/2) + pow(xpd + ypd + zmd, 0.5) )));
				Gzz_1d_cmplx[idx] = ((atan( (x-dx/2) * (y-dy/2)/((z-dz/2) * pow(xmd + ymd + zmd, 0.5)) ) + 
					                  atan( (x+dx/2) * (y+dy/2)/((z-dz/2) * pow(xpd + ypd + zmd, 0.5)) )) - 
					                 (atan( (x-dx/2) * (y+dy/2)/((z-dz/2) * pow(xmd + ypd + zmd, 0.5)) ) + 
					                  atan( (x+dx/2) * (y-dy/2)/((z-dz/2) * pow(xpd + ymd + zmd, 0.5)) ))) - 
					                ((atan( (x-dx/2) * (y-dy/2)/((z+dz/2) * pow(xmd + ymd + zpd, 0.5)) ) + 
					                  atan( (x+dx/2) * (y+dy/2)/((z+dz/2) * pow(xpd + ypd + zpd, 0.5)) )) - 
					                 (atan( (x-dx/2) * (y+dy/2)/((z+dz/2) * pow(xmd + ypd + zpd, 0.5)) ) + 
					                  atan( (x+dx/2) * (y-dy/2)/((z+dz/2) * pow(xpd + ymd + zpd, 0.5)) )));
			}
		}
	}
	return true;
}

int G_tensor(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad)
{
	
	
	// Demag Tensor
	if (!G_matrix(lx_zero_pad, ly_zero_pad, lz_zero_pad)) 
	{ printf("G_matrix() failed!\n"); }
    

	//---------- Set to Device ----------//
	cudaMemcpy(dev_Gxx_cufft, Gxx_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gxy_cufft, Gxy_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gxz_cufft, Gxz_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gyx_cufft, Gyx_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gyy_cufft, Gyy_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gyz_cufft, Gyz_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gzx_cufft, Gzx_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gzy_cufft, Gzy_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Gzz_cufft, Gzz_1d_cmplx, (lx_zero_pad)*(ly_zero_pad)*(lz_zero_pad)*sizeof(complex<float>), cudaMemcpyHostToDevice);
	
	

	//---------- CUFFT ----------//
	cufftHandle plan;
	if (cufftPlan3d(&plan, lz_zero_pad, ly_zero_pad, lx_zero_pad, CUFFT_C2C) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: Plan creation failed");};
	if (cufftExecC2C(plan, dev_Gxx_cufft, dev_Gxx_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gxy_cufft, dev_Gxy_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gxz_cufft, dev_Gxz_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gyx_cufft, dev_Gyx_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gyy_cufft, dev_Gyy_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gyz_cufft, dev_Gyz_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gzx_cufft, dev_Gzx_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gzy_cufft, dev_Gzy_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cufftExecC2C(plan, dev_Gzz_cufft, dev_Gzz_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");};
	if (cudaThreadSynchronize() != cudaSuccess){
		fprintf(stderr, "CUDA error: Failed to synchronize\n");
		return 1;}
	cufftDestroy(plan);
     
    
	

	return 1;

}

int Hms(int lx_zero_pad, int ly_zero_pad, int lz_zero_pad,
		dim3 grids2, dim3 blocks2)
{
	
	//----------CUFFT----------//
	cufftHandle plan1;
	if (cufftPlan3d(&plan1, lz_zero_pad, ly_zero_pad, lx_zero_pad, CUFFT_C2C) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: Plan creation failed");}
	if (cufftExecC2C(plan1, dev_Mx_cufft, dev_Mx_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");}
	if (cufftExecC2C(plan1, dev_My_cufft, dev_My_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");}
	if (cufftExecC2C(plan1, dev_Mz_cufft, dev_Mz_cufft, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");}
	if (cudaThreadSynchronize() != cudaSuccess){
		fprintf(stderr, "CUDA error: Failed to synchronize\n");
		return 1;
	}
	cufftDestroy(plan1);

	Kernel_CUFFT_M_times_G<<<grids2, blocks2>>>(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, BLK_SZ_Z, lx_zero_pad, ly_zero_pad,
		                                        dev_Gxx_cufft,  dev_Gxy_cufft,  dev_Gxz_cufft,
											    dev_Gyx_cufft,  dev_Gyy_cufft,  dev_Gyz_cufft,
											    dev_Gzx_cufft,  dev_Gzy_cufft,  dev_Gzz_cufft,
											    dev_Mx_cufft,   dev_My_cufft,   dev_Mz_cufft,
											    dev_Hd_x_cufft, dev_Hd_y_cufft, dev_Hd_z_cufft);

	//////---------- CUIFFT ----------//
    cufftHandle plan2;
	if (cufftPlan3d(&plan2, lz_zero_pad, ly_zero_pad, lx_zero_pad, CUFFT_C2C) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: Plan creation failed");}
	if (cufftExecC2C(plan2, dev_Hd_x_cufft, dev_Hd_x_cufft, CUFFT_INVERSE) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");}
	if (cufftExecC2C(plan2, dev_Hd_y_cufft, dev_Hd_y_cufft, CUFFT_INVERSE) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");}
	if (cufftExecC2C(plan2, dev_Hd_z_cufft, dev_Hd_z_cufft, CUFFT_INVERSE) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecC2C failed");}
	if (cudaThreadSynchronize() != cudaSuccess){
		fprintf(stderr, "CUDA error: Failed to synchronize\n");
		return 1;
	}
	cufftDestroy(plan2);

	cudaMemcpy(Hd_x_1d_cmplx, dev_Hd_x_cufft, (lx_zero_pad*ly_zero_pad*lz_zero_pad)*sizeof(complex<float>), cudaMemcpyDeviceToHost);
	cudaMemcpy(Hd_y_1d_cmplx, dev_Hd_y_cufft, (lx_zero_pad*ly_zero_pad*lz_zero_pad)*sizeof(complex<float>), cudaMemcpyDeviceToHost);
	cudaMemcpy(Hd_z_1d_cmplx, dev_Hd_z_cufft, (lx_zero_pad*ly_zero_pad*lz_zero_pad)*sizeof(complex<float>), cudaMemcpyDeviceToHost);

	// Normalize
	for (int i = 0; i < lx_zero_pad*ly_zero_pad*lz_zero_pad; i++){
		Hd_x_1d[i] = (double)real(Hd_x_1d_cmplx[i]) / (lx_zero_pad*ly_zero_pad*lz_zero_pad);
		Hd_y_1d[i] = (double)real(Hd_y_1d_cmplx[i]) / (lx_zero_pad*ly_zero_pad*lz_zero_pad);
		Hd_z_1d[i] = (double)real(Hd_z_1d_cmplx[i]) / (lx_zero_pad*ly_zero_pad*lz_zero_pad);
	}

	/////// Watch /////
	/*cudaMemcpy(Gxx_1d_cmplx, dev_Gxx_cufft, lx_zero_pad*ly_zero_pad*lz_zero_pad * sizeof(complex<float>), cudaMemcpyDeviceToHost);*/
	/*FILE *pfile;
	pfile = fopen("watch.out", "w");
	for (int idx=0; idx<lx_zero_pad*ly_zero_pad*lz_zero_pad; ++idx){
		if (idx%lx_zero_pad==0 && idx!=0) fprintf(pfile, "\n");
		if (idx%(lx_zero_pad*ly_zero_pad)==0 && idx!=0) fprintf(pfile, "\n\n");
		fprintf(pfile, "%15.5lf", Hd_z_1d[idx]);
	}
	fclose(pfile);*/
	/////// Watch /////

	// Shift the Hms field
	int idxx = 0;
	for (int k = int(lz_zero_pad/2+0.5)-1; k < int(lz_zero_pad/2+0.5)-1+mNz; k++){
		for (int j = int(ly_zero_pad/2+0.5)-1; j < int(ly_zero_pad/2+0.5)-1+mNy; j++){
			for (int i = int(lx_zero_pad/2+0.5)-1; i < int(lx_zero_pad/2+0.5)-1+mNx; i++){
				int idx = i + j*lx_zero_pad + k*lx_zero_pad*ly_zero_pad;
				Hd_x_1d_shift[idxx] = DmagFAC*Hd_x_1d[idx];
				Hd_y_1d_shift[idxx] = DmagFAC*Hd_y_1d[idx];
				Hd_z_1d_shift[idxx] = DmagFAC*Hd_z_1d[idx];
				idxx++;		
			}
		}
	}
	/////// Watch /////
	/*cudaMemcpy(Gxx_1d_cmplx, dev_Gxx_cufft, lx_zero_pad*ly_zero_pad*lz_zero_pad * sizeof(complex<float>), cudaMemcpyDeviceToHost);*/
	/*FILE *pfile;
	pfile = fopen("watch.out", "w");
	for (int idx=0; idx<mNx*mNy*mNz; ++idx){
		if (idx%mNx==0 && idx!=0) fprintf(pfile, "\n");
		if (idx%(mNx*mNy)==0 && idx!=0) fprintf(pfile, "\n\n");
		fprintf(pfile, "%15.5lf", Hd_z_1d_shift[idx]);
	}
	fclose(pfile);*/
	/////// Watch /////

	return 1;
}



