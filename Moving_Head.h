#include "Parameters.h"
#include "Parameters_input.h"
#include "LLG_kernel.cu"
#include "Outputs.h"
#include <cstdio>
#include <cstdlib>



static void ExtractField_from_FP(int mNz_i, int mNz_f, int f_x00){
	int idx;
	for (int k = mNz_i; k < mNz_f; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = f_x00; i < mNx; i++){
				idx = i + j*mNx + k*mNx*mNy;
				Happl_x_temp[idx] = FP[0][(i-f_x00)+j*mNx+k*mNx*mNy];
				Happl_y_temp[idx] = FP[1][(i-f_x00)+j*mNx+k*mNx*mNy];
				Happl_z_temp[idx] = FP[2][(i-f_x00)+j*mNx+k*mNx*mNy];
			}
			for (int i = 0; i < f_x00; i++){
				idx = i + j*mNx + k*mNx*mNy;
				Happl_x_temp[idx] = FP_trail[0][((mNx-f_x00)+i)+j*mNx+k*mNx*mNy];
				Happl_y_temp[idx] = FP_trail[1][((mNx-f_x00)+i)+j*mNx+k*mNx*mNy];
				Happl_z_temp[idx] = FP_trail[2][((mNx-f_x00)+i)+j*mNx+k*mNx*mNy];
			}
		}
	}
	for (int k = mNz_i; k < mNz_f; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				idx = i + j*mNx + k*mNx*mNy;
				memTemp0[idx] = abs(Happl_x_temp[idx]/10);
				memTemp1[idx] = abs(Happl_y_temp[idx]/10);
				memTemp2[idx] = abs(Happl_z_temp[idx]/10);
			}
		}
	}
}

static void FieldInching(int mNz_i, int mNz_f){
	int idx;
	for (int k = mNz_i; k < mNz_f; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){ 
				idx = i + j*mNx + k*mNx*mNy;
				Happl_x_temp[idx] = abs(Happl_x_temp[idx]) - memTemp0[idx];
				Happl_y_temp[idx] = abs(Happl_y_temp[idx]) - memTemp1[idx];
				Happl_z_temp[idx] = abs(Happl_z_temp[idx]) - memTemp2[idx];
				if (i < mNx-1){
					Happl_x_temp[(i+1)+j*mNx+k*mNx*mNy] = abs(Happl_x_temp[(i+1)+j*mNx+k*mNx*mNy]) + memTemp0[i+j*(mNx)+k*(mNx)*(mNy)];
					Happl_y_temp[(i+1)+j*mNx+k*mNx*mNy] = abs(Happl_y_temp[(i+1)+j*mNx+k*mNx*mNy]) + memTemp1[i+j*(mNx)+k*(mNx)*(mNy)];
					Happl_z_temp[(i+1)+j*mNx+k*mNx*mNy] = abs(Happl_z_temp[(i+1)+j*mNx+k*mNx*mNy]) + memTemp2[i+j*(mNx)+k*(mNx)*(mNy)];
				}
				/*if (Happl_z_temp[idx] > 0.){
					Happl_x_temp[idx] = abs(Happl_x_temp[idx]) - memTemp0[idx];
					Happl_y_temp[idx] = abs(Happl_y_temp[idx]) - memTemp1[idx];
					Happl_z_temp[idx] = abs(Happl_z_temp[idx]) - memTemp2[idx];
				}
				else { 
					Happl_x_temp[idx] = 0.;
					Happl_y_temp[idx] = 0.;
					Happl_z_temp[idx] = 0.;
				}*/
			}
		}
	}
}

bool Moving_Head(int t)
{
	
	if (EXT_FP_PROFILE){
		
		// Head Field Flying (mimicing recording process)
		#ifndef __FIELD_SPEC__
		if (t%(dxNt) == 0){
			if (f_x0 < mNx){
				ExtractField_from_FP(0,      mNz_1, f_x0);
				ExtractField_from_FP(mNz_12, mNz_2, f_x0);
				ExtractField_from_FP(mNz_23, mNz_3, f_x0);
				ExtractField_from_FP(mNz_34, mNz_4, f_x0);
				ExtractField_from_FP(mNz_45, mNz_5, f_x0);
				ExtractField_from_FP(mNz_56, mNz_6, f_x0);
				f_x0++;
			}
		}
		if (f_x0 < mNx){
			if (t%(dxNt/10) == 0){
				FieldInching(0,      mNz_1);
				FieldInching(mNz_12, mNz_2);
				FieldInching(mNz_23, mNz_3);
				FieldInching(mNz_34, mNz_4);
				FieldInching(mNz_45, mNz_5);
				FieldInching(mNz_56, mNz_6);
			}
		}
		if (t%sfcNt == 0){ 
			if (*idx_f < sfNc-1 ){
				*idx_f = *idx_f + 1;
			}
			else {
				*idx_f = sfNc-1;
			}
		}
		if (f_x0 < mNx){
			for (int k = 0; k < mNz; k++){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){ 
						int idx = i+j*mNx+k*mNx*mNy;
						Happl_x[idx] = Happl_x_temp[idx];
						Happl_y[idx] = Happl_y_temp[idx];
						Happl_z[idx] = Happl_z_temp[idx]*fSeq[*idx_f];
					}
				}
			}
		}
		else {
			for (int k = 0; k < mNz; k++){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){ 
						int idx = i+j*mNx+k*mNx*mNy;
						Happl_x[idx] = 0.;
						Happl_y[idx] = 0.;
						Happl_z[idx] = 0.;
					}
				}
			}
		}

		// For Output_1D_Happl_xyz_atTrans_XXX.dat
		if (t%sfcNt == 0){
			char buf1[3], buf2[80], buf3[4];
			sprintf(buf1, "%d", *idx_f);
			strcpy(buf2, "Output_1D_Happl_xyz_atTrans_");
			strcpy(buf3, ".out");
			strcat(buf2, buf1);
			strcat(buf2, buf3);
			if (!Output_Float_1D_Format_6col(mNx, mNy, mNz, Happl_x, Happl_y, Happl_z, buf2)) { printf("Output_Float_1D_Format_6col() fails!\n"); }
			strcpy(buf2, "Output_1D_Mxyz_atTrans_");
			strcpy(buf3, ".out");
			strcat(buf2, buf1);
			strcat(buf2, buf3);
			if (!Output_Float_1D_Format_6col(mNx, mNy, mNz, Mx, My, Mz, buf2)) { printf("Output_Float_1D_Format_6col() fails!\n"); }
			*Seq_f = 0;
		}
		if (t%(sfcNt/10) == 0){
			char buf4[4], buf5[80], buf6[4];
			sprintf(buf4, "%d.%d", *idx_f, *Seq_f);
			strcpy(buf5, "Output_1D_Happl_xyz_Seq_");
			strcpy(buf6, ".out");
			strcat(buf5, buf4);
			strcat(buf5, buf6);
			if (!Output_Float_1D_Format_6col(mNx, mNy, mNz, Happl_x, Happl_y, Happl_z, buf5)) { printf("Output_Float_1D_Format_6col() fails!\n"); }
			strcpy(buf5, "Output_1D_Mxyz_Seq_");
			strcpy(buf6, ".out");
			strcat(buf5, buf4);
			strcat(buf5, buf6);
			if (!Output_Float_1D_Format_6col(mNx, mNy, mNz, Mx, My, Mz, buf5)) { printf("Output_Float_1D_Format_6col() fails!\n"); }
			*Seq_f = *Seq_f + 1;
		}
		#endif

	}
		
	if (EXT_FP_PROFILE){
		
		#ifndef __ATHERMAL__
		if (t == 0 && !THERMAL){

			// Temperature Profile
			for (int i = 0; i < mNx; i++){
				for (int j = 0; j < mNy; j++){
					for (int k = 0; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						T[idx] = 0.0;
					}
				}
			}

			// Set magnetic parameters for each cell
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					for (int k = 0; k < mNz_1; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L1_Ms;
						host_Ku[idx]    = L1_Ku*(1+std_Ku[idx]*DEL_Hk*dHk1_scale);  
						host_Aex[idx]   = L1_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_1; k < mNz_12; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL12_Ms;
						host_Ku[idx]    = 0.0;  
						host_Aex[idx]   = L1_Aex*((double)BL12*(1+std_Aex[idx]*dBL12_scale));
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_12; k < mNz_2; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L2_Ms;
						host_Ku[idx]    = L2_Ku*(1+std_Ku[idx]*DEL_Hk*dHk2_scale);  
						host_Aex[idx]   = L2_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_2; k < mNz_23; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL23_Ms;
						host_Ku[idx]    = 0.0;  
						host_Aex[idx]   = L2_Aex*((double)BL23*(1+std_Aex[idx]*dBL23_scale));
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_23; k < mNz_3; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L3_Ms;
						host_Ku[idx]    = L3_Ku*(1+std_Ku[idx]*DEL_Hk*dHk3_scale);  
						host_Aex[idx]   = L3_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_3; k < mNz_34; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL34_Ms;
						host_Ku[idx]    = 0.0;  
						host_Aex[idx]   = L3_Aex*((double)BL34*(1+std_Aex[idx]*dBL34_scale));
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_34; k < mNz_4; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L4_Ms;
						host_Ku[idx]    = L4_Ku*(1+std_Ku[idx]*DEL_Hk*dHk4_scale);  
						host_Aex[idx]   = L4_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_4; k < mNz_45; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL45_Ms;
						host_Ku[idx]    = 0.0;  
						host_Aex[idx]   = L4_Aex*((double)BL45*(1+std_Aex[idx]*dBL45_scale));
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_45; k < mNz_5; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L5_Ms;
						host_Ku[idx]    = L5_Ku*(1+std_Ku[idx]*DEL_Hk*dHk5_scale);
						host_Aex[idx]   = L5_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_5; k < mNz_56; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL56_Ms;
						host_Ku[idx]    = 0.0;
						host_Aex[idx]   = L5_Aex*((double)BL56*(1+std_Aex[idx]*dBL56_scale));
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_56; k < mNz_6; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L6_Ms;
						host_Ku[idx]    = L6_Ku*(1+std_Ku[idx]*DEL_Hk*dHk6_scale);
						host_Aex[idx]   = L6_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L6_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_6; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = 1.0;
						host_Ku[idx]    = 0.0;
						host_Aex[idx]   = 0.0;
						host_alpha[idx] = 1.0;
						D[idx] = 1.0;
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
				}
			}

			// Set defects in CGC layer
			if (CGC_DEF){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){
						switch (CGC_label){
							case 1:
								for (int k = 0; k < mNz_1; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 2:
								for (int k = mNz_12; k < mNz_2; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 3:
								for (int k = mNz_23; k < mNz_3; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 4:
								for (int k = mNz_34; k < mNz_4; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 5:
								for (int k = mNz_45; k < mNz_5; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 6:
								for (int k = mNz_56; k < mNz_6; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
						}
					}
				}

				// Write to file
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ku, "host_Ku.out")) { printf("Output_Float_3D_Format() failed!\n"); }
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Aex, "host_Aex.out")) { printf("Output_Float_3D_Format() failed!\n"); }
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ms, "host_Ms.out")) { printf("Output_Float_3D_Format() failed!\n"); }
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_alpha, "host_alpha.out")) { printf("Output_Float_3D_Format() failed!\n"); }
			}
		}
		#endif

	}


	if (EXT_FP_PROFILE){
		
		#ifndef __ATHERMAL__
		if (t == 0 && !THERMAL){
			cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		}
		#endif

		#ifndef __THERMAL_NON_HAMR__
		if (THERMAL){
			cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		}
		#endif

		cudaMemcpy(dev_Happl_x, Happl_x,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_y, Happl_y,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_z, Happl_z,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	}

	return true;
}

bool NonMoving_Head(int t, double SweepFieldStep1){
	
	if (MH_LOOP){
		#ifndef __FIELD_SPEC__
		// External Field Sweeping (mimicing metrology)
		if ((t > 0) && (t%FieldSweepTimeStep == 0)){
			IniScalFact = IniScalFact - SweepFieldStep1;
			Happl = field1*IniScalFact;
		}
		
		// Set happl for each layer
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				for (int k = 0; k < mNz_1; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field1*IniScalFact*sin(angle1/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field1*IniScalFact*cos(angle1/180.0*PI);
				}
				for (int k = mNz_1; k < mNz_12; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field12*IniScalFact*sin(angle12/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field12*IniScalFact*cos(angle12/180.0*PI);
				}
				for (int k = mNz_12; k < mNz_2; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field2*IniScalFact*sin(angle2/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field2*IniScalFact*cos(angle2/180.0*PI);
				}
				for (int k = mNz_2; k < mNz_23; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field23*IniScalFact*sin(angle23/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field23*IniScalFact*cos(angle23/180.0*PI);
				}
				for (int k = mNz_23; k < mNz_3; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field3*IniScalFact*sin(angle3/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field3*IniScalFact*cos(angle3/180.0*PI);
				}
				for (int k = mNz_3; k < mNz_34; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field34*IniScalFact*sin(angle34/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field34*IniScalFact*cos(angle34/180.0*PI);
				}
				for (int k = mNz_34; k < mNz_4; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field4*IniScalFact*sin(angle4/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field4*IniScalFact*cos(angle4/180.0*PI);
				}
				for (int k = mNz_4; k < mNz_45; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field45*IniScalFact*sin(angle45/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field45*IniScalFact*cos(angle45/180.0*PI);
				}
				for (int k = mNz_45; k < mNz_5; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field5*IniScalFact*sin(angle5/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field5*IniScalFact*cos(angle5/180.0*PI);
				}
				for (int k = mNz_5; k < mNz_56; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field56*IniScalFact*sin(angle56/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field56*IniScalFact*cos(angle56/180.0*PI);
				}
				for (int k = mNz_56; k < mNz_6; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					Happl_x[idx] = field6*IniScalFact*sin(angle6/180.0*PI);
					Happl_y[idx] = 0.;
					Happl_z[idx] = field6*IniScalFact*cos(angle6/180.0*PI);
				}
			}
		}
		#endif
	}

	if (MH_LOOP){
		
		#ifndef __ATHERMAL__
		if (t == 0 && !THERMAL){
			
			// Temperature Profile
			for (int i = 0; i < mNx; i++){
				for (int j = 0; j < mNy; j++){
					for (int k = 0; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						T[idx] = 0.0;
					}
				}
			}

			// Set magnetic parameters for each cell
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					for (int k = 0; k < mNz_1; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L1_Ms;
						host_Ku[idx]    = L1_Ku*(1+std_Ku[idx]*DEL_Hk*dHk1_scale);  
						host_Aex[idx]   = L1_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_1; k < mNz_12; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL12_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L1_Aex*((double)BL12*(1+std_Aex[idx]*dBL12_scale));
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_12; k < mNz_2; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L2_Ms;
						host_Ku[idx]    = L2_Ku*(1+std_Ku[idx]*DEL_Hk*dHk2_scale);  
						host_Aex[idx]   = L2_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_2; k < mNz_23; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL23_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L2_Aex*((double)BL23*(1+std_Aex[idx]*dBL23_scale));
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_23; k < mNz_3; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L3_Ms;
						host_Ku[idx]    = L3_Ku*(1+std_Ku[idx]*DEL_Hk*dHk3_scale);  
						host_Aex[idx]   = L3_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_3; k < mNz_34; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL34_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L3_Aex*((double)BL34*(1+std_Aex[idx]*dBL34_scale));
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_34; k < mNz_4; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L4_Ms;
						host_Ku[idx]    = L4_Ku*(1+std_Ku[idx]*DEL_Hk*dHk4_scale);  
						host_Aex[idx]   = L4_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_4; k < mNz_45; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL45_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L4_Aex*((double)BL45*(1+std_Aex[idx]*dBL45_scale));
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_45; k < mNz_5; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L5_Ms;
						host_Ku[idx]    = L5_Ku*(1+std_Ku[idx]*DEL_Hk*dHk5_scale);
						host_Aex[idx]   = L5_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_5; k < mNz_56; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL56_Ms;
						host_Ku[idx]    = 0.;
						host_Aex[idx]   = L5_Aex*((double)BL56*(1+std_Aex[idx]*dBL56_scale));
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_56; k < mNz_6; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L6_Ms;
						host_Ku[idx]    = L6_Ku*(1+std_Ku[idx]*DEL_Hk*dHk6_scale);
						host_Aex[idx]   = L6_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L6_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_6; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = 1.0;
						host_Ku[idx]    = 0.0;
						host_Aex[idx]   = 0.0;
						host_alpha[idx] = 1.0;
						D[idx] = 0.0;
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
				}
			}

			// Set defects in CGC layer
			if (CGC_DEF){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){
						switch (CGC_label){
							case 1:
								for (int k = 0; k < mNz_1; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 2:
								for (int k = mNz_12; k < mNz_2; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 3:
								for (int k = mNz_23; k < mNz_3; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 4:
								for (int k = mNz_34; k < mNz_4; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 5:
								for (int k = mNz_45; k < mNz_5; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 6:
								for (int k = mNz_56; k < mNz_6; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
						}
					}
				}

				// Write to file
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ku, "host_Ku.out")) { printf("Output_Float_3D_Format() failed!\n"); }
			}
		}
		#endif

		

		#ifndef __THERMAL_NON_HAMR__
		if (t == 0 && THERMAL){
			
			// Temperature Profile
			for (int i = 0; i < mNx; i++){
				for (int j = 0; j < mNy; j++){
					for (int k = 0; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						T[idx] = Temperature;
					}
				}
			}

			// Set magnetic parameters for each cell
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					for (int k = 0; k < mNz_1; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L1_Ms;
						host_Ku[idx]    = L1_Ku*(1+std_Ku[idx]*DEL_Hk*dHk1_scale);  
						host_Aex[idx]   = L1_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_1; k < mNz_12; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL12_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L1_Aex*((double)BL12*(1+std_Aex[idx]*dBL12_scale));
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_12; k < mNz_2; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L2_Ms;
						host_Ku[idx]    = L2_Ku*(1+std_Ku[idx]*DEL_Hk*dHk2_scale);  
						host_Aex[idx]   = L2_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_2; k < mNz_23; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL23_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L2_Aex*((double)BL23*(1+std_Aex[idx]*dBL23_scale));
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_23; k < mNz_3; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L3_Ms;
						host_Ku[idx]    = L3_Ku*(1+std_Ku[idx]*DEL_Hk*dHk3_scale);  
						host_Aex[idx]   = L3_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_3; k < mNz_34; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL34_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L3_Aex*((double)BL34*(1+std_Aex[idx]*dBL34_scale));
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_34; k < mNz_4; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L4_Ms;
						host_Ku[idx]    = L4_Ku*(1+std_Ku[idx]*DEL_Hk*dHk4_scale);  
						host_Aex[idx]   = L4_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_4; k < mNz_45; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL45_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L4_Aex*((double)BL45*(1+std_Aex[idx]*dBL45_scale));
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_45; k < mNz_5; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L5_Ms;
						host_Ku[idx]    = L5_Ku*(1+std_Ku[idx]*DEL_Hk*dHk5_scale);
						host_Aex[idx]   = L5_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_5; k < mNz_56; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL56_Ms;
						host_Ku[idx]    = 0.;
						host_Aex[idx]   = L5_Aex*((double)BL56*(1+std_Aex[idx]*dBL56_scale));
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_56; k < mNz_6; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L6_Ms;
						host_Ku[idx]    = L6_Ku*(1+std_Ku[idx]*DEL_Hk*dHk6_scale);
						host_Aex[idx]   = L6_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L6_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_6; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = 1.0;
						host_Ku[idx]    = 0.0;
						host_Aex[idx]   = 0.0;
						host_alpha[idx] = 1.0;
						D[idx] = 0.0;
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
				}
			}

			// Set defects in CGC layer
			if (CGC_DEF){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){
						switch (CGC_label){
							case 1:
								for (int k = 0; k < mNz_1; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 2:
								for (int k = mNz_12; k < mNz_2; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 3:
								for (int k = mNz_23; k < mNz_3; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 4:
								for (int k = mNz_34; k < mNz_4; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 5:
								for (int k = mNz_45; k < mNz_5; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 6:
								for (int k = mNz_56; k < mNz_6; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
						}
					}
				}

				// Write to file
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ku, "host_Ku.out")) { printf("Output_Float_3D_Format() failed!\n"); }
			}
		}
		#endif

		#ifdef __HAMR__
		if (1){
			
			// Temperature Profile
			for (int i = 0; i < mNx; i++){
				for (int j = 0; j < mNy; j++){
					for (int k = 0; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						T[idx] = Temperature;
					}
				}
			}

			// Set magnetic parameters for each cell
			// Temperature-dependent input parameters for 1.5nm cells (h: high temperature, l: low temperature)
			double M1h1_alpha = 0.005112, M1h2_alpha = 0.003613,  M1h3_alpha = 3.536e-18, M1h4_alpha = 0.05414,
				   M1l1_alpha = 2.2e-5,   M1l2_alpha = 0.02,
				   M2h1_alpha = 0.005112, M2h2_alpha = 0.003613,  M2h3_alpha = 3.536e-18, M2h4_alpha = 0.05414,
				   M2l1_alpha = 2.2e-5,   M2l2_alpha = 0.02,
				   M3h1_alpha = 0.005112, M3h2_alpha = 0.003613,  M3h3_alpha = 3.536e-18, M3h4_alpha = 0.05414,
				   M3l1_alpha = 2.2e-5,   M3l2_alpha = 0.02;
			double T, A_mean, A_ratio;
	
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					for (int k = 0; k < mNz_1; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						T = uGrain[idx].getT()*(1 + uGrain[idx].getIndicator6()*dTc1_scale);  //dTc_scale is applied here for implementing Tc distribution.
						if (uGrain[idx].getIndicator1() != 0){
					uGrain[idx].setMs((136*M1_Ms_scale)*pow(739.4+Tc1_shift-T,0.3182)); //from Yipeng Jiao
					uGrain[idx].setKu(((-7.143e7/(763.1655+Tc1_shift))*T + (7.143e7*M1_Ku_scale))*(1+uGrain[idx].getIndicator4()*dHk1_scale)); //from Yipeng Jiao
					if (T >= 500 && T < 710+Tc1_shift){
						uGrain[idx].setAlpha(((M1h1_alpha*exp(M1h2_alpha*(T-Tc1_shift)) + M1h3_alpha*exp(M1h4_alpha*(T-Tc1_shift))) + M1_Alpha_shift)*Alpha_CALI); //set renormalized damping constant for high temperature range
					}
					else if (T >= 250 && T < 500){
						uGrain[idx].setAlpha( (M1l1_alpha*T + M1l2_alpha)*Alpha_CALI ); //set renormalized damping constant for low temperature range
					}
					else {
						uGrain[idx].setAlpha(1.);
					}
					if (T >= 250 && T <= 710+Tc1_shift){
						A_mean = -9.478e-17*pow((T-Tc1_shift),4) + (1.777e-13)*pow((T-Tc1_shift),3) + (-1.247e-10)*pow((T-Tc1_shift),2) +
                                 (3.631e-8)*pow((T-Tc1_shift),1) + (-2.218e-6);
						if (Tc1_shift > 0 && T <= 300+Tc1_shift) {
							A_mean = 1.482e-6;
						}
						A_ratio = 0.2408*exp(-pow(((T-Tc1_shift)-657.6)/20.38,2)) + 1.301*exp(-pow(((T-Tc1_shift)-300)/1.985e6,2)) +
                                 (-0.001452)*exp(-pow(((T-Tc1_shift)-370)/7.748,2));
					}
					else if (T == 0){
						uGrain[idx].setMs(M1_Ms);
						uGrain[idx].setKu(M1_Ku);
						uGrain[idx].setAex_xy(M1_Aex_xy);
						uGrain[idx].setAex_z(M1_Aex_z);
						uGrain[idx].setAex_mean(M1_Aex_xy*2/3 + M1_Aex_z*1/3);
						uGrain[idx].setAlpha(M1_alpha);
					}
					else {
						A_mean = 0;
						A_ratio = 1.0;
						uGrain[idx].setMs(100);
						uGrain[idx].setKu(0.);
						uGrain[idx].setAex_xy(0.);
						uGrain[idx].setAex_z(0.);
						uGrain[idx].setAex_mean(0.);
					}
					if (uGrain[idx].getAex_mean() >= 0){
						uGrain[idx].setAex_mean( A_mean );
					}
					else {
						uGrain[idx].setAex_mean( 0. );
					}
					uGrain[idx].setAex_z( 3*A_mean/(1+2*A_ratio) );
					uGrain[idx].setAex_xy( A_ratio*uGrain[idx].getAex_z() );
				}
				else if ( uGrain[idx].getIndicator1() == 0){
					uGrain[idx].setMs(1.);
					uGrain[idx].setKu(0.);
					uGrain[idx].setAex_xy(0.);
					uGrain[idx].setAex_z(0.);
					uGrain[idx].setAex_mean(0.);
					uGrain[idx].setAlpha(1.);
				}
				uGrain[idx].setD((2*kB*uGrain[idx].getT()*uGrain[idx].getAlpha())/(uGrain[idx].getMs()*delta_t*GAMMA*(delta_x*delta_y*delta_z)));
				uGrain[idx].setGamma(GAMMA/(1+pow(uGrain[idx].getAlpha(), 2.0)));
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					for (int k = 0; k < mNz_1; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L1_Ms;
						host_Ku[idx]    = L1_Ku*(1+std_Ku[idx]*DEL_Hk*dHk1_scale);  
						host_Aex[idx]   = L1_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_1; k < mNz_12; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL12_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L1_Aex*((double)BL12);
						host_alpha[idx] = L1_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_12; k < mNz_2; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L2_Ms;
						host_Ku[idx]    = L2_Ku*(1+std_Ku[idx]*DEL_Hk*dHk2_scale);  
						host_Aex[idx]   = L2_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_2; k < mNz_23; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL23_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L2_Aex*((double)BL23);
						host_alpha[idx] = L2_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_23; k < mNz_3; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L3_Ms;
						host_Ku[idx]    = L3_Ku*(1+std_Ku[idx]*DEL_Hk*dHk3_scale);  
						host_Aex[idx]   = L3_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_3; k < mNz_34; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL34_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L3_Aex*((double)BL34);
						host_alpha[idx] = L3_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_34; k < mNz_4; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L4_Ms;
						host_Ku[idx]    = L4_Ku*(1+std_Ku[idx]*DEL_Hk*dHk4_scale);  
						host_Aex[idx]   = L4_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_4; k < mNz_45; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL45_Ms;
						host_Ku[idx]    = 0.;  
						host_Aex[idx]   = L4_Aex*((double)BL45);
						host_alpha[idx] = L4_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_45; k < mNz_5; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L5_Ms;
						host_Ku[idx]    = L5_Ku*(1+std_Ku[idx]*DEL_Hk*dHk5_scale);
						host_Aex[idx]   = L5_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_5; k < mNz_56; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = BL56_Ms;
						host_Ku[idx]    = 0.;
						host_Aex[idx]   = L5_Aex*((double)BL56);
						host_alpha[idx] = L5_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_56; k < mNz_6; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = L6_Ms;
						host_Ku[idx]    = L6_Ku*(1+std_Ku[idx]*DEL_Hk*dHk6_scale);
						host_Aex[idx]   = L6_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
						host_alpha[idx] = L6_alpha;
						D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
					for (int k = mNz_6; k < mNz; k++){
						int idx = i + j*mNx + k*mNx*mNy;
						host_Ms[idx]    = 1.0;
						host_Ku[idx]    = 0.0;
						host_Aex[idx]   = 0.0;
						host_alpha[idx] = 1.0;
						D[idx] = 0.0;
						host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
					}
				}
			}

			// Set defects in CGC layer
			if (CGC_DEF){
				for (int j = 0; j < mNy; j++){
					for (int i = 0; i < mNx; i++){
						switch (CGC_label){
							case 1:
								for (int k = 0; k < mNz_1; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 2:
								for (int k = mNz_12; k < mNz_2; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 3:
								for (int k = mNz_23; k < mNz_3; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 4:
								for (int k = mNz_34; k < mNz_4; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 5:
								for (int k = mNz_45; k < mNz_5; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
							case 6:
								for (int k = mNz_56; k < mNz_6; k++){
									int idx = i + j*mNx + k*mNx*mNy;
									if (indicator7[idx] == 1) host_Ku[idx] = host_Ku[idx]*0.1;
								}
								break;
						}
					}
				}

				// Write to file
				if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ku, "host_Ku.out")) { printf("Output_Float_3D_Format() failed!\n"); }
			}
		}
		#endif

	}

	if (MH_LOOP){
		
		#ifndef __ATHERMAL__
		if (t == 0 && !THERMAL){
			cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		}
		#endif

		#ifndef __THERMAL_NON_HAMR__
		if (t == 0 && THERMAL){
			cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		}
		#endif
		
		#ifndef __HAMR__
		if (MODEL==1){
			cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		}
		#endif
		
		cudaMemcpy(dev_Happl_x, Happl_x,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_y, Happl_y,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_z, Happl_z,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	}

	return true;
}

bool DT_Analysis_Head(int t){
	
	
	if (t == 0){

		// Temperature Profile
		for (int i = 0; i < mNx; i++){
			for (int j = 0; j < mNy; j++){
				for (int k = 0; k < mNz; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					T[idx] = Temperature;
				}
			}
		}
		// Temperature-dependent input parameters for 2nm cells (h: high temperature, l: low temperature)
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				for (int k = 0; k < mNz_1; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L1_Ms;
					host_Ku[idx]    = L1_Ku*(1+std_Ku[idx]*DEL_Hk*dHk1_scale);  
					host_Aex[idx]   = L1_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L1_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_1; k < mNz_12; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL12_Ms;
					host_Ku[idx]    = 0.0;  
					host_Aex[idx]   = L1_Aex*((double)BL12*(1+std_Aex[idx]*dBL12_scale));
					host_alpha[idx] = L1_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_12; k < mNz_2; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L2_Ms;
					host_Ku[idx]    = L2_Ku*(1+std_Ku[idx]*DEL_Hk*dHk2_scale);  
					host_Aex[idx]   = L2_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L2_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_2; k < mNz_23; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL23_Ms;
					host_Ku[idx]    = 0.0;  
					host_Aex[idx]   = L2_Aex*((double)BL23*(1+std_Aex[idx]*dBL23_scale));
					host_alpha[idx] = L2_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_23; k < mNz_3; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L3_Ms;
					host_Ku[idx]    = L3_Ku*(1+std_Ku[idx]*DEL_Hk*dHk3_scale);  
					host_Aex[idx]   = L3_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L3_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_3; k < mNz_34; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL34_Ms;
					host_Ku[idx]    = 0.0;  
					host_Aex[idx]   = L3_Aex*((double)BL34*(1+std_Aex[idx]*dBL34_scale));
					host_alpha[idx] = L3_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_34; k < mNz_4; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L4_Ms;
					host_Ku[idx]    = L4_Ku*(1+std_Ku[idx]*DEL_Hk*dHk4_scale);  
					host_Aex[idx]   = L4_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L4_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_4; k < mNz_45; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL45_Ms;
					host_Ku[idx]    = 0.0;  
					host_Aex[idx]   = L4_Aex*((double)BL45*(1+std_Aex[idx]*dBL45_scale));
					host_alpha[idx] = L4_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_45; k < mNz_5; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L5_Ms;
					host_Ku[idx]    = L5_Ku*(1+std_Ku[idx]*DEL_Hk*dHk5_scale);
					host_Aex[idx]   = L5_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L5_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_5; k < mNz_56; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL56_Ms;
					host_Ku[idx]    = 0.0;
					host_Aex[idx]   = L5_Aex*((double)BL56*(1+std_Aex[idx]*dBL56_scale));
					host_alpha[idx] = L5_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_56; k < mNz_6; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L6_Ms;
					host_Ku[idx]    = L6_Ku*(1+std_Ku[idx]*DEL_Hk*dHk6_scale);
					host_Aex[idx]   = L6_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L6_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_6; k < mNz; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = 1.0;
					host_Ku[idx]    = 0.0;
					host_Aex[idx]   = 0.0;
					host_alpha[idx] = 1.0;
					D[idx] = 0.0;
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
			}
		}


		if (CGC_DEF){
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					switch (CGC_label){
					case 1:
						for (int k = 0; k < mNz_1; k++){
							int idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 2:
						for (int k = mNz_12; k < mNz_2; k++){
							int idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 3:
						for (int k = mNz_23; k < mNz_3; k++){
							int idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 4:
						for (int k = mNz_34; k < mNz_4; k++){
							int idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 5:
						for (int k = mNz_45; k < mNz_5; k++){
							int idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 6:
						for (int k = mNz_56; k < mNz_6; k++){
							int idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					}
				}
			}
		}
		if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ku, "host_Ku.out")) printf("Output_Float_3D_Format() failed!\n");

	}

	if (t == 0){
		cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_x, Happl_x,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_y, Happl_y,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_z, Happl_z,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	}

	/*if (THERMAL){
		cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	}*/


	return true;
}

bool CT_Analysis_Head(int t){
	
	int idx;
	if (t == 0){

		// Temperature Profile
		for (int i = 0; i < mNx; i++){
			for (int j = 0; j < mNy; j++){
				for (int k = 0; k < mNz; k++){
					idx = i + j*mNx + k*mNx*mNy;
					T[idx] = Temperature;
				}
			}
		}
		// Temperature-dependent input parameters for 2nm cells (h: high temperature, l: low temperature)
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				for (int k = 0; k < mNz_1; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L1_Ms;
					host_Ku[idx]    = L1_Ku*(1+std_Ku[idx]*DEL_Hk*dHk1_scale);  
					host_Aex[idx]   = L1_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L1_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_1; k < mNz_12; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL12_Ms;
					host_Ku[idx]    = 0.;  
					host_Aex[idx]   = L1_Aex*((double)BL12*(1+std_Aex[idx]*dBL12_scale));
					host_alpha[idx] = L1_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_12; k < mNz_2; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L2_Ms;
					host_Ku[idx]    = L2_Ku*(1+std_Ku[idx]*DEL_Hk*dHk2_scale);  
					host_Aex[idx]   = L2_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L2_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_2; k < mNz_23; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL23_Ms;
					host_Ku[idx]    = 0.;  
					host_Aex[idx]   = L2_Aex*((double)BL23*(1+std_Aex[idx]*dBL23_scale));
					host_alpha[idx] = L2_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_23; k < mNz_3; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L3_Ms;
					host_Ku[idx]    = L3_Ku*(1+std_Ku[idx]*DEL_Hk*dHk3_scale);  
					host_Aex[idx]   = L3_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L3_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_3; k < mNz_34; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL34_Ms;
					host_Ku[idx]    = 0.;  
					host_Aex[idx]   = L3_Aex*((double)BL34*(1+std_Aex[idx]*dBL34_scale));
					host_alpha[idx] = L3_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_34; k < mNz_4; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L4_Ms;
					host_Ku[idx]    = L4_Ku*(1+std_Ku[idx]*DEL_Hk*dHk4_scale);  
					host_Aex[idx]   = L4_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L4_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_4; k < mNz_45; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL45_Ms;
					host_Ku[idx]    = 0.;  
					host_Aex[idx]   = L4_Aex*((double)BL45*(1+std_Aex[idx]*dBL45_scale));
					host_alpha[idx] = L4_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_45; k < mNz_5; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L5_Ms;
					host_Ku[idx]    = L5_Ku*(1+std_Ku[idx]*DEL_Hk*dHk5_scale);
					host_Aex[idx]   = L5_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L5_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_5; k < mNz_56; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = BL56_Ms;
					host_Ku[idx]    = 0.;
					host_Aex[idx]   = L5_Aex*((double)BL56*(1+std_Aex[idx]*dBL56_scale));
					host_alpha[idx] = L5_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_56; k < mNz_6; k++){
					idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = L6_Ms;
					host_Ku[idx]    = L6_Ku*(1+std_Ku[idx]*DEL_Hk*dHk6_scale);
					host_Aex[idx]   = L6_Aex*(1+std_Aex[idx]*DEL_Aex)*1;
					host_alpha[idx] = L6_alpha;
					D[idx] = (2*kb*T[idx]*host_alpha[idx])/(host_Ms[idx]*delta_t*GAMMA*(delta_x*delta_y*delta_z));
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
				for (int k = mNz_6; k < mNz; k++){
					int idx = i + j*mNx + k*mNx*mNy;
					host_Ms[idx]    = 1.0;
					host_Ku[idx]    = 0.0;
					host_Aex[idx]   = 0.0;
					host_alpha[idx] = 1.0;
					D[idx] = 0.0;
					host_gamma[idx] = (1.76e7)/(1+pow(host_alpha[idx], 2.0));
				}
			}
		}

		if (CGC_DEF){
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					switch (CGC_label){
					case 1:
						for (int k = 0; k < mNz_1; k++){
							idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 2:
						for (int k = mNz_12; k < mNz_2; k++){
							idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 3:
						for (int k = mNz_23; k < mNz_3; k++){
							idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 4:
						for (int k = mNz_34; k < mNz_4; k++){
							idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 5:
						for (int k = mNz_45; k < mNz_5; k++){
							idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					case 6:
						for (int k = mNz_56; k < mNz_6; k++){
							idx = i + j*mNx + k*mNx*mNy;
							if (indicator7[idx] == 1){
								host_Ku[idx]    = host_Ku[idx]*0.1;
							}
						}
						break;
					}
				}
			}
		}
		if (!Output_Float_3D_Format(mNx, mNy, mNz, host_Ku, "host_Ku.out")) { printf("Output_Float_3D_Format() failed!\n"); }

	}

	if (t == 0){
		cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_x, Happl_x,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_y, Happl_y,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Happl_z, Happl_z,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	}

	/*if (THERMAL){
		cudaMemcpy(dev_Ms, host_Ms,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Ku, host_Ku,       (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_Aex, host_Aex,     (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_alpha, host_alpha, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_gamma, host_gamma, (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_D, D,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_T, T,			  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	}*/


	return true;
}
