#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <algorithm>
#include <vector>
#include "MiscFunctions.h"
using namespace std;

bool CalcMagLayerAve(long int tt, int count, double* Mz_CGC, vector<vector<coor2d_t> >& grain_coor, int NumOfGrains);
bool External_HeadField(int EXT_FP_PROFILE);



bool CalcMagLayerAve(long int tt, int count, double* Mz_CGC, vector<vector<coor2d_t> >& grain_coor, int NumOfGrains){
	
	#ifndef __UNIFORM_GRAIN_MEDIA__
	if (!VORO_GRAIN){
		int GrainNumAveFactor1, GrainNumAveFactor2, GrainNumAveFactor3, GrainNumAveFactor4, GrainNumAveFactor5, GrainNumAveFactor6;
		int idx;
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){
				for (int k = 0; k < mNz_1; k++){
					if (CGC_label==1) { GrainNumAveFactor1 = count; } else { GrainNumAveFactor1 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar1[tt] = Mx_bar1[tt] + Mx[idx]/(GrainNumAveFactor1*mNz_1);
						My_bar1[tt] = My_bar1[tt] + My[idx]/(GrainNumAveFactor1*mNz_1);
						Mz_bar1[tt] = Mz_bar1[tt] + Mz[idx]/(GrainNumAveFactor1*mNz_1);
					}
				}
				for (int k = mNz_12; k < mNz_2; k++){
					if (CGC_label==2) { GrainNumAveFactor2 = count; } else { GrainNumAveFactor2 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar2[tt] = Mx_bar2[tt] + Mx[idx]/(GrainNumAveFactor2*(mNz_2-mNz_12));
						My_bar2[tt] = My_bar2[tt] + My[idx]/(GrainNumAveFactor2*(mNz_2-mNz_12));
						Mz_bar2[tt] = Mz_bar2[tt] + Mz[idx]/(GrainNumAveFactor2*(mNz_2-mNz_12));
					}
				}
				for (int k = mNz_23; k < mNz_3; k++){
					if (CGC_label==3) { GrainNumAveFactor3 = count; } else { GrainNumAveFactor3 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar3[tt] = Mx_bar3[tt] + Mx[idx]/(GrainNumAveFactor3*(mNz_3-mNz_23));
						My_bar3[tt] = My_bar3[tt] + My[idx]/(GrainNumAveFactor3*(mNz_3-mNz_23));
						Mz_bar3[tt] = Mz_bar3[tt] + Mz[idx]/(GrainNumAveFactor3*(mNz_3-mNz_23));
					}
				}
				for (int k = mNz_34; k < mNz_4; k++){
					if (CGC_label==4) { GrainNumAveFactor4 = count; } else { GrainNumAveFactor4 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar4[tt] = Mx_bar4[tt] + Mx[idx]/(GrainNumAveFactor4*(mNz_4-mNz_34));
						My_bar4[tt] = My_bar4[tt] + My[idx]/(GrainNumAveFactor4*(mNz_4-mNz_34));
						Mz_bar4[tt] = Mz_bar4[tt] + Mz[idx]/(GrainNumAveFactor4*(mNz_4-mNz_34));
					}
				}
				for (int k = mNz_45; k < mNz_5; k++){
					if (CGC_label==5) { GrainNumAveFactor5 = count; } else { GrainNumAveFactor5 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar5[tt] = Mx_bar5[tt] + Mx[idx]/(GrainNumAveFactor5*(mNz_5-mNz_45));
						My_bar5[tt] = My_bar5[tt] + My[idx]/(GrainNumAveFactor5*(mNz_5-mNz_45));
						Mz_bar5[tt] = Mz_bar5[tt] + Mz[idx]/(GrainNumAveFactor5*(mNz_5-mNz_45));
					}
				}
				for (int k = mNz_56; k < mNz_6; k++){
					if (CGC_label==6) { GrainNumAveFactor6 = count; } else { GrainNumAveFactor6 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar6[tt] = Mx_bar6[tt] + Mx[idx]/(GrainNumAveFactor6*(mNz_6-mNz_56));
						My_bar6[tt] = My_bar6[tt] + My[idx]/(GrainNumAveFactor6*(mNz_6-mNz_56));
						Mz_bar6[tt] = Mz_bar6[tt] + Mz[idx]/(GrainNumAveFactor6*(mNz_6-mNz_56));
					}
				}
				for (int k = 0; k < mNz_1; k++){   //important
					if (CGC_label==1) { GrainNumAveFactor1 = count; } else { GrainNumAveFactor1 = mNx*mNy; }
					idx = i + j*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar[tt] = Mx_bar[tt] + Mx[idx]/(GrainNumAveFactor1*(mNz_1));  //important, because we are more interested in sigHs
						My_bar[tt] = My_bar[tt] + My[idx]/(GrainNumAveFactor1*(mNz_1));  //important, because we are more interested in sigHs
						Mz_bar[tt] = Mz_bar[tt] + Mz[idx]/(GrainNumAveFactor1*(mNz_1));  //important, because we are more interested in sigHs
					}
				}
			}
		}
	}
	#endif

	#ifndef __VORO_GRAIN_MEDIA__
	if (VORO_GRAIN){
		int GrainNumAveFactor;
		int idx;
		
		for (int k = 0; k < mNz_1; k++){
			if (CGC_label==1) GrainNumAveFactor = count; else GrainNumAveFactor = NumOfGrains;
			for (int i=0; i<NumOfGrains; i++){
				double Mx_bar_tmp = 0.0;
				double My_bar_tmp = 0.0;
				double Mz_bar_tmp = 0.0;
				for (int ii=0; ii<grain_coor[i].size(); ++ii){
					idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar_tmp = Mx_bar_tmp + Mx[idx]/(grain_coor[i].size());
						My_bar_tmp = My_bar_tmp + My[idx]/(grain_coor[i].size());
						Mz_bar_tmp = Mz_bar_tmp + Mz[idx]/(grain_coor[i].size());
					}
				}
				Mx_bar1[tt] = Mx_bar1[tt] + Mx_bar_tmp/(GrainNumAveFactor*mNz_1);
				My_bar1[tt] = My_bar1[tt] + My_bar_tmp/(GrainNumAveFactor*mNz_1);
				Mz_bar1[tt] = Mz_bar1[tt] + Mz_bar_tmp/(GrainNumAveFactor*mNz_1);
			}
		}

		for (int k = mNz_12; k < mNz_2; k++){
			if (CGC_label==2) GrainNumAveFactor = count; else GrainNumAveFactor = NumOfGrains;
			for (int i=0; i<NumOfGrains; i++){
				double Mx_bar_tmp = 0.0;
				double My_bar_tmp = 0.0;
				double Mz_bar_tmp = 0.0;
				for (int ii=0; ii<grain_coor[i].size(); ++ii){
					idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar_tmp = Mx_bar_tmp + Mx[idx]/(grain_coor[i].size());
						My_bar_tmp = My_bar_tmp + My[idx]/(grain_coor[i].size());
						Mz_bar_tmp = Mz_bar_tmp + Mz[idx]/(grain_coor[i].size());
					}
				}
				Mx_bar2[tt] = Mx_bar2[tt] + Mx_bar_tmp/(GrainNumAveFactor*(mNz_2-mNz_12));
				My_bar2[tt] = My_bar2[tt] + My_bar_tmp/(GrainNumAveFactor*(mNz_2-mNz_12));
				Mz_bar2[tt] = Mz_bar2[tt] + Mz_bar_tmp/(GrainNumAveFactor*(mNz_2-mNz_12));
			}
		}

		for (int k = mNz_23; k < mNz_3; k++){
			if (CGC_label==3) GrainNumAveFactor = count; else GrainNumAveFactor = NumOfGrains;
			for (int i=0; i<NumOfGrains; i++){
				double Mx_bar_tmp = 0.0;
				double My_bar_tmp = 0.0;
				double Mz_bar_tmp = 0.0;
				for (int ii=0; ii<grain_coor[i].size(); ++ii){
					idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar_tmp = Mx_bar_tmp + Mx[idx]/(grain_coor[i].size());
						My_bar_tmp = My_bar_tmp + My[idx]/(grain_coor[i].size());
						Mz_bar_tmp = Mz_bar_tmp + Mz[idx]/(grain_coor[i].size());
					}
				}
				Mx_bar3[tt] = Mx_bar3[tt] + Mx_bar_tmp/(GrainNumAveFactor*(mNz_3-mNz_23));
				My_bar3[tt] = My_bar3[tt] + My_bar_tmp/(GrainNumAveFactor*(mNz_3-mNz_23));
				Mz_bar3[tt] = Mz_bar3[tt] + Mz_bar_tmp/(GrainNumAveFactor*(mNz_3-mNz_23));
			}
		}
						
		for (int k = mNz_34; k < mNz_4; k++){
			if (CGC_label==4) GrainNumAveFactor = count; else GrainNumAveFactor = NumOfGrains;
			for (int i=0; i<NumOfGrains; i++){
				double Mx_bar_tmp = 0.0;
				double My_bar_tmp = 0.0;
				double Mz_bar_tmp = 0.0;
				for (int ii=0; ii<grain_coor[i].size(); ++ii){
					idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar_tmp = Mx_bar_tmp + Mx[idx]/(grain_coor[i].size());
						My_bar_tmp = My_bar_tmp + My[idx]/(grain_coor[i].size());
						Mz_bar_tmp = Mz_bar_tmp + Mz[idx]/(grain_coor[i].size());
					}
				}
				Mx_bar4[tt] = Mx_bar4[tt] + Mx_bar_tmp/(GrainNumAveFactor*(mNz_4-mNz_34));
				My_bar4[tt] = My_bar4[tt] + My_bar_tmp/(GrainNumAveFactor*(mNz_4-mNz_34));
				Mz_bar4[tt] = Mz_bar4[tt] + Mz_bar_tmp/(GrainNumAveFactor*(mNz_4-mNz_34));
			}
		}
			
		for (int k = mNz_45; k < mNz_5; k++){
			if (CGC_label==5) GrainNumAveFactor = count; else GrainNumAveFactor = NumOfGrains;
			for (int i=0; i<NumOfGrains; i++){
				double Mx_bar_tmp = 0.0;
				double My_bar_tmp = 0.0;
				double Mz_bar_tmp = 0.0;
				for (int ii=0; ii<grain_coor[i].size(); ++ii){
					idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar_tmp = Mx_bar_tmp + Mx[idx]/(grain_coor[i].size());
						My_bar_tmp = My_bar_tmp + My[idx]/(grain_coor[i].size());
						Mz_bar_tmp = Mz_bar_tmp + Mz[idx]/(grain_coor[i].size());
					}
				}
				Mx_bar5[tt] = Mx_bar5[tt] + Mx_bar_tmp/(GrainNumAveFactor*(mNz_5-mNz_45));
				My_bar5[tt] = My_bar5[tt] + My_bar_tmp/(GrainNumAveFactor*(mNz_5-mNz_45));
				Mz_bar5[tt] = Mz_bar5[tt] + Mz_bar_tmp/(GrainNumAveFactor*(mNz_5-mNz_45));
			}
		}

		for (int k = mNz_56; k < mNz_6; k++){
			if (CGC_label==6) GrainNumAveFactor = count; else GrainNumAveFactor = NumOfGrains;
			for (int i=0; i<NumOfGrains; i++){
				double Mx_bar_tmp = 0.0;
				double My_bar_tmp = 0.0;
				double Mz_bar_tmp = 0.0;
				for (int ii=0; ii<grain_coor[i].size(); ++ii){
					idx = grain_coor[i][ii].x + grain_coor[i][ii].y*mNx + k*mNx*mNy;
					if (indicator7[idx] == 0){
						Mx_bar_tmp = Mx_bar_tmp + Mx[idx]/(grain_coor[i].size());
						My_bar_tmp = My_bar_tmp + My[idx]/(grain_coor[i].size());
						Mz_bar_tmp = Mz_bar_tmp + Mz[idx]/(grain_coor[i].size());
					}
				}
				Mx_bar6[tt] = Mx_bar6[tt] + Mx_bar_tmp/(GrainNumAveFactor*(mNz_6-mNz_56));
				My_bar6[tt] = My_bar6[tt] + My_bar_tmp/(GrainNumAveFactor*(mNz_6-mNz_56));
				Mz_bar6[tt] = Mz_bar6[tt] + Mz_bar_tmp/(GrainNumAveFactor*(mNz_6-mNz_56));
			}
		}

		Mx_bar[tt] = Mx_bar1[tt];
		My_bar[tt] = My_bar1[tt];
		Mz_bar[tt] = Mz_bar1[tt];

	}
	#endif

	switch (CGC_label){
		case 1: Mz_CGC[tt-1] = Mz_bar1[tt-1]; break;
		case 2: Mz_CGC[tt-1] = Mz_bar2[tt-1]; break;
		case 3: Mz_CGC[tt-1] = Mz_bar3[tt-1]; break;
		case 4: Mz_CGC[tt-1] = Mz_bar4[tt-1]; break;
		case 5: Mz_CGC[tt-1] = Mz_bar5[tt-1]; break;
		case 6: Mz_CGC[tt-1] = Mz_bar6[tt-1]; break;
	}
	return true;
}

bool External_HeadField(int EXT_FP_PROFILE){
	int fNx_end, fNy_end, fNz_end, idx;
	if (EXT_FP_PROFILE){
		if (fNx0+fNx > mNx){fNx_end = mNx;}	else {fNx_end = fNx0+fNx;}
		if (fNy0+fNy > mNy){fNy_end = mNy;}	else {fNy_end = fNy0+fNy;}
		if (fNz0+fNz > mNz){fNz_end = mNz;}	else {fNz_end = fNz0+fNz;}
		for (int k = fNz0; k < (fNz0+fNz_end); k++){
			for (int j = fNy0; j < (fNy0+fNy_end); j++){
				for (int i = fNx0; i < (fNx0+fNx_end); i++){
					FP[0][(i-fNx0)+(j-fNy0)*mNx+(k-fNz0)*mNx*mNy] = FP_inp[0][(i)+(j)*fNx+(k)*fNx*fNy];
					FP[1][(i-fNx0)+(j-fNy0)*mNx+(k-fNz0)*mNx*mNy] = FP_inp[1][(i)+(j)*fNx+(k)*fNx*fNy];
					FP[2][(i-fNx0)+(j-fNy0)*mNx+(k-fNz0)*mNx*mNy] = FP_inp[2][(i)+(j)*fNx+(k)*fNx*fNy];
					FP_trail[0][(i-fNx0)+(j-fNy0)*mNx+(k-fNz0)*mNx*mNy] = FP_inp[0][(i-mNx)+(j)*fNx+(k)*fNx*fNy];
					FP_trail[1][(i-fNx0)+(j-fNy0)*mNx+(k-fNz0)*mNx*mNy] = FP_inp[1][(i-mNx)+(j)*fNx+(k)*fNx*fNy];
					FP_trail[2][(i-fNx0)+(j-fNy0)*mNx+(k-fNz0)*mNx*mNy] = FP_inp[2][(i-mNx)+(j)*fNx+(k)*fNx*fNy];
				}
			}
		}
		for (int k = 0; k < mNz; k++){
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){ 
					idx = i+j*mNx+k*mNx*mNy;
					Happl_x[idx] = 0.;
					Happl_y[idx] = 0.;
					Happl_z[idx] = 0.;
					Happl_x_temp[idx] = 0.;
					Happl_y_temp[idx] = 0.;
					Happl_z_temp[idx] = 0.;
				}
			}
		}
	}
	return true;
}

bool Set_delta_K_Aex_Tc(double *std_Ku, double *std_Aex, double *std_Tc, double *indicator4, double *indicator5, double *indicator6){
	
	for (int idx = 0; idx < mNx*mNy*mNz; idx++){
		std_Ku[idx] = indicator4[idx];
		/*fprintf(wfile5, "%12.3lf \n", std_Ku[idx]);*/
	}
	// Define delta Aex
	for (int idx = 0; idx < mNx*mNy*mNz; idx++){
		std_Aex[idx] = indicator5[idx];
		/*fprintf(wfile5, "%12.3lf \n", std_Aex[idx]);*/
	}
	// Define delta Tc
	for (int idx = 0; idx < mNx*mNy*mNz; idx++){
		std_Tc[idx] = indicator6[idx];
		/*fprintf(wfile5, "%12.3lf \n", std_Tc[idx]);*/
	}

	return true;
}

bool Set_Happl_DT(){
	int idx, idxx;
	for (int k = 0; k < mNz; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){ 
				idx = i+j*mNx+k*mNx*mNy;
				idxx = (DT_x0+i)+((DT_y0)*fNx+k*fNx*fNy);
				Happl_x[idx] = -FP_inp[0][idxx];
				Happl_y[idx] = -FP_inp[1][idxx];
				Happl_z[idx] = -FP_inp[2][idxx];
			}
		}
	}

	return true;
}

bool Set_Happl_CT(){
	int idx, idxx;
	for (int k = 0; k < mNz; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){ 
				idx = i+j*mNx+k*mNx*mNy;
				idxx = (CT_x0)+((CT_y0+j)*fNx+k*fNx*fNy);
				Happl_x[idx] = -FP_inp[0][idxx];
				Happl_y[idx] = -FP_inp[1][idxx];
				Happl_z[idx] = -FP_inp[2][idxx];
			}
		}
	}

	return true;
}

bool Reset_Happl(){
	int idx;
	for (int k = 0; k < mNz; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){ 
				idx = i+j*mNx+k*mNx*mNy;
				Happl_x[idx] = 0.;
				Happl_y[idx] = 0.;
				Happl_z[idx] = 0.;
			}
		}
	}

	return true;
}

bool Jitter_Calc(){
	
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);

	int idx, i;
	int *M1_z_DT = NULL; 
	double M1_z_DT_mean = 0, Mz_zDT_std = 0;
	//double *Mz_z_test = NULL;
	M1_z_DT = (int*) calloc(mNy, sizeof(int));
	//Mz_z_test = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	
	//Input_Float_1D_Format_1col(mNx, mNy, mNz, Mz_z_test, "Mz.test");
	// Looking at M1 magnetization (k = 0)
	for (int j = 0; j < mNy; j++){
		i = 0;
		while (i < mNx){
			if (Mz[i + j*mNx]*Mz[(i+1) + j*mNx] < 0){
				M1_z_DT[j] = i;
				i = mNx;
			}
			i++;
		}
	}
	Output_Int_1D_Format_1col(1, mNy, 1, M1_z_DT, "M1_z_DT.out");
	
	for (int j = 0; j < mNy; j++){
		M1_z_DT_mean = M1_z_DT_mean + double(M1_z_DT[j])/(mNy);
	}
	for (int j = 0; j < mNy; j++){
		Mz_zDT_std = Mz_zDT_std + pow(double(M1_z_DT[j])-M1_z_DT_mean, 2);
	}
	Mz_zDT_std = pow(Mz_zDT_std/mNy, 0.5)/M1_z_DT_mean;
	
	printf("\n---------------------------\n");
	printf("Trans Pos= %12.5f [au]\n", M1_z_DT_mean);
	printf("jitter=    %12.5f [au]\n", Mz_zDT_std);
	printf("---------------------------\n");
	printf("BL1=     %12.5f\n", BL12);
	printf("BL2=     %12.5f\n", BL23);
	printf("BL3=     %12.5f\n", BL34);
	printf("BL4=     %12.5f\n", BL45);
	printf("BL5=     %12.5f\n", BL56);
	printf("M1-Hk=   %12.5f [kOe]\n", L1_Hk);
	printf("M2-Hk=   %12.5f [kOe]\n", L2_Hk);
	printf("M3-Hk=   %12.5f [kOe]\n", L3_Hk);
	printf("M4-Hk=   %12.5f [kOe]\n", L4_Hk);
	printf("M5-Hk=   %12.5f [kOe]\n", L5_Hk);
	printf("M6-Hk=   %12.5f [kOe]\n", L6_Hk);
	printf("M1-thk=  %12d [nm]\n", mNz_1);
	printf("M2-thk=  %12d [nm]\n", mNz_2-mNz_12);
	printf("M3-thk=  %12d [nm]\n", mNz_3-mNz_23);
	printf("M4-thk=  %12d [nm]\n", mNz_4-mNz_34);
	printf("M5-thk=  %12d [nm]\n", mNz_5-mNz_45);
	printf("M5-thk=  %12d [nm]\n", mNz_6-mNz_56);
	printf("M1-Ms=   %12.5f [emu/cc]\n", L1_Ms);
	printf("M2-Ms=   %12.5f [emu/cc]\n", L2_Ms);
	printf("M3-Ms=   %12.5f [emu/cc]\n", L3_Ms);
	printf("M4-Ms=   %12.5f [emu/cc]\n", L4_Ms);
	printf("M5-Ms=   %12.5f [emu/cc]\n", L5_Ms);
	printf("M6-Ms=   %12.5f [emu/cc]\n", L6_Ms);
	printf("L1_Hex_l=%12.5f\n", L1_Hex_l);
	printf("L2_Hex_l=%12.5f\n", L2_Hex_l);
	printf("L3_Hex_l=%12.5f\n", L3_Hex_l);
	printf("L4_Hex_l=%12.5f\n", L4_Hex_l);
	printf("L5_Hex_l=%12.5f\n", L5_Hex_l);
	printf("L6_Hex_l=%12.5f\n", L6_Hex_l);
	printf("---------------------------\n");

	FILE *output1;
	timeinfo = localtime(&rawtime);
	output1 = fopen("output_Jitter.out", "w");
	fprintf(output1,"%s", asctime(timeinfo));
	fprintf(output1,"\n---------------------------\n");
	fprintf(output1,"Trans Pos= %12.5f [au]\n", M1_z_DT_mean);
	fprintf(output1,"jitter=    %12.5f [au]\n", Mz_zDT_std);
	fprintf(output1,"---------------------------\n");
	fprintf(output1,"BL1=     %12.5f\n", BL12);
	fprintf(output1,"BL2=     %12.5f\n", BL23);
	fprintf(output1,"BL3=     %12.5f\n", BL34);
	fprintf(output1,"BL4=     %12.5f\n", BL45);
	fprintf(output1,"BL5=     %12.5f\n", BL56);
	fprintf(output1,"M1-Hk=   %12.5f [kOe]\n", L1_Hk);
	fprintf(output1,"M2-Hk=   %12.5f [kOe]\n", L2_Hk);
	fprintf(output1,"M3-Hk=   %12.5f [kOe]\n", L3_Hk);
	fprintf(output1,"M4-Hk=   %12.5f [kOe]\n", L4_Hk);
	fprintf(output1,"M5-Hk=   %12.5f [kOe]\n", L5_Hk);
	fprintf(output1,"M6-Hk=   %12.5f [kOe]\n", L6_Hk);
	fprintf(output1,"M1-thk=  %12d [nm]\n", mNz_1);
	fprintf(output1,"M2-thk=  %12d [nm]\n", mNz_2-mNz_12);
	fprintf(output1,"M3-thk=  %12d [nm]\n", mNz_3-mNz_23);
	fprintf(output1,"M4-thk=  %12d [nm]\n", mNz_4-mNz_34);
	fprintf(output1,"M5-thk=  %12d [nm]\n", mNz_5-mNz_45);
	fprintf(output1,"M5-thk=  %12d [nm]\n", mNz_6-mNz_56);
	fprintf(output1,"M1-Ms=   %12.5f [emu/cc]\n", L1_Ms);
	fprintf(output1,"M2-Ms=   %12.5f [emu/cc]\n", L2_Ms);
	fprintf(output1,"M3-Ms=   %12.5f [emu/cc]\n", L3_Ms);
	fprintf(output1,"M4-Ms=   %12.5f [emu/cc]\n", L4_Ms);
	fprintf(output1,"M5-Ms=   %12.5f [emu/cc]\n", L5_Ms);
	fprintf(output1,"M6-Ms=   %12.5f [emu/cc]\n", L6_Ms);
	fprintf(output1,"L1_Hex_l=%12.5f\n", L1_Hex_l);
	fprintf(output1,"L2_Hex_l=%12.5f\n", L2_Hex_l);
	fprintf(output1,"L3_Hex_l=%12.5f\n", L3_Hex_l);
	fprintf(output1,"L4_Hex_l=%12.5f\n", L4_Hex_l);
	fprintf(output1,"L5_Hex_l=%12.5f\n", L5_Hex_l);
	fprintf(output1,"L6_Hex_l=%12.5f\n", L6_Hex_l);
	fprintf(output1,"---------------------------\n");
	
	fclose(output1);

	return true;
}

bool WPE_Calc(){
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	
	int idx, i;
	double *M1_z_CT = NULL; 
	double M1_z_CT_mean = 0, Mz_zCT_std = 0, WPE = 0;
	M1_z_CT = (double*) calloc(mNy, sizeof(double));
	

	// Looking at M1 magnetization (k = 0)
	for (int j = 0; j < mNy; j++){
		for (int i = 0; i < mNx; i++){
			M1_z_CT[j] = M1_z_CT[j] + Mz[i + j*mNx]/double(mNx);
		}
	}
	Output_Float_1D_Format_1col(1, mNy, 1, M1_z_CT, "M1_z_CT.out");
	
	i = 0;
	while (!(M1_z_CT[i]*M1_z_CT[i+1] <= 0)){
		i++;
	}
	
	
	WPE = double(i) + (-M1_z_CT[i])*(-1.)/(M1_z_CT[i]-M1_z_CT[i+1]);
	printf("\n---------------------------\n");
	printf("WPE=     %12.5lf [au]\n", WPE);
	printf("---------------------------\n");
	printf("BL1=     %12.5f\n", BL12);
	printf("BL2=     %12.5f\n", BL23);
	printf("BL3=     %12.5f\n", BL34);
	printf("BL4=     %12.5f\n", BL45);
	printf("BL5=     %12.5f\n", BL56);
	printf("M1-Hk=   %12.5f [kOe]\n", L1_Hk);
	printf("M2-Hk=   %12.5f [kOe]\n", L2_Hk);
	printf("M3-Hk=   %12.5f [kOe]\n", L3_Hk);
	printf("M4-Hk=   %12.5f [kOe]\n", L4_Hk);
	printf("M5-Hk=   %12.5f [kOe]\n", L5_Hk);
	printf("M6-Hk=   %12.5f [kOe]\n", L6_Hk);
	printf("M1-thk=  %12d [nm]\n", mNz_1);
	printf("M2-thk=  %12d [nm]\n", mNz_2-mNz_12);
	printf("M3-thk=  %12d [nm]\n", mNz_3-mNz_23);
	printf("M4-thk=  %12d [nm]\n", mNz_4-mNz_34);
	printf("M5-thk=  %12d [nm]\n", mNz_5-mNz_45);
	printf("M5-thk=  %12d [nm]\n", mNz_6-mNz_56);
	printf("M1-Ms=   %12.5f [emu/cc]\n", L1_Ms);
	printf("M2-Ms=   %12.5f [emu/cc]\n", L2_Ms);
	printf("M3-Ms=   %12.5f [emu/cc]\n", L3_Ms);
	printf("M4-Ms=   %12.5f [emu/cc]\n", L4_Ms);
	printf("M5-Ms=   %12.5f [emu/cc]\n", L5_Ms);
	printf("M6-Ms=   %12.5f [emu/cc]\n", L6_Ms);
	printf("L1_Hex_l=%12.5f\n", L1_Hex_l);
	printf("L2_Hex_l=%12.5f\n", L2_Hex_l);
	printf("L3_Hex_l=%12.5f\n", L3_Hex_l);
	printf("L4_Hex_l=%12.5f\n", L4_Hex_l);
	printf("L5_Hex_l=%12.5f\n", L5_Hex_l);
	printf("L6_Hex_l=%12.5f\n", L6_Hex_l);
	printf("---------------------------\n");

	FILE *output1;
	timeinfo = localtime(&rawtime);
	output1 = fopen("output_WPE.out", "w");
	fprintf(output1, "%s", asctime(timeinfo));
	fprintf(output1, "\n---------------------------\n");
	fprintf(output1, "WPE=     %12.5lf [au]\n", WPE);
	fprintf(output1, "---------------------------\n");
	fprintf(output1,"BL1=     %12.5f\n", BL12);
	fprintf(output1,"BL2=     %12.5f\n", BL23);
	fprintf(output1,"BL3=     %12.5f\n", BL34);
	fprintf(output1,"BL4=     %12.5f\n", BL45);
	fprintf(output1,"BL5=     %12.5f\n", BL56);
	fprintf(output1,"M1-Hk=   %12.5f [kOe]\n", L1_Hk);
	fprintf(output1,"M2-Hk=   %12.5f [kOe]\n", L2_Hk);
	fprintf(output1,"M3-Hk=   %12.5f [kOe]\n", L3_Hk);
	fprintf(output1,"M4-Hk=   %12.5f [kOe]\n", L4_Hk);
	fprintf(output1,"M5-Hk=   %12.5f [kOe]\n", L5_Hk);
	fprintf(output1,"M6-Hk=   %12.5f [kOe]\n", L6_Hk);
	fprintf(output1,"M1-thk=  %12d [nm]\n", mNz_1);
	fprintf(output1,"M2-thk=  %12d [nm]\n", mNz_2-mNz_12);
	fprintf(output1,"M3-thk=  %12d [nm]\n", mNz_3-mNz_23);
	fprintf(output1,"M4-thk=  %12d [nm]\n", mNz_4-mNz_34);
	fprintf(output1,"M5-thk=  %12d [nm]\n", mNz_5-mNz_45);
	fprintf(output1,"M5-thk=  %12d [nm]\n", mNz_6-mNz_56);
	fprintf(output1,"M1-Ms=   %12.5f [emu/cc]\n", L1_Ms);
	fprintf(output1,"M2-Ms=   %12.5f [emu/cc]\n", L2_Ms);
	fprintf(output1,"M3-Ms=   %12.5f [emu/cc]\n", L3_Ms);
	fprintf(output1,"M4-Ms=   %12.5f [emu/cc]\n", L4_Ms);
	fprintf(output1,"M5-Ms=   %12.5f [emu/cc]\n", L5_Ms);
	fprintf(output1,"M6-Ms=   %12.5f [emu/cc]\n", L6_Ms);
	fprintf(output1,"L1_Hex_l=%12.5f\n", L1_Hex_l);
	fprintf(output1,"L2_Hex_l=%12.5f\n", L2_Hex_l);
	fprintf(output1,"L3_Hex_l=%12.5f\n", L3_Hex_l);
	fprintf(output1,"L4_Hex_l=%12.5f\n", L4_Hex_l);
	fprintf(output1,"L5_Hex_l=%12.5f\n", L5_Hex_l);
	fprintf(output1,"L6_Hex_l=%12.5f\n", L6_Hex_l);
	fprintf(output1,"---------------------------\n");
	fclose(output1);


	return true;
}

bool sigHc_Calc(vector <vector<double> > & M1_z_SingleGrain_field, double* Happl1_sweep, int count, int NumOfGrains, vector<vector<coor2d_t> > & grain_coor){
	
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);

	#ifndef __UNIFORM_GRAIN_MEDIA__
	if (!VORO_GRAIN){
		
		// Find Hn of the bottom-most layer
		double *M1_Hn_gbyg = (double*) calloc(mNx*mNy, sizeof(double));
		double M1_Hn_mean = 0.0;
	
		// Find Hc of the bottom-most layer
		std::vector<double> M1_Hc_gbyg(mNx*mNy, 0.0);

		for (int j=0; j<mNy; j++){
			for (int i=0; i<mNx; i++){
				int idxx = i + j*mNx;
				int itt1 = 1;
				while (M1_z_SingleGrain_field[idxx][itt1]*M1_z_SingleGrain_field[idxx][itt1-1] > 0 && itt1 < ceil(TOTAL_TIME/FieldSweepTimeStep)+1){
					M1_Hc_gbyg[idxx] = Happl1_sweep[itt1-1];
					itt1++;
				}
				int itt2 = 100;
			
				// Calculate Hn
				/*while (abs(M1_z_SingleGrain_field[idxx][itt2]-M1_z_SingleGrain_field[idxx][itt2-100]) > L1_Ms*0.01 && itt2 < ceil(TOTAL_TIME/FieldSweepTimeStep)+1){
					M1_Hn_gbyg[idxx] = Happl1_sweep[itt2-100];
					itt2++;
				}*/
			}
		}
		std::ofstream vfile1;
		vfile1.open("M1_Hc_gbyg.dat");
		for (int i=0; i<mNx*mNy; i++) vfile1 << M1_Hc_gbyg[i] << std::endl;
	
		// Hc mean
		double M1_Hc_mean = 0;
		for (int j=0; j < mNy; j++){
			for (int i=0; i < mNx; i++){
				if (indicator7[i+j*mNx+(mNz_6-1)*mNx*mNy] == 0){
					M1_Hc_mean = M1_Hc_mean + abs(M1_Hc_gbyg[i+j*mNx])/(count);
					//M1_Hn_mean = M1_Hn_mean + abs(M1_Hn_gbyg[i+j*mNx])/(count);
				}
			}
		}

		// Hc variance
		double M1_Hc_var = 0;
		for (int j=0; j < mNy; j++){
			for (int i=0; i < mNx; i++){
				if (indicator7[i+j*mNx+(mNz_6-1)*mNx*mNy] == 0){
					M1_Hc_var = M1_Hc_var + pow(abs(M1_Hc_gbyg[i+j*mNx])-M1_Hc_mean, 2);
				}
			}
		}

		// Hc std
		double M1_Hc_std = 0;
		M1_Hc_std = pow(M1_Hc_var/(count-1), 0.5);


		//--------------------------------------------------
		// Cluster size calculation
		//--------------------------------------------------
	
		// Locate the time where M1_Hc=0
		int t_Mz_is_0;
		for (int tt=FieldSweepTimeStep; tt<TOTAL_TIME; tt++){
			if (Mz_bar1[tt]*Mz_bar1[tt-FieldSweepTimeStep]<0) {
				t_Mz_is_0 = tt;
				break;
			}
		}

		// Set the center coordinates for the 2D grid
		std::vector<double> xcenter_SFD(mNx, 0.);
		std::vector<double> ycenter_SFD(mNy, 0.);

		for (int i_SFD=0; i_SFD<mNx; ++i_SFD) xcenter_SFD[i_SFD] = 0.5 + i_SFD;
		for (int j_SFD=0; j_SFD<mNy; ++j_SFD) ycenter_SFD[j_SFD] = 0.5 + j_SFD;
	
		// Estimate pair correlation length
		std::vector<double> corr_ij(mNx, 0.);

		double r = 0.0;
		int tott_SFD = 0;
		int tk = 0;

		for (tott_SFD=0; tott_SFD<mNx; tott_SFD++)
		{
			r = r + 1;
			for (int i_SFD=0; i_SFD<mNx; i_SFD++){
				for (int j_SFD=0; j_SFD<mNy; j_SFD++){
					if (indicator7[i_SFD+j_SFD*mNx+(mNz-1)*mNx*mNy] == 0){
						tk = 0;
					
						for (int ii_SFD=0; ii_SFD<mNx; ii_SFD++){
							for (int jj_SFD=0; jj_SFD<mNy; jj_SFD++){
								if (indicator7[ii_SFD+jj_SFD*mNx+(mNz-1)*mNx*mNy] == 0){
									double rij = pow( (pow(xcenter_SFD[ii_SFD]-xcenter_SFD[i_SFD],2.0) + pow(ycenter_SFD[jj_SFD]-ycenter_SFD[j_SFD],2.0)), 0.5 );
									if (rij>=r-1 && rij<=r){
										tk = tk + 1;
										corr_ij[tott_SFD] = corr_ij[tott_SFD] + Mz_tt[ii_SFD+jj_SFD*mNx][t_Mz_is_0]*Mz_tt[i_SFD+j_SFD*mNx][t_Mz_is_0];
									}
								}
							}
						}
	
					}
				}
			}
			corr_ij[tott_SFD]=corr_ij[tott_SFD]/(tk*2);
			rrt_cor[tott_SFD]=r;
		}

		// Write to file
		std::ofstream vfile2;
		vfile2.open("corr_func.dat");
		for (int i=0; i<mNx; i++) vfile2 << rrt_cor[i] << "\t" << corr_ij[i] <<std::endl;



		double cluster_size=0.0;
		int jsign;
		double max_NOR = corr_ij[0];
	
		for (int i_SFD = 0; i_SFD < mNx; i_SFD++)	corr_ij[i_SFD] = corr_ij[i_SFD]/max_NOR;
	
		jsign = 1;
		for(int i_SFD=0; i_SFD<mNx; i_SFD++){	
			if(corr_ij[i_SFD]>=0.3679 && corr_ij[i_SFD+1]<0.3679){
				double h_left=corr_ij[i_SFD]-0.3679;
				double h_right=0.3679-corr_ij[i_SFD+1];
				cluster_size = rrt_cor[i_SFD]+h_left/(h_left+h_right);
				jsign=2;
			}
			if(jsign == 1)	cluster_size = 0.0;
		}
		//--------------------------------------------//

		// Print to screen
		timeinfo = localtime(&rawtime);
		printf("\nM1_SFD%=%7.4f\n", M1_Hc_std/M1_Hc_mean);
		printf("M1_Hc_mean=%9.2f\n", M1_Hc_mean);
		printf("M1_Hc_std=%8.2f\n", M1_Hc_std);
		printf("M1_Hn_mean=%8.2f\n", M1_Hn_mean);
		printf("M1_CLUSTER_SIZE=%7.4f\n", cluster_size);
		printf("BL1=%9.5f\n", BL12);
		printf("BL2=%9.5f\n", BL23);
		printf("BL3=%9.5f\n", BL34);
		printf("BL4=%9.5f\n", BL45);
		printf("BL5=%9.5f\n", BL56);
		if (CGC_DEF){
			switch (CGC_label){
			case 1:
				printf("M1_Hex_l=%9.5f\n", L1_Hex_l);
				break;
			case 2:
				printf("M2_Hex_l=%9.5f\n", L2_Hex_l);
				break;
			case 3:
				printf("M3_Hex_l=%9.5f\n", L3_Hex_l);
				break;
			case 4:
				printf("M4_Hex_l=%9.5f\n", L4_Hex_l);
				break;
			case 5:
				printf("M5_Hex_l=%9.5f\n", L5_Hex_l);
				break;
			case 6:
				printf("M6_Hex_l=%9.5f\n", L6_Hex_l);
				break;
			}
		}
	
		cout << asctime(timeinfo) << std::endl;
	
		// Write to file
		FILE *output1;
		output1 = fopen("output.out", "w");
		fprintf(output1,"%s", asctime(timeinfo));
		fprintf(output1,"------------------------\n");
		fprintf(output1,"\nM1_SFD%=%7.4f\n", M1_Hc_std/M1_Hc_mean);
		fprintf(output1,"M1_Hc_mean=%9.2f\n", M1_Hc_mean);
		fprintf(output1,"M1_Hc_std=%8.2f\n", M1_Hc_std);
		fprintf(output1,"M1_Hn_mean=%8.2f\n", M1_Hn_mean);
		fprintf(output1, "M1_CLUSTER_SIZE=%7.4f\n", cluster_size);
		fprintf(output1,"BL1=%9.5f\n", BL12);
		fprintf(output1,"BL2=%9.5f\n", BL23);
		fprintf(output1,"BL3=%9.5f\n", BL34);
		fprintf(output1,"BL4=%9.5f\n", BL45);
		fprintf(output1,"BL5=%9.5f\n", BL56);
		if (CGC_DEF){
			switch (CGC_label){
			case 1:
				fprintf(output1,"M1_Hex_l=%9.5f\n", L1_Hex_l);
				break;
			case 2:
				fprintf(output1,"M2_Hex_l=%9.5f\n", L2_Hex_l);
				break;
			case 3:
				fprintf(output1,"M3_Hex_l=%9.5f\n", L3_Hex_l);
				break;
			case 4:
				fprintf(output1,"M4_Hex_l=%9.5f\n", L4_Hex_l);
			break;
			case 5:
				fprintf(output1,"M5_Hex_l=%9.5f\n", L5_Hex_l);
				break;
			case 6:
				fprintf(output1,"M6_Hex_l=%9.5f\n", L6_Hex_l);
				break;
			}
		}

		fclose(output1);
	}
	#endif

	#ifndef __VORO_GRAIN_MEDIA__
	if (VORO_GRAIN){
		
		// Find Hc of the bottom-most layer
		std::vector<double> M1_Hc_gbyg(NumOfGrains, 0.0);
		
		for (int idx=0; idx<NumOfGrains; ++idx){
			int itt1 = 1;
			while (M1_z_SingleGrain_field[idx][itt1]*M1_z_SingleGrain_field[idx][itt1-1] > 0 && itt1 < ceil(TOTAL_TIME/FieldSweepTimeStep)+1){
				M1_Hc_gbyg[idx] = Happl1_sweep[itt1-1];
				itt1++;
			}
			int itt2 = 100;
		}
		std::ofstream vfile1;
		vfile1.open("M1_Hc_gbyg.dat");
		for (int i=0; i<NumOfGrains; i++) vfile1 << M1_Hc_gbyg[i] << std::endl;

		// Hc mean
		double M1_Hc_mean = 0.0;
		for (int i=0; i<NumOfGrains; i++){
			if (indicator7[grain_coor[i][0].x + grain_coor[i][0].y*mNx + (mNz_6-1)*mNx*mNy] == 0) M1_Hc_mean = M1_Hc_mean + abs(M1_Hc_gbyg[i])/(count);
		}

		// Hc variance
		double M1_Hc_var = 0.0;
		for (int i=0; i<NumOfGrains; i++){
			if (indicator7[grain_coor[i][0].x + grain_coor[i][0].y*mNx + (mNz_6-1)*mNx*mNy] == 0) M1_Hc_var = M1_Hc_var + pow(abs(M1_Hc_gbyg[i])-M1_Hc_mean, 2);
		}

		// Hc std
		double M1_Hc_std = 0.0;
		M1_Hc_std = pow(M1_Hc_var/(count-1), 0.5);


		// Print to screen
		timeinfo = localtime(&rawtime);
		printf("\nM1_SFD%=%7.4f\n", M1_Hc_std/M1_Hc_mean);
		printf("M1_Hc_mean=%9.2f\n", M1_Hc_mean);
		printf("M1_Hc_std=%8.2f\n", M1_Hc_std);
	//	printf("M1_Hn_mean=%8.2f\n", M1_Hn_mean);
	//	printf("M1_CLUSTER_SIZE=%7.4f\n", cluster_size);
		printf("BL1=%9.5f\n", BL12);
		printf("BL2=%9.5f\n", BL23);
		printf("BL3=%9.5f\n", BL34);
		printf("BL4=%9.5f\n", BL45);
		printf("BL5=%9.5f\n", BL56);
		if (CGC_DEF){
			switch (CGC_label){
			case 1:
				printf("M1_Hex_l=%9.5f\n", L1_Hex_l);
				break;
			case 2:
				printf("M2_Hex_l=%9.5f\n", L2_Hex_l);
				break;
			case 3:
				printf("M3_Hex_l=%9.5f\n", L3_Hex_l);
				break;
			case 4:
				printf("M4_Hex_l=%9.5f\n", L4_Hex_l);
				break;
			case 5:
				printf("M5_Hex_l=%9.5f\n", L5_Hex_l);
				break;
			case 6:
				printf("M6_Hex_l=%9.5f\n", L6_Hex_l);
				break;
			}
		}
		cout << "" << endl;	
		cout << asctime(timeinfo) << std::endl;
	
		// Write to file
		FILE *output1;
		output1 = fopen("output.out", "w");
		fprintf(output1,"%s", asctime(timeinfo));
		fprintf(output1,"------------------------\n");
		fprintf(output1,"\nM1_SFD%=%7.4f\n", M1_Hc_std/M1_Hc_mean);
		fprintf(output1,"M1_Hc_mean=%9.2f\n", M1_Hc_mean);
		fprintf(output1,"M1_Hc_std=%8.2f\n", M1_Hc_std);
	//	fprintf(output1,"M1_Hn_mean=%8.2f\n", M1_Hn_mean);
	//	fprintf(output1, "M1_CLUSTER_SIZE=%7.4f\n", cluster_size);
		fprintf(output1,"BL1=%9.5f\n", BL12);
		fprintf(output1,"BL2=%9.5f\n", BL23);
		fprintf(output1,"BL3=%9.5f\n", BL34);
		fprintf(output1,"BL4=%9.5f\n", BL45);
		fprintf(output1,"BL5=%9.5f\n", BL56);
		if (CGC_DEF){
			switch (CGC_label){
			case 1:
				fprintf(output1,"M1_Hex_l=%9.5f\n", L1_Hex_l);
				break;
			case 2:
				fprintf(output1,"M2_Hex_l=%9.5f\n", L2_Hex_l);
				break;
			case 3:
				fprintf(output1,"M3_Hex_l=%9.5f\n", L3_Hex_l);
				break;
			case 4:
				fprintf(output1,"M4_Hex_l=%9.5f\n", L4_Hex_l);
			break;
			case 5:
				fprintf(output1,"M5_Hex_l=%9.5f\n", L5_Hex_l);
				break;
			case 6:
				fprintf(output1,"M6_Hex_l=%9.5f\n", L6_Hex_l);
				break;
			}
		}

		fclose(output1);



	//	//--------------------------------------------------
	//	// Cluster size calculation
	//	//--------------------------------------------------
	//
	//	// Locate the time where M1_Hc=0
	//	int t_Mz_is_0;
	//	for (int tt=FieldSweepTimeStep; tt<TOTAL_TIME; tt++){
	//		if (Mz_bar1[tt]*Mz_bar1[tt-FieldSweepTimeStep]<0) {
	//			t_Mz_is_0 = tt;
	//			break;
	//		}
	//	}

	//	// Set the center coordinates for the 2D grid
	//	std::vector<double> xcenter_SFD(mNx, 0.);
	//	std::vector<double> ycenter_SFD(mNy, 0.);

	//	for (int i_SFD=0; i_SFD<mNx; ++i_SFD) xcenter_SFD[i_SFD] = 0.5 + i_SFD;
	//	for (int j_SFD=0; j_SFD<mNy; ++j_SFD) ycenter_SFD[j_SFD] = 0.5 + j_SFD;
	//
	//	// Estimate pair correlation length
	//	std::vector<double> corr_ij(mNx, 0.);

	//	double r = 0.0;
	//	int tott_SFD = 0;
	//	int tk = 0;

	//	for (tott_SFD=0; tott_SFD<mNx; tott_SFD++)
	//	{
	//		r = r + 1;
	//		for (int i_SFD=0; i_SFD<mNx; i_SFD++){
	//			for (int j_SFD=0; j_SFD<mNy; j_SFD++){
	//				if (indicator7[i_SFD+j_SFD*mNx+(mNz-1)*mNx*mNy] == 0){
	//					tk = 0;
	//				
	//					for (int ii_SFD=0; ii_SFD<mNx; ii_SFD++){
	//						for (int jj_SFD=0; jj_SFD<mNy; jj_SFD++){
	//							if (indicator7[ii_SFD+jj_SFD*mNx+(mNz-1)*mNx*mNy] == 0){
	//								double rij = pow( (pow(xcenter_SFD[ii_SFD]-xcenter_SFD[i_SFD],2.0) + pow(ycenter_SFD[jj_SFD]-ycenter_SFD[j_SFD],2.0)), 0.5 );
	//								if (rij>=r-1 && rij<=r){
	//									tk = tk + 1;
	//									corr_ij[tott_SFD] = corr_ij[tott_SFD] + Mz_tt[ii_SFD+jj_SFD*mNx][t_Mz_is_0]*Mz_tt[i_SFD+j_SFD*mNx][t_Mz_is_0];
	//								}
	//							}
	//						}
	//					}
	//
	//				}
	//			}
	//		}
	//		corr_ij[tott_SFD]=corr_ij[tott_SFD]/(tk*2);
	//		rrt_cor[tott_SFD]=r;
	//	}

	//	// Write to file
	//	std::ofstream vfile2;
	//	vfile2.open("corr_func.dat");
	//	for (int i=0; i<mNx; i++) vfile2 << rrt_cor[i] << "\t" << corr_ij[i] <<std::endl;



	//	double cluster_size=0.0;
	//	int jsign;
	//	double max_NOR = corr_ij[0];
	//
	//	for (int i_SFD = 0; i_SFD < mNx; i_SFD++)	corr_ij[i_SFD] = corr_ij[i_SFD]/max_NOR;
	//
	//	jsign = 1;
	//	for(int i_SFD=0; i_SFD<mNx; i_SFD++){	
	//		if(corr_ij[i_SFD]>=0.3679 && corr_ij[i_SFD+1]<0.3679){
	//			double h_left=corr_ij[i_SFD]-0.3679;
	//			double h_right=0.3679-corr_ij[i_SFD+1];
	//			cluster_size = rrt_cor[i_SFD]+h_left/(h_left+h_right);
	//			jsign=2;
	//		}
	//		if(jsign == 1)	cluster_size = 0.0;
	//	}
	//	//--------------------------------------------//

	//	// Print to screen
	//	timeinfo = localtime(&rawtime);
	//	printf("\nM1_SFD%=%7.4f\n", M1_Hc_std/M1_Hc_mean);
	//	printf("M1_Hc_mean=%9.2f\n", M1_Hc_mean);
	//	printf("M1_Hc_std=%8.2f\n", M1_Hc_std);
	//	printf("M1_Hn_mean=%8.2f\n", M1_Hn_mean);
	//	printf("M1_CLUSTER_SIZE=%7.4f\n", cluster_size);
	//	printf("BL1=%9.5f\n", BL12);
	//	printf("BL2=%9.5f\n", BL23);
	//	printf("BL3=%9.5f\n", BL34);
	//	printf("BL4=%9.5f\n", BL45);
	//	printf("BL5=%9.5f\n", BL56);
	//	if (CGC_DEF){
	//		switch (CGC_label){
	//		case 1:
	//			printf("M1_Hex_l=%9.5f\n", L1_Hex_l);
	//			break;
	//		case 2:
	//			printf("M2_Hex_l=%9.5f\n", L2_Hex_l);
	//			break;
	//		case 3:
	//			printf("M3_Hex_l=%9.5f\n", L3_Hex_l);
	//			break;
	//		case 4:
	//			printf("M4_Hex_l=%9.5f\n", L4_Hex_l);
	//			break;
	//		case 5:
	//			printf("M5_Hex_l=%9.5f\n", L5_Hex_l);
	//			break;
	//		case 6:
	//			printf("M6_Hex_l=%9.5f\n", L6_Hex_l);
	//			break;
	//		}
	//	}
	//
	//	cout << asctime(timeinfo) << std::endl;
	//
	//	// Write to file
	//	FILE *output1;
	//	output1 = fopen("output.out", "w");
	//	fprintf(output1,"%s", asctime(timeinfo));
	//	fprintf(output1,"------------------------\n");
	//	fprintf(output1,"\nM1_SFD%=%7.4f\n", M1_Hc_std/M1_Hc_mean);
	//	fprintf(output1,"M1_Hc_mean=%9.2f\n", M1_Hc_mean);
	//	fprintf(output1,"M1_Hc_std=%8.2f\n", M1_Hc_std);
	//	fprintf(output1,"M1_Hn_mean=%8.2f\n", M1_Hn_mean);
	//	fprintf(output1, "M1_CLUSTER_SIZE=%7.4f\n", cluster_size);
	//	fprintf(output1,"BL1=%9.5f\n", BL12);
	//	fprintf(output1,"BL2=%9.5f\n", BL23);
	//	fprintf(output1,"BL3=%9.5f\n", BL34);
	//	fprintf(output1,"BL4=%9.5f\n", BL45);
	//	fprintf(output1,"BL5=%9.5f\n", BL56);
	//	if (CGC_DEF){
	//		switch (CGC_label){
	//		case 1:
	//			fprintf(output1,"M1_Hex_l=%9.5f\n", L1_Hex_l);
	//			break;
	//		case 2:
	//			fprintf(output1,"M2_Hex_l=%9.5f\n", L2_Hex_l);
	//			break;
	//		case 3:
	//			fprintf(output1,"M3_Hex_l=%9.5f\n", L3_Hex_l);
	//			break;
	//		case 4:
	//			fprintf(output1,"M4_Hex_l=%9.5f\n", L4_Hex_l);
	//		break;
	//		case 5:
	//			fprintf(output1,"M5_Hex_l=%9.5f\n", L5_Hex_l);
	//			break;
	//		case 6:
	//			fprintf(output1,"M6_Hex_l=%9.5f\n", L6_Hex_l);
	//			break;
	//		}
	//	}

	//	fclose(output1);
	}
	#endif
	
	return true;
}

bool FP_inp_theta_n_Magnitude(){
	int idx;
	for (int k = 0; k < fNz; k++){
		for (int j = 0; j < fNy; j++){
			for (int i = 0; i < fNx; i++){ 
				idx = i+j*fNx+k*fNx*fNy;
				FP_mag[idx] = pow(pow(FP_inp[0][idx],2)+pow(FP_inp[1][idx],2)+pow(FP_inp[2][idx],2), 0.5);
				FP_theta[idx] = acos(abs(FP_inp[2][idx])/FP_mag[idx])/PI*180;
			}
		}
	}

	return true;
}

bool Set_devMxyz(double* host_Ms, double* dev_theta, double* dev_phi, double* dev_Mx, double* dev_My, double* dev_Mz){
	double *theta = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	double *phi = (double*) calloc(mNx*mNy*mNz, sizeof(double));
	cudaMemcpy(theta, dev_theta, mNx*mNy*mNz*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(phi, dev_phi, mNx*mNy*mNz*sizeof(double),cudaMemcpyDeviceToHost);
	for (int k=0; k<mNz; ++k){
		for (int j=0; j<mNy; ++j){
			for (int i=0; i<mNx; ++i){
				int idx = i+j*mNx+k*mNx*mNy;
				Mx[idx] = host_Ms[idx]*sin(theta[idx])*cos(phi[idx]);
				My[idx] = host_Ms[idx]*sin(theta[idx])*sin(phi[idx]);
				Mz[idx] = host_Ms[idx]*cos(theta[idx]);
			}
		}
	}
	cudaMemcpy(dev_Mx, Mx, mNx*mNy*mNz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_My, My, mNx*mNy*mNz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Mz, Mz, mNx*mNy*mNz*sizeof(double),cudaMemcpyHostToDevice);
	
	return true;
}

bool FieldReset(){
	
	for (int k = 0; k < mNz; k++){
		for (int j = 0; j < mNy; j++){
			for (int i = 0; i < mNx; i++){ 
				int idx = i+j*mNx+k*mNx*mNy;
				Happl_x[idx] = 0.0;
				Happl_y[idx] = 0.0;
				Happl_z[idx] = 0.0;
			}
		}
	}

	cudaMemcpy(dev_Happl_x, Happl_x,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Happl_y, Happl_y,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Happl_z, Happl_z,  (mNx)*(mNy)*(mNz)*sizeof(double), cudaMemcpyHostToDevice);

	return true;
}
