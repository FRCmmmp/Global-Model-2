#include "Outputs.h"
#include <fstream>


bool Output_Float_3D_Format(int Nx, int Ny, int Nz, double *varArray, char filename[]){
	FILE *o1;
	o1 = fopen(filename, "w");
	for (int k = 0; k < (Nz); k++){
		for (int j = 0; j < (Ny); j++){
			for (int i = 0; i < (Nx); i++){
				fprintf(o1, "%16.5e", varArray[i+j*(Nx)+k*(Nx)*(Ny)]);
			}
			fprintf(o1, "\n");
		}
		fprintf(o1, "\n\n");
	}
	fclose(o1);
	
	return true;
}

bool Output_Int_3D_Format(int Nx, int Ny, int Nz, int *varArray, char filename[]){
	FILE *o1;
	o1 = fopen(filename, "w");
	for (int k = 0; k < (Nz); k++){
		for (int j = 0; j < (Ny); j++){
			for (int i = 0; i < (Nx); i++){
				fprintf(o1, "%5d", varArray[i+j*(Nx)+k*(Nx)*(Ny)]);
			}
			fprintf(o1, "\n");
		}
		fprintf(o1, "\n\n");
	}
	fclose(o1);
	
	return true;
}

bool Output_Float_1D_Format_6col(int Nx, int Ny, int Nz, double *varArray1, double *varArray2, double *varArray3, char filename[]){
	FILE *o1;
	o1 = fopen(filename, "w");
	for (int k = 0; k < (Nz); k++){
		for (int j = 0; j < (Ny); j++){
			for (int i = 0; i < (Nx); i++){
				fprintf(o1, "%5d%5d%5d%16.3lf%16.3lf%16.3lf\n", i, j, k, varArray1[i+j*(Nx)+k*(Nx)*(Ny)], varArray2[i+j*(Nx)+k*(Nx)*(Ny)], varArray3[i+j*(Nx)+k*(Nx)*(Ny)]);
			}
		}
	}
	fclose(o1);
	
	return true;
}

bool Output_Float_1D_Format_1col(int Nx, int Ny, int Nz, double *varArray1, char filename[]){
	FILE *o1;
	o1 = fopen(filename, "w");
	for (int k = 0; k < (Nz); k++){
		for (int j = 0; j < (Ny); j++){
			for (int i = 0; i < (Nx); i++){
				fprintf(o1, "%16.3lf\n", varArray1[i+j*(Nx)+k*(Nx)*(Ny)]);
			}
		}
	}
	fclose(o1);
	
	return true;
}

bool Output_Int_1D_Format_1col(int Nx, int Ny, int Nz, int *varArray1, char filename[]){
	FILE *o1;
	o1 = fopen(filename, "w");
	for (int k = 0; k < (Nz); k++){
		for (int j = 0; j < (Ny); j++){
			for (int i = 0; i < (Nx); i++){
				fprintf(o1, "%5d\n", varArray1[i+j*(Nx)+k*(Nx)*(Ny)]);
			}
		}
	}
	fclose(o1);
	
	return true;
}

bool Input_Float_1D_Format_1col(int mNx, int mNy, int mNz, double *varArray1, char filename[]){
	FILE *rfile1;
	rfile1 = fopen(filename, "r");
	for (int i = 0; i < mNx*mNy*mNz; i++){
		fscanf(rfile1, "%lf", &varArray1[i]);
	}
	fclose(rfile1);
	return true;
}

bool Input_Indicator_1D_Format(int mNx, int mNy, int mNz, int *varArray1, int *varArray2, int *varArray3, double *varArray4, double *varArray5, double *varArray6, int *varArray7, char filename[]){
	
	std::ifstream rfile1(filename);
	int idx = 0;
	while (rfile1 >> varArray1[idx] >> varArray2[idx] >> varArray3[idx] >> varArray4[idx] >> varArray5[idx] >> varArray6[idx] >> varArray7[idx]){
		idx++;
	}

	/*FILE *rfile1;
	rfile1 = fopen(filename, "r");
	for (int i = 0; i < mNx*mNy*mNz; i++){
		fscanf(rfile1, "%d %d %d %lf %lf %lf %d", &varArray1[i], &varArray2[i], &varArray3[i], &varArray4[i], &varArray5[i], &varArray6[i], &varArray7[i]);
	}
	fclose(rfile1);*/
	
	return true;
}