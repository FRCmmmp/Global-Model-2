#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Parameters_input.h"
#include "Parameters.h"


//***********************************************************************************************************//
// setting global variables
//***********************************************************************************************************//
#ifdef __6_LAYER__
int deviceID;
int BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, BLK_SZ_Z;
int mNx, mNy, mNz;
int FieldSweepTimeStep; 
double SweepFieldStep;
double ScalFactFieldStep;
unsigned long long int iseed;
int mNz_1, mNz_2, mNz_3, mNz_4, mNz_5, mNz_6;
int nz_12, nz_23, nz_34, nz_45, nz_56;
int mNz_12, mNz_23, mNz_34, mNz_45, mNz_56;
long int TOTAL_TIME, EQUI_START_TIME;
double GrainDx, GrainDy, delta_x, delta_y, delta_z, delta_t, Temperature,
	   field1, field2, field3, field4, field5, field6,
	   angle1, angle2, angle3, angle4, angle5, angle6,
	   field12, field23, field34, field45, field56,
	   angle12, angle23, angle34, angle45, angle56,
	   IniScalFact;
double L1_Ms, L1_Hk, L1_Ku,  L1_Aex, L1_alpha, L1_Hex_L, 
       L2_Ms, L2_Hk, L2_Ku,  L2_Aex, L2_alpha, L2_Hex_L,
       L3_Ms, L3_Hk, L3_Ku,  L3_Aex, L3_alpha, L3_Hex_L,
	   L4_Ms, L4_Hk, L4_Ku,  L4_Aex, L4_alpha, L4_Hex_L,
	   L5_Ms, L5_Hk, L5_Ku,  L5_Aex, L5_alpha, L5_Hex_L,
	   L6_Ms, L6_Hk, L6_Ku,  L6_Aex, L6_alpha, L6_Hex_L,
	   BL12_Ms, BL23_Ms, BL34_Ms, BL45_Ms, BL56_Ms,
	   dHk1_scale, dHk2_scale, dHk3_scale, dHk4_scale, dHk5_scale, dHk6_scale;
float  L1_Hex_l, BL12, BL12_Hex_l, dBL12_scale,
       L2_Hex_l, BL23, BL23_Hex_l, dBL23_scale,
       L3_Hex_l, BL34, BL34_Hex_l, dBL34_scale,
	   L4_Hex_l, BL45, BL45_Hex_l, dBL45_scale,
	   L5_Hex_l, BL56, BL56_Hex_l, dBL56_scale,
	   L6_Hex_l;
int MODEL, DEMAG, THERMAL, DEL_Hk, DEL_Aex, DEL_Tc, VORO_GRAIN, CGC_DEF, MH_LOOP, AFC, EXT_FP_PROFILE, PRE_DEFINED_PATTERN, FULL_REC;
int AF_layer_label;
double def_perc;
int CGC_label;
int AC_DC_erased;
int DT_Rec_Analysis, CT_Rec_Analysis;
int DT_x0, DT_y0, CT_x0, CT_y0;

// Demag time interval
int DmagINT;

// Demag factor
double DmagFAC;

// Recording Simulation
int fNx, fNy, fNz, fNx0, fNy0, fNz0, f_x0; //head field profile dimensions
double sBL; //shortest bit length
double v; //head velocity [cm/s]
int sfNc; //number of the shortest field cycle (sfCNt)

// Varialbe derived from input data
int dxNt;  // Number of time steps for head to travel a distance delta_x at head speed v (dxNt, derived)
int sfcNt; //number of time steps of the shortest field cycle (sfcNt, derived)

// Varialbe derived from input data
int DEG_FREEDOM;


// Voronoi
int media_type;
int num_vgrains;
double minimum_sdistance;
double grain_bnd;
int num_vdomains;
double rgrain_size;
double rgrain_scaling;
double domain_bnd_qfactor;




#endif
//***********************************************************************************************************//




//***********************************************************************************************************//
// declaration of function
//***********************************************************************************************************//
int is_file_exist(const char *fileName);
int Input(void);
int Input_FieldProfile(int EXT_FP_PROFILE);
int Input_FieldSeq(int EXT_FP_PROFILE);
int Input_Mag0(int PRE_DEFINED_PATTERN);
//***********************************************************************************************************//








//***********************************************************************************************************//
/*
 * Check if a file exist using fopen() function
 * return 1 if the file exist otherwise return 0
*/
int is_file_exist(const char* filename){
    struct stat buffer;
    int exist = stat(filename,&buffer);
    if(exist == 0)
        return 1;
    else 
        return 0;
}



#ifdef __6_LAYER__
int Input(void)
{
	FILE *file1, *file2;
	char buff[MAX_CHAR_SZ], deli[MAX_CHAR_SZ];
	file2 = fopen("Config.dat", "w");
	file1 = fopen("Input_6layer.txt", "r");
	if (!file1){
		printf("open Input_6layer.txt Error! \n");
	}
	
	// Headerline
	fgets (buff, MAX_CHAR_SZ, file1);
	printf("%s", buff);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Set GPU device
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &deviceID);
	printf("%d \n", deviceID);
	fprintf(file2, "%s", buff);
	// BLOCK_SIZE_XYZ (Kernel cnofiguration)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &BLOCK_SIZE_X);
	printf("%d \n", BLOCK_SIZE_X);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &BLOCK_SIZE_Y);
	printf("%d \n", BLOCK_SIZE_Y);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &BLOCK_SIZE_Z);
	printf("%d \n", BLOCK_SIZE_Z);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &BLK_SZ_Z);
	printf("%d \n", BLK_SZ_Z);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Media Dimensions
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNx);
	printf("%d \n", mNx);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNy);
	printf("%d \n", mNy);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz);
	printf("%d \n", mNz);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz_1);
	printf("%d \n", mNz_1);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz_2);
	printf("%d \n", mNz_2);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz_3);
	printf("%d \n", mNz_3);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz_4);
	printf("%d \n", mNz_4);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz_5);
	printf("%d \n", mNz_5);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &mNz_6);
	printf("%d \n", mNz_6);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &nz_12);
	printf("%d \n", nz_12);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &nz_23);
	printf("%d \n", nz_23);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &nz_34);
	printf("%d \n", nz_34);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &nz_45);
	printf("%d \n", nz_45);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &nz_56);
	printf("%d \n", nz_56);
	fprintf(file2, "%s", buff);
	
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	
	// media type: 0=normal;  1=short-ranged ordering (sro)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &media_type);
	printf("%d \n", media_type);
	fprintf(file2, "%s", buff);
	// number of Voronoi grains       (for normal media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &num_vgrains);
	printf("%d \n", num_vgrains);
	fprintf(file2, "%s", buff);
	// minimum distance between seeds (for normal media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &minimum_sdistance);
	printf("%lf \n", minimum_sdistance);
	fprintf(file2, "%s", buff);
	// grain boundary                 (for normal media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &grain_bnd);
	printf("%lf \n", grain_bnd);
	fprintf(file2, "%s", buff);
	// number of domains              (for sro media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &num_vdomains);
	printf("%d \n", num_vdomains);
	fprintf(file2, "%s", buff);
	// grain size                     (for sro media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &rgrain_size);
	printf("%lf \n", rgrain_size);
	fprintf(file2, "%s", buff);
	// grain scaling                  (for sro media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &rgrain_scaling);
	printf("%lf \n", rgrain_scaling);
	fprintf(file2, "%s", buff);
	// domain-domain bnd control       (for sro media type)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &domain_bnd_qfactor);
	printf("%lf \n", domain_bnd_qfactor);
	fprintf(file2, "%s", buff);
	
	// Deliminator
	fgets (deli, MAX_CHAR_SZ, file1);
	printf("%s", deli);
	fprintf(file2, "%s", deli);
	
	// Integration Time
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%ld", &TOTAL_TIME);
	printf("%ld \n", TOTAL_TIME);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%ld", &EQUI_START_TIME);
	printf("%ld \n", EQUI_START_TIME);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Boolean for Demag
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DEMAG);
	printf("%d \n", DEMAG);
	fprintf(file2, "%s", buff);
	// Demag update time interval (DmagINT)
    fgets (buff, MAX_CHAR_SZ, file1);
    sscanf(buff, "%d", &DmagINT);
    printf("%d \n", DmagINT);
    fprintf(file2, "%s", buff);
    // Demag factor (DmagFAC)
    fgets (buff, MAX_CHAR_SZ, file1);
    sscanf(buff, "%lf", &DmagFAC);
    printf("%lf \n", DmagFAC);
    fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// delta_xyzt, grain diameter
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &GrainDx);
	printf("%e \n", GrainDx);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &GrainDy);
	printf("%e \n", GrainDy);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &delta_x);
	printf("%e \n", delta_x);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &delta_y);
	printf("%e \n", delta_y);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &delta_z);
	printf("%e \n", delta_z);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &delta_t);
	printf("%e \n", delta_t);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Random Seed for Hth
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &iseed);
	printf("%d \n", iseed);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &THERMAL);
	printf("%d \n", THERMAL);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &Temperature);
	printf("%lf \n", Temperature);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Initial field and field angle for major/minor loops
	// Layer1
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field1);
	printf("%e \n", field1);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle1);
	printf("%e \n", angle1);
	fprintf(file2, "%s", buff);
	// Layer12
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field12);
	printf("%e \n", field12);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle12);
	printf("%e \n", angle12);
	fprintf(file2, "%s", buff);
	// Layer2
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field2);
	printf("%e \n", field2);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle2);
	printf("%e \n", angle2);
	fprintf(file2, "%s", buff);
	// Layer23
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field23);
	printf("%e \n", field23);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle23);
	printf("%e \n", angle23);
	fprintf(file2, "%s", buff);
	// Layer3
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field3);
	printf("%e \n", field3);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle3);
	printf("%e \n", angle3);
	fprintf(file2, "%s", buff);
	// Layer34
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field34);
	printf("%e \n", field34);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle34);
	printf("%e \n", angle34);
	fprintf(file2, "%s", buff);
	// Layer4
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field4);
	printf("%e \n", field4);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle4);
	printf("%e \n", angle4);
	fprintf(file2, "%s", buff);
	// Layer45
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field45);
	printf("%e \n", field45);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle45);
	printf("%e \n", angle45);
	fprintf(file2, "%s", buff);
	// Layer5
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field5);
	printf("%e \n", field5);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle5);
	printf("%e \n", angle5);
	fprintf(file2, "%s", buff);
	// Layer45
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field56);
	printf("%e \n", field56);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle56);
	printf("%e \n", angle56);
	fprintf(file2, "%s", buff);
	// Layer6
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &field6);
	printf("%e \n", field6);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &angle6);
	printf("%e \n", angle6);
	fprintf(file2, "%s", buff);
	// Initial Scaling Factor
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &IniScalFact);
	printf("%e \n", IniScalFact);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Layer 1
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L1_Ms);
	printf("%e \n", L1_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L1_Ku);
	printf("%e \n", L1_Ku);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L1_Aex);
	printf("%e \n", L1_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &L1_Hex_l);
	printf("%e \n", L1_Hex_l);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &dHk1_scale);
	printf("%e \n", dHk1_scale);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L1_alpha);
	printf("%e \n", L1_alpha);
	fprintf(file2, "%s", buff);
	// BL12
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL12);
	printf("%e \n", BL12);
	fprintf(file2, "%s", buff);
	// dBL12_scale
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &dBL12_scale);
	printf("%e \n", dBL12_scale);
	fprintf(file2, "%s", buff);
	// BL12_Ms
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &BL12_Ms);
	printf("%e \n", (double)BL12_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL12_Hex_l);
	printf("%e \n", BL12_Hex_l);
	fprintf(file2, "%s", buff);

	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);

	// Layer 2
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L2_Ms);
	printf("%e \n", L2_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L2_Ku);
	printf("%e \n", L2_Ku);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L2_Aex);
	printf("%e \n", L2_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &L2_Hex_l);
	printf("%e \n", L2_Hex_l);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &dHk2_scale);
	printf("%e \n", dHk2_scale);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L2_alpha);
	printf("%e \n", L2_alpha);
	fprintf(file2, "%s", buff);
	// BL23
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL23);
	printf("%e \n", BL23);
	fprintf(file2, "%s", buff);
	// dBL23_scale
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &dBL23_scale);
	printf("%e \n", dBL23_scale);
	fprintf(file2, "%s", buff);
	// BL23_Ms
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &BL23_Ms);
	printf("%e \n", BL23_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL23_Hex_l);
	printf("%e \n", BL23_Hex_l);
	fprintf(file2, "%s", buff);

	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);

	// Layer 3
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L3_Ms);
	printf("%e \n", L3_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L3_Ku);
	printf("%e \n", L3_Ku);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L3_Aex);
	printf("%e \n", L3_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &L3_Hex_l);
	printf("%e \n", L3_Hex_l);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &dHk3_scale);
	printf("%e \n", dHk3_scale);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L3_alpha);
	printf("%e \n", L3_alpha);
	fprintf(file2, "%s", buff);
	// BL34
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL34);
	printf("%e \n", BL34);
	fprintf(file2, "%s", buff);
	// dBL34_scale
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &dBL34_scale);
	printf("%e \n", dBL34_scale);
	fprintf(file2, "%s", buff);
	// BL34_Ms
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &BL34_Ms);
	printf("%e \n", BL34_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL34_Hex_l);
	printf("%e \n", BL34_Hex_l);
	fprintf(file2, "%s", buff);

	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);

	// Layer 4
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L4_Ms);
	printf("%e \n", L4_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L4_Ku);
	printf("%e \n", L4_Ku);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L4_Aex);
	printf("%e \n", L4_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &L4_Hex_l);
	printf("%e \n", L4_Hex_l);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &dHk4_scale);
	printf("%e \n", dHk4_scale);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L4_alpha);
	printf("%e \n", L4_alpha);
	fprintf(file2, "%s", buff);
	// BL45
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL45);
	printf("%e \n", BL45);
	fprintf(file2, "%s", buff);
	// dBL45_scale
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &dBL45_scale);
	printf("%e \n", dBL45_scale);
	fprintf(file2, "%s", buff);
	// BL45_Ms
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &BL45_Ms);
	printf("%e \n", BL45_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL45_Hex_l);
	printf("%e \n", BL45_Hex_l);
	fprintf(file2, "%s", buff);

	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);

	// Layer 5
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L5_Ms);
	printf("%e \n", L5_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L5_Ku);
	printf("%e \n", L5_Ku);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L5_Aex);
	printf("%e \n", L5_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &L5_Hex_l);
	printf("%e \n", L5_Hex_l);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &dHk5_scale);
	printf("%e \n", dHk5_scale);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L5_alpha);
	printf("%e \n", L5_alpha);
	fprintf(file2, "%s", buff);
	// BL56
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL56);
	printf("%e \n", BL56);
	fprintf(file2, "%s", buff);
	// dBL56_scale
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &dBL56_scale);
	printf("%e \n", dBL56_scale);
	fprintf(file2, "%s", buff);
	// BL56_Ms
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &BL56_Ms);
	printf("%e \n", BL56_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &BL56_Hex_l);
	printf("%e \n", BL56_Hex_l);
	fprintf(file2, "%s", buff);

	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);

	// Layer 6
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L6_Ms);
	printf("%e \n", L6_Ms);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L6_Ku);
	printf("%e \n", L6_Ku);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L6_Aex);
	printf("%e \n", L6_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%f", &L6_Hex_l);
	printf("%e \n", L6_Hex_l);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &dHk6_scale);
	printf("%e \n", dHk6_scale);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &L6_alpha);
	printf("%e \n", L6_alpha);
	fprintf(file2, "%s", buff);
	
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	
	// MODELLING METHOD
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &MODEL);
	printf("%d\n", MODEL);
	fprintf(file2, "%s", buff);

	// delta Hk, delta Aex, delta Tc --> turn on/off
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DEL_Hk);
	printf("%d \n", DEL_Hk);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DEL_Aex);
	printf("%d \n", DEL_Aex);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DEL_Tc);
	printf("%d \n", DEL_Tc);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &VORO_GRAIN);
	printf("%d \n", VORO_GRAIN);
	fprintf(file2, "%s", buff);
	// Hysteresis Loop Mode?
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &MH_LOOP);
	printf("%d \n", MH_LOOP);
	fprintf(file2, "%s", buff);
	// AFC Mode?
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &AFC);
	printf("%d \n", AFC);
	fprintf(file2, "%s", buff);
	// AF layer label?
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &AF_layer_label);
	printf("%d \n", AF_layer_label);
	fprintf(file2, "%s", buff);
	// 3D Recording Analysis
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &FULL_REC);
	printf("%d \n", FULL_REC);
	fprintf(file2, "%s", buff);
	// Input external field profile?
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &EXT_FP_PROFILE);
	printf("%d \n", EXT_FP_PROFILE);
	fprintf(file2, "%s", buff);
	// Pre-defined pattern?
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &PRE_DEFINED_PATTERN);
	printf("%d \n", PRE_DEFINED_PATTERN);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &CGC_DEF);
	printf("%d \n", CGC_DEF);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &CGC_label);
	printf("%d \n", CGC_label);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &def_perc);
	printf("%f \n", def_perc);
	fprintf(file2, "%s", buff);
	// AC- or DC-erased
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &AC_DC_erased);
	printf("%d \n", AC_DC_erased);
	fprintf(file2, "%s", buff);
	// DT Recording Analysis
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DT_Rec_Analysis);
	printf("%d \n", DT_Rec_Analysis);
	fprintf(file2, "%s", buff);
	// DT_x0 (for DT Recording Analysis)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DT_x0);
	printf("%d \n", DT_x0);
	fprintf(file2, "%s", buff);
	// DT_y0 (for DT Recording Analysis)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &DT_y0);
	printf("%d \n", DT_y0);
	fprintf(file2, "%s", buff);
	// CT Recording Analysis
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &CT_Rec_Analysis);
	printf("%d \n", CT_Rec_Analysis);
	fprintf(file2, "%s", buff);
	// CT_x0 (for CT Recording Analysis)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &CT_x0);
	printf("%d \n", CT_x0);
	fprintf(file2, "%s", buff);
	// CT_y0 (for CT Recording Analysis)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &CT_y0);
	printf("%d \n", CT_y0);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	// Sweeping Rate for Loops
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &FieldSweepTimeStep);
	printf("%d \n", FieldSweepTimeStep);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &SweepFieldStep);
	printf("%lf \n", SweepFieldStep);
	fprintf(file2, "%s", buff);
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &ScalFactFieldStep);
	printf("%lf \n", ScalFactFieldStep);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%s", &deli);
	printf("%s \n", deli);
	fprintf(file2, "%s", buff);
	/* Head Field */
	// x
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &fNx);
	printf("%d \n", fNx);
	fprintf(file2, "%s", buff);
	// y
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &fNy);
	printf("%d \n", fNy);
	fprintf(file2, "%s", buff);
	// z
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &fNz);
	printf("%d \n", fNz);
	fprintf(file2, "%s", buff);
	// x0
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &fNx0);
	printf("%d \n", fNx0);
	fprintf(file2, "%s", buff);
	// y0
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &fNy0);
	printf("%d \n", fNy0);
	fprintf(file2, "%s", buff);
	// z0
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &fNz0);
	printf("%d \n", fNz0);
	fprintf(file2, "%s", buff);
	// f_x0
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &f_x0);
	printf("%d \n", f_x0);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (deli, MAX_CHAR_SZ, file1);
	printf("%s", deli);
	fprintf(file2, "%s", deli);
	/* Recording Process */
	// Shortest Bit Length
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &sBL);
	printf("%7.3e \n", sBL);
	fprintf(file2, "%s", buff);
	// Head Velocity
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%lf", &v);
	printf("%lf \n", v);
	fprintf(file2, "%s", buff);
	// number of the shortest field cycles (sfNc)
	fgets (buff, MAX_CHAR_SZ, file1);
	sscanf(buff, "%d", &sfNc);
	printf("%d \n", sfNc);
	fprintf(file2, "%s", buff);
	// Deliminator
	fgets (deli, MAX_CHAR_SZ, file1);
	printf("%s\n", deli);
	fprintf(file2, "%s\n", deli);





	// Varialbe derived from input data
	int nz_1, nz_2, nz_3, nz_4, nz_5, nz_6, n_BLs;
	
	// Make mNz be multiple of 2 (for CUDA)
	n_BLs = nz_12 + nz_23 + nz_34 + nz_45 + nz_56;
	if ((n_BLs+mNz_6)%2 == 0) mNz = mNz_6 + n_BLs;
	else mNz = mNz_6 + n_BLs + 1;

	// Make mNx and mNy be multiple of 2 (for CUDA)
	if (mNx%2 == 0) mNx = mNx;
	else mNx = mNx + 1;
	if (mNy%2 == 0) mNy = mNy;
	else mNy = mNy + 1;

	// #'s monolayers for mag layers
	nz_1 = mNz_1;
	nz_2 = mNz_2-mNz_1;
	nz_3 = mNz_3-mNz_2;
	nz_4 = mNz_4-mNz_3;
	nz_5 = mNz_5-mNz_4;
	nz_6 = mNz_6-mNz_5;
	
	// z coordinates for all layers
	mNz_1 = nz_1;
	mNz_2 = mNz_1 + nz_2 + nz_12;
	mNz_3 = mNz_2 + nz_3 + nz_23;
	mNz_4 = mNz_3 + nz_4 + nz_34;
	mNz_5 = mNz_4 + nz_5 + nz_45;
	mNz_6 = mNz_5 + nz_6 + nz_56;
	mNz_12 = mNz_1 + nz_12;
	mNz_23 = mNz_2 + nz_23;
	mNz_34 = mNz_3 + nz_34;
	mNz_45 = mNz_4 + nz_45;
	mNz_56 = mNz_5 + nz_56;

	// Write to screen
	std::cout << "------------mNx and mNy------------"  <<std::endl;
	std::cout << setw(7) << "mNx    ="  << setw(6) << mNx  <<std::endl;
	std::cout << setw(7) << "mNy    ="  << setw(6) << mNy  <<std::endl;
	std::cout << "------------#'s monolayers for mag layers------------"  <<std::endl;
	std::cout << setw(7) << "nNz_1 ="  << setw(6) << nz_1  <<std::endl;
	std::cout << setw(7) << "nNz_2 ="  << setw(6) << nz_2  <<std::endl;
	std::cout << setw(7) << "nNz_3 ="  << setw(6) << nz_3  <<std::endl;
	std::cout << setw(7) << "nNz_4 ="  << setw(6) << nz_4  <<std::endl;
	std::cout << setw(7) << "nNz_5 ="  << setw(6) << nz_5  <<std::endl;
	std::cout << setw(7) << "nNz_6 ="  << setw(6) << nz_6  <<std::endl;
	std::cout << "------------#'s monolayers for break layers------------"  <<std::endl;
	std::cout << setw(7) << "nNz_12 ="  << setw(6) << nz_12  <<std::endl;
	std::cout << setw(7) << "nNz_23 ="  << setw(6) << nz_23  <<std::endl;
	std::cout << setw(7) << "nNz_34 ="  << setw(6) << nz_34  <<std::endl;
	std::cout << setw(7) << "nNz_45 ="  << setw(6) << nz_45  <<std::endl;
	std::cout << setw(7) << "nNz_56 ="  << setw(6) << nz_56  <<std::endl;
	std::cout << "------------z coordinates for all layers------------"  <<std::endl;
	std::cout << setw(8) << "mNz_1  =" << setw(6) << mNz_1  <<std::endl;
	std::cout << setw(8) << "mNz_12 =" << setw(6) << mNz_12  <<std::endl;
	std::cout << setw(8) << "mNz_2  ="  << setw(6) << mNz_2  <<std::endl;
	std::cout << setw(8) << "mNz_23 =" << setw(6) << mNz_23  <<std::endl;
	std::cout << setw(8) << "mNz_3  ="  << setw(6) << mNz_3  <<std::endl;
	std::cout << setw(8) << "mNz_34 =" << setw(6) << mNz_34  <<std::endl;
	std::cout << setw(8) << "mNz_4  ="  << setw(6) << mNz_4  <<std::endl;
	std::cout << setw(8) << "mNz_45 =" << setw(6) << mNz_45  <<std::endl;
	std::cout << setw(8) << "mNz_5  ="  << setw(6) << mNz_5  <<std::endl;
	std::cout << setw(8) << "mNz_56 =" << setw(6) << mNz_56  <<std::endl;
	std::cout << setw(8) << "mNz_6  ="  << setw(6) << mNz_6  <<std::endl;
	std::cout << setw(8) << "mNz    ="  << setw(6) << mNz  <<std::endl;
	std::cout << "----------------------------------------------------"  <<std::endl;
	

	fprintf(file2, "///////////////////////////////////////////////////////// \n");
	fprintf(file2, "%d    // mNz_1\n", mNz_1);
	fprintf(file2, "%d    // mNz_12\n", mNz_12);
	fprintf(file2, "%d    // mNz_2\n", mNz_2);
	fprintf(file2, "%d    // mNz_23\n", mNz_23);
	fprintf(file2, "%d    // mNz_3\n", mNz_3);
	fprintf(file2, "%d    // mNz_34\n", mNz_34);
	fprintf(file2, "%d    // mNz_4\n", mNz_4);
	fprintf(file2, "%d    // mNz_45\n", mNz_45);
	fprintf(file2, "%d    // mNz_5\n", mNz_5);
	fprintf(file2, "%d    // mNz_56\n", mNz_56);
	fprintf(file2, "%d    // mNz_6\n", mNz_6);
	fprintf(file2, "%d    // mNz\n", mNz);






	DEG_FREEDOM = mNx*mNy*mNz*3;
    L1_Hk = L1_Ku/L1_Ms*2;
    L2_Hk = L2_Ku/L2_Ms*2;
    L3_Hk = L3_Ku/L3_Ms*2;
    L4_Hk = L4_Ku/L4_Ms*2;
    L5_Hk = L5_Ku/L5_Ms*2;
    L6_Hk = L6_Ku/L6_Ms*2;
	L1_Hex_L = (double)2*L1_Hex_l*L1_Aex/L1_Ms/pow(GrainDx, 2);   // Lateral exchange field
	L2_Hex_L = (double)2*L2_Hex_l*L2_Aex/L2_Ms/pow(GrainDx, 2);   // Lateral exchange field
	L3_Hex_L = (double)2*L3_Hex_l*L3_Aex/L3_Ms/pow(GrainDx, 2);   // Lateral exchange field
	L4_Hex_L = (double)2*L4_Hex_l*L4_Aex/L4_Ms/pow(GrainDx, 2);   // Lateral exchange field
	L5_Hex_L = (double)2*L5_Hex_l*L5_Aex/L5_Ms/pow(GrainDx, 2);   // Lateral exchange field
	L6_Hex_L = (double)2*L6_Hex_l*L6_Aex/L6_Ms/pow(GrainDx, 2);   // Lateral exchange field
	fprintf(file2, "///////////////////////////////////////////////////////// \n");
	fprintf(file2, "%12.4f    // L1_Hk Oe\n", L1_Hk); 
	fprintf(file2, "%12.4f    // L2_Hk Oe\n", L2_Hk); 
	fprintf(file2, "%12.4f    // L3_Hk Oe\n", L3_Hk); 
	fprintf(file2, "%12.4f    // L4_Hk Oe\n", L4_Hk); 
	fprintf(file2, "%12.4f    // L5_Hk Oe\n", L5_Hk); 
	fprintf(file2, "%12.4f    // L6_Hk Oe\n", L6_Hk); 
	fprintf(file2, "///////////////////////////////////////////////////////// \n");
	fprintf(file2, "%6.3e    // L1_Hex_L (Oe)\n", L1_Hex_L); 
	fprintf(file2, "%6.3e    // L2_Hex_L (Oe)\n", L2_Hex_L); 
	fprintf(file2, "%6.3e    // L3_Hex_L (Oe)\n", L3_Hex_L); 
	fprintf(file2, "%6.3e    // L4_Hex_L (Oe)\n", L4_Hex_L); 
	fprintf(file2, "%6.3e    // L5_Hex_L (Oe)\n", L5_Hex_L); 
	fprintf(file2, "%6.3e    // L6_Hex_L (Oe)\n", L6_Hex_L);

	// Number of time steps for head to travel a distance delta_x at head speed v (dxNt, derived)
	if (v==0.0) dxNt = 1000000000;
	else dxNt = (int)(ceil)(delta_x/v/delta_t);
	printf("%d \n",  dxNt);
	fprintf(file2, "%d      // Number of time steps for head to travel a distance delta_x at head speed v \n", dxNt);
	// Number of time steps of the shortest field cycle (sfcNt, derived)
	if (v==0.0) sfcNt = 1000000000;
	else sfcNt = (int)ceil(sBL/v/delta_t);
	printf("%d \n", sfcNt);
	fprintf(file2, "%d     // Number of time steps of the shortest field cycle (derived) \n", sfcNt);
	fprintf(file2, "///////////////////////////////////////////////////////// \n");



	
	fclose(file1);
	fclose(file2);
	return 1;
} 

int Input_FieldProfile(int EXT_FP_PROFILE){
	if (EXT_FP_PROFILE){
		FILE *rfile1;
		rfile1 = fopen("FieldProfile.inp", "r");
		for (int i = 0; i < fNx*fNy*fNz; i++){
			fscanf(rfile1, "%lf %lf %lf %lf %lf %lf", &fGrid[0][i], &fGrid[1][i], &fGrid[2][i], &FP_inp[0][i], &FP_inp[1][i], &FP_inp[2][i]);
		}
		fclose(rfile1);
		printf("Read-in FieldProfile.inp Succeeded! \n");
	}
	return 1;
}

int Input_FieldSeq(int EXT_FP_PROFILE){
	if (EXT_FP_PROFILE && !MH_LOOP){
		FILE *rfile1;
		rfile1 = fopen("FieldTimeSequence.inp", "r");
		for (int i = 0; i < sfNc; i++){
			fscanf(rfile1, "%lf", &fSeq[i]);
		}
		fclose(rfile1);
		printf("Read-in FieldTimeSequence.inp succeeded!! \n");
	}
	return 1;
}

int Input_Mag0(int PRE_DEFINED_PATTERN){
	char filename[30] = "Mag0.dat";
	if (PRE_DEFINED_PATTERN && is_file_exist(filename)){
		FILE *rfile1 = fopen("Mag0.dat", "r");
		if (rfile1 != NULL){
			for (int i = 0; i < mNx*mNy*mNz; i++){
				fscanf(rfile1, "%lf %lf %lf %lf %lf %lf", &mGrid[0][i], &mGrid[1][i], &mGrid[2][i], &Mag0[0][i], &Mag0[1][i], &Mag0[2][i]);
			}
			fclose(rfile1);	
		}
		printf("Read-in Mag0.dat succeeded!! \n");

		for (int i = 0; i < mNx*mNy*mNz; i++){
			theta0[i] = acos(Mag0[2][i]/pow(pow(Mag0[0][i],2)+pow(Mag0[1][i],2)+pow(Mag0[2][i],2), 0.5));
		}
	}
	return 1;
}


#endif




