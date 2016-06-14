#ifndef __PARAMETER_INPUT__
#define __PARAMETER_INPUT__

#define __6_LAYER__


//***********************************************************************************************************//
// input parameter variable declartion in header file
//***********************************************************************************************************//
#ifdef __6_LAYER__
#define   MAX_CHAR_SZ   200
extern int deviceID;
extern int BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z, BLK_SZ_Z;
extern int mNx, mNy, mNz;
extern int FieldSweepTimeStep; 
extern double SweepFieldStep;
extern double ScalFactFieldStep;
extern unsigned long long int iseed;
extern int mNz_1, mNz_2, mNz_3, mNz_4, mNz_5, mNz_6;
extern int nz_12, nz_23, nz_34, nz_45, nz_56;
extern int mNz_12, mNz_23, mNz_34, mNz_45, mNz_56;
extern long int TOTAL_TIME, EQUI_START_TIME;
extern double GrainDx, GrainDy, delta_x, delta_y, delta_z, delta_t, Temperature,
	          field1, field2, field3, field4, field5, field6,
			  angle1, angle2, angle3, angle4, angle5, angle6,
			  field12, field23, field34, field45, field56,
			  angle12, angle23, angle34, angle45, angle56,
			  IniScalFact;
extern double L1_Ms, L1_Hk, L1_Ku,  L1_Aex, L1_alpha, L1_Hex_L, 
              L2_Ms, L2_Hk, L2_Ku,  L2_Aex, L2_alpha, L2_Hex_L,
              L3_Ms, L3_Hk, L3_Ku,  L3_Aex, L3_alpha, L3_Hex_L,
	          L4_Ms, L4_Hk, L4_Ku,  L4_Aex, L4_alpha, L4_Hex_L,
	          L5_Ms, L5_Hk, L5_Ku,  L5_Aex, L5_alpha, L5_Hex_L,
	          L6_Ms, L6_Hk, L6_Ku,  L6_Aex, L6_alpha, L6_Hex_L,
			  BL12_Ms, BL23_Ms, BL34_Ms, BL45_Ms, BL56_Ms,
			  dHk1_scale, dHk2_scale, dHk3_scale, dHk4_scale, dHk5_scale, dHk6_scale;
extern float  L1_Hex_l, BL12, BL12_Hex_l, dBL12_scale,
              L2_Hex_l, BL23, BL23_Hex_l, dBL23_scale,
              L3_Hex_l, BL34, BL34_Hex_l, dBL34_scale,
	          L4_Hex_l, BL45, BL45_Hex_l, dBL45_scale,
	          L5_Hex_l, BL56, BL56_Hex_l, dBL56_scale,
	          L6_Hex_l;
extern int MODEL, DEMAG, THERMAL, DEL_Hk, DEL_Aex, DEL_Tc, VORO_GRAIN, CGC_DEF, MH_LOOP, AFC, EXT_FP_PROFILE, PRE_DEFINED_PATTERN, FULL_REC;
extern int AF_layer_label;
extern double def_perc;
extern int CGC_label;
extern int AC_DC_erased;
extern int DT_Rec_Analysis, CT_Rec_Analysis;
extern int DT_x0, DT_y0, CT_x0, CT_y0;

// Demag time interval
extern int DmagINT;

// Demag factor (DmagFAC)
extern double DmagFAC;

// Recording Simulation
extern int fNx, fNy, fNz, fNx0, fNy0, fNz0, f_x0; //head field profile dimensions
extern double sBL; //shortest bit length
extern double v; //head velocity [cm/s]
extern int sfNc; //number of the shortest field cycle (sfCNt)

// Varialbe derived from input data
extern int dxNt;  // Number of time steps for head to travel a distance delta_x at head speed v (dxNt, derived)
extern int sfcNt; //number of time steps of the shortest field cycle (sfcNt, derived)


// Varialbe derived from input data
extern int DEG_FREEDOM;

// Voronoi
extern int media_type;
extern int num_vgrains;
extern double minimum_sdistance;
extern double grain_bnd;
extern int num_vdomains;
extern double rgrain_size;
extern double rgrain_scaling;
extern double domain_bnd_qfactor;



#endif
//***********************************************************************************************************//



//***********************************************************************************************************//
// input parameter variable declartion in header file
//***********************************************************************************************************//
extern int Input(void);
extern int Input_FieldProfile(int);
extern int Input_FieldSeq(int);
extern int Input_Mag0(int);
//***********************************************************************************************************//



#endif
