#ifndef __PARAMETER__
#define __PARAMETER__

#include <math.h>
#include <complex>
#include <cufft.h>
#include "Parameters_input.h"

using namespace std;

// Constant
#define  PI              (double)acos(-1.0)
#define  ALPHA           (double)0.1
#define  GAMMA           (1.76e7) // /(1+pow(ALPHA, 2.0))
#define  ZERO            0.0
#define  Ini_THETA_Up    0.01
#define  Ini_THETA_Down  PI-0.01

const double kb = 1.38e-16;

// Variables
double static Aex = 1.00e-6, Ku = 4.14e7, Ms = 925.07, Ms_soft = 1400; 
long static idum0;
int static flag = 0;

// Array for device
static double *dev_Aex = NULL, *dev_Ku = NULL, *dev_Ms = NULL, *dev_alpha = NULL, *dev_gamma = NULL;
static double *dev_Aex_temp = NULL, *dev_Ku_temp = NULL, *dev_Ms_temp = NULL, *dev_alpha_temp = NULL;

static double theta = 0, phi = 0;

static double *dev_GasArray = NULL;

static double *dev_theta = NULL, *dev_phi = NULL;

static double *dev_Hth_x = NULL, *dev_Hth_y = NULL, *dev_Hth_z = NULL;
static double *dev_Hk_x  = NULL, *dev_Hk_y  = NULL, *dev_Hk_z  = NULL;
static double *dev_Ha_x  = NULL, *dev_Ha_y  = NULL, *dev_Ha_z  = NULL;
static double *dev_Hal_x  = NULL, *dev_Hal_y  = NULL, *dev_Hal_z  = NULL;
static double *dev_Hd_x  = NULL, *dev_Hd_y  = NULL, *dev_Hd_z  = NULL;
static double *dev_Hint_x  = NULL, *dev_Hint_y  = NULL, *dev_Hint_z  = NULL;

static double *dev_Happl_x = NULL, *dev_Happl_y = NULL, *dev_Happl_z = NULL;
static double *dev_D = NULL, *dev_T = NULL;

static double *dev_a_theta = NULL, *dev_b_theta = NULL, *dev_c_theta = NULL, *dev_d_theta = NULL,
			  *dev_a_phi = NULL, *dev_b_phi = NULL, *dev_c_phi = NULL, *dev_d_phi = NULL;




static double *dev_Aex_XP = NULL, *dev_Aex_XM = NULL, *dev_Aex_YP = NULL, *dev_Aex_YM = NULL, *dev_Aex_ZP = NULL, *dev_Aex_ZM = NULL;  //zyliu
static double *dev_Ms_XP = NULL, *dev_Ms_XM = NULL, *dev_Ms_YP = NULL, *dev_Ms_YM = NULL, *dev_Ms_ZP = NULL, *dev_Ms_ZM = NULL; //zyliu
static double *dev_H1_x = NULL, *dev_H1_y = NULL, *dev_H1_z = NULL, *dev_H2_x = NULL, *dev_H2_y = NULL, *dev_H2_z = NULL, //zyliu
			  *dev_H3_x = NULL, *dev_H3_y = NULL, *dev_H3_z = NULL, *dev_H4_x = NULL, *dev_H4_y = NULL, *dev_H4_z = NULL; //zyliu





static double *dev_d_theta_d_t = NULL, *dev_d_phi_d_t = NULL;
static double *dev_Mx = NULL, *dev_My = NULL, *dev_Mz = NULL;
static double *dev_M_temp_x = NULL, *dev_M_temp_y = NULL, *dev_M_temp_z = NULL;

static int    *dev_indicator1 = NULL, *dev_indicator1_temp = NULL, 
			  *dev_indicator2 = NULL, *dev_indicator2_temp = NULL,  // Voronoi-cell inner point indicator 
              *dev_indicator3 = NULL;


static double *dev_watch2 = NULL;


// Array for host
static double *host_Aex = NULL, *host_Ku = NULL, *host_Ms = NULL, *host_alpha = NULL, *host_gamma = NULL;

extern int     *indicator1, *indicator2, *indicator3, *indicator7;  // indicator1 = is_in_polygon; indicator2 = which_polygon
static double  *indicator4 = NULL, *indicator5 = NULL, *indicator6 = NULL;

static double  *T = NULL, *D = NULL;
extern double  *Happl_x, *Happl_y, *Happl_z;
static double  *Hint_x = NULL, *Hint_y = NULL, *Hint_z = NULL; 
static double  Happl;

extern double  *Mx, *My, *Mz; 
static double  *watch100 = NULL, *watch2 = NULL, *watch3 = NULL, *watch4 = NULL, *watch5 = NULL;

static float   *Mx_float = NULL, *My_float = NULL, *Mz_float = NULL;
static double  *Mx_t_bar = NULL, *My_t_bar = NULL, *Mz_t_bar = NULL;
static int     *watch100_int = NULL, *watch2_int = NULL, *watch3_int = NULL;


extern double  *Mx_bar1, *My_bar1, *Mz_bar1, 
               *Mx_bar2, *My_bar2, *Mz_bar2,
			   *Mx_bar3, *My_bar3, *Mz_bar3,
			   *Mx_bar4, *My_bar4, *Mz_bar4,
			   *Mx_bar5, *My_bar5, *Mz_bar5,
			   *Mx_bar6, *My_bar6, *Mz_bar6,
			   *Mx_bar,  *My_bar,  *Mz_bar;
static double  *M_bar = NULL,
			   *Mz_LayerLayer = NULL, *M1_z_SingleGrain = NULL, 
			   *Hint12_x_SingleGrain = NULL, *Hint12_y_SingleGrain = NULL, *Hint12_z_SingleGrain = NULL;

static double  *std_Ku = NULL, *std_Aex = NULL, *std_Tc = NULL;
static double  *gasarray = NULL;


/*static double x[(TOTAL_TIME>>1)+1], y[(TOTAL_TIME>>1)+1];
static double a, b, xbar, ybar, xsqr_sum, xy_sum, std_torq, std_Keff, Keff,
			  sin_2theta_bar, torque_pp_bar, sin_2theta_bar_temp, torque_pp_bar_temp,	
			  Ku_bar_t, torque_pp_bar_t, sin_2theta_bar_t;*/
static double M_bar_t, *Mx_bar_t, My_bar_t, Mz_bar_t;

// Head moving variable; head_dist: head's moving distance; NextBit: used to flip the applied field.
static int HeadDist = -1, NextBit = 0;    





// FFT on Hms
// Arrays
/*static float   *Mx_1d   = NULL,  *My_1d   = NULL,  *Mz_1d   = NULL,
*Gxx_1d  = NULL,  *Gxy_1d  = NULL,  *Gxz_1d  = NULL,
*Gyx_1d  = NULL,  *Gyy_1d  = NULL,  *Gyz_1d  = NULL,
*Gzx_1d  = NULL,  *Gzy_1d  = NULL,  *Gzz_1d  = NULL;*/


extern complex<float>  *Hd_x_1d_cmplx, *Hd_y_1d_cmplx, *Hd_z_1d_cmplx,
					   *Mx_1d_cmplx,   *My_1d_cmplx,   *Mz_1d_cmplx,
                       *Gxx_1d_cmplx,  *Gxy_1d_cmplx,  *Gxz_1d_cmplx, 
                       *Gyx_1d_cmplx,  *Gyy_1d_cmplx,  *Gyz_1d_cmplx, 
                       *Gzx_1d_cmplx,  *Gzy_1d_cmplx,  *Gzz_1d_cmplx;

extern double *Hd_x_1d, *Hd_y_1d, *Hd_z_1d;

extern double *Hd_x_1d_shift, *Hd_y_1d_shift, *Hd_z_1d_shift;


extern cufftComplex *dev_Gxx_cufft,  *dev_Gxy_cufft,  *dev_Gxz_cufft,
                    *dev_Gyx_cufft,  *dev_Gyy_cufft,  *dev_Gyz_cufft,
                    *dev_Gzx_cufft,  *dev_Gzy_cufft,  *dev_Gzz_cufft,
                    *dev_Mx_cufft,   *dev_My_cufft,   *dev_Mz_cufft,
                    *dev_Hd_x_cufft, *dev_Hd_y_cufft, *dev_Hd_z_cufft;



// Variables defined for recording simulation
static double *dev_Mag0_z = NULL;
extern double *memTemp0, *memTemp1, *memTemp2, *memTemp3,
	          *Happl_x_temp, *Happl_y_temp, *Happl_z_temp;
extern double *theta0;
extern double *mGrid[];
extern double *fGrid[];
extern double *FP_inp[], *FP_trail[], *FP[];
extern double  *FP_theta, *FP_mag;
extern double *Mag0[];
extern double *fSeq;
extern int *idx_f, *idx_t; //t_x0 is a reference x coordinate that indicates the starting trailing edge of the temperture profile;

extern int *Seq_f;

// Variables defined for cluster size calculation
extern double *rrt_cor;
extern double **Mz_tt;



#endif
