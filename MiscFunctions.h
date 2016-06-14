#ifndef __MiscFunctions__
#define __MiscFunctions__

#include "Parameters.h"
#include "Parameters_input.h"
#include "Outputs.h"
#include <time.h>
#include <vector>
#include <cstdlib>


struct vec2d_t{
  double x;
  double y;
};

struct coor2d_t{
	int x;
	int y;
};

extern bool CalcMagLayerAve(long int, int, double*, std::vector<std::vector<coor2d_t> > &, int);
extern bool External_HeadField(int);
extern bool Set_delta_K_Aex_Tc(double*, double*, double*, double*, double*, double*);
extern bool Set_Happl_DT();
extern bool Set_Happl_CT();
extern bool Reset_Happl();
extern bool Jitter_Calc();
extern bool WPE_Calc();
extern bool sigHc_Calc(std::vector<std::vector<double> >&, double*, int, int, std::vector<std::vector<coor2d_t> > & );
extern bool FP_inp_theta_n_Magnitude();
extern bool Set_devMxyz(double* host_Ms, double* dev_theta, double* dev_phi, double* dev_Mx, double* dev_My, double* dev_Mz);
extern bool FieldReset();

#endif
