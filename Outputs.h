#ifndef __OUTPUTS__
#define __OUTPUTS__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>




// Functions
extern bool Output_Float_3D_Format(int, int, int, double*, char*);
extern bool Output_Int_3D_Format(int, int, int, int*, char*);
extern bool Output_Float_1D_Format_6col(int, int, int, double*, double*, double*, char*);
extern bool Output_Float_1D_Format_1col(int, int, int, double*, char*);
extern bool Output_Int_1D_Format_1col(int, int, int, int*, char*);
extern bool Input_Float_1D_Format_1col(int, int, int, double*, char*);
extern bool Input_Indicator_1D_Format(int, int, int, int*, int*, int*, double*, double*, double*, int*, char*);

#endif