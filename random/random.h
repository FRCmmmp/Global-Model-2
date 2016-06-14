
//
#ifndef RANDOM_H_
#define RANDOM_H_
#include "mtrand.h"
namespace mtrandom
//==========================================================
// Namespace mtrandom
//==========================================================
{
	extern MTRand grnd; /// single sequence of random numbers
	extern double gaussian();
	extern double gaussianc(MTRand&);

}


#endif /*RANDOM_H_*/
