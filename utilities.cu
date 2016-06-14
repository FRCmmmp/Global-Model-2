  //------------------------------------------------------------------------------
  //
  //         A script that contains some useful functions for general purposes
  //
  //         (c) P W Huang, Seagate Technology (2016). All rights reserved.
  //
  //------------------------------------------------------------------------------


#include <cstdlib>

// Only valid for 32-bit integer variables
// Ideal to be "unsigned int"
int Rnd_upto_pow2(int v){
    if (v >= 0){
        v--;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        v++;
        return v;
    }
    else return -1;
}