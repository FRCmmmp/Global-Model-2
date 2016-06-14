///========================================================================================================
///                                  Function to decide if point is within polygon
///                                                      Version 2.0
///                                                  R F Evans 11/09/2012
///========================================================================================================
#include "internal.h"

bool point_in_polygon(double x, double y, std::vector<vec2d_t>& poly){

        x+=1e-10; // Added tiny amount to coordinates to include points at 0,0,0
        y+=1e-10;

        int j=poly.size()-1 ;
        bool oddNodes=false;

        for (int i=0; i<poly.size(); i++) {
                if ((poly[i].y<y && poly[j].y>=y) || (poly[j].y<y && poly[i].y>=y)) {
                        if (poly[i].x+(y-poly[i].y)/(poly[j].y-poly[i].y)*(poly[j].x-poly[i].x)<x) {
                                oddNodes=!oddNodes;
                        }
                }
                j=i;
        }

  return oddNodes;

}