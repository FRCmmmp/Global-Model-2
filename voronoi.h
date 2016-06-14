//------------------------------------------------------------------------------
//
//         Header file for voronoi media module
//
//         (c) P W Huang, Seagate Technology (2016). All rights reserved.
//
//------------------------------------------------------------------------------
#include <vector>
struct vec2d_t{
  double x;
  double y;
};

void granular_media(vec2d_t media_size, int num_grains, double min_distance, double shrink_dist, double shrink_factor, std::vector<vec2d_t>& seed_array, std::vector <std::vector <vec2d_t> > &  vertex_array);
void sro_media(vec2d_t media_size, int num_large_grains, double grain_size, double grain_scaling, double domain_bnd_qfactor, std::vector<vec2d_t>& seed_array, std::vector <std::vector <vec2d_t> > &  vertex_array);