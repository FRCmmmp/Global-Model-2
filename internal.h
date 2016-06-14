//------------------------------------------------------------------------------
//
//         Header file for voronoi media internal functions
//
//         (c) P W Huang, Seagate Technology (2016). All rights reserved.
//
//------------------------------------------------------------------------------
#include <vector>

#include "voronoi.h"

// internal functions only
void populate_vertex_points(std::vector <vec2d_t> & grain_coord_array, std::vector <std::vector < vec2d_t> > &  grain_vertices_array);
bool point_in_polygon(double x, double y, std::vector<vec2d_t>& poly);
bool convert_to_old_voro_output_format(std::vector<vec2d_t>& poly);
bool polygon_area(std::vector <std::vector < vec2d_t> > &  grain_vertices_array, std::vector<double> &   grain_area_array);