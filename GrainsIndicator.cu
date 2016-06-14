#include "GrainsIndicator.h"
#include "voronoi.h"

using namespace std;






int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}
//Arguments:
//nvert: Number of vertices in the polygon. Whether to repeat the first vertex at the end has been discussed in the article referred above.
// vertx, verty: Arrays containing the x- and y-coordinates of the polygon's vertices.
// testx, testy: X- and y-coordinate of the test point.
//It's both short and efficient and works both for convex and concave polygons. As suggested before, you should check the bounding rectangle 
// first and treat polygon holes separately.
//The idea behind this is pretty simple. The author describes it as follows:
//I run a semi-infinite ray horizontally (increasing x, fixed y) out from the test point, and count how many edges it crosses. At each crossing, 
// the ray switches between inside and outside. This is called the Jordan curve theorem.


void randperm(int n,int perm[])
{
	int i, j, t;

	for(i=0; i<n; i++)
		perm[i] = i;
	for(i=0; i<n; i++) {
		j = rand()%(n-i)+i;
		t = perm[j];
		perm[j] = perm[i];
		perm[i] = t;
	}
}
/* 
* randperm.c
*
* Copyright (c) 1990 Michael E. Hohmeyer,
*       hohmeyer@icemcfd.com
* Permission is granted to modify and re-distribute this code in any manner
* as long as this notice is preserved.  All standard disclaimers apply.
*       
*/


enum voronoi_t {normal=0, sro=1};
bool GrainsIndicator(){
	
	#ifndef __VORONOI__
	if (VORO_GRAIN){
		double x_seed, y_seed, x_vertex, y_vertex;
		int NumOfSide_voro;
		double voro_vertex_x, voro_vertex_y;
		int *NumOfSide_voro_array = NULL;
		double *voro_vertex_array[2];
		double grid_size_x, grid_size_y, grid_size_z;
		double medium_x, medium_y, medium_z, medium_z1, medium_z12, medium_z2, medium_z23, medium_z3, medium_z34, medium_z4, medium_z45, medium_z5, medium_z56, medium_z6;
		int grid_x, grid_y, grid_z, grid_z1, grid_z12, grid_z2, grid_z23, grid_z3, grid_z34, grid_z4, grid_z45, grid_z5, grid_z56, grid_z6;
		double *vertex_x = NULL, *vertex_y = NULL, *vertex_z = NULL;
		double *grid_vertex[2], *grid_center[4];  // grid_vertex[2]'s (x,y) coordinates; grid_center[4]'s (x,y,logic,poly_lable)
		int idx;
		int memSize_out2, memSize_out3;
		int *logic_idx = NULL;
		int *poly_lable = NULL;
		int start_point;
		double *x_poly = NULL, *y_poly = NULL;
		int idx_temp;
		int *AC_erased = NULL;
		int *DEF_array = NULL;
		double *indicator[7];
		double *R1 = NULL, *R2 = NULL, *R2_cap = NULL, *R3 = NULL;
		
		// PRN Generator (Linux)
		MTRand grnd;
		double *indicator3D[7];
		double AC_erased_ratio;

		



		//------------------------------------------------------------------------------
		//
		//         Function to generate random grain structure by Richard Evans
		//
		//         (c) P W Huang, Seagate Technology (2016). All rights reserved.
		//
		//------------------------------------------------------------------------------
		// define voronoi type
		voronoi_t voronoi_type;
		switch (media_type){
		case 0:
			voronoi_type = normal;
			break;
		case 1:
			voronoi_type = sro;
			break;
		}

		// define structures for voronoi data
		std::vector <vec2d_t> seed_array;
		std::vector <std::vector <vec2d_t> > vertex_array;

		// define media size and number of grains
		vec2d_t media_size;
		media_size.x = mNx*delta_x*1e7;
		media_size.y = mNy*delta_y*1e7;
		const int num_grains = num_vgrains;                             //number of "grains" for granular_media() 
		const int num_domains = num_vdomains;                           //number of "domains" for sro_media()
		const double grain_size = rgrain_size*delta_x*1e7;              //grain size for sro_media()
		const double minimum_distance = minimum_sdistance*delta_x*1e7;  //for granular_media()
		const double shrink_dist = -grain_bnd/2*delta_x*1e7;            //for grain boundaries distance
	    const double shrink_factor = 0.85;                              //for grain boundaries factor

		// generate seed points
		switch(voronoi_type){
			case normal:
				granular_media(media_size, num_grains, minimum_distance, shrink_dist, shrink_factor, seed_array, vertex_array);
				break;

			case sro:
				sro_media(media_size, num_domains, grain_size, rgrain_scaling, domain_bnd_qfactor, seed_array, vertex_array);
				break;

			default:
				granular_media(media_size, num_grains, minimum_distance, shrink_dist, shrink_factor, seed_array, vertex_array);
				break;
		}

		// Set voro_vertex_array[];
		ifstream fin2, fin3;
		fin2.open("out2.out",ios::in);
		assert(!fin2.fail());
		fin3.open("out3.out",ios::in);
		assert(!fin3.fail());
		memSize_out2 = -1;
		while (!fin2.eof()){
			fin2 >> NumOfSide_voro;
			memSize_out2++;
		}
		fin2.close();

		// Set NumOfSide_voro_array[]
		NumOfSide_voro_array = (int*) calloc(memSize_out2, sizeof(int));
		memSize_out3 = -1;
		while (!fin3.eof()){
			fin3 >> voro_vertex_x >> voro_vertex_y;
			memSize_out3++;
		}
		fin3.close();

		voro_vertex_array[0] = (double*) calloc(memSize_out3, sizeof(double));
		voro_vertex_array[1] = (double*) calloc(memSize_out3, sizeof(double));

		// Read-in out2 and out3 files for vertices and NumOfSide arrays
		fin2.open("out2.out",ios::in);
		assert(!fin2.fail());
		fin3.open("out3.out",ios::in);
		assert(!fin3.fail());
		idx = -1;
		while (!fin2.eof()){
			idx++;
			fin2 >> NumOfSide_voro_array[idx];
		}
		fin2.close();
		printf("Loaded the files 'out2.out, out3.out' successfully \n");

		idx = 0;
		while (fin3 >> voro_vertex_array[0][idx] >> voro_vertex_array[1][idx]){
			voro_vertex_array[0][idx] = voro_vertex_array[0][idx];
			voro_vertex_array[1][idx] = voro_vertex_array[1][idx];
			idx++;
		}

									
		//-------- Create 2D grids --------//
		medium_x   = mNx*delta_x    * 1e7; // in [nm]
		medium_y   = mNy*delta_y    * 1e7; // in [nm]
		medium_z   = mNz*delta_z    * 1e7; // in [nm]
		medium_z1  = mNz_1*delta_z  * 1e7; // in [nm]
		medium_z12 = mNz_12*delta_z * 1e7; // in [nm]
		medium_z2  = mNz_2*delta_z  * 1e7; // in [nm]
		medium_z23 = mNz_23*delta_z * 1e7; // in [nm]
		medium_z3  = mNz_3*delta_z  * 1e7; // in [nm]
		medium_z34 = mNz_34*delta_z * 1e7; // in [nm]
		medium_z4  = mNz_4*delta_z  * 1e7; // in [nm]
		medium_z45 = mNz_45*delta_z * 1e7; // in [nm]
		medium_z5  = mNz_5*delta_z  * 1e7; // in [nm]
		medium_z56 = mNz_56*delta_z * 1e7; // in [nm]
		medium_z6  = mNz_6*delta_z  * 1e7; // in [nm]
		medium_z   = mNz*delta_z    * 1e7; // in [nm]

		grid_size_x = delta_x * 1e7; // in [nm]
		grid_size_y = delta_y * 1e7; // in [nm]
		grid_size_z = delta_z * 1e7; // in [nm]

		// "plus 0.5" is for truncating issue in cpp compiler
		grid_x   = (int)(medium_x/grid_size_x    + 0.5);
		grid_y   = (int)(medium_y/grid_size_y    + 0.5);
		grid_z   = (int)(medium_z/grid_size_z    + 0.5);
		grid_z1  = (int)(medium_z1/grid_size_z   + 0.5);
		grid_z12 = (int)(medium_z12/grid_size_z  + 0.5);
		grid_z2  = (int)(medium_z2/grid_size_z   + 0.5);
		grid_z23 = (int)(medium_z23/grid_size_z  + 0.5);
		grid_z3  = (int)(medium_z3/grid_size_z   + 0.5);
		grid_z34 = (int)(medium_z34/grid_size_z  + 0.5);
		grid_z4  = (int)(medium_z4/grid_size_z   + 0.5);
		grid_z45 = (int)(medium_z45/grid_size_z  + 0.5);
		grid_z5  = (int)(medium_z5/grid_size_z   + 0.5);
		grid_z56 = (int)(medium_z56/grid_size_z  + 0.5);
		grid_z6  = (int)(medium_z6/grid_size_z   + 0.5);
		
		vertex_x = (double*) calloc(grid_x+1, sizeof(double));
		vertex_y = (double*) calloc(grid_y+1, sizeof(double));
		vertex_z = (double*) calloc(grid_z+1, sizeof(double));
		for (int i = 0; i < grid_x+1; i++){ vertex_x[i] = i*grid_size_x; }
		for (int i = 0; i < grid_y+1; i++){ vertex_y[i] = i*grid_size_y; }
		for (int i = 0; i < grid_z+1; i++){ vertex_z[i] = i*grid_size_z; }

		grid_vertex[0] = (double*) calloc((grid_x+1)*(grid_y+1), sizeof(double));
		grid_vertex[1] = (double*) calloc((grid_x+1)*(grid_y+1), sizeof(double));
		grid_center[0] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		grid_center[1] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		grid_center[2] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		grid_center[3] = (double*) calloc((grid_x)*(grid_y), sizeof(double));

		idx = -1;
		for (int i = 0; i < (grid_x+1); i++){
			for (int j = 0; j < (grid_y+1); j++){
				idx = idx + 1;
				grid_vertex[0][idx] = i*grid_size_x;
				grid_vertex[1][idx] = j*grid_size_y;
			}
		}
		idx = -1;
		for (int i = 0; i < (grid_x); i++){
			for (int j = 0; j < (grid_y); j++){
				idx = idx + 1;
				grid_center[0][idx] = grid_size_x/2 + i*grid_size_x;
				grid_center[1][idx] = grid_size_y/2 + j*grid_size_y;
			}
		}



		///--------------- Mapping with 2D grids-------------------///
		logic_idx = (int*) calloc(grid_x*grid_y, sizeof(int));
		poly_lable = (int*) calloc(grid_x*grid_y, sizeof(int));
		start_point = 0;
		
		// i_poly: ith polygon
		for (int i_poly = 0; i_poly < memSize_out2; i_poly++){
			x_poly = (double*) calloc(NumOfSide_voro_array[i_poly], sizeof(double));
			y_poly = (double*) calloc(NumOfSide_voro_array[i_poly], sizeof(double));
			idx_temp = 0;

			// v_poly: vertices of ith polygon
			for (int v_poly = start_point; v_poly < (start_point+NumOfSide_voro_array[i_poly]); v_poly++){
				x_poly[idx_temp] = voro_vertex_array[0][v_poly];
				y_poly[idx_temp] = voro_vertex_array[1][v_poly];
				idx_temp++;
			}
			for (int i = 0; i < grid_x*grid_y; i++){
				if (pnpoly(NumOfSide_voro_array[i_poly], x_poly, y_poly, grid_center[0][i], grid_center[1][i])){ 
					logic_idx[i] = 1;
					poly_lable[i] = i_poly+1;
				}
			}
			start_point = start_point+NumOfSide_voro_array[i_poly];
		}
		for (int i = 0; i < grid_x*grid_y; i++){
			grid_center[2][i] = logic_idx[i];
			grid_center[3][i] = poly_lable[i];
		}

		AC_erased = (int*) calloc(memSize_out2, sizeof(int));
		idx = 0;
		for (int i = 0; i < memSize_out2; i++){
			idx++;
			AC_erased[i] = idx;
		}
		randperm(memSize_out2, AC_erased);
		if (AC_DC_erased == 1) AC_erased_ratio = 0.5; else AC_erased_ratio = 1.0;
		for (int i = 0; i < memSize_out2; i++){
			if (AC_erased[i] <= floor(memSize_out2*AC_erased_ratio)){  // Define AC-erased ratio. E.g., 0.5-->50% up and 50% down, 0.2-->20% up and 80% down
				AC_erased[i] = 1;
			}
			else {
				AC_erased[i] = 2;
			}
		}
		for (int i = 0; i < 7; i++){
			indicator[i] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		}
		idx = -1;
		for (int k = 0; k < grid_y; k++){
			for (int j = k; j < grid_x*grid_y; j = j + grid_y){
				idx++;
				indicator[0][idx] = grid_center[2][j];
				indicator[1][idx] = grid_center[3][j];
				if (indicator[1][idx] != 0){
					indicator[2][idx] = AC_erased[(int)indicator[1][idx]-1];
				}
			}
		}
	
		R1 = (double*) calloc(memSize_out2, sizeof(double));
		R2 = (double*) calloc(memSize_out2, sizeof(double));
		R2_cap = (double*) calloc(memSize_out2, sizeof(double));
		R3 = (double*) calloc(memSize_out2, sizeof(double));
	
		grnd.seed(iseed);
		for (int i = 0; i < memSize_out2; i++){
			R1[i] = mtrandom::gaussian();
			R2[i] = mtrandom::gaussian();
			R2_cap[i] = mtrandom::gaussian();
			R3[i] = mtrandom::gaussian();
		}
		idx = -1;
		for (int k = 0; k < grid_y; k++){
			for (int j = k; j < grid_x*grid_y; j = j + grid_y){
				idx++;
				if (indicator[1][idx] != 0){
					indicator[3][idx] = R1[(int)indicator[1][idx]-1];
					indicator[4][idx] = R2[(int)indicator[1][idx]-1];
					indicator[5][idx] = R3[(int)indicator[1][idx]-1];
				}
			}
		}
	
		for (int i = 0; i < 7; i++){
			indicator3D[i] = (double*) calloc((grid_x)*(grid_y)*(grid_z), sizeof(double));
		}
		
		idx = -1;
		// Layer 1
		for (int k = 0; k < grid_z1; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*1;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 12
		for (int k = grid_z1; k < grid_z12; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*12;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 2
		for (int k = grid_z12; k < grid_z2; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*2;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 23
		for (int k = grid_z2; k < grid_z23; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*23;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 3
		for (int k = grid_z23; k < grid_z3; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*3;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 34
		for (int k = grid_z3; k < grid_z34; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*34;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 4
		for (int k = grid_z34; k < grid_z4; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*4;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 45
		for (int k = grid_z4; k < grid_z45; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*45;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 5
		for (int k = grid_z45; k < grid_z5; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*5;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 56
		for (int k = grid_z5; k < grid_z56; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*56;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 6
		for (int k = grid_z56; k < grid_z6; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*6;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer zero-padding
		for (int k = grid_z6; k < grid_z; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = 7;
				indicator3D[1][idx] = 7000;
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}

		// Setup random defects in CGC (not useful here)
		if (CGC_DEF == 1){
			DEF_array = (int*) calloc(memSize_out2, sizeof(int));
			idx = 0;
			for (int i = 0; i < memSize_out2; i++){
				idx++;
				DEF_array[i] = idx;
			}
			randperm(memSize_out2, DEF_array);
			for (int i = 0; i < memSize_out2; i++){
				if (DEF_array[i] <= int(memSize_out2*def_perc + 0.5)){  
					DEF_array[i] = 1;
				}
				else {
					DEF_array[i] = 0;
				}
			}
			idx = -1;
			for (int k = 0; k < grid_y; k++){
				for (int j = k; j < grid_x*grid_y; j = j + grid_y){
					idx++;
					indicator[0][idx] = grid_center[2][j];
					indicator[1][idx] = grid_center[3][j];
					if (indicator[1][idx] != 0){
						indicator[6][idx] = DEF_array[(int)indicator[1][idx]-1];
					}
				}
			}
			
			// Set random defects in CGC
			switch (CGC_label){
			case 1:
				for (int k=0; k<grid_z1; k++){
					for (int j=0; j<grid_y; j++){
						for (int i=0; i<grid_x; i++){
							indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = indicator[6][i+j*grid_x];
						}
					}
				}
				break;
			case 2:
				for (int k=grid_z12; k<grid_z2; k++){
					for (int j=0; j<grid_y; j++){
						for (int i=0; i<grid_x; i++){
							indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = indicator[6][i+j*grid_x];
						}
					}
				}
				break;
			case 3:
				for (int k=grid_z23; k<grid_z3; k++){
					for (int j=0; j<grid_y; j++){
						for (int i=0; i<grid_x; i++){
							indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = indicator[6][i+j*grid_x];
						}
					}
				}
				break;
			case 4:
				for (int k=grid_z3; k<grid_z34; k++){
					for (int j=0; j<grid_y; j++){
						for (int i=0; i<grid_x; i++){
							indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = indicator[6][i+j*grid_x];
						}
					}
				}
				break;
			case 5:
				for (int k=grid_z45; k<grid_z5; k++){
					for (int j=0; j<grid_y; j++){
						for (int i=0; i<grid_x; i++){
							indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = indicator[6][i+j*grid_x];
						}
					}
				}
				break;
			case 6:
				for (int k=grid_z56; k<grid_z6; k++){
					for (int j=0; j<grid_y; j++){
						for (int i=0; i<grid_x; i++){
							indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = indicator[6][i+j*grid_x];
						}
					}
				}
				break;
			}
		}

		// Write to file "indicator.inp"
		FILE *p1;
		p1 = fopen("indicator.inp", "w");
		for (int i = 0; i < grid_x*grid_y*grid_z; i++){
			fprintf(p1, "%5d %5d %5d %10.2lf %10.4lf %10.4lf %5d\n", (int)indicator3D[0][i], (int)indicator3D[1][i], (int)indicator3D[2][i], indicator3D[3][i], indicator3D[4][i], indicator3D[5][i], (int)indicator3D[6][i]);
		}
		fclose(p1);
		printf("Read-in out2, out3, and indicator.inp generation succeeded!\n");
		printf("out2 ----> number of sides for each grain\n");
		printf("out3 ----> vertices for each grain\n");

	}
	#endif

	#ifndef __UNIFORM__
	else {
		
		int NumOfSide_uni = 4;
		int n_X, n_Y, n_Z;  // number of grains in x- and y-direction;
		double x0, y0;
		int idx;
		int memSize_out2, memSize_out3;
		int *NumOfSide_uni_array = NULL;
		double *uni_vertex_array[2];
		double grid_size_x, grid_size_y, grid_size_z;
		double medium_x, medium_y, medium_z, medium_z1, medium_z12, medium_z2, medium_z23, medium_z3, medium_z34, medium_z4, medium_z45, medium_z5, medium_z56, medium_z6;
		int grid_x, grid_y, grid_z, grid_z1, grid_z12, grid_z2, grid_z23, grid_z3, grid_z34, grid_z4, grid_z45, grid_z5, grid_z56, grid_z6;
		double *vertex_x = NULL, *vertex_y = NULL, *vertex_z = NULL;
		double *grid_vertex[2], *grid_center[4];  // grid_vertex[2]'s (x,y) coordinates; grid_center[4]'s (x,y,logic,poly_lable)
		int *logic_idx = NULL;
		int *poly_lable = NULL;
		int start_point;
		double *x_poly = NULL, *y_poly = NULL;
		int idx_temp;
		int *AC_erased = NULL;
		int *DEF_array = NULL;
		double *indicator[7];
		double *R1 = NULL, *R2 = NULL, *R2_cap = NULL, *R3 = NULL;
		MTRand grnd;
		double *indicator3D[7];
		double AC_erased_ratio;


		//-------- Create 2D Media Cells --------//
		medium_x  = mNx*delta_x * 1e7; // in [nm]
		medium_y  = mNy*delta_y * 1e7; // in [nm]
		medium_z  = mNz*delta_z * 1e7; // in [nm]
		medium_z1  = mNz_1*delta_z * 1e7; // in [nm]
		medium_z12 = mNz_12*delta_z * 1e7; // in [nm]
		medium_z2 = mNz_2*delta_z * 1e7; // in [nm]
		medium_z23 = mNz_23*delta_z * 1e7; // in [nm]
		medium_z3 = mNz_3*delta_z * 1e7; // in [nm]
		medium_z34 = mNz_34*delta_z * 1e7; // in [nm]
		medium_z4 = mNz_4*delta_z * 1e7; // in [nm]
		medium_z45 = mNz_45*delta_z * 1e7; // in [nm]
		medium_z5 = mNz_5*delta_z * 1e7; // in [nm]
		medium_z56 = mNz_56*delta_z * 1e7; // in [nm]
		medium_z6 = mNz_6*delta_z * 1e7; // in [nm]

		
		
		n_X = (int)(medium_x*1e-7/GrainDx + 0.5);
		n_Y = (int)(medium_y*1e-7/GrainDy + 0.5);
		//n_Z = (int)(medium_z*1e-7/GrainDz);
		memSize_out2 = n_X*n_Y;
		memSize_out3 = n_X*n_Y*4;

		
		uni_vertex_array[0] = (double*) calloc(memSize_out3, sizeof(double)); //vertices coordinates
		uni_vertex_array[1] = (double*) calloc(memSize_out3, sizeof(double)); //vertices coordinates
		NumOfSide_uni_array = (int*) calloc(memSize_out2, sizeof(int));

		idx = 0;
		for (int j = 0; j < n_Y; j++){
			y0 = j*GrainDy*1e7;
			for (int i = 0; i < n_X; i++){
				x0 = i*GrainDx*1e7;
				uni_vertex_array[0][idx]   = x0 + 0*GrainDx*1e7;
				uni_vertex_array[0][idx+1] = x0 + 1*GrainDx*1e7;
				uni_vertex_array[0][idx+2] = x0 + 1*GrainDx*1e7;
				uni_vertex_array[0][idx+3] = x0 + 0*GrainDx*1e7;
				uni_vertex_array[1][idx]   = y0 + 0*GrainDy*1e7;
				uni_vertex_array[1][idx+1] = y0 + 0*GrainDy*1e7;
				uni_vertex_array[1][idx+2] = y0 + 1*GrainDy*1e7;
				uni_vertex_array[1][idx+3] = y0 + 1*GrainDy*1e7;
				idx = idx + 4;
			}
		}
		for (int i = 0; i < memSize_out2; i++){
			NumOfSide_uni_array[i] = 4;
		}
		


		
		
		
		//-------- Create 2D grids --------//
		grid_size_x = delta_x * 1e7; // in [nm]
		grid_size_y = delta_y * 1e7; // in [nm]
		grid_size_z = delta_z * 1e7; // in [nm]

		grid_x   = (int)(medium_x/grid_size_x    + 0.5); // "plus 0.5" is for truncating issue in cpp compiler
		grid_y   = (int)(medium_y/grid_size_y    + 0.5);
		grid_z   = (int)(medium_z/grid_size_z    + 0.5);
		grid_z1  = (int)(medium_z1/grid_size_z   + 0.5);
		grid_z12 = (int)(medium_z12/grid_size_z  + 0.5);
		grid_z2  = (int)(medium_z2/grid_size_z   + 0.5);
		grid_z23 = (int)(medium_z23/grid_size_z  + 0.5);
		grid_z3  = (int)(medium_z3/grid_size_z   + 0.5);
		grid_z34 = (int)(medium_z34/grid_size_z  + 0.5);
		grid_z4  = (int)(medium_z4/grid_size_z   + 0.5);
		grid_z45 = (int)(medium_z45/grid_size_z  + 0.5);
		grid_z5  = (int)(medium_z5/grid_size_z   + 0.5);
		grid_z56 = (int)(medium_z56/grid_size_z  + 0.5);
		grid_z6  = (int)(medium_z6/grid_size_z   + 0.5);
		
		vertex_x = (double*) calloc(grid_x+1, sizeof(double));
		vertex_y = (double*) calloc(grid_y+1, sizeof(double));
		vertex_z = (double*) calloc(grid_z+1, sizeof(double));
		for (int i = 0; i < grid_x+1; i++){ vertex_x[i] = i*grid_size_x; }
		for (int i = 0; i < grid_y+1; i++){ vertex_y[i] = i*grid_size_y; }
		for (int i = 0; i < grid_z+1; i++){ vertex_z[i] = i*grid_size_z; }

		grid_vertex[0] = (double*) calloc((grid_x+1)*(grid_y+1), sizeof(double));
		grid_vertex[1] = (double*) calloc((grid_x+1)*(grid_y+1), sizeof(double));
		grid_center[0] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		grid_center[1] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		grid_center[2] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		grid_center[3] = (double*) calloc((grid_x)*(grid_y), sizeof(double));

		idx = -1;
		for (int i = 0; i < (grid_x+1); i++){
			for (int j = 0; j < (grid_y+1); j++){
				idx = idx + 1;
				grid_vertex[0][idx] = i*grid_size_x;
				grid_vertex[1][idx] = j*grid_size_y;
			}
		}
		idx = -1;
		for (int i = 0; i < (grid_x); i++){
			for (int j = 0; j < (grid_y); j++){
				idx = idx + 1;
				grid_center[0][idx] = grid_size_x/2 + i*grid_size_x;
				grid_center[1][idx] = grid_size_y/2 + j*grid_size_y;
			}
		}



		///--------------- Mapping with 2D grids-------------------///
		logic_idx = (int*) calloc(grid_x*grid_y, sizeof(int));
		poly_lable = (int*) calloc(grid_x*grid_y, sizeof(int));
		start_point = 0;
		for (int i_poly = 0; i_poly < memSize_out2; i_poly++){                                  // i_poly: ith polygon
			x_poly = (double*) calloc(NumOfSide_uni_array[i_poly], sizeof(double));
			y_poly = (double*) calloc(NumOfSide_uni_array[i_poly], sizeof(double));
			idx_temp = 0;
			for (int v_poly = start_point; v_poly < (start_point+NumOfSide_uni_array[i_poly]); v_poly++){    // v_poly: vertices of ith polygon
				x_poly[idx_temp] = uni_vertex_array[0][v_poly];
				y_poly[idx_temp] = uni_vertex_array[1][v_poly];
				idx_temp++;
			}
			for (int i = 0; i < grid_x*grid_y; i++){
				if (pnpoly(NumOfSide_uni_array[i_poly], x_poly, y_poly, grid_center[0][i], grid_center[1][i])){ 
					logic_idx[i] = 1;
					poly_lable[i] = i_poly+1;
				}
			}
			start_point = start_point+NumOfSide_uni_array[i_poly];
		}
		for (int i = 0; i < grid_x*grid_y; i++){
			grid_center[2][i] = logic_idx[i];
			grid_center[3][i] = poly_lable[i];
		}

		AC_erased = (int*) calloc(/*(grid_x)*(grid_y)*/memSize_out2, sizeof(int));
		idx = 0;
		for (int i = 0; i < /*grid_x*grid_y*/memSize_out2; i++){
			idx++;
			AC_erased[i] = idx;
		}
		randperm(memSize_out2, AC_erased);
		if (AC_DC_erased == 1){ AC_erased_ratio = 0.5; } else { AC_erased_ratio = 1.0; }
		for (int i = 0; i < memSize_out2; i++){
			if (AC_erased[i] <= floor(memSize_out2*AC_erased_ratio)){  
				AC_erased[i] = 1;
			}
			else {
				AC_erased[i] = 2;
			}
		}
		for (int i = 0; i < 7; i++){
			indicator[i] = (double*) calloc((grid_x)*(grid_y), sizeof(double));
		}
		idx = -1;
		for (int k = 0; k < grid_y; k++){
			for (int j = k; j < grid_x*grid_y; j = j + grid_y){
				idx++;
				indicator[0][idx] = grid_center[2][j];
				indicator[1][idx] = grid_center[3][j];
				if (indicator[1][idx] != 0){
					indicator[2][idx] = AC_erased[(int)indicator[1][idx]-1];
				}
			}
		}
	
	
	
		R1 = (double*) calloc(memSize_out2, sizeof(double));
		R2 = (double*) calloc(memSize_out2, sizeof(double));
		R2_cap = (double*) calloc(memSize_out2, sizeof(double));
		R3 = (double*) calloc(memSize_out2, sizeof(double));
	
		grnd.seed(iseed);
		for (int i = 0; i < memSize_out2; i++){
			R1[i]     = mtrandom::gaussian();
			R2[i]     = mtrandom::gaussian();
			R2_cap[i] = mtrandom::gaussian();
			R3[i]     = mtrandom::gaussian();
		}
		idx = -1;
		for (int k = 0; k < grid_y; k++){
			for (int j = k; j < grid_x*grid_y; j = j + grid_y){
				idx++;
				if (indicator[1][idx] != 0){
					indicator[3][idx] = R1[(int)indicator[1][idx]-1];
					indicator[4][idx] = R2[(int)indicator[1][idx]-1];
					indicator[5][idx] = R3[(int)indicator[1][idx]-1];
				}
			}
		}

		if (DT_Rec_Analysis){
			idx = 0;
			for (int j = 0; j < mNy; j++){
				for (int i = 0; i < mNx; i++){
					indicator[3][idx] = R1[j];
					indicator[4][idx] = R2[j];
					indicator[5][idx] = R3[j];
					idx++;
				}
			}
		}
		/*if (CT_Rec_Analysis){
			idx = 0;
			for (int i = 0; i < mNx; i++){
				for (int j = 0; j < mNy; j++){
					indicator[3][idx] = R1[i];
					indicator[4][idx] = R2[i];
					indicator[5][idx] = R3[i];
					idx++;
				}
			}
		}*/

	
		for (int i = 0; i < 7; i++){
			indicator3D[i] = (double*) calloc((grid_x)*(grid_y)*(grid_z), sizeof(double));
		}
		idx = -1;
		// Layer 1
		for (int k = 0; k < grid_z1; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*1;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 12
		for (int k = grid_z1; k < grid_z12; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*12;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 2
		for (int k = grid_z12; k < grid_z2; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*2;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 23
		for (int k = grid_z2; k < grid_z23; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*23;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 3
		for (int k = grid_z23; k < grid_z3; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*3;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 34
		for (int k = grid_z3; k < grid_z34; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*34;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 4
		for (int k = grid_z34; k < grid_z4; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*4;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 45
		for (int k = grid_z4; k < grid_z45; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*45;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 5
		for (int k = grid_z45; k < grid_z5; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*5;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 56
		for (int k = grid_z5; k < grid_z56; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*56;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer 6
		for (int k = grid_z56; k < grid_z6; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = indicator[0][j]*6;
				indicator3D[1][idx] = indicator[1][j];
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}
		// Layer zero-padding
		for (int k = grid_z6; k < grid_z; k++){
			for (int j = 0; j < grid_x*grid_y; j++){
				idx++;
				indicator3D[0][idx] = 7;
				indicator3D[1][idx] = 7000;
				indicator3D[2][idx] = indicator[2][j];
				indicator3D[3][idx] = indicator[3][j];
				indicator3D[4][idx] = indicator[4][j];
				indicator3D[5][idx] = indicator[5][j];
			}
		}

		// Setup random defects in CGC (not useful here)
		if (CGC_DEF == 1){
			DEF_array = (int*) calloc(memSize_out2, sizeof(int));
			idx = 0;
			for (int i = 0; i < memSize_out2; i++){
				idx++;
				DEF_array[i] = idx;
			}
			randperm(memSize_out2, DEF_array);
			for (int i = 0; i < memSize_out2; i++){
				if (DEF_array[i] <= int(memSize_out2*def_perc + 0.5)){  
					DEF_array[i] = 3;
				}
				else {
					DEF_array[i] = 4;
				}
			}
			idx = -1;
			for (int k = 0; k < grid_y; k++){
				for (int j = k; j < grid_x*grid_y; j = j + grid_y){
					idx++;
					indicator[0][idx] = grid_center[2][j];
					indicator[1][idx] = grid_center[3][j];
					if (indicator[1][idx] != 0){
						indicator[6][idx] = DEF_array[(int)indicator[1][idx]-1];
					}
				}
			}
			
			// Set uniform defects in CGC
			switch (CGC_label){
			case 1:
				for (int k = 0; k < grid_z1; k++){
					for (int j = 0; j < grid_y; j++){
						for (int i = 0; i < grid_x; i++){
							if (i%3 == 0 && j%3 == 0){
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 1;
							}
							else {
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 0;
							}
						}
					}
				}
				break;
			case 2:
				for (int k = grid_z12; k < grid_z2; k++){
					for (int j = 0; j < grid_y; j++){
						for (int i = 0; i < grid_x; i++){
							if (i%3 == 0 && j%3 == 0){
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 1;
							}
							else {
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 0;
							}
						}
					}
				}
				break;
			case 3:
				for (int k = grid_z23; k < grid_z3; k++){
					for (int j = 0; j < grid_y; j++){
						for (int i = 0; i < grid_x; i++){
							if (i%3 == 0 && j%3 == 0){
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 1;
							}
							else {
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 0;
							}
						}
					}
				}
				break;
			case 4:
				for (int k = grid_z34; k < grid_z4; k++){
					for (int j = 0; j < grid_y; j++){
						for (int i = 0; i < grid_x; i++){
							if (i%3 == 0 && j%3 == 0){
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 1;
							}
							else {
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 0;
							}
						}
					}
				}
				break;
			case 5:
				for (int k = grid_z45; k < grid_z5; k++){
					for (int j = 0; j < grid_y; j++){
						for (int i = 0; i < grid_x; i++){
							if (i%3 == 0 && j%3 == 0){
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 1;
							}
							else {
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 0;
							}
						}
					}
				}
				break;
			case 6:
				for (int k = grid_z56; k < grid_z6; k++){
					for (int j = 0; j < grid_y; j++){
						for (int i = 0; i < grid_x; i++){
							if (i%3 == 0 && j%3 == 0){
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 1;
							}
							else {
								indicator3D[6][i+j*grid_x+k*grid_x*grid_y] = 0;
							}
						}
					}
				}
				break;
			}
		}
		

		// Write to file "indicator.inp"
		FILE *p1;
		p1 = fopen("indicator.inp", "w");
		for (int i = 0; i < grid_x*grid_y*grid_z; i++){
			fprintf(p1, "%5d %5d %5d %10.2lf %10.4lf %10.4lf %5d\n", (int)indicator3D[0][i], (int)indicator3D[1][i], (int)indicator3D[2][i], indicator3D[3][i], indicator3D[4][i], indicator3D[5][i], (int)indicator3D[6][i]);
		}
		fclose(p1);
	}
	#endif

	return true;
 
 }