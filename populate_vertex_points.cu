//-----------------------------------------------------------------------------
//
// (C) P Huang Seagate Technology 2016
//
// ----------------------------------------------------------------------------
//
#include "qvoronoi.h"
#include "internal.h"
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>

void populate_vertex_points(std::vector <vec2d_t> & grain_coord_array, std::vector <std::vector < vec2d_t> > &  grain_vertices_array){
	//========================================================================================================
	//		 				Function to populate voronoi vertices for grains using qhull
	//
	//														Version 1.0
	//
	//												R F Evans 15/07/2009
	//
	//========================================================================================================
	//		Locally allocated variables: vertex_array
	//========================================================================================================

	const int num_grains=grain_coord_array.size();
	std::stringstream grain_file_sstr;
	std::stringstream voronoi_file_sstr;

	grain_file_sstr << "grains.out";
	voronoi_file_sstr << "voronoi.out";

	std::string grain_file = grain_file_sstr.str();
	std::string voronoi_file = voronoi_file_sstr.str();
	const char* grain_filec = grain_file.c_str();
	const char* voronoi_filec = voronoi_file.c_str();

    // erase data and resize grain vertices array to num_grains
	for(int i=0; i<grain_vertices_array.size(); ++i) grain_vertices_array[i].resize(0);
	grain_vertices_array.resize(0);
    grain_vertices_array.resize(num_grains);

	//--------------------------------------
	// Scale grain coordinates
	//--------------------------------------

	double scale_factor=0.0;

	//--------------------------------------------------------------
	// Calculate max grain coordinate
	//--------------------------------------------------------------
	for(int i=0;i<num_grains;i++){
		if(grain_coord_array[i].x>scale_factor){
			scale_factor=grain_coord_array[i].x;
		}
		if(grain_coord_array[i].y>scale_factor){
			scale_factor=grain_coord_array[i].y;
		}
	}

	//--------------------------------------------------------------
	// Output grain coordindates and Scale to be unit length (0:1)
	//--------------------------------------------------------------

	std::ofstream grain_coord_file;
  	grain_coord_file.open (grain_filec);

	grain_coord_file << "#============================================" << std::endl;
	grain_coord_file << "# Grain Coordinate file for input into qhull" << std::endl;
	grain_coord_file << "#============================================" << std::endl;
	grain_coord_file << "# Grain Scaling Factor" << std::endl;
	grain_coord_file << "#\t" << scale_factor << std::endl;
	grain_coord_file << 2 << std::endl;
	grain_coord_file << "# Number of Grains" << std::endl;
	grain_coord_file << num_grains << std::endl;
	grain_coord_file << "# Grains Coordinates" << std::endl;

	for(int i=0;i<num_grains;i++){
		grain_coord_file << grain_coord_array[i].x/scale_factor-0.5 << "\t";
 		grain_coord_file << grain_coord_array[i].y/scale_factor-0.5 << std::endl;
	}

	grain_coord_file.close();

   //--------------------------------------------------------
   // Use voronoi library creating an import and export temporary files
   //--------------------------------------------------------
   FILE *inputqv, *outputqv;
   int qargc=3;
   const char *qargv[3]={"qvoronoi", "-o", "-Fv"};
   inputqv=fopen(grain_file.c_str(),"r");
   outputqv=fopen(voronoi_file.c_str(),"w");
   qvoronoi(qargc, const_cast<char**>(qargv), inputqv, outputqv);
   fclose(outputqv);
   fclose(inputqv);

   //--------------------------------------------------------
   // Read in number of Voronoi vertices
   //--------------------------------------------------------

	int dimensions,num_vertices,num_points,one;

	std::ifstream vertices_file;
	vertices_file.open (voronoi_filec);
	vertices_file >> dimensions;
	vertices_file >> num_vertices >> num_points >> one;


	//----------------------------------------------------------
	// Allocate vertex_array
	//----------------------------------------------------------
	std::vector<vec2d_t> vertex_array(num_vertices);

	//--------------------------------------
	// Read in Voronoi vertices and rescale
	//--------------------------------------

	for(int i=0;i<num_vertices;i++){
		vertices_file >> vertex_array[i].x;
		vertices_file >> vertex_array[i].y;
	}

	for(int i=0;i<num_vertices;i++){
		vertex_array[i].x=(vertex_array[i].x+0.5)*scale_factor;
		vertex_array[i].y=(vertex_array[i].y+0.5)*scale_factor;
	}

   //--------------------------------------
   // Read in Voronoi vertex associations
   //--------------------------------------
	 vec2d_t tmp;
	 for(int i=0;i<num_grains;i++){
      int num_assoc_vertices;	// Number of vertices associated with point i
      int vertex_number;		// temporary vertex number
      vertices_file >> num_assoc_vertices;
      bool inf=false;
      for(int j=0;j<num_assoc_vertices;j++){
         vertices_file >> vertex_number;
				 tmp.x = vertex_array[vertex_number].x;
				 tmp.y = vertex_array[vertex_number].y;
         grain_vertices_array[i].push_back(tmp);
         // check for unbounded grains
         if(vertex_number==0) inf=true;
         // check for bounded grains with vertices outside bounding box
         //if((vertex_array[vertex_number].x<0.0) || (vertex_array[vertex_number].x>cs::system_dimensions.x)) inf=true;
         //if((vertex_array[vertex_number].y<0.0) || (vertex_array[vertex_number].y>cs::system_dimensions.y)) inf=true;
      }

      //-------------------------------------------------------------------
      // Set unbounded grains to zero vertices for later removal
      //-------------------------------------------------------------------
      if(inf==true){
         grain_vertices_array[i].resize(0);
      }
   }
   vertices_file.close();

   //-------------------------------------------------------------------
   // Remove Temporary voronoi files
   //-------------------------------------------------------------------
   //{
   //   std::stringstream rmfsstr;
   //   #ifdef WIN_COMPILE
   //      //rmfsstr << "del /f " << grain_file;
   //   #else
   //      rmfsstr << "rm -f " << grain_file;
   //   #endif
   //   std::string rmfstr = rmfsstr.str();
   //   const char* rmfcstr = rmfstr.c_str();
   //   int sysstat=system(rmfcstr);
   //   if(sysstat!=0) {
		 // std::cerr << "Error removing temporary grain files" << std::endl;
	  //}
   //}
   //{
   //   std::stringstream rmfsstr;
   //   #ifdef WIN_COMPILE
   //      rmfsstr << "del /f " << voronoi_file;
   //   #else
   //      rmfsstr << "rm -f " << voronoi_file;
   //   #endif
   //   std::string rmfstr = rmfsstr.str();
   //   const char* rmfcstr = rmfstr.c_str();
   //   int sysstat=system(rmfcstr);
   //   if(sysstat!=0) {
		 // std::cerr << "Error removing temporary voronoi files" << std::endl;
	  //}
   //}

   return;

}
