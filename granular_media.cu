//------------------------------------------------------------------------------
//
//         Function to generate random grain structure
//
//         (c) P W Huang, Seagate Technology (2016). All rights reserved.
//
//------------------------------------------------------------------------------

#include "internal.h"
#include "random.h"
#include "clipper.h"
#include "Parameters_input.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#define M_PI  3.14159265358979323846
#include <cmath>

enum grain_bnd_t {same_dist=0, same_scaling=1};

void granular_media(vec2d_t media_size, int num_grains, double min_distance, double shrink_dist, double shrink_factor, std::vector<vec2d_t>& seed_array, std::vector <std::vector < vec2d_t> > &  vertex_array){

  const double md2 = min_distance*min_distance;
  // reserve() may not be necessary here
  seed_array.reserve(num_grains);
  
  
  // initialze random seed
  mtrandom::grnd.seed(iseed);  // PRN Generator (Linux)
  
  vec2d_t start;
  start.x = mtrandom::grnd()*media_size.x;
  start.y = mtrandom::grnd()*media_size.y;
  seed_array.push_back(start);

  int max_iterations = 10*num_grains;
  int counter=0;

  // generate random points
  while(seed_array.size()<num_grains && counter<max_iterations){
    vec2d_t tmp;
    tmp.x = mtrandom::grnd()*media_size.x;
    tmp.y = mtrandom::grnd()*media_size.y;

    // loop over all existing points and check grain is not too close
    bool overlapping=false;
    for(int j=0; j<seed_array.size(); ++j){
      double rij2 = (tmp.x-seed_array[j].x)*(tmp.x-seed_array[j].x) + (tmp.y-seed_array[j].y)*(tmp.y-seed_array[j].y);
      if(rij2<md2){
        overlapping=true;
        // break out of loop
        break;
      }
    }
    // Only accept non-overlapping grains
    if(overlapping==false) seed_array.push_back(tmp);

    // increment attempt counter
    counter++;

  }

  if(counter==max_iterations) std::cerr << "Warning: Unable to find " << num_grains << " non-touching grains in sample. " << seed_array.size() << "grains generated" << std::endl;

  // translate voronoi seeds to centre on zero
  for(int g=0; g<seed_array.size();++g){
    seed_array[g].x -= 0.5*media_size.x;
    seed_array[g].y -= 0.5*media_size.y;
  }

  // Generate list of vertices for grains
  populate_vertex_points(seed_array, vertex_array);

  // remove edge vertices
  std::vector <std::vector < vec2d_t> > tmp_vertex_array(0);
  std::vector<vec2d_t> tmp_seed_array(0);

  for(int g=0; g<seed_array.size();++g){
    bool out_of_bounds = false;
    for(int v=0; v<vertex_array[g].size(); ++v){
      if((vertex_array[g][v].x < -0.5*media_size.x || vertex_array[g][v].x > 0.5*media_size.x) ||
         (vertex_array[g][v].y < -0.5*media_size.y || vertex_array[g][v].y > 0.5*media_size.y) ) out_of_bounds = true;
    }
    // if grain is well bounded, keep it
    if(out_of_bounds == false){
      tmp_seed_array.push_back(seed_array[g]);
      tmp_vertex_array.push_back(vertex_array[g]);
    }
  }

  seed_array = tmp_seed_array;
  vertex_array = tmp_vertex_array;

  
  /*// grain rounding
  const double area_cutoff = 0.7; // Needs to be an input variable
  double grain_size = 10.0;       // Needs to be an input variable
  for(int g=0; g<seed_array.size();++g){

    // get number of vertices for each grain
    const int nv = vertex_array[g].size();

      // allocate temporary array for area calculation
      std::vector <vec2d_t> rnd;
      rnd.resize(48);

      double area_frac=area_cutoff;
      double deltar=0.02*grain_size; // nm
      double radius=0.0; // Starting radius
      double area=0.0; // Starting area

      std::vector < vec2d_t> tmp_vertex_array = vertex_array[g];

      //----------------------------------
      // precalculate voronoi area
      //----------------------------------
      double varea=0.0;
      for(int r=0;r<100;r++){
        radius+=deltar;

        for(int i=0;i<48;i++){
          double theta = 2.0*M_PI*double(i)/48.0;
          double x = radius*cos(theta);
          double y = radius*sin(theta);

          // Check to see if site is within polygon
          if(point_in_polygon(x+seed_array[g].x,y+seed_array[g].y,tmp_vertex_array)==true){
            rnd[i].x=x;
            rnd[i].y=y;
          }
        }

        //update area
        varea=0.0;
        for(int i=0;i<48;i++){
          int nvi = i+1;
          if(nvi>=48) nvi=0;
          varea+=0.5*sqrt((rnd[nvi].x-rnd[i].x)*(rnd[nvi].x-rnd[i].x))*sqrt((rnd[nvi].y-rnd[i].y)*(rnd[nvi].y-rnd[i].y));
        }
      }

      //-----------------------------------------
      // calculate rounding
      //-----------------------------------------
      // reset polygon positions and radius
      for(int idx=0; idx<48;idx++){
        rnd[idx].x=0.0;
        rnd[idx].y=0.0;
      }
      radius=0.0;
      // expand polygon
      for(int r=0;r<100;r++){
      if(area<area_frac*varea){
        radius+=deltar;

        //loop over coordinates
        for(int i=0;i<48;i++){
          double theta = 2.0*M_PI*double(i)/48.0;
          double x = radius*cos(theta);
          double y = radius*sin(theta);

          // Check to see if site is within polygon
          if(point_in_polygon(x+seed_array[g].x,y+seed_array[g].y,tmp_vertex_array)==true){
            rnd[i].x=x;
            rnd[i].y=y;
          }
        }

        //update area
        area=0.0;
        for(int i=0;i<48;i++){
          int nvi = i+1;
          if(nvi>=48) nvi=0;
          area+=0.5*sqrt((rnd[nvi].x-rnd[i].x)*(rnd[nvi].x-rnd[i].x))*sqrt((rnd[nvi].y-rnd[i].y)*(rnd[nvi].y-rnd[i].y));
        }
      }
    }

    // set new polygon points
    vertex_array[g]=rnd;
    // translate new vertices
    for(int v=0; v<vertex_array[g].size(); ++v){
      vertex_array[g][v].x += seed_array[g].x;
      vertex_array[g][v].y += seed_array[g].y;
    }
  } // end of grain loop*/

  
  
  //-------------------------------------
  // shrink grains for grain boundaries
  //-------------------------------------
  // define grain boundary type
  grain_bnd_t grain_bnd_type = same_dist;
  //grain_bnd_t grain_bnd_type = same_scaling;
  std::vector <std::vector <vec2d_t> > shrinked_vertex_array;
  switch(grain_bnd_type){
  // grain boundaries are of the same width
  case same_dist:
	{

		// find middle point between adjacent vertices for uniformly shrinking distance
		#ifdef __CALC_MIDDLE_POINT__
		std::vector <std::vector < vec2d_t> > tmp_mid_vertex_array = vertex_array;
		for(int g=0; g<seed_array.size();++g){
			for(int v=0; v<vertex_array[g].size(); ++v){
				if (v+1 == vertex_array[g].size()){
					tmp_mid_vertex_array[g][v].x = (vertex_array[g][v].x + vertex_array[g][0].x)/2;
					tmp_mid_vertex_array[g][v].y = (vertex_array[g][v].y + vertex_array[g][0].y)/2;
				}
				else {
					tmp_mid_vertex_array[g][v].x = (vertex_array[g][v].x + vertex_array[g][v+1].x)/2;
					tmp_mid_vertex_array[g][v].y = (vertex_array[g][v].y + vertex_array[g][v+1].y)/2;
				}
			}
		}
		#endif

		// employ Clipper Offset  (from Clipper Library) to offset polygons;
		// the ClipperLib of this version only handles points of integer coordinates;
		// need to do upscale (when converting double to int) and downscale (when converting int to double) steps;
		// random boundary widths can be implemented if adding distributions to the scaling factors.
		ClipperLib::Path subj;
		ClipperLib::Paths solution;
		ClipperLib::ClipperOffset co;
		std::vector <std::vector < vec2d_t> > tmp_shrinked_vertex_array;
		for (int g=0; g<seed_array.size();++g){
			ClipperLib::Path subj;
			ClipperLib::Paths solution;
			ClipperLib::ClipperOffset co;
			for(int v=0; v<vertex_array[g].size(); ++v){

				// 1000000 as the up- or down-scaling factor should be large enough
				long long int tmp_int_x = (int)((vertex_array[g][v].x)*1000000),    
							  tmp_int_y = (int)((vertex_array[g][v].y)*1000000);
				subj << ClipperLib::IntPoint( tmp_int_x, tmp_int_y);
			}
			co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
			co.Execute(solution, shrink_dist*1000000);
			for (int g=0; g<solution.size(); ++g){
				std::vector < vec2d_t> tmp_shrinked_single_grain_vertex_array;
				for (int v=0; v<solution[g].size(); ++v){
					vec2d_t tmp_double;
					tmp_double.x = (double)(solution[g][v].X)/1000000;
					tmp_double.y = (double)(solution[g][v].Y)/1000000;
					tmp_shrinked_single_grain_vertex_array.push_back(tmp_double);
				}
				tmp_shrinked_vertex_array.push_back(tmp_shrinked_single_grain_vertex_array);
			}
		}
		shrinked_vertex_array = tmp_shrinked_vertex_array;
		break;
	}

	// grain boundaries are of different widths;
	// random boundary widths can be implemented if adding distributions to the scaling factors.
	case same_scaling:
	{
		for(int g=0; g<seed_array.size();++g){
			for(int v=0; v<vertex_array[g].size(); ++v){

				vec2d_t tmp;
				// calculate reduced vertex coordinate
				tmp.x = (vertex_array[g][v].x - seed_array[g].x)*0.85+seed_array[g].x;
				tmp.y = (vertex_array[g][v].y - seed_array[g].y)*shrink_factor+seed_array[g].y;
				vertex_array[g][v] = tmp;
			}
		}
		break;
	}

	default:
	{
		// find middle point between adjacent vertices for uniformly shrinking distance
		#ifdef __CALC_MIDDLE_POINT__
		std::vector <std::vector < vec2d_t> > tmp_mid_vertex_array = vertex_array;
		for(int g=0; g<seed_array.size();++g){
			for(int v=0; v<vertex_array[g].size(); ++v){
				if (v+1 == vertex_array[g].size()){
					tmp_mid_vertex_array[g][v].x = (vertex_array[g][v].x + vertex_array[g][0].x)/2;
					tmp_mid_vertex_array[g][v].y = (vertex_array[g][v].y + vertex_array[g][0].y)/2;
				}
				else {
					tmp_mid_vertex_array[g][v].x = (vertex_array[g][v].x + vertex_array[g][v+1].x)/2;
					tmp_mid_vertex_array[g][v].y = (vertex_array[g][v].y + vertex_array[g][v+1].y)/2;
				}
			}
		}
		#endif

		// employ Clipper Offset  (from Clipper Library) to offset polygons;
		// the ClipperLib of this version only handles points of integer coordinates;
		// need to do upscale (when converting double to int) and downscale (when converting int to double) steps;
		// random boundary widths can be implemented if adding distributions to the scaling factors.
		ClipperLib::Path subj;
		ClipperLib::Paths solution;
		ClipperLib::ClipperOffset co;
		std::vector <std::vector < vec2d_t> > tmp_shrinked_vertex_array;
		for (int g=0; g<seed_array.size();++g){
			ClipperLib::Path subj;
			ClipperLib::Paths solution;
			ClipperLib::ClipperOffset co;
			for(int v=0; v<vertex_array[g].size(); ++v){

				// 1000000 as the up- or down-scaling factor should be large enough
				long long int tmp_int_x = (int)((vertex_array[g][v].x)*1000000), 
							  tmp_int_y = (int)((vertex_array[g][v].y)*1000000);
				subj << ClipperLib::IntPoint( tmp_int_x, tmp_int_y);

			}
			co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
			co.Execute(solution, shrink_dist*1000000);
			for (int g=0; g<solution.size(); ++g){
				std::vector < vec2d_t> tmp_shrinked_single_grain_vertex_array;
				for (int v=0; v<solution[g].size(); ++v){

					vec2d_t tmp_double;
					tmp_double.x = (double)(solution[g][v].X)/1000000;
					tmp_double.y = (double)(solution[g][v].Y)/1000000;
					tmp_shrinked_single_grain_vertex_array.push_back(tmp_double);

				}
				tmp_shrinked_vertex_array.push_back(tmp_shrinked_single_grain_vertex_array);
			}
		}
		shrinked_vertex_array = tmp_shrinked_vertex_array;
		break;
	}
  }

  // translate voronoi vertices to beginning point on zero
  for (int g=0; g<shrinked_vertex_array.size(); ++g){
	  for (int v=0; v<shrinked_vertex_array[g].size(); ++v){
		  shrinked_vertex_array[g][v].x += 0.5*media_size.x;
		  shrinked_vertex_array[g][v].y += 0.5*media_size.y;
	  }
  }

  // remove those seeds with zero grain area
  std::vector<vec2d_t>::iterator it;
  std::vector <std::vector <vec2d_t> >::iterator iit;
  std::vector <vec2d_t> tmp_seed_array1 = seed_array;
  std::vector <std::vector <vec2d_t> > tmp_vertex_array1 = shrinked_vertex_array;
  int g = 0;
  while(g<tmp_vertex_array1.size()){	
	  //std::cout << tmp_vertex_array1[g].size() <<std::endl;
	  if (tmp_vertex_array1[g].size()==0){
		iit = tmp_vertex_array1.begin() + g;
		tmp_vertex_array1.erase(iit);
		//std::cout << tmp_vertex_array1[g].size() <<std::endl;
	  }
		else g++;
  }

  // save grain vertices to file 
  std::ofstream vfile1;
  vfile1.open("out3.out");
  for(int g=0; g<tmp_vertex_array1.size();++g){
	  for(int v=0; v<tmp_vertex_array1[g].size(); ++v){
		  vfile1 << tmp_vertex_array1[g][v].x << "\t" << tmp_vertex_array1[g][v].y << std::endl;
	  }
  }
  std::ofstream vfile2;
  vfile2.open("out2.out");
  for(int g=0; g<tmp_vertex_array1.size();++g){
	  vfile2 << tmp_vertex_array1[g].size() << std::endl;
  }




  //-------------------------------------
  // Estimate grains' diamteres for mean
  // and standard deviation.
  //-------------------------------------

  // calculate area for each grain
  std::vector<double> grain_area_array;
  polygon_area(shrinked_vertex_array, grain_area_array);
 
  std::vector<double> grain_dia_array;
  for (int g=0; g<shrinked_vertex_array.size(); ++g){
	  grain_dia_array.push_back( sqrt(grain_area_array[g]*4/M_PI) );
  }
  
  // calculate mean diameter
  double mean_dia = 0.0f;
  for (int g=0; g<grain_dia_array.size(); ++g){	mean_dia = mean_dia + grain_dia_array[g]/grain_dia_array.size(); }

  // calcualte std diameter
  double std_dia = 0.0f;
  for (int g=0; g<grain_dia_array.size(); ++g){	
	  std_dia = std_dia + (grain_dia_array[g]-mean_dia)*(grain_dia_array[g]-mean_dia)/grain_dia_array.size(); 
  }
  std_dia = sqrt(std_dia);
  
  
  // Save Voronoi grain properties to file
  std::ofstream vfile3;
  vfile3.open("media_voronoi_info.out");
  vfile3<<std::setw(15)<<"diameter std= "<<std::setw(10)<<std_dia<<" [nm]"<<std::endl
	    <<std::setw(15)<<"diameter mean="<<std::setw(10)<<mean_dia<<" [nm]"<<std::endl
		<<std::setw(10)<<"std/mean="<<std::setw(15)<<std_dia/mean_dia*100<<" [%]"<<std::endl
		<<std::setw(15)<<"# grains=     "<<std::setw(10)<<grain_dia_array.size()<<""<<std::endl;
  std::cout <<std::setw(15)<<"diameter std= "<<std::setw(10)<<std_dia<<" [nm]"<<std::endl
	        <<std::setw(15)<<"diameter mean="<<std::setw(10)<<mean_dia<<" [nm]"<<std::endl
		    <<std::setw(10)<<"std/mean="<<std::setw(15)<<std_dia/mean_dia*100<<" [%]"<<std::endl;

  // Save vertices to file in gnuplot format
  /*std::ofstream vfile;
  vfile.open("vertices.out");
  for(int g=0; g<seed_array.size();++g){
    for(int v=0; v<vertex_array[g].size(); ++v){
      vfile << vertex_array[g][v].x << "\t" << vertex_array[g][v].y << std::endl;
    }
    vfile << vertex_array[g][0].x << "\t" << vertex_array[g][0].y << std::endl << std::endl;
  }*/

  return;

}