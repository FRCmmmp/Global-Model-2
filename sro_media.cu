//------------------------------------------------------------------------------
//
//         Function to generate random grain structure with hexagonal
//         sub ordering (short-range ordered (sro) media)
//
//         (c) P W Huang, Seagate Technology (2016). All rights reserved.
//
//------------------------------------------------------------------------------

#include "internal.h"
#include "random.h"
#include "Parameters_input.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#define M_PI  3.14159265358979323846
#include <cmath>

void sro_media(vec2d_t media_size, int num_large_grains, double grain_size, double grain_scaling, double domain_bnd_qfactor, std::vector<vec2d_t>& seed_array, std::vector <std::vector < vec2d_t> > &  vertex_array){

  const double media_area = media_size.x*media_size.y;
  const double grain_area = media_area/num_large_grains;
  const double md2 = grain_area;
  const double large_grain_size = sqrt(md2);


  // initialze random seed
  mtrandom::grnd.seed(iseed);  // PRN Generator (Linux)

  // calculate large grain voronoi but for 3x3 area
  // reserve() may not be necessary here
  seed_array.reserve(9*num_large_grains);
  vec2d_t start;
  start.x = mtrandom::grnd()*(3*media_size.x);
  start.y = mtrandom::grnd()*(3*media_size.y);
  seed_array.push_back(start);

  int max_iterations = 100*9*num_large_grains;
  int counter=0;

  
  // generate random points
  while(seed_array.size()<num_large_grains*9 && counter<max_iterations){
    vec2d_t tmp;
    tmp.x = mtrandom::grnd()*3*media_size.x;
    tmp.y = mtrandom::grnd()*3*media_size.y;

	// loop over all existing points and check grain is not too close
    bool overlapping=false;
    for(int j=0; j<seed_array.size(); ++j){
      double rij2 = (tmp.x-seed_array[j].x)*(tmp.x-seed_array[j].x) + (tmp.y-seed_array[j].y)*(tmp.y-seed_array[j].y);
	  if(rij2<0.5*md2){ // 0.5 need to be set as input
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

  
  
  if(counter==max_iterations) std::cerr << "Warning: Unable to find " << num_large_grains << " non-touching grains in sample. " << "\n"<<seed_array.size() << " grains generated" << std::endl;

  // translate voronoi seeds to centre on zero
  for(int g=0; g<seed_array.size();++g){
	  seed_array[g].x -= 1.5*media_size.x;
      seed_array[g].y -= 1.5*media_size.y;
  }

  // Generate list of vertices for large grains
  populate_vertex_points(seed_array, vertex_array);

  // remove edge vertices
  std::vector <std::vector < vec2d_t> > tmp_vertex_array(0);
  std::vector<vec2d_t> tmp_seed_array(0);
  
  for(int g=0; g<seed_array.size();++g){
	  bool out_of_bounds = false;
      for(int v=0; v<vertex_array[g].size(); ++v){
        if((vertex_array[g][v].x < -1.5*media_size.x || vertex_array[g][v].x > 1.5*media_size.x) ||
           (vertex_array[g][v].y < -1.5*media_size.y || vertex_array[g][v].y > 1.5*media_size.y) ) out_of_bounds = true;
      }
      // if grain is well bounded, keep it
      if(out_of_bounds == false){
        tmp_seed_array.push_back(seed_array[g]);
        tmp_vertex_array.push_back(vertex_array[g]);
      }
  }

    seed_array = tmp_seed_array;
    vertex_array = tmp_vertex_array;



  // Save vertices to file in gnuplot format
  /*std::ofstream vfile;
  vfile.open("vertices.out");
  for(int g=0; g<seed_array.size();++g){
    for(int v=0; v<vertex_array[g].size(); ++v){
      vfile << vertex_array[g][v].x << "\t" << vertex_array[g][v].y << std::endl;
    }
    vfile << vertex_array[g][0].x << "\t" << vertex_array[g][0].y << std::endl << std::endl;
  }*/

  //watch
  /*std::cout << seed_array.size() << "\t"<< vertex_array.size()<<std::endl;
  std::ofstream vfile1;
  vfile1.open("out3_new.out");
  for(int g=0; g<vertex_array.size();++g){
	for(int v=0; v<vertex_array[g].size(); ++v){
		vfile1 << vertex_array[g][v].x << "\t" << vertex_array[g][v].y << std::endl;
	}
  }
  std::ofstream vfile2;
  vfile2.open("out2_new.out");
  for(int g=0; g<vertex_array.size();++g){
	vfile2 << vertex_array[g].size() << std::endl;
  }*/




    // Declare list of sub-grain seed points
    std::vector<vec2d_t> sub_seed_array;

    //create media vertices
    std::vector<vec2d_t> media_vertex_array(4);
    media_vertex_array[0].x = -media_size.x*0.5;
    media_vertex_array[0].y = -media_size.y*0.5;
    media_vertex_array[1].x = -media_size.x*0.5;
    media_vertex_array[1].y = +media_size.y*0.5;
    media_vertex_array[2].x = +media_size.x*0.5;
    media_vertex_array[2].y = +media_size.y*0.5;
    media_vertex_array[3].x = +media_size.x*0.5;
    media_vertex_array[3].y = -media_size.y*0.5;

    // Loop over all large grains
    for(int g=0; g<seed_array.size(); ++g){

	  // generate a random angle and factors for rotational matrix
	  const double theta = 2.0*M_PI*mtrandom::grnd();
	  const double cost = cos(theta);
	  const double sint = sin(theta);

	  // calculate number of points in x and y directions
	  const int npx = large_grain_size*3;
	  const int npy = int(large_grain_size*3.0*sqrt(3.0)/3.0);

	  // calculate threshold for grain-grain separation
	  const double gs2 = grain_size*grain_size*domain_bnd_qfactor; // 0.999 need to be set as input

	  // Loop over all points in hexagonal lattice and check which ones are within grain
	  for (int i=0; i<npx; ++i){
		  for (int j=0; j<npy; ++j){
			  for (int parity=0;parity<2;++parity){
				  const double x = double(i)*grain_size + parity*grain_size*0.5 - npx*grain_size*0.5;
				  const double y = double(j)*grain_size*2.0*sqrt(3.0)/2.0+parity*grain_size*sqrt(3.0)/2.0 - npy*grain_size*0.5;
				  vec2d_t tmp;
				  // determine points rotated about origin and translated to grain centre
				  tmp.x = (x*cost - y*sint) + seed_array[g].x;
				  tmp.y = (x*sint + y*cost) + seed_array[g].y;
				  // check that translated point is in polygon and in media area
				  if(point_in_polygon(tmp.x, tmp.y, vertex_array[g]) && point_in_polygon(tmp.x, tmp.y, media_vertex_array)){
					  // check for overlapping grains
					  bool overlapping=false;
					  for(int j=0; j<sub_seed_array.size(); ++j){
						  double rij2 = (tmp.x-sub_seed_array[j].x)*(tmp.x-sub_seed_array[j].x) + (tmp.y-sub_seed_array[j].y)*(tmp.y-sub_seed_array[j].y);
						  if(rij2<gs2){
							  overlapping = true;
							  // break out of loop
							  break;
						  }
					  }
					  if(overlapping==false || sub_seed_array.size()==0) sub_seed_array.push_back(tmp);
				  }
			  }
		  }
	  }
	}

    // save lattice points to file
    std::ofstream pfile;
    pfile.open("grains.out");
    for(int i=0; i<sub_seed_array.size(); i++){
		pfile << sub_seed_array[i].x << "\t" << sub_seed_array[i].y << std::endl;
	}

    // save seed points to main seed file
    seed_array = sub_seed_array;

    // Regenerate list of vertices for small grains
    //populate_vertex_points(seed_array, vertex_array);

    // Generate circular vertices for small grains
    std::vector <vec2d_t> circle;
	for(int i=0; i<20; ++i){
		vec2d_t tmp;
		tmp.x = grain_size*(grain_scaling/2)*cos(double(i)*2.0*M_PI/20.0); //input
		tmp.y = grain_size*(grain_scaling/2)*sin(double(i)*2.0*M_PI/20.0); //input
		circle.push_back(tmp);
	}

    // erase data and resize grain vertices array to num_grains
	for(int i=0; i<vertex_array.size(); ++i) vertex_array[i].resize(0);
	vertex_array.resize(0);
	vertex_array.resize(sub_seed_array.size());
	for(int i=0; i<sub_seed_array.size(); i++){
		for(int v=0; v<circle.size(); ++v){
			vec2d_t tmp;
			// translate vertices
			tmp.x = circle[v].x+sub_seed_array[i].x;
			tmp.y = circle[v].y+sub_seed_array[i].y;
			vertex_array[i].push_back(tmp);
		}
	}

	// translate voronoi vertices to beginning point on zero
	for (int g=0; g<vertex_array.size(); ++g){
		for (int v=0; v<vertex_array[g].size(); ++v){
			vertex_array[g][v].x += 0.5*media_size.x;
			vertex_array[g][v].y += 0.5*media_size.y;
		}
	}


	// save grain vertices to file
    std::ofstream vfile1;
    vfile1.open("out3.out");
	for(int g=0; g<vertex_array.size();++g){
		for(int v=0; v<vertex_array[g].size(); ++v){
			vfile1 << vertex_array[g][v].x << "\t" << vertex_array[g][v].y << std::endl;
		}
	}
    std::ofstream vfile2;
    vfile2.open("out2.out");
    for(int g=0; g<vertex_array.size();++g){
		vfile2 << vertex_array[g].size() << std::endl;
	}
	
  return;

}
