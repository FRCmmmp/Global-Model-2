#include "internal.h"
#include <cstdlib>
#include <iostream>


bool polygon_area(std::vector <std::vector < vec2d_t> > &  grain_vertices_array, std::vector<double> &   grain_area_array){
	
	//std::cout<<grain_vertices_array.size()<<std::endl;
	for (int g=0; g<grain_vertices_array.size(); ++g){
		double area = 0;
		int tmp_v = grain_vertices_array[g].size()-1;
		for (int v=0; v<grain_vertices_array[g].size(); ++v){
			area += (grain_vertices_array[g][tmp_v].x + grain_vertices_array[g][v].x)*
				    (grain_vertices_array[g][tmp_v].y - grain_vertices_array[g][v].y);
			tmp_v = v;
		}
		area = abs(area)*0.5;
		grain_area_array.push_back(area);
	}
	std::cout<<"# grains= "<<grain_area_array.size()<<std::endl;

	return true; 
}