//------------------------------------------------------------------------------
//
//         Header file for useful functions for general purposes
//
//         (c) P W Huang, Seagate Technology (2016). All rights reserved.
//
//------------------------------------------------------------------------------

#include <vector>
#include <numeric>
#include <algorithm>

int Rnd_upto_pow2(int v);

//template <typename T, typename Compare>
//std::vector<std::size_t> sort_permutation(const std::vector<T> & vec, Compare& compare){
//	std::vector<std::size_t> p(vec.size());
//	std::iota(p.begin(), p.end(), 0);
//	std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
//	return p;
//}

//template <typename T>
//std::vector<T> apply_permutation(const std::vector<T> & vec, const std::vector<std::size_t> & p){
//	
//	std::vector<T> sorted_vec(p.size());
//	std::transform(p.begin(), p.end(), sorted_vec.begin(), [&](std::size_t i) { return vec[i]; });
//	return sorted_vec;
//}
