#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include <vector>

class neighbor
{
	// * creating instance
	// base_grid d_base_grid;

public:
	// // TODO: generate neighborhood using link-list algorithm
	// std::vector<std::vector<int>> link_list(const int np, const std::vector<double> &sp,
	// 										const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &zp, int neighbor_scale);
	// ! -- obsolete --
	//// void link_list(const int np, const std::vector<double> &sp, const std::vector<double> &xp, const std::vector<double> &yp,
	////     std::vector<int> &pair_i, std::vector<int> &pair_j, int neighbor_scale);
	// TODO: generate neighborhood using direct searching
	std::vector<std::vector<int>> direct_find(const int np, const std::vector<double> &sp,
											  const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &zp, const int neighbor_scale);
	// ! -- obsolete --
	//// void direct_find(const int np, const std::vector<double> &sp, const std::vector<double> &xp, const std::vector<double> &yp,
	//// 				 std::vector<int> &pair_i, std::vector<int> &pair_j, const int neighbor_scale);
};

#endif