#pragma once
#ifndef SLAM_H
#define SLAM_H

#include "init_visual.h"
#include "Euclid_Bound.h"

struct bound_entry
{
	int id;
	double bound_value;
	bool is_LB;

	bool operator < (const bound_entry& another) const {
		return bound_value < another.bound_value;
	}
};

void envelope_point_set(statistics& stat, vector<int>& E_k);
void bound_list(statistics& stat, vector<int>& E_k, vector<bound_entry>& List);
void SLAM_visual(statistics& stat);

#endif