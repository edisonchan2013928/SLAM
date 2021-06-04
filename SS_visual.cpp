#include "SS_visual.h"

double SCAN(double*q, statistics& stat)
{
	double temp_value;
	double sq_dist_value;
	double dist_value;
	double incr_value = 0;

	for (int i = 0; i < stat.n; i++)
	{
		sq_dist_value = sq_euclid_dist(q, stat.featureVector[i], stat.dim);

		if (sq_dist_value > stat.b_value*stat.b_value)
			continue;

		if (stat.kernel_type == 0) //Uniform kernel
			incr_value += 1.0 - (1.0 / stat.b_value);
		if (stat.kernel_type == 1) //Epanechnikov kernel
			incr_value += 1.0 - (1.0 / (stat.b_value*stat.b_value))*sq_dist_value;
		if (stat.kernel_type == 2) //Quartic kernel
		{
			temp_value = 1.0 - (1.0 / (stat.b_value*stat.b_value))*sq_dist_value;
			incr_value += temp_value * temp_value;
		}
	}

	return incr_value;
}

void KDE_visual(statistics& stat)
{
	double KDE_value;

	for (int r = 0; r < stat.n_row; r++)
	{
		for (int c = 0; c < stat.n_col; c++)
		{
			KDE_value = SCAN(stat.queryVector[r*stat.n_col + c], stat);
			stat.out_visual[r][c] = KDE_value;
		}
	}
}

double refinement(Node*curNode, statistics& stat)
{
	double f_cur;
	double temp_value;
	double sq_dist_value;
	double dist_value;
	int id;
	double*q = stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c];
	f_cur = 0;

	for (int i = 0; i < (int)curNode->idList.size(); i++)
	{
		id = curNode->idList[i];
		sq_dist_value = sq_euclid_dist(q, stat.featureVector[id], stat.dim);

		if (sq_dist_value > stat.b_value*stat.b_value)
			continue;

		if (stat.kernel_type == 0) //Uniform kernel
			f_cur += 1.0 / stat.b_value;

		if (stat.kernel_type == 1) //Epanechnikov kernel
			f_cur += 1.0 - (1.0 / (stat.b_value*stat.b_value))*sq_dist_value;
		if (stat.kernel_type == 2) //Quartic kernel
		{
			temp_value = 1.0 - (1.0 / (stat.b_value*stat.b_value))*sq_dist_value;
			f_cur += temp_value * temp_value;
		}
	}

	return f_cur;
}