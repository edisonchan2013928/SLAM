#include "SLAM.h"

void envelope_point_set(statistics& stat, vector<int>& E_k)
{
	for (int i = 0; i < stat.n; i++)
	{
		if (fabs(stat.featureVector[i][stat.static_coord] - stat.k) < stat.b_value)
			E_k.push_back(i);
	}
}

void bound_list(statistics& stat, vector<int>& E_k, vector<bound_entry>& List)
{
	double pm_value;
	int id;

	bound_entry b_entry_ell;
	bound_entry b_entry_u;

	for (int i = 0; i < (int)E_k.size(); i++)
	{
		id = E_k[i];
		pm_value = sqrt(stat.b_value*stat.b_value - (stat.k - stat.featureVector[id][stat.static_coord])*(stat.k - stat.featureVector[id][stat.static_coord]));
		
		b_entry_ell.id = id;
		b_entry_ell.bound_value = stat.featureVector[id][stat.dynamic_coord] - pm_value;
		b_entry_ell.is_LB = true;

		b_entry_u.id = id;
		b_entry_u.bound_value = stat.featureVector[id][stat.dynamic_coord] + pm_value;
		b_entry_u.is_LB = false;

		List.push_back(b_entry_ell);
		List.push_back(b_entry_u);
	}
}

void init_SLAM(statistics& stat)
{
	stat.A_L_ell = new double[stat.dim];
	stat.A_U_ell = new double[stat.dim];
	stat.S_L_ell = 0;
	stat.S_U_ell = 0;
	stat.R_q_card = 0;
	stat.A_R_q = new double[stat.dim];
	stat.S_R_q = 0;

	for (int q_id = 0; q_id < stat.dynamic_pixel_size; q_id++)
	{
		double*query = new double[stat.dim];
		stat.query_list.push_back(query);
		stat.result_list.push_back(0);
	}

	for (int d = 0; d < stat.dim; d++)
	{
		stat.A_L_ell[d] = 0;
		stat.A_U_ell[d] = 0;
		stat.A_R_q[d] = 0;
	}

	if (stat.method == 6 || stat.method == 8) //SLAM_BUCKET
	{
		vector<int> idList;
		for (int b = 0; b <= stat.dynamic_pixel_size; b++)
		{
			stat.B_L.push_back(idList);
			stat.B_U.push_back(idList);
		}
	}
}

void clear_SLAM(statistics& stat)
{
	for (int d = 0; d < stat.dim; d++)
	{
		stat.A_L_ell[d] = 0;
		stat.A_U_ell[d] = 0;
	}
	stat.S_L_ell = 0;
	stat.S_U_ell = 0;
	stat.L_ell.clear();
	stat.U_ell.clear();

	if (stat.method == 6 || stat.method == 8) //SLAM_BUCKET
	{
		for (int b = 0; b <= stat.dynamic_pixel_size; b++)
		{
			stat.B_L[b].clear();
			stat.B_U[b].clear();
		}
	}
}

void SLAM_SORT(statistics& stat)
{
	vector<int> E_k;
	vector<bound_entry> List;
	int i_q, i_b;
	double ip;
	double p_SquareNorm;
	bool p_finish = false;

	envelope_point_set(stat, E_k);
	bound_list(stat, E_k, List);

	sort(List.begin(), List.end());
	i_q = 0;
	i_b = 0;
	
	if (E_k.size() == 0) //degenerate case 1
	{
		for (int i_q = 0; i_q < stat.dynamic_pixel_size; i_q++)
			stat.result_list[i_q] = 0;
		return;
	}
	
	if (stat.dynamic_pixel_size == 0) //degenerate case 2
		return;

	//while(i_q < stat.dynamic_pixel_size)
	//while (i_q + i_b <= 2 * E_k.size() + stat.dynamic_pixel_size - 2)
	while (i_q < stat.dynamic_pixel_size)
	{
		if (p_finish == true || stat.query_list[i_q][stat.dynamic_coord] <= List[i_b].bound_value)
		{
			stat.q_SquareNorm = computeSqNorm(stat.query_list[i_q], stat.dim);
			stat.R_q_card = stat.L_ell.size() - stat.U_ell.size();
			for (int d = 0; d < stat.dim; d++)
				stat.A_R_q[d] = stat.A_L_ell[d] - stat.A_U_ell[d];
			stat.S_R_q = stat.S_L_ell - stat.S_U_ell;

			ip = inner_product(stat.query_list[i_q], stat.A_R_q, stat.dim);
			stat.result_list[i_q] = stat.R_q_card - (1.0 / (stat.b_value*stat.b_value))
				*(stat.R_q_card*stat.q_SquareNorm - 2 * ip + stat.S_R_q);

			i_q++;
			//i_q = min(i_q + 1, stat.dynamic_pixel_size - 1);
		}
		else
		{
			p_SquareNorm = 0;
			if (List[i_b].is_LB == true)
			{
				stat.L_ell.push_back(List[i_b].id);
				for (int d = 0; d < stat.dim; d++)
				{
					stat.A_L_ell[d] += stat.featureVector[List[i_b].id][d];
					p_SquareNorm += stat.featureVector[List[i_b].id][d] * stat.featureVector[List[i_b].id][d];
				}
					
				stat.S_L_ell += p_SquareNorm;
			}
			else
			{
				stat.U_ell.push_back(List[i_b].id);
				for (int d = 0; d < stat.dim; d++)
				{
					stat.A_U_ell[d] += stat.featureVector[List[i_b].id][d];
					p_SquareNorm += stat.featureVector[List[i_b].id][d] * stat.featureVector[List[i_b].id][d];
				}

				stat.S_U_ell += p_SquareNorm;
			}

			i_b++;
			if (i_b >= 2 * (int)E_k.size())
				p_finish = true;
			//i_b = min(i_b + 1, (int)E_k.size() - 1);
		}
	}
}

void SLAM_BUCKET(statistics& stat)
{
	vector<int> E_k;
	vector<bound_entry> List;
	double ip;
	double p_L_SquareNorm, p_U_SquareNorm;
	int id;
	int i_l, i_u;

	envelope_point_set(stat, E_k);
	bound_list(stat, E_k, List);

	for (int i = 0; i < (int)List.size(); i++)
	{
		id = List[i].id;
		if (List[i].is_LB == true)
		{
			i_l = (int)max(ceil((List[i].bound_value - stat.query_list[0][stat.dynamic_coord]) / stat.gap), 0.0);
			stat.B_L[i_l].push_back(id);
		}
		if (List[i].is_LB == false)
		{
			i_u = (int)min(ceil((List[i].bound_value - stat.query_list[0][stat.dynamic_coord]) / stat.gap), (double)stat.dynamic_pixel_size);
			stat.B_U[i_u].push_back(id);
		}
	}

	for (int i_q = 0; i_q < stat.dynamic_pixel_size; i_q++)
	{
		for (int e = 0; e < stat.B_L[i_q].size(); e++)
		{
			p_L_SquareNorm = 0;
			stat.L_ell.push_back(stat.B_L[i_q][e]);
			for (int d = 0; d < stat.dim; d++)
			{
				stat.A_L_ell[d] += stat.featureVector[stat.B_L[i_q][e]][d];
				p_L_SquareNorm += stat.featureVector[stat.B_L[i_q][e]][d]
					* stat.featureVector[stat.B_L[i_q][e]][d];
			}
			
			stat.S_L_ell += p_L_SquareNorm;
		}
		
		for (int e = 0; e < stat.B_U[i_q].size(); e++)
		{
			p_U_SquareNorm = 0;
			stat.U_ell.push_back(stat.B_U[i_q][e]);
			for (int d = 0; d < stat.dim; d++)
			{
				stat.A_U_ell[d] += stat.featureVector[stat.B_U[i_q][e]][d];
				p_U_SquareNorm += stat.featureVector[stat.B_U[i_q][e]][d]
					* stat.featureVector[stat.B_U[i_q][e]][d];
			}

			stat.S_U_ell += p_U_SquareNorm;
		}

		stat.q_SquareNorm = computeSqNorm(stat.query_list[i_q], stat.dim);
		stat.R_q_card = stat.L_ell.size() - stat.U_ell.size();
		for (int d = 0; d < stat.dim; d++)
			stat.A_R_q[d] = stat.A_L_ell[d] - stat.A_U_ell[d];
		stat.S_R_q = stat.S_L_ell - stat.S_U_ell;

		ip = inner_product(stat.query_list[i_q], stat.A_R_q, stat.dim);
		stat.result_list[i_q] = stat.R_q_card - (1.0 / (stat.b_value*stat.b_value))
			*(stat.R_q_card*stat.q_SquareNorm - 2 * ip + stat.S_R_q);
	}
}

void SLAM_scan_x(statistics& stat, bool is_bucket)
{
	for (int c = 0; c < stat.n_col; c++)
	{
		for (int r = 0; r < stat.n_row; r++)
		{
			stat.query_list[r][0] = stat.queryVector[r*stat.n_col + c][0];
			stat.query_list[r][1] = stat.queryVector[r*stat.n_col + c][1];
		}

		stat.k = stat.query_list[0][1]; //set the parameter k

		if (is_bucket == false)
			SLAM_SORT(stat);
		if (is_bucket == true)
		{
			stat.gap = stat.query_list[1][0] - stat.query_list[0][0];
			SLAM_BUCKET(stat);
		}

		for (int r = 0; r < stat.n_row; r++)
			stat.out_visual[r][c] = stat.result_list[r];
		clear_SLAM(stat);
	}
}

void SLAM_scan_y(statistics& stat, bool is_bucket)
{
	for (int r = 0; r < stat.n_row; r++)
	{
		for (int c = 0; c < stat.n_col; c++)
		{
			stat.query_list[c][0] = stat.queryVector[r*stat.n_col + c][0];
			stat.query_list[c][1] = stat.queryVector[r*stat.n_col + c][1];
		}

		stat.k = stat.query_list[0][0]; //set the parameter k

		if (is_bucket == false)
			SLAM_SORT(stat);
		if (is_bucket == true)
		{
			stat.gap = stat.query_list[1][1] - stat.query_list[0][1];
			SLAM_BUCKET(stat);
		}

		for (int c = 0; c < stat.n_col; c++)
			stat.out_visual[r][c] = stat.result_list[c];
		clear_SLAM(stat);
	}
}

void SLAM_visual(statistics& stat)
{
	stat.dynamic_pixel_size = stat.n_row;
	stat.static_pixel_size = stat.n_col;
	stat.static_coord = 1;
	stat.dynamic_coord = 0;

	if (stat.method == 7 || stat.method == 8)
	{
		if (stat.n_col > stat.n_row)
		{
			stat.dynamic_pixel_size = stat.n_col;
			stat.static_pixel_size = stat.n_row;
			stat.static_coord = 0;
			stat.dynamic_coord = 1;
		}
	}

	init_SLAM(stat);

	if (stat.method == 5)
		SLAM_scan_x(stat, false);
	if (stat.method == 6)
		SLAM_scan_x(stat, true);

	if (stat.method == 7) //SLAM_SORT^{(RAO)}
	{
		if (stat.n_row > stat.n_col)
			SLAM_scan_x(stat, false);
		else
			SLAM_scan_y(stat, false);
	}
	if (stat.method == 8) //SLAM_BUCKET^{(RAO)}
	{
		if (stat.n_row > stat.n_col)
			SLAM_scan_x(stat, true);
		else
			SLAM_scan_y(stat, true);
	}
}

