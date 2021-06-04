#pragma once
#ifndef INIT_VISUAL_H
#define INIT_VISUAL_H

#include "Library.h"

const double inf = 999999999999;
const double small_epsilon = 0.000000000000001;
const double pi = 3.14159265358979323846;

struct statistics
{
	int n; //number of points in datasets
	double b_value; //bandwidth value
	double b_para; //Gan et al. SIGMOD2017
	int n_row; //number of discrete region in the row, e.g. 100
	int n_col; //number of discrete region in the col, e.g. 100
	char*outMatrixFileName; //output file name

	double row_L, row_U; //row region, e.g. [-2,2]
	double col_L, col_U; //col region e.g. [-2,2]
	double incr_row; //incremental row e.g. (2-(-2))/100=0.04
	double incr_col; //incremental col e.g. (2-(-2))/100=0.04
	double**featureVector; //feature vector of all data points
	double**queryVector; //query vector
	double**out_visual; //output the visualization with size n_row x n_col

	//Used in indexing framework
	int cur_r;
	int cur_c;
	double q_SquareNorm;
	double**query_boundary;

	//Used in SLAM
	//(If static_coord is 1 (y-coordinate), dynamic_coord is 0 (x-coordinate))
	//(If static_coord is 0 (x-coordinate), dynamic_coord is 1 (y-coordinate))
	int static_coord;
	int dynamic_coord;
	int static_pixel_size;
	int dynamic_pixel_size;
	vector<int> L_ell;
	vector<int> U_ell;
	double*A_L_ell;
	double*A_U_ell;
	double S_L_ell;
	double S_U_ell;
	double k;
	vector<double*> query_list;
	vector<double> result_list;
	double R_q_card;
	double*A_R_q;
	double S_R_q;
	//Used in SLAM_BUCKET
	vector< vector<int> > B_L;
	vector< vector<int> > B_U;
	double gap;

	//Used in Z-order method
	int ori_n;

	const int dim = 2; //dim = 2 (2d visualization)
	const int leafCapacity = 40; //leaf capacity for the P set
	int method; //chosen method
	double epsilon; //epsilon value for epsilon-KVQ
	int kernel_type; //kernel type //0: Uniform kernel, 1: Epanechnikov kernel 2: Quartic kernel
};

void initStat(int argc, char**argv, statistics& stat);
void updateRegion(statistics& stat);
void extract_FeatureVector(char*fileName, statistics& stat);
void find_bandwidth(statistics& stat);
void initQuery(statistics& stat);
void init_outVisual(statistics& stat);
void output_visual(statistics& stat);

#endif