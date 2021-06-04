#include "init_visual.h"

void initStat(int argc, char**argv, statistics& stat)
{
	char*dataFileName = (char*)argv[1];
	stat.outMatrixFileName = (char*)argv[2];
	stat.b_para = atof(argv[3]);
	stat.method = atoi(argv[4]);
	stat.n_row = atoi(argv[5]);
	stat.n_col = atoi(argv[6]);
	stat.kernel_type = atoi(argv[7]);
	stat.epsilon = atof(argv[8]); //Used in method aKDE and QUAD
	
	//debug
	/*char*dataFileName = (char*)"../../../Datasets/test/test";
	stat.outMatrixFileName = (char*)"./Results/test_M8";
	stat.b_para = 1;
	stat.method = 8;
	stat.n_row = 16;
	stat.n_col = 32;
	stat.kernel_type = 1;
	stat.epsilon = 0.05;*/

	extract_FeatureVector(dataFileName, stat);
	updateRegion(stat);
	initQuery(stat);
	init_outVisual(stat);
	find_bandwidth(stat);
}

void updateRegion(statistics& stat)
{
	stat.row_L = inf;
	stat.row_U = -inf;
	stat.col_L = inf;
	stat.col_U = -inf;

	for (int i = 0; i < stat.n; i++)
	{
		if (stat.featureVector[i][0] < stat.row_L)
			stat.row_L = stat.featureVector[i][0];
		if (stat.featureVector[i][0] > stat.row_U)
			stat.row_U = stat.featureVector[i][0];

		if (stat.featureVector[i][1] < stat.col_L)
			stat.col_L = stat.featureVector[i][1];
		if (stat.featureVector[i][1] > stat.col_U)
			stat.col_U = stat.featureVector[i][1];
	}
}

void extract_FeatureVector(char*fileName, statistics& stat)
{
	//load data to feature array
	fstream file;
	int n;
	int dim;

	file.open(fileName);
	if (file.is_open() == false)
	{
		cout << "Cannot Open File!" << endl;
		exit(1);
	}

	file >> n;
	file >> dim;

	stat.n = n;
	if (stat.dim != dim)
	{
		cout << "dimension is not matched!" << endl;
		exit(0);
	}

	stat.featureVector = new double*[n];
	for (int i = 0; i < n; i++)
		stat.featureVector[i] = new double[dim];

	for (int i = 0; i < n; i++)
		for (int d = 0; d < dim; d++)
			file >> stat.featureVector[i][d];

	file.close();
}

void find_bandwidth(statistics& stat)
{
	double sum_x, sum_y;
	double mean_x, mean_y;
	double sd_x, sd_y;
	double h_x, h_y;

	//Using Scott's rule to obtain the bandwidth for the kernel function

	//Compute mean value
	sum_x = 0; sum_y = 0;
	for (int i = 0; i < stat.n; i++)
	{
		sum_x += stat.featureVector[i][0];
		sum_y += stat.featureVector[i][1];
	}
	mean_x = sum_x / stat.n; mean_y = sum_y / stat.n;

	sd_x = 0; sd_y = 0;
	for (int i = 0; i < stat.n; i++)
	{
		sd_x += (stat.featureVector[i][0] - mean_x)*(stat.featureVector[i][0] - mean_x) / (stat.n - 1);
		sd_y += (stat.featureVector[i][1] - mean_y)*(stat.featureVector[i][1] - mean_y) / (stat.n - 1);
	}
	sd_x = sqrt(sd_x); sd_y = sqrt(sd_y);
	h_x = stat.b_para * pow((double)stat.n, -1.0 / 6.0)*sd_x; //d = 2
	h_y = stat.b_para * pow((double)stat.n, -1.0 / 6.0)*sd_y; //d = 2

	stat.b_value = sqrt(h_x*h_x + h_y * h_y);
}

void initQuery(statistics& stat)
{
	int total_q = stat.n_row*stat.n_col;
	double x_coord;
	double y_coord;
	stat.queryVector = new double*[stat.n_row*stat.n_col];

	if (stat.n_row != 1 || stat.n_col != 1)
	{
		stat.incr_row = (stat.row_U - stat.row_L) / (stat.n_row - 1);
		stat.incr_col = (stat.col_U - stat.col_L) / (stat.n_col - 1);
	}

	if (stat.n_row == 1)
		stat.incr_row = 0;
	if (stat.n_col == 1)
		stat.incr_col = 0;

	for (int q = 0; q < total_q; q++)
		stat.queryVector[q] = new double[stat.dim];

	for (int r = 0; r < stat.n_row; r++)
	{
		x_coord = stat.row_L + r * stat.incr_row;
		for (int c = 0; c < stat.n_col; c++)
		{
			y_coord = stat.col_L + c * stat.incr_col;
			stat.queryVector[r*stat.n_col + c][0] = x_coord;
			stat.queryVector[r*stat.n_col + c][1] = y_coord;
		}
	}
}

void init_outVisual(statistics& stat)
{
	stat.out_visual = new double*[stat.n_row];
	for (int r = 0; r < stat.n_row; r++)
		stat.out_visual[r] = new double[stat.n_col];
}

void output_visual(statistics& stat)
{
	fstream outMatrixFile;
	double x_coord, y_coord, KDE_value;

	outMatrixFile.open(stat.outMatrixFileName, ios::in | ios::out | ios::trunc);

	if (outMatrixFile.is_open() == false)
	{
		cout << "Cannot open output file!" << endl;
		exit(0);
	}

	for (int r = 0; r < stat.n_row; r++)
	{
		for (int c = 0; c < stat.n_col; c++)
		{
			x_coord = stat.queryVector[r*stat.n_col + c][0];
			y_coord = stat.queryVector[r*stat.n_col + c][1];
			KDE_value = stat.out_visual[r][c];

			outMatrixFile << x_coord << " " << y_coord << " " << KDE_value << endl;
		}
	}

	outMatrixFile.close();
}