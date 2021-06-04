#include "Euclid_Bound.h"

double ell_MBR(double*q,double**boundary,int dim)
{
	double ell=0;
	double temp;
	for(int d=0;d<dim;d++)
	{
		if(q[d]<boundary[d][0])
			temp=boundary[d][0]-q[d];
		else
		{
			if(q[d]>boundary[d][1])
				temp=q[d]-boundary[d][1];
			else
				temp=0;
		}
		
		ell+=temp*temp;
	}
	ell=sqrt(ell);

	return ell;
}

double u_MBR(double*q,double**boundary,int dim)
{
	double u=0;
	double temp;
	for(int d=0;d<dim;d++)
	{
		temp=max(fabs(q[d]-boundary[d][0]),fabs(q[d]-boundary[d][1]));
		u+=temp*temp;
	}
	u=sqrt(u);

	return u;
}

double u_tri(double*q,double*center,int dim,double radius,double& obt_dist)
{
	double u;
	obt_dist=0;

	for(int d=0;d<dim;d++)
		obt_dist+=(q[d]-center[d])*(q[d]-center[d]);

	obt_dist=sqrt(obt_dist);

	u=obt_dist+radius;

	return u;

}

double euclid_dist(double*q,double*p,int dim)
{
	double dist=0;
	for(int d=0;d<dim;d++)
		dist+=(q[d]-p[d])*(q[d]-p[d]);
	return sqrt(dist);
}

double sq_euclid_dist(double*q, double*p, int dim)
{
	double dist = 0;
	for (int d = 0; d < dim; d++)
		dist += (q[d] - p[d])*(q[d] - p[d]);
	return dist;
}

double computeSqNorm(double*q, int dim)
{
	double sqNorm = 0;
	for (int d = 0; d < dim; d++)
		sqNorm += q[d] * q[d];

	return sqNorm;
}

double inner_product(double*q, double*p, int dim)
{
	double ip = 0;

	for (int d = 0; d < dim; d++)
		ip += q[d] * p[d];

	return ip;
}

double ell_MBR_MBR(double**boundary_1, double**boundary_2, int dim)
{
	double MIN_Q_P;
	double ell;

	ell = 0;
	for (int d = 0; d < dim; d++)
	{
		if (boundary_1[d][0] > boundary_2[d][1])
			MIN_Q_P = boundary_1[d][0] - boundary_2[d][1];
		else
		{
			if (boundary_2[d][0] > boundary_1[d][1])
				MIN_Q_P = boundary_2[d][0] - boundary_1[d][1];
			else
				MIN_Q_P = 0;
		}
		ell += MIN_Q_P * MIN_Q_P;
	}

	return sqrt(ell);
}

double u_MBR_MBR(double**boundary_1, double**boundary_2, int dim)
{
	double MAX_Q_P;
	double u;

	u = 0;

	for (int d = 0; d < dim; d++)
	{
		MAX_Q_P = max(fabs(boundary_1[d][1] - boundary_2[d][0]), fabs(boundary_2[d][1] - boundary_1[d][0]));
		u += MAX_Q_P * MAX_Q_P;
	}

	return sqrt(u);
}


//Used in 2d-visualization
/*double corner_maxDist_sq(double*a_G_avg, double**boundary)
{
	double max_dist_sq = 0;
	double dist_sq_1 = 0;
	double dist_sq_2 = 0;
	double dist_sq_3 = 0;
	double dist_sq_4 = 0;

	dist_sq_1 = (a_G_avg[0] - boundary[0][0])*(a_G_avg[0] - boundary[0][0]) +
		(a_G_avg[1] - boundary[1][0])*(a_G_avg[1] - boundary[1][0]);
	dist_sq_2 = (a_G_avg[0] - boundary[0][1])*(a_G_avg[0] - boundary[0][1]) +
		(a_G_avg[1] - boundary[1][0])*(a_G_avg[1] - boundary[1][0]);
	dist_sq_3 = (a_G_avg[0] - boundary[0][0])*(a_G_avg[0] - boundary[0][0]) +
		(a_G_avg[1] - boundary[1][1])*(a_G_avg[1] - boundary[1][1]);
	dist_sq_4 = (a_G_avg[0] - boundary[0][1])*(a_G_avg[0] - boundary[0][1]) +
		(a_G_avg[1] - boundary[1][1])*(a_G_avg[1] - boundary[1][1]);

	max_dist_sq = max(max(dist_sq_1, dist_sq_2), max(dist_sq_3, dist_sq_4));

	return max_dist_sq;
}*/