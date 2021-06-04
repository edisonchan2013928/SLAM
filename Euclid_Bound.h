#ifndef EUCLID_BOUND_H
#define EUCLID_BOUND_H

#include "init_visual.h"

double ell_MBR(double*q,double**boundary,int dim);
double u_MBR(double*q,double**boundary,int dim);
double u_tri(double*q,double*center,int dim,double radius,double& obt_dist);
double euclid_dist(double*q,double*p,int dim);
double sq_euclid_dist(double*q, double*p, int dim);
double computeSqNorm(double*q, int dim);
double inner_product(double*q, double*p, int dim);

//Used in dual-group setting
double ell_MBR_MBR(double**boundary_1, double**boundary_2, int dim);
double u_MBR_MBR(double**boundary_1, double**boundary_2, int dim);

//corner max-dist
//double corner_maxDist_sq(double*a_G_avg, double**boundary);

#endif