#pragma once
#ifndef KD_TREE_H
#define KD_TREE_H

#include "Tree.h"
#include "SS_visual.h"

class kdNode : public Node  //Gan_SIGMOD17
{
public:
	double**boundary; //boundary

	double LB(double*q, statistics& stat);
	double UB(double*q, statistics& stat);

	int check_condition_boundaries(statistics& stat); //Three conditions: Cover (0), Intersect (1) and Not intersect (2)

	void update_Aug(Node*node, Tree*t) {}

	kdNode*createNode();
};

class kdQuadAugNode : public kdNode //QUAD
{
public:
	double*a_G;
	double S_G;

	//facilitates online (sharing) computation
	double gamma_sum;

	double LB(double*q, statistics& stat);
	double UB(double*q, statistics& stat);

	void update_AugInfo(kdQuadAugNode*node, Tree*t);

	void update_Aug(Node*node, Tree*t);
	kdQuadAugNode*createNode();
};

class kdTree : public Tree
{
public:
	//Member functions
	kdTree(int dim, double**dataMatrix, int leafCapacity);

	void getNode_Boundary(kdNode*node);
	double obtain_SplitValue(kdNode*node, int split_Dim);
	void KD_Tree_Recur(kdNode*node, int split_Dim);
	void build_kdTree(statistics& stat);

	void updateAugment(kdNode*node);

	//RQS (online phase)
	void init_RQS(statistics& stat);
	double RQS_Recur(kdNode*node, statistics& stat);
	void RQS(statistics& stat);
	void obtain_boundary(statistics& stat);
};

#endif