#pragma once
#ifndef BALL_TREE_H
#define BALL_TREE_H

#include "Tree.h"
#include "SS_visual.h"

class ballNode : public Node
{
public:
	double*center;
	double radius;

	double LB(double*q, statistics& stat);
	double UB(double*q, statistics& stat);
	void update_Augment(double**dataMatrix, int dim);
	void update_Aug(Node*node, Tree*t) {}
	ballNode*createNode();
};

class ballTree : public Tree
{
public:
	const int rand_num = 5;
	//offline
	ballTree(int dim, double**dataMatrix, int leafCapacity);
	void divide_node(ballNode*node, ballNode*leftNode, ballNode*rightNode);
	void ballTree_Recur(ballNode*node);
	void build_ballTree(statistics& stat);

	//RQS method
	double RQS_Recur(ballNode*node, statistics& stat);
	void RQS(statistics& stat);
};

#endif