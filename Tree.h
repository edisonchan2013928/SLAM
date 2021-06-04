#pragma once
#ifndef TREE_H
#define TREE_H

#include "init_visual.h"
#include "Euclid_Bound.h"

class Node;
class Tree;

class Node
{
public:
	vector<int> idList; //idxs
	vector<Node*> childVector;

	virtual double LB(double*q, statistics& stat) = 0;
	virtual double UB(double*q, statistics& stat) = 0;

	virtual Node*createNode() = 0;
	virtual void update_Aug(Node*node, Tree*t) = 0;
};

class Tree
{
public:
	int dim;
	double**dataMatrix;
	int leafCapacity; //Set to be 20 in tKDC
	Node*rootNode;

	virtual void RQS(statistics& stat) = 0;
};

#endif