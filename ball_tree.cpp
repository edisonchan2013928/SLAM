#include "ball_tree.h"

double ballNode::LB(double*q, statistics& stat)
{
	cout << "No bound function for ball tree!" << endl;
	exit(0);
	return 0;
}

double ballNode::UB(double*q, statistics& stat)
{
	cout << "No bound function for ball tree!" << endl;
	exit(0);
	return 0;
}

void ballNode::update_Augment(double**dataMatrix, int dim)
{
	int id;
	int num_of_points;
	double cur_radius;

	center = new double[dim];
	for (int d = 0; d < dim; d++)
		center[d] = 0;

	num_of_points = (int)idList.size();
	//obtain the center for this node
	for (int i = 0; i < num_of_points; i++)
	{
		id = idList[i];
		for (int d = 0; d < dim; d++)
			center[d] += dataMatrix[id][d] / num_of_points;
	}

	radius = 0;
	for (int i = 0; i < num_of_points; i++)
	{
		id = idList[i];
		cur_radius = euclid_dist(center, dataMatrix[id], dim);
		radius = max(radius, cur_radius);
	}

}

ballNode*ballNode::createNode()
{
	return new ballNode();
}

//ballTree
ballTree::ballTree(int dim, double**dataMatrix, int leafCapacity)
{
	this->dim = dim;
	this->dataMatrix = dataMatrix;
	this->leafCapacity = leafCapacity;
}

void ballTree::divide_node(ballNode*node, ballNode*leftNode, ballNode*rightNode)
{
	int index_1, index_2;
	int id, id_1, id_2;
	int best_id_1, best_id_2;
	double cur_dist;
	double max_dist = 0;
	int num_of_points;
	num_of_points = (int)node->idList.size();
	int half_size = (int)floor(num_of_points / 2.0);

	for (int r = 0; r < rand_num; r++)
	{
		index_1 = rand() % num_of_points; index_2 = rand() % num_of_points;
		id_1 = node->idList[index_1]; id_2 = node->idList[index_2];
		cur_dist = euclid_dist(dataMatrix[id_1], dataMatrix[id_2], dim);

		if (cur_dist > max_dist)
		{
			max_dist = cur_dist;
			best_id_1 = id_1;
			best_id_2 = id_2;
		}
	}

	for (int i = 0; i < num_of_points; i++)
	{
		id = node->idList[i];
		if ((int)leftNode->idList.size() >= half_size)
		{
			rightNode->idList.push_back(id);
			continue;
		}
		if ((int)rightNode->idList.size() >= half_size)
		{
			leftNode->idList.push_back(id);
			continue;
		}

		if (euclid_dist(dataMatrix[id], dataMatrix[id_1], dim) < euclid_dist(dataMatrix[id], dataMatrix[id_2], dim))
			leftNode->idList.push_back(id);
		else
			rightNode->idList.push_back(id);
	}
}

void ballTree::ballTree_Recur(ballNode*node)
{
	ballNode*leftNode;
	ballNode*rightNode;
	int num_of_points;

	num_of_points = (int)node->idList.size();
	//base case
	if (num_of_points <= leafCapacity)
		return;

	//create two children
	leftNode = node->createNode();
	rightNode = node->createNode();
	divide_node(node, leftNode, rightNode);
	leftNode->update_Augment(dataMatrix, dim);
	rightNode->update_Augment(dataMatrix, dim);

	ballTree_Recur(leftNode);
	ballTree_Recur(rightNode);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void ballTree::build_ballTree(statistics& stat)
{
	for (int i = 0; i < stat.n; i++)
		rootNode->idList.push_back(i);
	((ballNode*)rootNode)->update_Augment(dataMatrix, dim);

	ballTree_Recur((ballNode*)rootNode);
}

double ballTree::RQS_Recur(ballNode*node, statistics& stat)
{
	ballNode*childNode;
	double dist_value;
	double density_value;

	if ((int)node->childVector.size() == 0)
		return refinement(node, stat);

	density_value = 0;
	for (int c = 0; c < (int)node->childVector.size(); c++)
	{
		childNode = (ballNode*)node->childVector[c];
		dist_value = euclid_dist(stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c], childNode->center, dim);

		if (dist_value + childNode->radius <= stat.b_value)
		{
			density_value += refinement(childNode, stat);
			continue;
		}

		if (dist_value >= stat.b_value + childNode->radius)
			continue;

		density_value += RQS_Recur(childNode, stat);
	}

	return density_value;
}

void ballTree::RQS(statistics& stat)
{
	double dist_value = euclid_dist(stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c], ((ballNode*)rootNode)->center, dim);
	double density_value;

	if (dist_value + ((ballNode*)rootNode)->radius <= stat.b_value)
		density_value = refinement(rootNode, stat);
	else if (dist_value >= stat.b_value + ((ballNode*)rootNode)->radius)
		density_value = 0;
	else
		density_value = RQS_Recur((ballNode*)rootNode, stat);

	stat.out_visual[stat.cur_r][stat.cur_c] = density_value;
}