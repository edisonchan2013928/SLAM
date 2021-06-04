#include "kd_tree.h"

//kdNode
double kdNode::LB(double*q, statistics& stat)
{
	double ub;
	double temp_value;
	double LB_value;

	ub = u_MBR(q, boundary, stat.dim);
	if (ub > stat.b_value)
		return 0;

	if (stat.kernel_type == 0) //Uniform kernel
		LB_value = idList.size() / stat.b_value;
	if (stat.kernel_type == 1) //Epanechnikov kernel
		LB_value = idList.size()*(1.0 - (1.0 / (stat.b_value*stat.b_value))*ub*ub);
	if (stat.kernel_type == 2) //Quartic kernel
	{
		temp_value = 1.0 - (1.0 / (stat.b_value*stat.b_value))*ub*ub;
		LB_value = idList.size()*temp_value*temp_value;
	}

	return LB_value;
}

double kdNode::UB(double*q, statistics& stat)
{
	double lb;
	double temp_value;
	double UB_value;

	lb = ell_MBR(q, boundary, stat.dim);
	if (lb > stat.b_value)
		return 0;

	if (stat.kernel_type == 0) //Uniform kernel
		UB_value = idList.size() / stat.b_value;
	if (stat.kernel_type == 1) //Epanechnikov kernel
		UB_value = idList.size()*(1.0 - (1.0 / (stat.b_value*stat.b_value))*lb*lb);
	if (stat.kernel_type == 2) //Quartic kernel
	{
		temp_value = 1.0 - (1.0 / (stat.b_value*stat.b_value))*lb*lb;
		UB_value = idList.size()*temp_value*temp_value;
	}

	return UB_value;
}

int kdNode::check_condition_boundaries(statistics& stat)
{
	//initialize to the value 0 (Cover)
	int condition = 0;
	int cur_condition = 2; //This state is "not intersect"

	for (int d = 0; d < stat.dim; d++)
	{
		if (stat.query_boundary[d][0] <= boundary[d][0] && boundary[d][1] <= stat.query_boundary[d][1])
			cur_condition = 0;
		if (boundary[d][0] <= stat.query_boundary[d][0] && stat.query_boundary[d][1] <= boundary[d][1])
			cur_condition = 1;
		if (stat.query_boundary[d][0] <= boundary[d][0] && boundary[d][0] <= stat.query_boundary[d][1] && stat.query_boundary[d][1] <= boundary[d][1])
			cur_condition = 1;
		if (boundary[d][0] <= stat.query_boundary[d][0] && stat.query_boundary[d][0] <= boundary[d][1] && boundary[d][1] <= stat.query_boundary[d][1])
			cur_condition = 1;

		if (cur_condition == 2) //N + C => N, N + I => N, //N + N => N 
		{
			condition = 2;
			break;
		}

		if (cur_condition == 1 && (condition == 0 || condition == 1))
			condition = 1;
	}

	return condition;
}

kdNode*kdNode::createNode()
{
	return new kdNode();
}

//kdQuadAugNode
double kdQuadAugNode::LB(double*q, statistics& stat)
{
	int node_size;
	double ip = 0;
	double LB_value = 0;
	//double bandwidth = stat.b_value;
	double sq_bandwidth;

	for (int d = 0; d < stat.dim; d++)
		ip = ip + q[d] * a_G[d];

	if (stat.kernel_type == 0) //Uniform kernel
	{
		//code here
		return LB_value;
	}
	if (stat.kernel_type == 1) //Epanechnikov kernel
	{
		sq_bandwidth = stat.b_value * stat.b_value;
		node_size = this->idList.size();
		LB_value = max(node_size - (node_size / sq_bandwidth)*stat.q_SquareNorm
			+ (2.0 / sq_bandwidth) *ip - (1.0 / sq_bandwidth) *S_G, 0.0);
		return LB_value;
	}
	if (stat.kernel_type == 2) //Quartic kernel
	{
		//code here
		return LB_value;
	}

	cout << "Error! We cannot support this type of kernel functions!" << endl;
	exit(0);
}

double kdQuadAugNode::UB(double*q, statistics& stat)
{
	int node_size;
	double ip = 0;
	double UB_value = 0;
	double lb, ub;
	double x_min, x_max;
	double y_min, y_max;
	double m, c;
	//double bandwidth = stat.b_value;
	double sq_bandwidth;

	for (int d = 0; d < stat.dim; d++)
		ip = ip + q[d] * a_G[d];

	if (stat.kernel_type == 0) //Uniform kernel
	{
		//code here
		return UB_value;
	}
	if (stat.kernel_type == 1) //Epanechnikov kernel
	{
		sq_bandwidth = stat.b_value * stat.b_value;
		node_size = this->idList.size();
		ub = u_MBR(q, this->boundary, stat.dim);
		if (ub < stat.b_value)
			UB_value = node_size - (node_size / sq_bandwidth)*stat.q_SquareNorm
			+ (2.0 / sq_bandwidth) *ip - (1.0 / sq_bandwidth) *S_G;
		else
		{
			lb = ell_MBR(q, this->boundary, stat.dim);
			if (lb > stat.b_value)
				return 0;

			x_min = lb * lb;
			x_max = ub * ub;
			y_min = 1 - 1 / (sq_bandwidth)*x_min;
			y_max = 0;

			m = (y_max - y_min) / (x_max - x_min);
			c = (x_max*y_min - x_min * y_max) / (x_max - x_min);
			UB_value = m * (node_size*stat.q_SquareNorm - 2 * ip + S_G) + c * node_size;
		}
		return UB_value;
	}
	if (stat.kernel_type == 2) //Quartic kernel
	{
		//code here
		return UB_value;
	}

	cout << "Error! We cannot support this type of kernel functions!" << endl;
	exit(0);
}

void kdQuadAugNode::update_AugInfo(kdQuadAugNode*node, Tree*t)
{
	int id;
	double norm_square;
	//obtain vec a_G
	node->a_G = new double[t->dim];

	for (int d = 0; d < t->dim; d++)
		node->a_G[d] = 0;

	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		for (int d = 0; d < t->dim; d++)
			node->a_G[d] += t->dataMatrix[id][d];
	}

	//obtain S_G
	node->S_G = 0;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		norm_square = 0;
		for (int d = 0; d < t->dim; d++)
			norm_square = norm_square + t->dataMatrix[id][d] * t->dataMatrix[id][d];

		node->S_G += norm_square;
	}

	if ((int)node->idList.size() <= t->leafCapacity) //this is the leaf node
		return;

	update_AugInfo((kdQuadAugNode*)node->childVector[0], t);
	update_AugInfo((kdQuadAugNode*)node->childVector[1], t);
}

void kdQuadAugNode::update_Aug(Node*node, Tree*t)
{
	this->update_AugInfo((kdQuadAugNode*)node, t);
}

kdQuadAugNode*kdQuadAugNode::createNode()
{
	return new kdQuadAugNode();
}

//Construct kd-tree
kdTree::kdTree(int dim, double**dataMatrix, int leafCapacity)
{
	this->dim = dim;
	this->dataMatrix = dataMatrix;
	this->leafCapacity = leafCapacity;
}

void kdTree::getNode_Boundary(kdNode*node)
{
	//double sum_alpha;
	int id;
	node->boundary = new double*[dim];
	for (int d = 0; d < dim; d++)
		node->boundary[d] = new double[2];

	for (int d = 0; d < dim; d++)
	{
		node->boundary[d][0] = inf;
		node->boundary[d][1] = -inf;
	}

	for (int d = 0; d < dim; d++)
	{
		for (int i = 0; i < (int)node->idList.size(); i++)
		{
			id = node->idList[i];
			if (dataMatrix[id][d] < node->boundary[d][0])
				node->boundary[d][0] = dataMatrix[id][d];

			if (dataMatrix[id][d] > node->boundary[d][1])
				node->boundary[d][1] = dataMatrix[id][d];
		}
	}
}

double kdTree::obtain_SplitValue(kdNode*node, int split_Dim)
{
	vector<double> tempVector;
	int id;
	int middle_left, middle_right, middle;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		tempVector.push_back(dataMatrix[id][split_Dim]);
	}

	sort(tempVector.begin(), tempVector.end());

	if ((int)tempVector.size() % 2 == 0)//even number
	{
		middle_right = (int)tempVector.size() / 2;
		middle_left = middle_right - 1;

		return ((tempVector[middle_left] + tempVector[middle_right]) / 2.0);
	}
	else
	{
		middle = ((int)tempVector.size() - 1) / 2;
		return tempVector[middle];
	}

	tempVector.clear();
}

void kdTree::KD_Tree_Recur(kdNode*node, int split_Dim)
{
	int id;
	int counter;
	//base case
	if ((int)node->idList.size() <= leafCapacity)
		return;

	//added for testing
	//split_Dim=select_split_Dim(node);
	double splitValue = obtain_SplitValue(node, split_Dim); //code here

	//create two children
	kdNode*leftNode;
	kdNode*rightNode;

	leftNode = node->createNode();
	rightNode = node->createNode();

	counter = 0;
	int halfSize = ((int)node->idList.size()) / 2;
	for (int i = 0; i < (int)node->idList.size(); i++)
	{
		id = node->idList[i];
		if (dataMatrix[id][split_Dim] <= splitValue && counter <= halfSize)
		{
			leftNode->idList.push_back(id);
			counter++;
		}
		else
			rightNode->idList.push_back(id);
	}

	getNode_Boundary(leftNode);
	getNode_Boundary(rightNode);

	KD_Tree_Recur(leftNode, (split_Dim + 1) % dim);
	KD_Tree_Recur(rightNode, (split_Dim + 1) % dim);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void kdTree::build_kdTree(statistics& stat)
{
	for (int i = 0; i < stat.n; i++)
		rootNode->idList.push_back(i);

	getNode_Boundary((kdNode*)rootNode);
	KD_Tree_Recur((kdNode*)rootNode, 0);
}

void kdTree::updateAugment(kdNode*node)
{
	node->update_Aug(node, this);
}


void kdTree::init_RQS(statistics& stat)
{
	stat.query_boundary = new double*[stat.dim];
	for (int d = 0; d < stat.dim; d++)
		stat.query_boundary[d] = new double[2];
}

//Follows the algorithm in the book (Computational Geometry: Algorithms and Applications (Second Edition) p.103)
double kdTree::RQS_Recur(kdNode*node, statistics& stat)
{
	int condition;
	kdNode*child_node;
	double density_value;

	if (node->childVector.size() == 0) //child node
		return refinement(node, stat);
	else
	{
		density_value = 0;
		for (int c = 0; c < (int)node->childVector.size(); c++)
		{
			child_node = (kdNode*)node->childVector[c];
			condition = child_node->check_condition_boundaries(stat);
			if (condition == 0) //Cover
				density_value += refinement(child_node, stat);
			else if (condition == 1) //Intersect
				density_value += RQS_Recur(child_node, stat);
		}

		return density_value;
	}
}

void kdTree::RQS(statistics& stat)
{
	int condition;
	double density_value = 0;

	obtain_boundary(stat);
	condition = ((kdNode*)rootNode)->check_condition_boundaries(stat);
	if (condition == 0)
		density_value = refinement(rootNode, stat);
	if (condition == 1)
		density_value = RQS_Recur((kdNode*)rootNode, stat);

	stat.out_visual[stat.cur_r][stat.cur_c] = density_value;
}

void kdTree::obtain_boundary(statistics& stat)
{
	double*cur_q = stat.queryVector[stat.cur_r*stat.n_col + stat.cur_c];
	stat.query_boundary[0][0] = cur_q[0] - stat.b_value;
	stat.query_boundary[0][1] = cur_q[0] + stat.b_value;
	stat.query_boundary[1][0] = cur_q[1] - stat.b_value;
	stat.query_boundary[1][0] = cur_q[1] + stat.b_value;
}