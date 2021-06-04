#include "alg_visual.h"

int main(int argc, char**argv)
{
	statistics stat;
	initStat(argc, argv, stat);
	visual_Algorithm(stat);
}