#pragma once
#ifndef SS_VISUAL_H
#define SS_VISUAL_H

#include "init_visual.h"
#include "Euclid_Bound.h"
#include "kd_tree.h"

double SCAN(double*q, statistics& stat);
void KDE_visual(statistics& stat);
double refinement(Node*curNode, statistics& stat); //refinement in the leaf node

#endif