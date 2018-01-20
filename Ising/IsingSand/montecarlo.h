
#pragma once
#include "model.h"
#include <iostream>

void monteCarloSimulate(model& m);
void print(model& m);
void exportSlope(model& m, FILE* fp, int i);
void exportSlope2D(model& m, FILE* fp, int i);