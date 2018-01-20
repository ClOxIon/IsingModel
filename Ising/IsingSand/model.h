#pragma once

#include "config2.h"
#include <vector>

class model {
	__int8 model3D[MODEL_SIZE_X][MODEL_SIZE_Y][MODEL_SIZE_Z];
	__int8 blankref[2];
public:
	
	double unpushed_parts;
	
	model(const int mode);
	void push(double parts);
	__int8& getdata(int x, int y, int z);
};
//std::vector<model> modelMemory;