
#include "model.h"
#include <random>
model::model(const int mode) {
#ifndef M_PI
#	define M_PI 3.14159265358979323846
#endif
	//model3D.reserve(pow(MODEL_SIZE,3));
	
	/**
	mode  determines initial condition.
	0 : random( 50% filled particle)
	1 : random( density determied by mu)
	2 : ball at the middle of the box(size determined by mu)
    3 : cone at the bottom. 




	*/
	
	std::random_device rd;  //Will be used	 to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(0, 32767);
	blankref[0] = 0;
	blankref[1] = 1;
	double r = 0;
	switch (mode)
	{
			
	case 0:
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_X; j++)
				for (int k = 0; k < MODEL_SIZE_X; k++)
					model3D[i][j][k] = dis(gen) >= 16384;

		break;
	case 1:
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_X; j++)
				for (int k = 0; k < MODEL_SIZE_X; k++)
					model3D[i][j][k] = dis(gen) < 32768.0*N/MODEL_SIZE_X / MODEL_SIZE_X / MODEL_SIZE_X;
		break;
	case 2:
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_X; j++)
				for (int k = 0; k < MODEL_SIZE_X; k++)
					model3D[i][j][k] = ((i-MODEL_SIZE_X/2.0)*(i - MODEL_SIZE_X / 2.0)+ (j - MODEL_SIZE_X / 2.0)*(j - MODEL_SIZE_X / 2.0)+ (k - MODEL_SIZE_X / 2.0)*(k - MODEL_SIZE_X / 2.0)<pow(3.0/4*N/ M_PI,1/3))&& dis(gen) >= 10000 ? 1 : 0;
		break;
	case 3:
		r = pow(3.0 * N / 3.141592, 1.0 / 3);
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_X; j++)
				for (int k = 0; k < MODEL_SIZE_X; k++)
					model3D[i][j][k] =( (i - MODEL_SIZE_X / 2.0)*(i - MODEL_SIZE_X / 2.0) + (j - MODEL_SIZE_X / 2.0)*(j - MODEL_SIZE_X / 2.0)<r*r*(1-k/r)*(1-k/r)&&k<r) ? 1 : 0;
		break;
	case 4:
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_X; j++)
				for (int k = 0; k < MODEL_SIZE_X; k++)
					model3D[i][j][k] = (k<N/MODEL_SIZE_X/MODEL_SIZE_X&&dis(gen) >= 16384);

		break;
	case 5:
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_Y; j++)
				for (int k = 0; k < MODEL_SIZE_Z; k++)
					model3D[i][j][k] = ((i - MODEL_SIZE_X / 2.0)*(i - MODEL_SIZE_X / 2.0) + (j - MODEL_SIZE_Y / 2.0)*(j - MODEL_SIZE_Y / 2.0)<2*N/MODEL_SIZE_Z/M_PI&&k>MODEL_SIZE_Z/2);

		break;
	case 6:
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < MODEL_SIZE_Y; j++)
				for (int k = 0; k < MODEL_SIZE_Z; k++)
					model3D[i][j][k] = 0;
		for (int i = 0; i < MODEL_SIZE_X; i++)
			for (int j = 0; j < N/MODEL_SIZE_X/MODEL_SIZE_Z; j++)
				for (int k = 0; k < MODEL_SIZE_Z; k++)
					model3D[i][j][k] = 1;
		break;
	default:
		break;
	}
}

__int8& model::getdata(int x, int y, int z) {
	blankref[1] = 1;
	blankref[0] = 0;
	if (z < 0)return blankref[1];
	if (z >= MODEL_SIZE_Z || x < 0 || x >= MODEL_SIZE_X || y < 0 || y >= MODEL_SIZE_Y)return blankref[1];
	return std::ref(model3D[x][y][z]);
}