#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <random>
#include <Windows.h>

#include "configuration.h"
//tetrahedral model		
_inline int bv(bool b) {
	return (static_cast<int>(b) << 1) - 1;
}
//qwetqwe
int main() {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(0, 32767);
	//std::binomial_distribution<> i(1,p);
		int row, col, dE, tE = 0; bool t, b, r, l;
		srand(time(0));
		bool model2D[MODEL_SIZE][MODEL_SIZE];
		for (int i = 0; i < MODEL_SIZE; i++) {
			for (int j = 0; j < MODEL_SIZE; j++)
				model2D[i][j] = dis(gen) >= 16384;
		}

		for (int i = 0; i < N; i++) {
			tE = 0;
			row = static_cast<int>(dis(gen) / 32767.0*MODEL_SIZE);
			col = static_cast<int>(dis(gen) / 32767.0*MODEL_SIZE);
			int s = bv(model2D[row][col]);
			t = (row == 0) ? model2D[MODEL_SIZE - 1][col] : model2D[row - 1][col];
			b = (row == MODEL_SIZE - 1) ? model2D[1][col] : model2D[row + 1][col];
			l = (col == 0) ? model2D[row][MODEL_SIZE] : model2D[row][col - 1];
			r = (col == MODEL_SIZE - 1) ? model2D[row][1] : model2D[row][col + 1];
			dE = s *(2*(bv(t) + bv(b) + bv(l) + bv(r))-u*T);//2*4*(2*upR-1)d
			model2D[row][col] ^= (dE < 0) | (dis(gen) < 32767 * exp(-dE / T));//chemical potential : -T(dS/dN)(U,V). beware of probability overflow!    dE_Max = 2*Coord.    exp(-u*s) ,  exp(-dE/T+u*s) / exp(u*s)

			for (int i = 0; i < MODEL_SIZE; i++) {
				for (int j = 0; j < MODEL_SIZE; j++) {
					t = (i == 0) ? model2D[MODEL_SIZE - 1][j] : model2D[i - 1][j];
					b = (i == MODEL_SIZE - 1) ? model2D[1][j] : model2D[i + 1][j];
					l = (j == 0) ? model2D[i][MODEL_SIZE] : model2D[i][j - 1];
					r = (j == MODEL_SIZE - 1) ? model2D[i][1] : model2D[i][j + 1];
					tE -= bv(model2D[i][j]) *(bv(t) + bv(b) + bv(l) + bv(r));
				}
			}
			if (i%CALC_TO_PRINT == 0) {
				int sum = 0;
			for (int i = 0; i < MODEL_SIZE; i++) {
				for (int j = 0; j < MODEL_SIZE; j++) {
					printf("%d", model2D[i][j]);
					sum += model2D[i][j];
				}
				printf("\n");

			}
			
				printf("AvgE = %f\n", tE / static_cast<double>(MODEL_SIZE*MODEL_SIZE));
				printf("Monte Carlo Calcs : %d\n", i);
				printf("Number of UP Particle : %d / %d\n\n", sum, MODEL_SIZE*MODEL_SIZE);
				Sleep(300);
			}
		}
}
double upDownRatioCalc() {
	double r[2] = { 0.5 ,0.5};
	for (int i = 0; abs(r[i % 2] / r[i % 2 + 1] - 1) < 1e-5; i++) {
		double q = r[i % 2];
		//r[i % 2 + 1] = ((pow(1 - q, 4) + q*pow(1 - q, 3))*exp(-u) + (1 - pow(1 - q, 4) - q*pow(1 - q, 3))*exp(-8 * (2 * q - 1) / T))/(
			//pow(q,4)+(1-q)*pow(q,3)+(1- pow(q, 4) - (1 - q)*pow(q, 3))*e();
	
	}
	return 0.1;
}