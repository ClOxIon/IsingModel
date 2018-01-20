#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "configuration.h"
_inline int bv(bool b) {
	return static_cast<int>(b) * 2 - 1;
}

int main() {
	int row, col,dE;bool t, b, r, l;
	srand(time(0));
	bool model2D[MODEL_SIZE][MODEL_SIZE];
	for (int i = 0; i < MODEL_SIZE; i++) {
		for (int j = 0; j < MODEL_SIZE; j++)
			model2D[i][j] = rand() >= 16384;
	}
	for (int i = 0; i < N; i++) {
		row = rand() / 32767.0*MODEL_SIZE;
		col = rand() / 32767.0*MODEL_SIZE;
		
		t = (row == 0) ? model2D[MODEL_SIZE - 1][col] : model2D[row-1][col];
		b = (row == MODEL_SIZE - 1) ? model2D[1][col] : model2D[row + 1][col];
		l = (col == 0) ? model2D[row][MODEL_SIZE] : model2D[row][col - 1];
		r = (col == MODEL_SIZE - 1) ? model2D[row][1] : model2D[row][col + 1];
		dE = 2 * bv(model2D[row][col]) *(bv(t) + bv(b) + bv(l) + bv(r));
		model2D[row][col] ^= (dE < 0|rand()<32767*exp(-dE/(k*T)));
	}
	for (int i = 0; i < MODEL_SIZE; i++) {
		for (int j = 0; j < MODEL_SIZE; j++)
			printf("%d", model2D[i][j]);
		printf("\n");

	}

	
}