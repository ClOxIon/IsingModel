/*
#include <random>
#include <vector>
#include "montecarlo.h"
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>
#ifndef M_PI
#	define M_PI 3.14159265358979323846
#endif

_inline int bv(bool b) {
	return b ? 1 : -1;
}
_inline int rv(bool b) {
	return b ? 0 : 1;
}
_inline bool getdata(int x, int y, int z, model m) {
	if (z < 0)return 1;
	if (z >= MODEL_SIZE||x < 0 || x >= MODEL_SIZE||y<0||y>=MODEL_SIZE)return 0;
	return m.model3D[x][y][z];

}
_inline int iv(bool b) {
	return b ? 1 : 0;
}
void monteCarloSimulate(model m) {
	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

	std::random_device rd; //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(0, 32767);


	gp << "set xrange [0:"<<MODEL_SIZE<<"]\n";
	gp << "set yrange [0:" << MODEL_SIZE << "]\n";
	gp << "set zrange [0:" << MODEL_SIZE << "]\n";
	gp << "set hidden3d nooffset\n";
	for (int i = 0; i < TEST_NUM; i++) {
		int tE = 0;
		int x = static_cast<int>(dis(gen) / 32767.0*MODEL_SIZE);
		int y = static_cast<int>(dis(gen) / 32767.0*MODEL_SIZE);
		int z = static_cast<int>(dis(gen) / 32767.0*MODEL_SIZE);

		bool t = getdata(x, y, z + 1, m);
		bool b = getdata(x, y, z - 1, m);
		bool u = getdata(x, y+1, z, m);
		bool d = getdata(x, y-1, z, m);
		bool lu, ld, ru, rd;
		if (x % 2 == 0) {
			lu = getdata(x-1, y, z, m);
			ru = getdata(x+1 , y , z, m);
			ld = getdata(x - 1, y -1 , z, m);
			rd = getdata(x+1 , y -1, z, m);
		}
		else {
			lu = getdata(x-1 , y + 1, z, m);
			ru = getdata(x+1, y + 1, z, m);
			ld = getdata(x-1, y  , z, m);
			rd = getdata(x+1, y , z, m);
		}
		double dE = bv(m.model3D[x][y][z])*((ip+ipb*iv(b))*(iv(u) + iv(d) + iv(lu) + iv(ru) + iv(ld) + iv(rd)) + it*iv(t) - ib*rv(b) - 0*g*(z - mu1)+mu2);//very important equation
		m.model3D[x][y][z] ^= (dE < 0) | (dis(gen) < 32767 * exp(-dE / T));//chemical potential : -T(dS/dN)(U,V). 
		for (int i = 0; i < MODEL_SIZE; i++)
			for (int j = 0; j < MODEL_SIZE; j++)
				for (int k = 0; k < MODEL_SIZE; k++) {

					tE += m.model3D[i][j][k] * k;//gravitational potential energy
				}

		if (i%CALC_TO_PRINT == 0) {
			//system("CLS");
			std::vector<boost::tuple<double, double, double> > pts;
			int numParticle[MODEL_SIZE] = { 0 };
			int sum = 0;
			for (int k = MODEL_SIZE - 1; k >= 0; k--) {
				for (int i = 0; i < MODEL_SIZE; i++) {
					for (int j = 0; j < MODEL_SIZE; j++) {
						
						if(m.model3D[i][j][k])pts.push_back(boost::make_tuple(i,(j+((i + 1) % 2)*0.5)*0.866,k));
						//m.model3D[i][j][k] ? printf("o ") : printf("  ");
						numParticle[k] += m.model3D[i][j][k];

					}
					
					//printf("\n");//
					//if (i % 2 == 0)printf(" ");//
					
				}
				//printf("\n");//
				sum += numParticle[k];
			}
			gp << "splot '-' with points\n"; gp.send1d(pts);
			printf("Distribution : \n");
			for (int k = MODEL_SIZE - 1; k >= 0; k--) {
				printf("%d\n", numParticle[k]);
			}
			printf("AvgE = %f\n", tE / static_cast<double>(sum));
			printf("Monte Carlo Calcs : %d\n", i);
			printf("Number of UP Particle : %d / %d\n\n", sum, MODEL_SIZE*MODEL_SIZE*MODEL_SIZE);
			//printf("%d\n", sum);


		}
	}

}
*/