#include <vector>
#include "montecarlo.h"
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>
#include <random>
#include <ctime>
void print(model& m) {

	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
	gp << "set xrange [0:" << MODEL_SIZE_Z << "]\n";
	gp << "set yrange [0:" << MODEL_SIZE_Z << "]\n";
	gp << "set zrange [0:" << MODEL_SIZE_Z << "]\n";
	gp << "set hidden3d nooffset\n";





	

		std::vector<boost::tuple<double, double, double> > pts;
		pts.reserve(2*N);
		int numParticle[MODEL_SIZE_Z] = { 0 };
		int nPsum = 0;
		double tE = 0;
		for (int k = MODEL_SIZE_Z - 1; k >= 0; k--) {
			for (int i = 0; i < MODEL_SIZE_X; i++) {
				for (int j = 0; j < MODEL_SIZE_Y; j++) {

					if (m.getdata(i,j,k)&&rand()<32768/HIDE_RATIO)pts.push_back(boost::make_tuple(i*0.8660254+0.5773503 * (k % 2), j + (i % 2)*0.5 , k*0.8164965));//(sqrt3)/2,(sqrt3)/3,sqrt(2/3)
					//m.model3D[i][j][k] ? printf("o ") : printf("  ");
					numParticle[k] += m.getdata(i,j,k);
					tE += m.getdata(i, j, k) * k * 0.8164965 * RADIUS * 2 *MASS* G;//gravitational potential energy
				}

				//printf("\n");//
				//if (i % 2 == 0)printf(" ");//

			}
			//printf("\n");//
			nPsum += numParticle[k];
		}
		//gp << "plot"<<gp.file1d(pts,"file.dat")<<"with points\n";//
		gp << "splot '-' with points\n"; gp.send1d(pts);
		/*
		printf("Distribution : \n");
		for (int k = MODEL_SIZE - 1; k >= 0; k--) {
		printf("%d\n", numParticle[k]);
		}
		*/
		printf("\nAvgE = %f\n", tE / static_cast<double>(nPsum));
		
		printf("Number of Particle : %d / %d\n", nPsum, MODEL_SIZE_X*MODEL_SIZE_Y*MODEL_SIZE_Z);
		//printf("%d\n", sum);
	

		
		if (nPsum < N / 5) {
			printf("particles are almost gone ; stopping the simulation...\n"); return;

		}
		

		
		system("pause");
}
void exportSlope(model& m, FILE* fp, int i) {
	

	fprintf(fp, "%d,", i);
	//divide model into circles :
	double avgHeight[15];
	for (int i = 1; i < STEP_SLOPE; i++) {
		int heightSum = 0;
		int cnt = 0;
		
		for (int x = 0; x<MODEL_SIZE_X; x++)
			for (int y = 0; y < MODEL_SIZE_Y; y++)
			{
				double rSquare = (x - MODEL_SIZE_X / 2.0)*(x - MODEL_SIZE_X / 2.0)*0.8660254*0.8660254 +
					(y - MODEL_SIZE_Y / 2.0)*(y - MODEL_SIZE_Y / 2.0);
				if (rSquare<MODEL_SIZE_X*MODEL_SIZE_X / pow(2*STEP_SLOPE,2) * i*i&&rSquare>MODEL_SIZE_X*MODEL_SIZE_X / pow(2*STEP_SLOPE, 2) * (i - 1)*(i - 1)) {
					int z = MODEL_SIZE_Z;
					while (!m.getdata(x, y, --z)); heightSum += z+1; cnt++;
				}
			}
		avgHeight[i - 1] = (cnt)?static_cast<double>(heightSum) / cnt:-1;
		fprintf_s(fp, "%f,", avgHeight[i - 1]);
	}
	fprintf_s(fp, ",");
	for (int i = 1; i < 16; i++)fprintf(fp,"%f,",(avgHeight[i-1]-avgHeight[i])/MODEL_SIZE_X*2*STEP_SLOPE);
	fprintf_s(fp, ",");
	int hsum = 0;
	int sum = 0;
	for (int x = 0; x < MODEL_SIZE_X; x++)
		for (int y = 0; y < MODEL_SIZE_Y; y++)
			for (int z = 0; z < MODEL_SIZE_Z; z++) { int e = m.getdata(x, y, z);  hsum += z * e; sum += e; }


	fprintf(fp,"%f", hsum * 0.8164965 * MASS*G / sum);
	fprintf_s(fp, "\n");

}
void exportSlope2D(model& m, FILE* fp, int i) {


	fprintf(fp, "%d,", i);
	//divide model into circles :
	double avgHeight[15];
	for (int i = 1; i < STEP_SLOPE; i++) {
		int heightSum = 0;
		int cnt = 0;

		for (int x = 0; x<MODEL_SIZE_X; x++)
			for (int y = MODEL_SIZE_Y*i/STEP_SLOPE; y < MODEL_SIZE_Y*(i+1)/STEP_SLOPE; y++)
			{
	
				
					int z = MODEL_SIZE_Z;
					while (!m.getdata(x, y, --z)); heightSum += z + 1; cnt++;
				
			}
		avgHeight[i - 1] = (cnt) ? static_cast<double>(heightSum)*0.8164965 / cnt : -1;
		fprintf_s(fp, "%f,", avgHeight[i - 1]);
	}
	fprintf_s(fp, ",");
	for (int i = 1; i < 16; i++)fprintf(fp, "%f,", (avgHeight[i - 1] - avgHeight[i]) / MODEL_SIZE_Y * STEP_SLOPE);
	fprintf_s(fp, ",");
	int hsum = 0;
	int sum = 0;
	for (int x = 0; x < MODEL_SIZE_X; x++)
		for (int y = 0; y < MODEL_SIZE_Y; y++)
			for (int z = 0; z < MODEL_SIZE_Z; z++) { int e = m.getdata(x, y, z);  hsum += z * e; sum += e; }


	fprintf(fp, "%f", hsum * 0.8164965 * MASS*G / sum);
	fprintf_s(fp, "\n");

}
