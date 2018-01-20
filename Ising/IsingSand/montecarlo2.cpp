//this file implements path-probability metropolis method#include <random>
// control k d
/*
frictional energy lost by
vertical disturbing particle when moving horizental =
frictional energy lost by 
vertical disturbing particle when moving semi-vertical = (when particle is in semi-stable state, it is not fully supported, so N remains same 1/root6 MG)
2 * MU * MASS * G * RADIUS * (1/root6)(=VERT_N_COEF) * root3 * ln3/2 

frictional energy lost by
horizental disturbing particle when moving horizental = 
2 * MU * MASS * G * RADIUS * 2root3 * (ln3-1)/2 * [(3/2)^(3/2)/2 or (2root2/3)^3/root3](average = ALPHA_K)

*/

#include "montecarlo.h"
#include <random>
double getPotential(model& m) {
	int hsum = 0;
	int sum = 0;
	for (int x = 0; x < MODEL_SIZE_X; x++)
		for (int y = 0; y < MODEL_SIZE_Y; y++)
			for (int z = 0; z < MODEL_SIZE_Z; z++) { int e = m.getdata(x, y, z);  hsum += z * e ; sum += e; }
		

	return hsum * 0.8164965 * MASS*G/sum;

}
void monteCarloSimulate(model& m) {


	std::random_device rd; //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> disx(0, MODEL_SIZE_X - 1);//360448 = 32768 * 11
	std::uniform_int_distribution<> disy(0, MODEL_SIZE_Y - 1);//360448 = 32768 * 11
	std::uniform_int_distribution<> disz(0, MODEL_SIZE_Z - 1);//360448 = 32768 * 11
	std::uniform_real_distribution<double> disR(0, 1);
	
	const int DIS_SIZE = 360447;
	double potInit = getPotential(m);
	double potFnl = 0;
	double k_prob = 0.01;
	int particleCnt = 0;
	for (int i = 0; ; i++) {

		particleCnt++;
	a:
		unsigned int x = static_cast<int>(disx(gen));
		unsigned int y = static_cast<int>(disy(gen));
		unsigned int z = static_cast<int>(disz(gen));
		__int8& e = m.getdata(x, y, z);//when a particle gets into pit and gone, it;s lost. forever.
		//while(!e&&z!=0)z--;
		double K = (disR(gen) < k_prob) ? T : 0; //initial kinetic energy
	b:
		
		
		if (!e || z < 0 || z >= MODEL_SIZE_Z || y < 0 || y >= MODEL_SIZE_Y || x < 0 || x >= MODEL_SIZE_X)goto a;
		//when x,y,z is out of scope, reroll.
		

		__int8& t = m.getdata(x, y, z + 1);//change these into integer
		__int8& b = m.getdata(x, y, z - 1);

		__int8& tu = m.getdata(x + z % 2 * 2 - 1, y + x % 2, z + 1);
		__int8& td = m.getdata(x + z % 2 * 2 - 1, y + x % 2 - 1, z + 1);
		__int8& bu = m.getdata(x + z % 2 * 2 - 1, y + x % 2, z - 1);
		__int8& bd = m.getdata(x + z % 2 * 2 - 1, y + x % 2 - 1, z - 1);

		__int8& u = m.getdata(x, y + 1, z);
		__int8& d = m.getdata(x, y - 1, z);

		__int8&	lu = m.getdata(x - 1, y + x % 2, z);//(x % 2) ? getdata(x - 1, y + 1, z, m) : getdata(x - 1, y, z, m); 
		__int8&	ru = m.getdata(x + 1, y + x % 2, z);//(x % 2) ? getdata(x + 1, y + 1, z, m) : getdata(x + 1, y, z, m);
		__int8& ld = m.getdata(x - 1, y + x % 2 - 1, z);//(x % 2) ? getdata(x - 1, y, z, m) : getdata(x - 1, y - 1, z, m);
		__int8&	rd = m.getdata(x + 1, y + x % 2 - 1, z);//(x % 2) ? getdata(x + 1, y, z, m) : getdata(x + 1, y - 1, z, m);

		
		if (t + b + tu + td + bu + bd + u + d + lu + ru + ld + rd == 12)goto a;
		/*
		double pt = (T - g) / (u + d + lu + ru + ld + rd) * 2 / mu;
		double pb = (T + g) / (u + d + lu + ru + ld + rd) * 2 / mu;
		double pu = (T) / (t/2 + b/2 + lu + ru) / mu;
		double pd = (T) / (t/2 + b/2 + ld + rd) / mu;
		double plu = (T) / (u + t/2 + b/2 + ld) / mu;
		double pru = (T) / (u + t/2 + b/2 + rd) / mu;
		double pld = (T) / (d + t/2 + b/2 + lu) / mu;
		double prd = (T) / (d + t/2 + b/2 + ru) / mu;
		*/
		//stablization is still needed
		//maybe we should only give K horizental direction.
		// bool addition is safe ; it is converted to int.
		//
		if (i == CALC_TO_PRINT)return;
		
		if (i!=0&&i%STEP_POT == 0) {
			potFnl = getPotential(m);
			k_prob = (particleCnt)?(potInit - potFnl)/particleCnt/T*EN_CONV:0;
			potInit = potFnl;
			particleCnt = 0;
		}
		
		int count1 = 3-(b + bd + bu);
		if (count1) {
			if (!b&&1.0 / count1-- > disR(gen)) {//disturbing particle:  bd, bu,(z%2)?(lu, ld):(ru,rd)
				double dE = (0.8164965 - VERT_N_COEF * MU*1.902852*(bd + bu+((z % 2 )? lu + ld : ru + rd)))*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e=0; b=1;
					if (K > 0) { i++; z--; goto b; }
					
				}
				continue;
			}
			if (!bd&&1.0 / count1-- > disR(gen)) {
				//disturbing particle:  b, bu,d,(z%2)?rd:ld}
				double dE = (0.8164965 - VERT_N_COEF * MU*1.902852*(b + bu + d + ((z % 2) ?rd :ld)))*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; bd= 1;
					if (K > 0) { i++; y += x % 2 - 1; x += z % 2 * 2 - 1; z--;  goto b; }//WARN : ORDER IS IMPORTANT

				}
				continue;
			}
			if (!bu&&1.0 / count1-- > disR(gen)) {
				//disturbing particle:  bd, b,u,(z%2)?ru:lu}
				double dE = (0.8164965 - VERT_N_COEF * MU*1.902852*(bd + b + u + ((z % 2) ? ru : lu)))*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; bu= 1;
					if (K > 0) { i++; y += x % 2 ; x += z % 2 * 2 - 1; z--;  goto b; }//WARN : ORDER IS IMPORTANT 

				}
				continue;
			}
		}
		
		//-----------------------------------------------------------------------------------------------------------------------------------
		int count2 = 9-(t + tu + td + u + d + lu + ld + ru + rd);
		if (count2) {
			if (!t&&1.0 / count2-- > disR(gen)) {//disturbing particle:  td, tu,(z%2)?(lu, ld):(ru,rd)
				double dE = (- 0.8164965 - VERT_N_COEF * MU*1.902852*(td + tu + ((z % 2) ? lu + ld : ru + rd)))*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; t= 1;
					if (K > 0) { i++; z--;  goto b; }

				}
				continue;
			}
			if (!td&&1.0 / count2-- > disR(gen)) {
				//disturbing particle:  t, tu,d,(z%2)?rd:ld}
				double dE = (- 0.8164965 - VERT_N_COEF * MU*1.902852*(t + tu + d + ((z % 2) ? rd : ld)))*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; td= 1;
					if (K > 0) { i++; y += x % 2 - 1; x += z % 2 * 2 - 1; z++;  goto b; }//WARN : ORDER IS IMPORTANT

				}
				continue;
			}
			if (!tu&&1.0 / count2-- > disR(gen)) {
				//disturbing particle:  td, t,u,(z%2)?ru:lu}
				double dE = (- 0.8164965 - VERT_N_COEF * MU*1.902852*(td + t + u + ((z % 2) ? ru : lu)))*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; tu= 1;
					if (K > 0) { i++; y += x % 2; x += z % 2 * 2 - 1; z++;  goto b; }//WARN : ORDER IS IMPORTANT 

				}
				continue;
			}
			//-----------------------------------------------------------------------------------------------------------------
			if (!u&&1.0 / count2-- > disR(gen)) {//disturbing particle:  bu, tu, lu, ru
				double dE = (- VERT_N_COEF *1.902852*(bu + tu) - ALPHA_K* 0.3416030*(lu+ru))*MU*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; u= 1;
					if (K > 0) { i++; y++;  goto b; }

				}
				continue;
			}
			if (!d&&1.0 / count2-- > disR(gen)) {//disturbing particle:  bd, td, ld, rd
				double dE = (-VERT_N_COEF * 1.902852*(bd + td) - ALPHA_K * 0.3416030*(ld + rd))*MU*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; d= 1;
					if (K > 0) { i++; y--;  goto b; }

				}
				continue;
			}
			if (!ld&&1.0 / count2-- > disR(gen)) {
				//disturbing particle:  lu, d,(z%2)?b+t:bd+td}
				double dE = (-VERT_N_COEF * 1.902852*(z%2?b+t:bd+td) - ALPHA_K * 0.3416030*(lu + d))*MU*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; ld= 1;
					if (K > 0) { i++; y += x % 2 - 1; x--;  goto b; }//WARN : ORDER IS IMPORTANT

				}
				continue;
			}
			if (!lu&&1.0 / count2-- > disR(gen)) {
				//disturbing particle:  ld, u,(z%2)?b+t:bu+tu}
				double dE = (-VERT_N_COEF * 1.902852*(z % 2 ? b + t : bu + tu) - ALPHA_K * 0.3416030*(ld + u))*MU*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; lu= 1;
					if (K > 0) { i++; y += x % 2; x--;  goto b; }//WARN : ORDER IS IMPORTANT 

				}
				continue;
			}
			if (!rd&&1.0 / count2-- > disR(gen)) {
				//disturbing particle:  ru, d,(z%2)?bd+td:b+t}
				double dE = (-VERT_N_COEF * 1.902852*(z % 2 ? bd + td : b + t) - ALPHA_K * 0.3416030*(ru + d))*MU*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; rd= 1;
					if (K > 0) { i++; y += x % 2 - 1; x++;  goto b; }//WARN : ORDER IS IMPORTANT

				}
				continue;
			}
			if (!ru&&1.0 / count2-- > disR(gen)) {
				//disturbing particle:  rd, u,(z%2)?bu+tu:b+t}
				double dE = (-VERT_N_COEF * 1.902852*(z % 2 ? bu + tu : b + t) - ALPHA_K * 0.3416030*(rd + u))*MU*MASS*G*RADIUS;
				if (K + dE / 2 > 0) {
					K += dE;
					e= 0; ru= 1;
					if (K > 0) { i++; y += x % 2; x++;  goto b; }//WARN : ORDER IS IMPORTANT 

				}
				continue;
			}
		}
		



		





	}

}
