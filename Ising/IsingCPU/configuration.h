#pragma once
const int MODEL_SIZE = 40;
const int N = 10000000;
const double T = 100;// e/k
const int CALC_TO_PRINT = 10000;
double u = 1; //(dS / dN)(U, V) , u/e 
//-uT = chemical potential
//-ln(n!)-ln((N-n)!)=-nlnn-lnn-(N-n)ln(N-n)-ln(N-n) => - lnn - 1/n + 1/N-n + ln(N-n)
const int coordination = 4;
const double upRatio = 0.2;