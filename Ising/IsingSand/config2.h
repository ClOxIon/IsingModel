#pragma once
#define MODEL_SIZE_X 100
#define MODEL_SIZE_Y 100
#define MODEL_SIZE_Z 100
extern const int TEST_NUM;
extern const int PRINT_NUM;
extern const double T;// thermal kinetic energy
extern const int CALC_TO_PRINT;
extern const double MU; //coef of friction * average load * unit length
extern const double G;//gravitational accel. * mass of particle * unit length
extern const int N;
extern const double RADIUS;//meter
extern const double MASS; //kg per particle
extern const double ALPHA_K;//coefficient of the d^2U/dx^2
extern const double VERT_N_COEF;//coefficient of the vertical N : 1/root6 when there is only one particle on top and slope is sufficiently low.
extern const double INIT_K_PROB;//probability of initial kinetic energy
extern const int HIDE_RATIO;
extern const int STEP_POT;
extern const double EN_CONV;
extern const int STEP_SLOPE;
extern const int STEP_PUSH;
extern const double FLOW_PER_STEP;