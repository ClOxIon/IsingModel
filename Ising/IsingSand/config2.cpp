#include "config2.h"



const int TEST_NUM = 100000000;

const int CALC_TO_PRINT = 10000000;
const int PRINT_NUM = 50;
const double MU = 0.2;//should not be zero : infinite loop
const int N = 10000;
const double G = 10;//gravitational accelation
const double RADIUS = 0.01;//meter
const double MASS = 0.10;//kg per particle
const double ALPHA_K = 0.8037388;//coefficient of the d^2U/dx^2
const double VERT_N_COEF = 0.4082482; //coefficient of the vertical N : 1/root6 when there is only one particle on top and slope is sufficiently low.
const double T = VERT_N_COEF * 1.902852*MU*MASS*G*RADIUS*2;//critical value
const double INIT_K_PROB = 0.1;
const int HIDE_RATIO = 10;
const int STEP_POT = 700;//
const double EN_CONV = 1000;
const int STEP_SLOPE = 20;
const int STEP_PUSH = 200;
const double FLOW_PER_STEP = 0.01;
//������ : �庮�� ���� ��ŭ T�� �ָ� ��� ������ ������. �������� 1ĭ�� �����̴� �ൿ�� ���̸� �� �������� ���.
//T�� �� �ָ� �庮�� �� �Ѿ �� �Ȱ��� �׿�������.