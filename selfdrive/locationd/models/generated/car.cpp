#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5629442633965842572) {
   out_5629442633965842572[0] = delta_x[0] + nom_x[0];
   out_5629442633965842572[1] = delta_x[1] + nom_x[1];
   out_5629442633965842572[2] = delta_x[2] + nom_x[2];
   out_5629442633965842572[3] = delta_x[3] + nom_x[3];
   out_5629442633965842572[4] = delta_x[4] + nom_x[4];
   out_5629442633965842572[5] = delta_x[5] + nom_x[5];
   out_5629442633965842572[6] = delta_x[6] + nom_x[6];
   out_5629442633965842572[7] = delta_x[7] + nom_x[7];
   out_5629442633965842572[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3598908807355077436) {
   out_3598908807355077436[0] = -nom_x[0] + true_x[0];
   out_3598908807355077436[1] = -nom_x[1] + true_x[1];
   out_3598908807355077436[2] = -nom_x[2] + true_x[2];
   out_3598908807355077436[3] = -nom_x[3] + true_x[3];
   out_3598908807355077436[4] = -nom_x[4] + true_x[4];
   out_3598908807355077436[5] = -nom_x[5] + true_x[5];
   out_3598908807355077436[6] = -nom_x[6] + true_x[6];
   out_3598908807355077436[7] = -nom_x[7] + true_x[7];
   out_3598908807355077436[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_53256178129770178) {
   out_53256178129770178[0] = 1.0;
   out_53256178129770178[1] = 0;
   out_53256178129770178[2] = 0;
   out_53256178129770178[3] = 0;
   out_53256178129770178[4] = 0;
   out_53256178129770178[5] = 0;
   out_53256178129770178[6] = 0;
   out_53256178129770178[7] = 0;
   out_53256178129770178[8] = 0;
   out_53256178129770178[9] = 0;
   out_53256178129770178[10] = 1.0;
   out_53256178129770178[11] = 0;
   out_53256178129770178[12] = 0;
   out_53256178129770178[13] = 0;
   out_53256178129770178[14] = 0;
   out_53256178129770178[15] = 0;
   out_53256178129770178[16] = 0;
   out_53256178129770178[17] = 0;
   out_53256178129770178[18] = 0;
   out_53256178129770178[19] = 0;
   out_53256178129770178[20] = 1.0;
   out_53256178129770178[21] = 0;
   out_53256178129770178[22] = 0;
   out_53256178129770178[23] = 0;
   out_53256178129770178[24] = 0;
   out_53256178129770178[25] = 0;
   out_53256178129770178[26] = 0;
   out_53256178129770178[27] = 0;
   out_53256178129770178[28] = 0;
   out_53256178129770178[29] = 0;
   out_53256178129770178[30] = 1.0;
   out_53256178129770178[31] = 0;
   out_53256178129770178[32] = 0;
   out_53256178129770178[33] = 0;
   out_53256178129770178[34] = 0;
   out_53256178129770178[35] = 0;
   out_53256178129770178[36] = 0;
   out_53256178129770178[37] = 0;
   out_53256178129770178[38] = 0;
   out_53256178129770178[39] = 0;
   out_53256178129770178[40] = 1.0;
   out_53256178129770178[41] = 0;
   out_53256178129770178[42] = 0;
   out_53256178129770178[43] = 0;
   out_53256178129770178[44] = 0;
   out_53256178129770178[45] = 0;
   out_53256178129770178[46] = 0;
   out_53256178129770178[47] = 0;
   out_53256178129770178[48] = 0;
   out_53256178129770178[49] = 0;
   out_53256178129770178[50] = 1.0;
   out_53256178129770178[51] = 0;
   out_53256178129770178[52] = 0;
   out_53256178129770178[53] = 0;
   out_53256178129770178[54] = 0;
   out_53256178129770178[55] = 0;
   out_53256178129770178[56] = 0;
   out_53256178129770178[57] = 0;
   out_53256178129770178[58] = 0;
   out_53256178129770178[59] = 0;
   out_53256178129770178[60] = 1.0;
   out_53256178129770178[61] = 0;
   out_53256178129770178[62] = 0;
   out_53256178129770178[63] = 0;
   out_53256178129770178[64] = 0;
   out_53256178129770178[65] = 0;
   out_53256178129770178[66] = 0;
   out_53256178129770178[67] = 0;
   out_53256178129770178[68] = 0;
   out_53256178129770178[69] = 0;
   out_53256178129770178[70] = 1.0;
   out_53256178129770178[71] = 0;
   out_53256178129770178[72] = 0;
   out_53256178129770178[73] = 0;
   out_53256178129770178[74] = 0;
   out_53256178129770178[75] = 0;
   out_53256178129770178[76] = 0;
   out_53256178129770178[77] = 0;
   out_53256178129770178[78] = 0;
   out_53256178129770178[79] = 0;
   out_53256178129770178[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2357487969533120980) {
   out_2357487969533120980[0] = state[0];
   out_2357487969533120980[1] = state[1];
   out_2357487969533120980[2] = state[2];
   out_2357487969533120980[3] = state[3];
   out_2357487969533120980[4] = state[4];
   out_2357487969533120980[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2357487969533120980[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2357487969533120980[7] = state[7];
   out_2357487969533120980[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8557443704727234356) {
   out_8557443704727234356[0] = 1;
   out_8557443704727234356[1] = 0;
   out_8557443704727234356[2] = 0;
   out_8557443704727234356[3] = 0;
   out_8557443704727234356[4] = 0;
   out_8557443704727234356[5] = 0;
   out_8557443704727234356[6] = 0;
   out_8557443704727234356[7] = 0;
   out_8557443704727234356[8] = 0;
   out_8557443704727234356[9] = 0;
   out_8557443704727234356[10] = 1;
   out_8557443704727234356[11] = 0;
   out_8557443704727234356[12] = 0;
   out_8557443704727234356[13] = 0;
   out_8557443704727234356[14] = 0;
   out_8557443704727234356[15] = 0;
   out_8557443704727234356[16] = 0;
   out_8557443704727234356[17] = 0;
   out_8557443704727234356[18] = 0;
   out_8557443704727234356[19] = 0;
   out_8557443704727234356[20] = 1;
   out_8557443704727234356[21] = 0;
   out_8557443704727234356[22] = 0;
   out_8557443704727234356[23] = 0;
   out_8557443704727234356[24] = 0;
   out_8557443704727234356[25] = 0;
   out_8557443704727234356[26] = 0;
   out_8557443704727234356[27] = 0;
   out_8557443704727234356[28] = 0;
   out_8557443704727234356[29] = 0;
   out_8557443704727234356[30] = 1;
   out_8557443704727234356[31] = 0;
   out_8557443704727234356[32] = 0;
   out_8557443704727234356[33] = 0;
   out_8557443704727234356[34] = 0;
   out_8557443704727234356[35] = 0;
   out_8557443704727234356[36] = 0;
   out_8557443704727234356[37] = 0;
   out_8557443704727234356[38] = 0;
   out_8557443704727234356[39] = 0;
   out_8557443704727234356[40] = 1;
   out_8557443704727234356[41] = 0;
   out_8557443704727234356[42] = 0;
   out_8557443704727234356[43] = 0;
   out_8557443704727234356[44] = 0;
   out_8557443704727234356[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8557443704727234356[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8557443704727234356[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8557443704727234356[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8557443704727234356[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8557443704727234356[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8557443704727234356[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8557443704727234356[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8557443704727234356[53] = -9.8000000000000007*dt;
   out_8557443704727234356[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8557443704727234356[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8557443704727234356[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8557443704727234356[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8557443704727234356[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8557443704727234356[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8557443704727234356[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8557443704727234356[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8557443704727234356[62] = 0;
   out_8557443704727234356[63] = 0;
   out_8557443704727234356[64] = 0;
   out_8557443704727234356[65] = 0;
   out_8557443704727234356[66] = 0;
   out_8557443704727234356[67] = 0;
   out_8557443704727234356[68] = 0;
   out_8557443704727234356[69] = 0;
   out_8557443704727234356[70] = 1;
   out_8557443704727234356[71] = 0;
   out_8557443704727234356[72] = 0;
   out_8557443704727234356[73] = 0;
   out_8557443704727234356[74] = 0;
   out_8557443704727234356[75] = 0;
   out_8557443704727234356[76] = 0;
   out_8557443704727234356[77] = 0;
   out_8557443704727234356[78] = 0;
   out_8557443704727234356[79] = 0;
   out_8557443704727234356[80] = 1;
}
void h_25(double *state, double *unused, double *out_1115798248010839447) {
   out_1115798248010839447[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5831744492098185822) {
   out_5831744492098185822[0] = 0;
   out_5831744492098185822[1] = 0;
   out_5831744492098185822[2] = 0;
   out_5831744492098185822[3] = 0;
   out_5831744492098185822[4] = 0;
   out_5831744492098185822[5] = 0;
   out_5831744492098185822[6] = 1;
   out_5831744492098185822[7] = 0;
   out_5831744492098185822[8] = 0;
}
void h_24(double *state, double *unused, double *out_3529007286521349464) {
   out_3529007286521349464[0] = state[4];
   out_3529007286521349464[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3659094893092686256) {
   out_3659094893092686256[0] = 0;
   out_3659094893092686256[1] = 0;
   out_3659094893092686256[2] = 0;
   out_3659094893092686256[3] = 0;
   out_3659094893092686256[4] = 1;
   out_3659094893092686256[5] = 0;
   out_3659094893092686256[6] = 0;
   out_3659094893092686256[7] = 0;
   out_3659094893092686256[8] = 0;
   out_3659094893092686256[9] = 0;
   out_3659094893092686256[10] = 0;
   out_3659094893092686256[11] = 0;
   out_3659094893092686256[12] = 0;
   out_3659094893092686256[13] = 0;
   out_3659094893092686256[14] = 1;
   out_3659094893092686256[15] = 0;
   out_3659094893092686256[16] = 0;
   out_3659094893092686256[17] = 0;
}
void h_30(double *state, double *unused, double *out_3557753197258034570) {
   out_3557753197258034570[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1304048161970577624) {
   out_1304048161970577624[0] = 0;
   out_1304048161970577624[1] = 0;
   out_1304048161970577624[2] = 0;
   out_1304048161970577624[3] = 0;
   out_1304048161970577624[4] = 1;
   out_1304048161970577624[5] = 0;
   out_1304048161970577624[6] = 0;
   out_1304048161970577624[7] = 0;
   out_1304048161970577624[8] = 0;
}
void h_26(double *state, double *unused, double *out_4227657190368867009) {
   out_4227657190368867009[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2090241173224129598) {
   out_2090241173224129598[0] = 0;
   out_2090241173224129598[1] = 0;
   out_2090241173224129598[2] = 0;
   out_2090241173224129598[3] = 0;
   out_2090241173224129598[4] = 0;
   out_2090241173224129598[5] = 0;
   out_2090241173224129598[6] = 0;
   out_2090241173224129598[7] = 1;
   out_2090241173224129598[8] = 0;
}
void h_27(double *state, double *unused, double *out_2700364704332613072) {
   out_2700364704332613072[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3527642233154520841) {
   out_3527642233154520841[0] = 0;
   out_3527642233154520841[1] = 0;
   out_3527642233154520841[2] = 0;
   out_3527642233154520841[3] = 1;
   out_3527642233154520841[4] = 0;
   out_3527642233154520841[5] = 0;
   out_3527642233154520841[6] = 0;
   out_3527642233154520841[7] = 0;
   out_3527642233154520841[8] = 0;
}
void h_29(double *state, double *unused, double *out_4161069081934512087) {
   out_4161069081934512087[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1814279506284969808) {
   out_1814279506284969808[0] = 0;
   out_1814279506284969808[1] = 1;
   out_1814279506284969808[2] = 0;
   out_1814279506284969808[3] = 0;
   out_1814279506284969808[4] = 0;
   out_1814279506284969808[5] = 0;
   out_1814279506284969808[6] = 0;
   out_1814279506284969808[7] = 0;
   out_1814279506284969808[8] = 0;
}
void h_28(double *state, double *unused, double *out_8458349618326351028) {
   out_8458349618326351028[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3777909777850296059) {
   out_3777909777850296059[0] = 1;
   out_3777909777850296059[1] = 0;
   out_3777909777850296059[2] = 0;
   out_3777909777850296059[3] = 0;
   out_3777909777850296059[4] = 0;
   out_3777909777850296059[5] = 0;
   out_3777909777850296059[6] = 0;
   out_3777909777850296059[7] = 0;
   out_3777909777850296059[8] = 0;
}
void h_31(double *state, double *unused, double *out_354890393280493669) {
   out_354890393280493669[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5862390453975146250) {
   out_5862390453975146250[0] = 0;
   out_5862390453975146250[1] = 0;
   out_5862390453975146250[2] = 0;
   out_5862390453975146250[3] = 0;
   out_5862390453975146250[4] = 0;
   out_5862390453975146250[5] = 0;
   out_5862390453975146250[6] = 0;
   out_5862390453975146250[7] = 0;
   out_5862390453975146250[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5629442633965842572) {
  err_fun(nom_x, delta_x, out_5629442633965842572);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3598908807355077436) {
  inv_err_fun(nom_x, true_x, out_3598908807355077436);
}
void car_H_mod_fun(double *state, double *out_53256178129770178) {
  H_mod_fun(state, out_53256178129770178);
}
void car_f_fun(double *state, double dt, double *out_2357487969533120980) {
  f_fun(state,  dt, out_2357487969533120980);
}
void car_F_fun(double *state, double dt, double *out_8557443704727234356) {
  F_fun(state,  dt, out_8557443704727234356);
}
void car_h_25(double *state, double *unused, double *out_1115798248010839447) {
  h_25(state, unused, out_1115798248010839447);
}
void car_H_25(double *state, double *unused, double *out_5831744492098185822) {
  H_25(state, unused, out_5831744492098185822);
}
void car_h_24(double *state, double *unused, double *out_3529007286521349464) {
  h_24(state, unused, out_3529007286521349464);
}
void car_H_24(double *state, double *unused, double *out_3659094893092686256) {
  H_24(state, unused, out_3659094893092686256);
}
void car_h_30(double *state, double *unused, double *out_3557753197258034570) {
  h_30(state, unused, out_3557753197258034570);
}
void car_H_30(double *state, double *unused, double *out_1304048161970577624) {
  H_30(state, unused, out_1304048161970577624);
}
void car_h_26(double *state, double *unused, double *out_4227657190368867009) {
  h_26(state, unused, out_4227657190368867009);
}
void car_H_26(double *state, double *unused, double *out_2090241173224129598) {
  H_26(state, unused, out_2090241173224129598);
}
void car_h_27(double *state, double *unused, double *out_2700364704332613072) {
  h_27(state, unused, out_2700364704332613072);
}
void car_H_27(double *state, double *unused, double *out_3527642233154520841) {
  H_27(state, unused, out_3527642233154520841);
}
void car_h_29(double *state, double *unused, double *out_4161069081934512087) {
  h_29(state, unused, out_4161069081934512087);
}
void car_H_29(double *state, double *unused, double *out_1814279506284969808) {
  H_29(state, unused, out_1814279506284969808);
}
void car_h_28(double *state, double *unused, double *out_8458349618326351028) {
  h_28(state, unused, out_8458349618326351028);
}
void car_H_28(double *state, double *unused, double *out_3777909777850296059) {
  H_28(state, unused, out_3777909777850296059);
}
void car_h_31(double *state, double *unused, double *out_354890393280493669) {
  h_31(state, unused, out_354890393280493669);
}
void car_H_31(double *state, double *unused, double *out_5862390453975146250) {
  H_31(state, unused, out_5862390453975146250);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
