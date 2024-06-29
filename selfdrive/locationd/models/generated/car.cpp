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
void err_fun(double *nom_x, double *delta_x, double *out_227391381542670811) {
   out_227391381542670811[0] = delta_x[0] + nom_x[0];
   out_227391381542670811[1] = delta_x[1] + nom_x[1];
   out_227391381542670811[2] = delta_x[2] + nom_x[2];
   out_227391381542670811[3] = delta_x[3] + nom_x[3];
   out_227391381542670811[4] = delta_x[4] + nom_x[4];
   out_227391381542670811[5] = delta_x[5] + nom_x[5];
   out_227391381542670811[6] = delta_x[6] + nom_x[6];
   out_227391381542670811[7] = delta_x[7] + nom_x[7];
   out_227391381542670811[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4332500849444111319) {
   out_4332500849444111319[0] = -nom_x[0] + true_x[0];
   out_4332500849444111319[1] = -nom_x[1] + true_x[1];
   out_4332500849444111319[2] = -nom_x[2] + true_x[2];
   out_4332500849444111319[3] = -nom_x[3] + true_x[3];
   out_4332500849444111319[4] = -nom_x[4] + true_x[4];
   out_4332500849444111319[5] = -nom_x[5] + true_x[5];
   out_4332500849444111319[6] = -nom_x[6] + true_x[6];
   out_4332500849444111319[7] = -nom_x[7] + true_x[7];
   out_4332500849444111319[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4301467965206663409) {
   out_4301467965206663409[0] = 1.0;
   out_4301467965206663409[1] = 0;
   out_4301467965206663409[2] = 0;
   out_4301467965206663409[3] = 0;
   out_4301467965206663409[4] = 0;
   out_4301467965206663409[5] = 0;
   out_4301467965206663409[6] = 0;
   out_4301467965206663409[7] = 0;
   out_4301467965206663409[8] = 0;
   out_4301467965206663409[9] = 0;
   out_4301467965206663409[10] = 1.0;
   out_4301467965206663409[11] = 0;
   out_4301467965206663409[12] = 0;
   out_4301467965206663409[13] = 0;
   out_4301467965206663409[14] = 0;
   out_4301467965206663409[15] = 0;
   out_4301467965206663409[16] = 0;
   out_4301467965206663409[17] = 0;
   out_4301467965206663409[18] = 0;
   out_4301467965206663409[19] = 0;
   out_4301467965206663409[20] = 1.0;
   out_4301467965206663409[21] = 0;
   out_4301467965206663409[22] = 0;
   out_4301467965206663409[23] = 0;
   out_4301467965206663409[24] = 0;
   out_4301467965206663409[25] = 0;
   out_4301467965206663409[26] = 0;
   out_4301467965206663409[27] = 0;
   out_4301467965206663409[28] = 0;
   out_4301467965206663409[29] = 0;
   out_4301467965206663409[30] = 1.0;
   out_4301467965206663409[31] = 0;
   out_4301467965206663409[32] = 0;
   out_4301467965206663409[33] = 0;
   out_4301467965206663409[34] = 0;
   out_4301467965206663409[35] = 0;
   out_4301467965206663409[36] = 0;
   out_4301467965206663409[37] = 0;
   out_4301467965206663409[38] = 0;
   out_4301467965206663409[39] = 0;
   out_4301467965206663409[40] = 1.0;
   out_4301467965206663409[41] = 0;
   out_4301467965206663409[42] = 0;
   out_4301467965206663409[43] = 0;
   out_4301467965206663409[44] = 0;
   out_4301467965206663409[45] = 0;
   out_4301467965206663409[46] = 0;
   out_4301467965206663409[47] = 0;
   out_4301467965206663409[48] = 0;
   out_4301467965206663409[49] = 0;
   out_4301467965206663409[50] = 1.0;
   out_4301467965206663409[51] = 0;
   out_4301467965206663409[52] = 0;
   out_4301467965206663409[53] = 0;
   out_4301467965206663409[54] = 0;
   out_4301467965206663409[55] = 0;
   out_4301467965206663409[56] = 0;
   out_4301467965206663409[57] = 0;
   out_4301467965206663409[58] = 0;
   out_4301467965206663409[59] = 0;
   out_4301467965206663409[60] = 1.0;
   out_4301467965206663409[61] = 0;
   out_4301467965206663409[62] = 0;
   out_4301467965206663409[63] = 0;
   out_4301467965206663409[64] = 0;
   out_4301467965206663409[65] = 0;
   out_4301467965206663409[66] = 0;
   out_4301467965206663409[67] = 0;
   out_4301467965206663409[68] = 0;
   out_4301467965206663409[69] = 0;
   out_4301467965206663409[70] = 1.0;
   out_4301467965206663409[71] = 0;
   out_4301467965206663409[72] = 0;
   out_4301467965206663409[73] = 0;
   out_4301467965206663409[74] = 0;
   out_4301467965206663409[75] = 0;
   out_4301467965206663409[76] = 0;
   out_4301467965206663409[77] = 0;
   out_4301467965206663409[78] = 0;
   out_4301467965206663409[79] = 0;
   out_4301467965206663409[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6090612947722609695) {
   out_6090612947722609695[0] = state[0];
   out_6090612947722609695[1] = state[1];
   out_6090612947722609695[2] = state[2];
   out_6090612947722609695[3] = state[3];
   out_6090612947722609695[4] = state[4];
   out_6090612947722609695[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6090612947722609695[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6090612947722609695[7] = state[7];
   out_6090612947722609695[8] = state[8];
}
void F_fun(double *state, double dt, double *out_5556904933626449349) {
   out_5556904933626449349[0] = 1;
   out_5556904933626449349[1] = 0;
   out_5556904933626449349[2] = 0;
   out_5556904933626449349[3] = 0;
   out_5556904933626449349[4] = 0;
   out_5556904933626449349[5] = 0;
   out_5556904933626449349[6] = 0;
   out_5556904933626449349[7] = 0;
   out_5556904933626449349[8] = 0;
   out_5556904933626449349[9] = 0;
   out_5556904933626449349[10] = 1;
   out_5556904933626449349[11] = 0;
   out_5556904933626449349[12] = 0;
   out_5556904933626449349[13] = 0;
   out_5556904933626449349[14] = 0;
   out_5556904933626449349[15] = 0;
   out_5556904933626449349[16] = 0;
   out_5556904933626449349[17] = 0;
   out_5556904933626449349[18] = 0;
   out_5556904933626449349[19] = 0;
   out_5556904933626449349[20] = 1;
   out_5556904933626449349[21] = 0;
   out_5556904933626449349[22] = 0;
   out_5556904933626449349[23] = 0;
   out_5556904933626449349[24] = 0;
   out_5556904933626449349[25] = 0;
   out_5556904933626449349[26] = 0;
   out_5556904933626449349[27] = 0;
   out_5556904933626449349[28] = 0;
   out_5556904933626449349[29] = 0;
   out_5556904933626449349[30] = 1;
   out_5556904933626449349[31] = 0;
   out_5556904933626449349[32] = 0;
   out_5556904933626449349[33] = 0;
   out_5556904933626449349[34] = 0;
   out_5556904933626449349[35] = 0;
   out_5556904933626449349[36] = 0;
   out_5556904933626449349[37] = 0;
   out_5556904933626449349[38] = 0;
   out_5556904933626449349[39] = 0;
   out_5556904933626449349[40] = 1;
   out_5556904933626449349[41] = 0;
   out_5556904933626449349[42] = 0;
   out_5556904933626449349[43] = 0;
   out_5556904933626449349[44] = 0;
   out_5556904933626449349[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5556904933626449349[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5556904933626449349[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5556904933626449349[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5556904933626449349[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5556904933626449349[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5556904933626449349[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5556904933626449349[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5556904933626449349[53] = -9.8000000000000007*dt;
   out_5556904933626449349[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5556904933626449349[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5556904933626449349[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5556904933626449349[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5556904933626449349[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5556904933626449349[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5556904933626449349[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5556904933626449349[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5556904933626449349[62] = 0;
   out_5556904933626449349[63] = 0;
   out_5556904933626449349[64] = 0;
   out_5556904933626449349[65] = 0;
   out_5556904933626449349[66] = 0;
   out_5556904933626449349[67] = 0;
   out_5556904933626449349[68] = 0;
   out_5556904933626449349[69] = 0;
   out_5556904933626449349[70] = 1;
   out_5556904933626449349[71] = 0;
   out_5556904933626449349[72] = 0;
   out_5556904933626449349[73] = 0;
   out_5556904933626449349[74] = 0;
   out_5556904933626449349[75] = 0;
   out_5556904933626449349[76] = 0;
   out_5556904933626449349[77] = 0;
   out_5556904933626449349[78] = 0;
   out_5556904933626449349[79] = 0;
   out_5556904933626449349[80] = 1;
}
void h_25(double *state, double *unused, double *out_2902874667874020476) {
   out_2902874667874020476[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4038467864605791628) {
   out_4038467864605791628[0] = 0;
   out_4038467864605791628[1] = 0;
   out_4038467864605791628[2] = 0;
   out_4038467864605791628[3] = 0;
   out_4038467864605791628[4] = 0;
   out_4038467864605791628[5] = 0;
   out_4038467864605791628[6] = 1;
   out_4038467864605791628[7] = 0;
   out_4038467864605791628[8] = 0;
}
void h_24(double *state, double *unused, double *out_7947337562859115535) {
   out_7947337562859115535[0] = state[4];
   out_7947337562859115535[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8784977191572300662) {
   out_8784977191572300662[0] = 0;
   out_8784977191572300662[1] = 0;
   out_8784977191572300662[2] = 0;
   out_8784977191572300662[3] = 0;
   out_8784977191572300662[4] = 1;
   out_8784977191572300662[5] = 0;
   out_8784977191572300662[6] = 0;
   out_8784977191572300662[7] = 0;
   out_8784977191572300662[8] = 0;
   out_8784977191572300662[9] = 0;
   out_8784977191572300662[10] = 0;
   out_8784977191572300662[11] = 0;
   out_8784977191572300662[12] = 0;
   out_8784977191572300662[13] = 0;
   out_8784977191572300662[14] = 1;
   out_8784977191572300662[15] = 0;
   out_8784977191572300662[16] = 0;
   out_8784977191572300662[17] = 0;
}
void h_30(double *state, double *unused, double *out_8581150449551952113) {
   out_8581150449551952113[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2878222476885825127) {
   out_2878222476885825127[0] = 0;
   out_2878222476885825127[1] = 0;
   out_2878222476885825127[2] = 0;
   out_2878222476885825127[3] = 0;
   out_2878222476885825127[4] = 1;
   out_2878222476885825127[5] = 0;
   out_2878222476885825127[6] = 0;
   out_2878222476885825127[7] = 0;
   out_2878222476885825127[8] = 0;
}
void h_26(double *state, double *unused, double *out_1459936804513372850) {
   out_1459936804513372850[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7779971183479847852) {
   out_7779971183479847852[0] = 0;
   out_7779971183479847852[1] = 0;
   out_7779971183479847852[2] = 0;
   out_7779971183479847852[3] = 0;
   out_7779971183479847852[4] = 0;
   out_7779971183479847852[5] = 0;
   out_7779971183479847852[6] = 0;
   out_7779971183479847852[7] = 1;
   out_7779971183479847852[8] = 0;
}
void h_27(double *state, double *unused, double *out_8387958699214852931) {
   out_8387958699214852931[0] = state[3];
}
void H_27(double *state, double *unused, double *out_703459165085400216) {
   out_703459165085400216[0] = 0;
   out_703459165085400216[1] = 0;
   out_703459165085400216[2] = 0;
   out_703459165085400216[3] = 1;
   out_703459165085400216[4] = 0;
   out_703459165085400216[5] = 0;
   out_703459165085400216[6] = 0;
   out_703459165085400216[7] = 0;
   out_703459165085400216[8] = 0;
}
void h_29(double *state, double *unused, double *out_2373992662071331058) {
   out_2373992662071331058[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3388453821200217311) {
   out_3388453821200217311[0] = 0;
   out_3388453821200217311[1] = 1;
   out_3388453821200217311[2] = 0;
   out_3388453821200217311[3] = 0;
   out_3388453821200217311[4] = 0;
   out_3388453821200217311[5] = 0;
   out_3388453821200217311[6] = 0;
   out_3388453821200217311[7] = 0;
   out_3388453821200217311[8] = 0;
}
void h_28(double *state, double *unused, double *out_1653128554850472032) {
   out_1653128554850472032[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6092302578853681391) {
   out_6092302578853681391[0] = 1;
   out_6092302578853681391[1] = 0;
   out_6092302578853681391[2] = 0;
   out_6092302578853681391[3] = 0;
   out_6092302578853681391[4] = 0;
   out_6092302578853681391[5] = 0;
   out_6092302578853681391[6] = 0;
   out_6092302578853681391[7] = 0;
   out_6092302578853681391[8] = 0;
}
void h_31(double *state, double *unused, double *out_6540324196128042826) {
   out_6540324196128042826[0] = state[8];
}
void H_31(double *state, double *unused, double *out_4007821902728831200) {
   out_4007821902728831200[0] = 0;
   out_4007821902728831200[1] = 0;
   out_4007821902728831200[2] = 0;
   out_4007821902728831200[3] = 0;
   out_4007821902728831200[4] = 0;
   out_4007821902728831200[5] = 0;
   out_4007821902728831200[6] = 0;
   out_4007821902728831200[7] = 0;
   out_4007821902728831200[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_227391381542670811) {
  err_fun(nom_x, delta_x, out_227391381542670811);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4332500849444111319) {
  inv_err_fun(nom_x, true_x, out_4332500849444111319);
}
void car_H_mod_fun(double *state, double *out_4301467965206663409) {
  H_mod_fun(state, out_4301467965206663409);
}
void car_f_fun(double *state, double dt, double *out_6090612947722609695) {
  f_fun(state,  dt, out_6090612947722609695);
}
void car_F_fun(double *state, double dt, double *out_5556904933626449349) {
  F_fun(state,  dt, out_5556904933626449349);
}
void car_h_25(double *state, double *unused, double *out_2902874667874020476) {
  h_25(state, unused, out_2902874667874020476);
}
void car_H_25(double *state, double *unused, double *out_4038467864605791628) {
  H_25(state, unused, out_4038467864605791628);
}
void car_h_24(double *state, double *unused, double *out_7947337562859115535) {
  h_24(state, unused, out_7947337562859115535);
}
void car_H_24(double *state, double *unused, double *out_8784977191572300662) {
  H_24(state, unused, out_8784977191572300662);
}
void car_h_30(double *state, double *unused, double *out_8581150449551952113) {
  h_30(state, unused, out_8581150449551952113);
}
void car_H_30(double *state, double *unused, double *out_2878222476885825127) {
  H_30(state, unused, out_2878222476885825127);
}
void car_h_26(double *state, double *unused, double *out_1459936804513372850) {
  h_26(state, unused, out_1459936804513372850);
}
void car_H_26(double *state, double *unused, double *out_7779971183479847852) {
  H_26(state, unused, out_7779971183479847852);
}
void car_h_27(double *state, double *unused, double *out_8387958699214852931) {
  h_27(state, unused, out_8387958699214852931);
}
void car_H_27(double *state, double *unused, double *out_703459165085400216) {
  H_27(state, unused, out_703459165085400216);
}
void car_h_29(double *state, double *unused, double *out_2373992662071331058) {
  h_29(state, unused, out_2373992662071331058);
}
void car_H_29(double *state, double *unused, double *out_3388453821200217311) {
  H_29(state, unused, out_3388453821200217311);
}
void car_h_28(double *state, double *unused, double *out_1653128554850472032) {
  h_28(state, unused, out_1653128554850472032);
}
void car_H_28(double *state, double *unused, double *out_6092302578853681391) {
  H_28(state, unused, out_6092302578853681391);
}
void car_h_31(double *state, double *unused, double *out_6540324196128042826) {
  h_31(state, unused, out_6540324196128042826);
}
void car_H_31(double *state, double *unused, double *out_4007821902728831200) {
  H_31(state, unused, out_4007821902728831200);
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
