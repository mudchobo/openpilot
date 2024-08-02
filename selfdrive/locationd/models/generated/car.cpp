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
void err_fun(double *nom_x, double *delta_x, double *out_7163063568989676842) {
   out_7163063568989676842[0] = delta_x[0] + nom_x[0];
   out_7163063568989676842[1] = delta_x[1] + nom_x[1];
   out_7163063568989676842[2] = delta_x[2] + nom_x[2];
   out_7163063568989676842[3] = delta_x[3] + nom_x[3];
   out_7163063568989676842[4] = delta_x[4] + nom_x[4];
   out_7163063568989676842[5] = delta_x[5] + nom_x[5];
   out_7163063568989676842[6] = delta_x[6] + nom_x[6];
   out_7163063568989676842[7] = delta_x[7] + nom_x[7];
   out_7163063568989676842[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3754260694415491349) {
   out_3754260694415491349[0] = -nom_x[0] + true_x[0];
   out_3754260694415491349[1] = -nom_x[1] + true_x[1];
   out_3754260694415491349[2] = -nom_x[2] + true_x[2];
   out_3754260694415491349[3] = -nom_x[3] + true_x[3];
   out_3754260694415491349[4] = -nom_x[4] + true_x[4];
   out_3754260694415491349[5] = -nom_x[5] + true_x[5];
   out_3754260694415491349[6] = -nom_x[6] + true_x[6];
   out_3754260694415491349[7] = -nom_x[7] + true_x[7];
   out_3754260694415491349[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8878594360494595028) {
   out_8878594360494595028[0] = 1.0;
   out_8878594360494595028[1] = 0;
   out_8878594360494595028[2] = 0;
   out_8878594360494595028[3] = 0;
   out_8878594360494595028[4] = 0;
   out_8878594360494595028[5] = 0;
   out_8878594360494595028[6] = 0;
   out_8878594360494595028[7] = 0;
   out_8878594360494595028[8] = 0;
   out_8878594360494595028[9] = 0;
   out_8878594360494595028[10] = 1.0;
   out_8878594360494595028[11] = 0;
   out_8878594360494595028[12] = 0;
   out_8878594360494595028[13] = 0;
   out_8878594360494595028[14] = 0;
   out_8878594360494595028[15] = 0;
   out_8878594360494595028[16] = 0;
   out_8878594360494595028[17] = 0;
   out_8878594360494595028[18] = 0;
   out_8878594360494595028[19] = 0;
   out_8878594360494595028[20] = 1.0;
   out_8878594360494595028[21] = 0;
   out_8878594360494595028[22] = 0;
   out_8878594360494595028[23] = 0;
   out_8878594360494595028[24] = 0;
   out_8878594360494595028[25] = 0;
   out_8878594360494595028[26] = 0;
   out_8878594360494595028[27] = 0;
   out_8878594360494595028[28] = 0;
   out_8878594360494595028[29] = 0;
   out_8878594360494595028[30] = 1.0;
   out_8878594360494595028[31] = 0;
   out_8878594360494595028[32] = 0;
   out_8878594360494595028[33] = 0;
   out_8878594360494595028[34] = 0;
   out_8878594360494595028[35] = 0;
   out_8878594360494595028[36] = 0;
   out_8878594360494595028[37] = 0;
   out_8878594360494595028[38] = 0;
   out_8878594360494595028[39] = 0;
   out_8878594360494595028[40] = 1.0;
   out_8878594360494595028[41] = 0;
   out_8878594360494595028[42] = 0;
   out_8878594360494595028[43] = 0;
   out_8878594360494595028[44] = 0;
   out_8878594360494595028[45] = 0;
   out_8878594360494595028[46] = 0;
   out_8878594360494595028[47] = 0;
   out_8878594360494595028[48] = 0;
   out_8878594360494595028[49] = 0;
   out_8878594360494595028[50] = 1.0;
   out_8878594360494595028[51] = 0;
   out_8878594360494595028[52] = 0;
   out_8878594360494595028[53] = 0;
   out_8878594360494595028[54] = 0;
   out_8878594360494595028[55] = 0;
   out_8878594360494595028[56] = 0;
   out_8878594360494595028[57] = 0;
   out_8878594360494595028[58] = 0;
   out_8878594360494595028[59] = 0;
   out_8878594360494595028[60] = 1.0;
   out_8878594360494595028[61] = 0;
   out_8878594360494595028[62] = 0;
   out_8878594360494595028[63] = 0;
   out_8878594360494595028[64] = 0;
   out_8878594360494595028[65] = 0;
   out_8878594360494595028[66] = 0;
   out_8878594360494595028[67] = 0;
   out_8878594360494595028[68] = 0;
   out_8878594360494595028[69] = 0;
   out_8878594360494595028[70] = 1.0;
   out_8878594360494595028[71] = 0;
   out_8878594360494595028[72] = 0;
   out_8878594360494595028[73] = 0;
   out_8878594360494595028[74] = 0;
   out_8878594360494595028[75] = 0;
   out_8878594360494595028[76] = 0;
   out_8878594360494595028[77] = 0;
   out_8878594360494595028[78] = 0;
   out_8878594360494595028[79] = 0;
   out_8878594360494595028[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_32589037743347182) {
   out_32589037743347182[0] = state[0];
   out_32589037743347182[1] = state[1];
   out_32589037743347182[2] = state[2];
   out_32589037743347182[3] = state[3];
   out_32589037743347182[4] = state[4];
   out_32589037743347182[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_32589037743347182[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_32589037743347182[7] = state[7];
   out_32589037743347182[8] = state[8];
}
void F_fun(double *state, double dt, double *out_2355328765180068083) {
   out_2355328765180068083[0] = 1;
   out_2355328765180068083[1] = 0;
   out_2355328765180068083[2] = 0;
   out_2355328765180068083[3] = 0;
   out_2355328765180068083[4] = 0;
   out_2355328765180068083[5] = 0;
   out_2355328765180068083[6] = 0;
   out_2355328765180068083[7] = 0;
   out_2355328765180068083[8] = 0;
   out_2355328765180068083[9] = 0;
   out_2355328765180068083[10] = 1;
   out_2355328765180068083[11] = 0;
   out_2355328765180068083[12] = 0;
   out_2355328765180068083[13] = 0;
   out_2355328765180068083[14] = 0;
   out_2355328765180068083[15] = 0;
   out_2355328765180068083[16] = 0;
   out_2355328765180068083[17] = 0;
   out_2355328765180068083[18] = 0;
   out_2355328765180068083[19] = 0;
   out_2355328765180068083[20] = 1;
   out_2355328765180068083[21] = 0;
   out_2355328765180068083[22] = 0;
   out_2355328765180068083[23] = 0;
   out_2355328765180068083[24] = 0;
   out_2355328765180068083[25] = 0;
   out_2355328765180068083[26] = 0;
   out_2355328765180068083[27] = 0;
   out_2355328765180068083[28] = 0;
   out_2355328765180068083[29] = 0;
   out_2355328765180068083[30] = 1;
   out_2355328765180068083[31] = 0;
   out_2355328765180068083[32] = 0;
   out_2355328765180068083[33] = 0;
   out_2355328765180068083[34] = 0;
   out_2355328765180068083[35] = 0;
   out_2355328765180068083[36] = 0;
   out_2355328765180068083[37] = 0;
   out_2355328765180068083[38] = 0;
   out_2355328765180068083[39] = 0;
   out_2355328765180068083[40] = 1;
   out_2355328765180068083[41] = 0;
   out_2355328765180068083[42] = 0;
   out_2355328765180068083[43] = 0;
   out_2355328765180068083[44] = 0;
   out_2355328765180068083[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2355328765180068083[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2355328765180068083[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2355328765180068083[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2355328765180068083[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2355328765180068083[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2355328765180068083[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2355328765180068083[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2355328765180068083[53] = -9.8000000000000007*dt;
   out_2355328765180068083[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2355328765180068083[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2355328765180068083[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2355328765180068083[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2355328765180068083[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2355328765180068083[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2355328765180068083[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2355328765180068083[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2355328765180068083[62] = 0;
   out_2355328765180068083[63] = 0;
   out_2355328765180068083[64] = 0;
   out_2355328765180068083[65] = 0;
   out_2355328765180068083[66] = 0;
   out_2355328765180068083[67] = 0;
   out_2355328765180068083[68] = 0;
   out_2355328765180068083[69] = 0;
   out_2355328765180068083[70] = 1;
   out_2355328765180068083[71] = 0;
   out_2355328765180068083[72] = 0;
   out_2355328765180068083[73] = 0;
   out_2355328765180068083[74] = 0;
   out_2355328765180068083[75] = 0;
   out_2355328765180068083[76] = 0;
   out_2355328765180068083[77] = 0;
   out_2355328765180068083[78] = 0;
   out_2355328765180068083[79] = 0;
   out_2355328765180068083[80] = 1;
}
void h_25(double *state, double *unused, double *out_1013193562225998629) {
   out_1013193562225998629[0] = state[6];
}
void H_25(double *state, double *unused, double *out_993794256298495203) {
   out_993794256298495203[0] = 0;
   out_993794256298495203[1] = 0;
   out_993794256298495203[2] = 0;
   out_993794256298495203[3] = 0;
   out_993794256298495203[4] = 0;
   out_993794256298495203[5] = 0;
   out_993794256298495203[6] = 1;
   out_993794256298495203[7] = 0;
   out_993794256298495203[8] = 0;
}
void h_24(double *state, double *unused, double *out_4596273345226045289) {
   out_4596273345226045289[0] = state[4];
   out_4596273345226045289[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4551719855100818339) {
   out_4551719855100818339[0] = 0;
   out_4551719855100818339[1] = 0;
   out_4551719855100818339[2] = 0;
   out_4551719855100818339[3] = 0;
   out_4551719855100818339[4] = 1;
   out_4551719855100818339[5] = 0;
   out_4551719855100818339[6] = 0;
   out_4551719855100818339[7] = 0;
   out_4551719855100818339[8] = 0;
   out_4551719855100818339[9] = 0;
   out_4551719855100818339[10] = 0;
   out_4551719855100818339[11] = 0;
   out_4551719855100818339[12] = 0;
   out_4551719855100818339[13] = 0;
   out_4551719855100818339[14] = 1;
   out_4551719855100818339[15] = 0;
   out_4551719855100818339[16] = 0;
   out_4551719855100818339[17] = 0;
}
void h_30(double *state, double *unused, double *out_7601021488702293792) {
   out_7601021488702293792[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3512127214805743830) {
   out_3512127214805743830[0] = 0;
   out_3512127214805743830[1] = 0;
   out_3512127214805743830[2] = 0;
   out_3512127214805743830[3] = 0;
   out_3512127214805743830[4] = 1;
   out_3512127214805743830[5] = 0;
   out_3512127214805743830[6] = 0;
   out_3512127214805743830[7] = 0;
   out_3512127214805743830[8] = 0;
}
void h_26(double *state, double *unused, double *out_8544531341603956227) {
   out_8544531341603956227[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2747709062575561021) {
   out_2747709062575561021[0] = 0;
   out_2747709062575561021[1] = 0;
   out_2747709062575561021[2] = 0;
   out_2747709062575561021[3] = 0;
   out_2747709062575561021[4] = 0;
   out_2747709062575561021[5] = 0;
   out_2747709062575561021[6] = 0;
   out_2747709062575561021[7] = 1;
   out_2747709062575561021[8] = 0;
}
void h_27(double *state, double *unused, double *out_1616509446902476146) {
   out_1616509446902476146[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5735721285989687047) {
   out_5735721285989687047[0] = 0;
   out_5735721285989687047[1] = 0;
   out_5735721285989687047[2] = 0;
   out_5735721285989687047[3] = 1;
   out_5735721285989687047[4] = 0;
   out_5735721285989687047[5] = 0;
   out_5735721285989687047[6] = 0;
   out_5735721285989687047[7] = 0;
   out_5735721285989687047[8] = 0;
}
void h_29(double *state, double *unused, double *out_4564235437257961107) {
   out_4564235437257961107[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4022358559120136014) {
   out_4022358559120136014[0] = 0;
   out_4022358559120136014[1] = 1;
   out_4022358559120136014[2] = 0;
   out_4022358559120136014[3] = 0;
   out_4022358559120136014[4] = 0;
   out_4022358559120136014[5] = 0;
   out_4022358559120136014[6] = 0;
   out_4022358559120136014[7] = 0;
   out_4022358559120136014[8] = 0;
}
void h_28(double *state, double *unused, double *out_987895873906221505) {
   out_987895873906221505[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1060040457949394560) {
   out_1060040457949394560[0] = 1;
   out_1060040457949394560[1] = 0;
   out_1060040457949394560[2] = 0;
   out_1060040457949394560[3] = 0;
   out_1060040457949394560[4] = 0;
   out_1060040457949394560[5] = 0;
   out_1060040457949394560[6] = 0;
   out_1060040457949394560[7] = 0;
   out_1060040457949394560[8] = 0;
}
void h_31(double *state, double *unused, double *out_3227977152407240354) {
   out_3227977152407240354[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1024440218175455631) {
   out_1024440218175455631[0] = 0;
   out_1024440218175455631[1] = 0;
   out_1024440218175455631[2] = 0;
   out_1024440218175455631[3] = 0;
   out_1024440218175455631[4] = 0;
   out_1024440218175455631[5] = 0;
   out_1024440218175455631[6] = 0;
   out_1024440218175455631[7] = 0;
   out_1024440218175455631[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_7163063568989676842) {
  err_fun(nom_x, delta_x, out_7163063568989676842);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3754260694415491349) {
  inv_err_fun(nom_x, true_x, out_3754260694415491349);
}
void car_H_mod_fun(double *state, double *out_8878594360494595028) {
  H_mod_fun(state, out_8878594360494595028);
}
void car_f_fun(double *state, double dt, double *out_32589037743347182) {
  f_fun(state,  dt, out_32589037743347182);
}
void car_F_fun(double *state, double dt, double *out_2355328765180068083) {
  F_fun(state,  dt, out_2355328765180068083);
}
void car_h_25(double *state, double *unused, double *out_1013193562225998629) {
  h_25(state, unused, out_1013193562225998629);
}
void car_H_25(double *state, double *unused, double *out_993794256298495203) {
  H_25(state, unused, out_993794256298495203);
}
void car_h_24(double *state, double *unused, double *out_4596273345226045289) {
  h_24(state, unused, out_4596273345226045289);
}
void car_H_24(double *state, double *unused, double *out_4551719855100818339) {
  H_24(state, unused, out_4551719855100818339);
}
void car_h_30(double *state, double *unused, double *out_7601021488702293792) {
  h_30(state, unused, out_7601021488702293792);
}
void car_H_30(double *state, double *unused, double *out_3512127214805743830) {
  H_30(state, unused, out_3512127214805743830);
}
void car_h_26(double *state, double *unused, double *out_8544531341603956227) {
  h_26(state, unused, out_8544531341603956227);
}
void car_H_26(double *state, double *unused, double *out_2747709062575561021) {
  H_26(state, unused, out_2747709062575561021);
}
void car_h_27(double *state, double *unused, double *out_1616509446902476146) {
  h_27(state, unused, out_1616509446902476146);
}
void car_H_27(double *state, double *unused, double *out_5735721285989687047) {
  H_27(state, unused, out_5735721285989687047);
}
void car_h_29(double *state, double *unused, double *out_4564235437257961107) {
  h_29(state, unused, out_4564235437257961107);
}
void car_H_29(double *state, double *unused, double *out_4022358559120136014) {
  H_29(state, unused, out_4022358559120136014);
}
void car_h_28(double *state, double *unused, double *out_987895873906221505) {
  h_28(state, unused, out_987895873906221505);
}
void car_H_28(double *state, double *unused, double *out_1060040457949394560) {
  H_28(state, unused, out_1060040457949394560);
}
void car_h_31(double *state, double *unused, double *out_3227977152407240354) {
  h_31(state, unused, out_3227977152407240354);
}
void car_H_31(double *state, double *unused, double *out_1024440218175455631) {
  H_31(state, unused, out_1024440218175455631);
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
