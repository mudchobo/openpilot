#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_227391381542670811);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4332500849444111319);
void car_H_mod_fun(double *state, double *out_4301467965206663409);
void car_f_fun(double *state, double dt, double *out_6090612947722609695);
void car_F_fun(double *state, double dt, double *out_5556904933626449349);
void car_h_25(double *state, double *unused, double *out_2902874667874020476);
void car_H_25(double *state, double *unused, double *out_4038467864605791628);
void car_h_24(double *state, double *unused, double *out_7947337562859115535);
void car_H_24(double *state, double *unused, double *out_8784977191572300662);
void car_h_30(double *state, double *unused, double *out_8581150449551952113);
void car_H_30(double *state, double *unused, double *out_2878222476885825127);
void car_h_26(double *state, double *unused, double *out_1459936804513372850);
void car_H_26(double *state, double *unused, double *out_7779971183479847852);
void car_h_27(double *state, double *unused, double *out_8387958699214852931);
void car_H_27(double *state, double *unused, double *out_703459165085400216);
void car_h_29(double *state, double *unused, double *out_2373992662071331058);
void car_H_29(double *state, double *unused, double *out_3388453821200217311);
void car_h_28(double *state, double *unused, double *out_1653128554850472032);
void car_H_28(double *state, double *unused, double *out_6092302578853681391);
void car_h_31(double *state, double *unused, double *out_6540324196128042826);
void car_H_31(double *state, double *unused, double *out_4007821902728831200);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}