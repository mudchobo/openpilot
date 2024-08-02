#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_4099822328291796061);
void live_err_fun(double *nom_x, double *delta_x, double *out_2280464617631615560);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_9046333686982459057);
void live_H_mod_fun(double *state, double *out_3010053490760277106);
void live_f_fun(double *state, double dt, double *out_57496166683870770);
void live_F_fun(double *state, double dt, double *out_4161128495582430796);
void live_h_4(double *state, double *unused, double *out_307609936800723994);
void live_H_4(double *state, double *unused, double *out_8214539338947963896);
void live_h_9(double *state, double *unused, double *out_9028607953761326231);
void live_H_9(double *state, double *unused, double *out_7973349692318373251);
void live_h_10(double *state, double *unused, double *out_9019171616117907999);
void live_H_10(double *state, double *unused, double *out_8919668588098911707);
void live_h_12(double *state, double *unused, double *out_703346885059506475);
void live_H_12(double *state, double *unused, double *out_7593440313900370229);
void live_h_35(double *state, double *unused, double *out_1633645418999018003);
void live_H_35(double *state, double *unused, double *out_4847877281575356520);
void live_h_32(double *state, double *unused, double *out_3891903989372884448);
void live_H_32(double *state, double *unused, double *out_4716712496634616051);
void live_h_13(double *state, double *unused, double *out_786724181952337315);
void live_H_13(double *state, double *unused, double *out_7488329551323887377);
void live_h_14(double *state, double *unused, double *out_9028607953761326231);
void live_H_14(double *state, double *unused, double *out_7973349692318373251);
void live_h_33(double *state, double *unused, double *out_1221708896165955769);
void live_H_33(double *state, double *unused, double *out_1697320276936498916);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}