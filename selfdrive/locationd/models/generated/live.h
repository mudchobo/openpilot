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
void live_H(double *in_vec, double *out_2165821824918486469);
void live_err_fun(double *nom_x, double *delta_x, double *out_3216996441859096695);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3507935633117718260);
void live_H_mod_fun(double *state, double *out_3893200104662780740);
void live_f_fun(double *state, double dt, double *out_3255786806568313974);
void live_F_fun(double *state, double dt, double *out_605085115782193670);
void live_h_4(double *state, double *unused, double *out_1720682332973874395);
void live_H_4(double *state, double *unused, double *out_7718129108773990032);
void live_h_9(double *state, double *unused, double *out_6654091283807724021);
void live_H_9(double *state, double *unused, double *out_430910173509542562);
void live_h_10(double *state, double *unused, double *out_8850294366417281417);
void live_H_10(double *state, double *unused, double *out_6402008631363930705);
void live_h_12(double *state, double *unused, double *out_763950336527271914);
void live_H_12(double *state, double *unused, double *out_2698672700742028237);
void live_h_35(double *state, double *unused, double *out_2733917337483254325);
void live_H_35(double *state, double *unused, double *out_4351467051401382656);
void live_h_32(double *state, double *unused, double *out_8397599966418024468);
void live_H_32(double *state, double *unused, double *out_2410628714068373855);
void live_h_13(double *state, double *unused, double *out_6024542805875939123);
void live_H_13(double *state, double *unused, double *out_3379336200210811815);
void live_h_14(double *state, double *unused, double *out_6654091283807724021);
void live_H_14(double *state, double *unused, double *out_430910173509542562);
void live_h_33(double *state, double *unused, double *out_82084314392906247);
void live_H_33(double *state, double *unused, double *out_1200910046762525052);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}