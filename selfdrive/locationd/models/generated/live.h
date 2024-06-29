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
void live_H(double *in_vec, double *out_1654711282601658745);
void live_err_fun(double *nom_x, double *delta_x, double *out_5734741546333750642);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_997272261207104900);
void live_H_mod_fun(double *state, double *out_6112760305448342788);
void live_f_fun(double *state, double dt, double *out_5727438172870790233);
void live_F_fun(double *state, double dt, double *out_7019330850193311240);
void live_h_4(double *state, double *unused, double *out_6704644497324442738);
void live_H_4(double *state, double *unused, double *out_457666873631447355);
void live_h_9(double *state, double *unused, double *out_6041008979035423841);
void live_H_9(double *state, double *unused, double *out_6829552061633000115);
void live_h_10(double *state, double *unused, double *out_510443638577489006);
void live_H_10(double *state, double *unused, double *out_2575019310555904300);
void live_h_12(double *state, double *unused, double *out_724301251006094885);
void live_H_12(double *state, double *unused, double *out_163432151416146312);
void live_h_35(double *state, double *unused, double *out_1094377590467689528);
void live_H_35(double *state, double *unused, double *out_2908995183741160021);
void live_h_32(double *state, double *unused, double *out_1396309647877035484);
void live_H_32(double *state, double *unused, double *out_3040159968681900490);
void live_h_13(double *state, double *unused, double *out_279510485795163836);
void live_H_13(double *state, double *unused, double *out_8153037128550556083);
void live_h_14(double *state, double *unused, double *out_6041008979035423841);
void live_H_14(double *state, double *unused, double *out_6829552061633000115);
void live_h_33(double *state, double *unused, double *out_721462086167171751);
void live_H_33(double *state, double *unused, double *out_6059552188380017625);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}