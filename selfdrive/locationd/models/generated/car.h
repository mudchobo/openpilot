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
void car_err_fun(double *nom_x, double *delta_x, double *out_5629442633965842572);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3598908807355077436);
void car_H_mod_fun(double *state, double *out_53256178129770178);
void car_f_fun(double *state, double dt, double *out_2357487969533120980);
void car_F_fun(double *state, double dt, double *out_8557443704727234356);
void car_h_25(double *state, double *unused, double *out_1115798248010839447);
void car_H_25(double *state, double *unused, double *out_5831744492098185822);
void car_h_24(double *state, double *unused, double *out_3529007286521349464);
void car_H_24(double *state, double *unused, double *out_3659094893092686256);
void car_h_30(double *state, double *unused, double *out_3557753197258034570);
void car_H_30(double *state, double *unused, double *out_1304048161970577624);
void car_h_26(double *state, double *unused, double *out_4227657190368867009);
void car_H_26(double *state, double *unused, double *out_2090241173224129598);
void car_h_27(double *state, double *unused, double *out_2700364704332613072);
void car_H_27(double *state, double *unused, double *out_3527642233154520841);
void car_h_29(double *state, double *unused, double *out_4161069081934512087);
void car_H_29(double *state, double *unused, double *out_1814279506284969808);
void car_h_28(double *state, double *unused, double *out_8458349618326351028);
void car_H_28(double *state, double *unused, double *out_3777909777850296059);
void car_h_31(double *state, double *unused, double *out_354890393280493669);
void car_H_31(double *state, double *unused, double *out_5862390453975146250);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}