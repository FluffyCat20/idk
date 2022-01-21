#ifndef DATA_2D_H
#define DATA_2D_H
#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

struct data_node_2d{
  std::vector<double> U;//rho, rho*u, rho*v, e
  std::vector<double> F, G;
  //F: rho*u, p + rho*u^2, rho*u*v, (e + p)*u
  //G: rho*v, rho*u*v, p + rho*v^2, (e + p)*v
  double rho, p, u, v, u_abs, e, a; //u_abs = velocity vector length
  double gamma;

  void calc_U_when_values_known() {
    U.resize(4);
    U[0] = rho;
    U[1] = rho*u;
    U[2] = rho*v;
    U[3] = e;
  }
  void calc_F_when_values_known() {
    F.resize(4);
    F[0] = rho*u;
    F[1] = p + rho*u*u;
    F[2] = rho*u*v;
    F[3] = (e + p)*u;
  }
  void calc_G_when_values_known() {
    G.resize(4);
    G[0] = rho*v;
    G[1] = rho*u*v;
    G[2] = p + rho*v*v;
    G[3] = (e + p)*v;
  }

  void calc_values_from_U() {
    rho = U[0];
    u = U[1]/rho;
    v = U[2]/rho;
    e = U[3];

    u_abs = std::sqrt(u*u + v*v);
    p = (gamma - 1)*(e - rho*u_abs*u_abs*0.5);
  }

  data_node_2d(double rho_, double p_, double u_, double v_, double gamma_)
    : rho(rho_), p(p_), u(u_), v(v_), gamma(gamma_) {
    u_abs = std::sqrt(u*u + v*v);
    e = p/(gamma - 1) + rho*u_abs*u_abs*0.5;
    a = std::sqrt(gamma*p/rho);
    calc_U_when_values_known();
    calc_F_when_values_known();
    calc_G_when_values_known();
  }

};

struct data_2d{

  std::string method_name = "";

  double x_left, x_right;
  double y_bottom, y_top;

  size_t size_x, size_y;
  double t_end;

  double delta_x, delta_y, delta_t;

  size_t x0; //initial (vertical!) discontinuity node number

  double rho_left, p_left, u_left, v_left,
         rho_right, p_right, u_right, v_right;
   //initial data for vertical discontinuity

  double courant_number;

  double gamma;

  double mach;
  double omega;

  std::vector<std::vector<data_node_2d>> mesh;

  bool stop_now = false;

  data_2d(std::ifstream& config_file_path) {

    get_data_from_config(config_file_path);

    shock_wave_initialization();

    delta_x = (x_right - x_left)/size_x;
    delta_y = (y_top - y_bottom)/size_y;
    x0 = size_x/2;
    mesh.resize(size_y);
    for (size_t i = 0; i < size_y; ++i) {
      mesh[i].reserve(size_x);
    }
    for (size_t i = 0; i < size_y/2; ++i) {
      for (size_t j = 0; j < x0; ++j) {
        mesh[i].emplace_back(omega*rho_left, p_left, u_left, v_left, gamma);
      }
      for (size_t j = x0; j < size_x; ++j) {
        mesh[i].emplace_back(rho_right, p_right, u_right, v_right, gamma);
      }
    }
    for (size_t i = size_y/2; i < size_y; ++i) {
      for (size_t j = 0; j < x0; ++j) {
        mesh[i].emplace_back(rho_left, p_left, u_left, v_left, gamma);
      }
      for (size_t j = x0; j < size_x; ++j) {
        mesh[i].emplace_back(rho_right, p_right, u_right, v_right, gamma);
      }
    }
  }

  void get_data_from_config(std::ifstream& input);
  double calc_delta_t();
  void mac_cormack(const data_2d& prev_grid);
  void mac_cormack_predictor_step(data_2d& predictor_grid)const;
  void mac_cormack_corrector_step(
      const data_2d& prev_grid, const data_2d& predictor_grid);
  void shock_wave_initialization();

  void boundary_conditions();

  void write_out_file(std::ofstream& outfile);

  data_2d& operator=(const data_2d& other);
};


#endif // DATA_2D_H
