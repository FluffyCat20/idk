#ifndef DATA_H
#define DATA_H
#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

const double eps = 0.000000001;
const double eps_for_rf_rf = 0.0001;

struct data_node{
  std::vector<double> U; //rho, rho*u, e
  std::vector<double> F; //rho*u, p + rho*u^2, (e+p)*u
  double rho, p, u, e, a;
  double gamma;

  data_node(double rho_, double p_, double u_, double gamma_)
    : rho(rho_), p(p_), u(u_), gamma(gamma_) {
    e = p/(gamma - 1) + rho*u*u*0.5;
    a = sqrt(gamma*p/rho);
    U.reserve(3);
    F.reserve(3);
    U.push_back(rho);
    U.push_back(rho*u);
    U.push_back(e);
    F.push_back(rho*u);
    F.push_back(p + rho*u*u);
    F.push_back((e + p)*u);
  }

  int calc_F_in_U();
};

struct data{

  std::string method_name = "";

  double x_left, x_right;
  size_t x0; //initial discontinuity coordinate, actually a node number
  double rho_left, p_left, u_left, rho_right, p_right, u_right;
  double courant_number;
  size_t size;
  double t_end;
  std::vector<data_node> mesh;
  double gamma = 1.4;
  double delta_x;
  double delta_t = -1.0;

  bool stop_now = false;

  std::vector<std::vector<double>> D; //artificial viscosity terms

  data(std::ifstream& config_file_path){
    get_data_from_config(config_file_path);    
    delta_x = (x_right - x_left)/size;
    x0 = size/2;
    mesh.reserve(size);
    for (size_t i = 0; i < x0; ++i){
      mesh.push_back(data_node(rho_left, p_left, u_left, gamma));
    }
    for (size_t i = x0; i < size; ++i){
      mesh.push_back(data_node(rho_right, p_right, u_right, gamma));
    }
    /*for (int i = 0; i < size; ++i){
      double x = double(i)/size*10.0 - 5.0;
      double expon = exp(-x*x*0.5);
      mesh.push_back(data_node(expon + 1.0, 1.0, 1.0, gamma));
    }*/

    if (method_name == "mac_cormack+davis"){
      D.resize(size-1);//D[j] = D(j+1/2)
      for (size_t i = 0; i < size - 1; ++i) {
        D[i].resize(3);
      }
      D[0] = {0, 0, 0};
      D[size-2] = {0, 0, 0};
    }
  }

  void get_data_from_config(std::ifstream& input);
  double calc_delta_t();
  void calc_new_U_with_lax_friedrichs(const data& prev_grid);
  void calc_D(); //Davis artificial viscosity terms
  void mac_cormack_predictor_step (const data& prev_grid);
  void mac_cormack_corrector_step(const data& prev_grid, const data& predictor_grid);
  void calc_new_U_with_mac_cormack(data& prev_grid);
  void calc_new_U_with_mac_cormack_davis(data& prev_grid);
  void calc_new_U_and_F_with_godunov(data& prev_grid, double current_t);
  void calc_F_in_U();
  void calc_F_btw_U(double& max_velocity);
  data& operator=(const data& other);

};

int Smooth_Davis_Maksimov(data& new_grid, const data& prev_grid);

#endif // DATA_H
