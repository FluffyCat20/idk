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
  int dimension;
  std::vector<double> U;
  std::vector<double> F, G; //fluxes
  std::vector<double> u; //velocity
  double u_abs;
  double rho, p, e, a;
  double gamma;

  void calc_F();
  void calc_G();
  void init_U_F_1D();
  void init_U_F_G_2D();

  data_node(int dim, double rho_, double p_, std::vector<double> u_, double gamma_)
    : dimension(dim), rho(rho_), p(p_), u(u_), gamma(gamma_) {
    a = sqrt(gamma*p/rho);
    switch (dimension) {
    case 1: {
      e = p/(gamma - 1) + rho*u[0]*u[0]*0.5;
      init_U_F_1D();
      break;
    }
    case 2: {
      u_abs = std::sqrt(u[0]*u[0] + u[1]*u[1]);
      e = p/(gamma - 1) + rho*u_abs*u_abs*0.5;
      init_U_F_G_2D();
      break;
    }
    default: {
      break;
    }
    }
  } //end of constructor

  int calc_F_in_node();
};

struct data{

  int dimension;

  std::string method_name = "";

  double x_left, x_right;
  double y_bottom, y_top;
  size_t x0; //initial discontinuity coordinate, actually a node number
  double rho_left, p_left, rho_right, p_right;
  std::vector<double> u_left, u_right;
  double courant_number;
  size_t size_x, size_y;
  double t_end;
  std::vector<data_node> mesh_1D;
  std::vector<std::vector<data_node>> mesh;
  double gamma;
  double delta_x, delta_y;
  double delta_t;

  bool stop_now = false;

  std::vector<std::vector<double>> D; //artificial viscosity terms

  void initialize_row(std::vector<data_node>& row) {
    row.reserve(size_x);
    for (size_t i = 0; i < x0; ++i){
      row.emplace_back(dimension, rho_left, p_left, u_left, gamma);
    }
    for (size_t i = x0; i < size_x; ++i){
      row.emplace_back(dimension, rho_right, p_right, u_right, gamma);
    }
  }

  data(std::ifstream& config_file_path) {
    get_data_from_config(config_file_path);
    delta_x = (x_right - x_left)/size_x;
    x0 = size_x/2;
    if (dimension == 1) {
      size_y = 3;
    }
    if (dimension == 1) {
      initialize_row(mesh_1D);
      /*for (int i = 0; i < size; ++i){
        double x = double(i)/size*10.0 - 5.0;
        double expon = exp(-x*x*0.5);
        mesh.push_back(data_node(expon + 1.0, 1.0, 1.0, gamma));
      }*/

      if (method_name == "mac_cormack+davis"){
        D.resize(size_x - 1);//D[j] = D(j+1/2)
        for (size_t i = 0; i < size_x - 1; ++i) {
          D[i].resize(3);
        }
        D[0] = {0, 0, 0};
        D[size_x-2] = {0, 0, 0};
      }
    }
    delta_y = (y_top - y_bottom)/size_y;
    mesh.resize(size_y);
    for (auto& row : mesh) {
      initialize_row(row);
    }
  }

  void get_data_from_config(std::ifstream& input);
  double calc_delta_t();
  void calc_new_U_with_lax_friedrichs(const data& prev_grid);
  void calc_D(); //Davis artificial viscosity terms
  void mac_cormack_predictor_step(const data& prev_grid);
  void mac_cormack_corrector_step(const data& prev_grid, const data& predictor_grid);
  void calc_new_U_with_mac_cormack(data& prev_grid);
  void calc_new_U_with_mac_cormack_davis(data& prev_grid);
  void calc_new_U_and_F_with_godunov(data& prev_grid, double current_t);
  void calc_F_in_nodes();
  void calc_F_btw_nodes(double& max_velocity);
  void calc_boundary_conditions();
  data& operator=(const data& other);

};

int Smooth_Davis_Maksimov(data& new_grid, const data& prev_grid);

#endif // DATA_H
