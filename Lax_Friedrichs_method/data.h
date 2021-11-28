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
    : rho(rho_), p(p_), u(u_), gamma(gamma_){
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

  int calc_F();
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
  void calc_F();
  data& operator=(const data& other);
  std::vector<std::vector<double> > calc_error_norm(double* R1, double* P1, double* U1, size_t size);

};

class exact_solution {
public:

  //initial values
  double x_left, x_right;

  double p1, rho1, u1,
    p2, rho2, u2;

  double a1, a2, v1, v2;       // v = 1/rho

  double gamma;

  size_t size;
  double t_end;

  //to calculate:
  double u3 = 0, p3,
    rho31, rho32,
    a31 = 0, a32 = 0,
    d1 = 0, d2 = 0;

  double delta_x;

  //left and right from contact discontinuity
  char left, right; //s - shock wave, r - riemann fan, v - vacuum
  int case_number;
  /* 0: Sh-Sh
   * 1: Sh-Rf
   * 2: Rf-Sh
   * 3: Rf-Rf
   * 4: Vacuum */

  double p_max = 1.0;

  double shift = 0.0; //for output: inside solver, interval must be symmetric

  std::vector<double> waves_velocities;

  exact_solution(double p1_, double rho1_, double u1_, double p2_, double rho2_, double u2_,
                 double gamma_, double t_end_, double x_left_, double x_right_, size_t size_)
    : p1(p1_), rho1(rho1_), u1(u1_), p2(p2_), rho2(rho2_), u2(u2_), gamma(gamma_),
      t_end(t_end_), x_left(x_left_), x_right(x_right_), size(size_) {
    //in this constructor we assume that interval is symmetric!
    v1 = 1.0/rho1;
    v2 = 1.0/rho2;
    a1 = std::sqrt(gamma*p1/rho1);
    a2 = std::sqrt(gamma*p2/rho2);
    delta_x = (x_right-x_left)/size;
    p_max = std::max(p1, p2);
  };

  exact_solution(double p1_, double rho1_, double u1_,
                 double p2_, double rho2_, double u2_,
                 double gamma_) :
    p1(p1_), rho1(rho1_), u1(u1_), p2(p2_), rho2(rho2_), u2(u2_), gamma(gamma_) {

    //in this constructor we assume that interval is symmetric!
    v1 = 1.0/rho1;
    v2 = 1.0/rho2;
    a1 = std::sqrt(gamma*p1/rho1);
    a2 = std::sqrt(gamma*p2/rho2);
    p_max = std::max(p1, p2);
    waves_velocities.reserve(5);
  }

  exact_solution(const data& dt) :
    rho1(dt.rho_left), p1(dt.p_left), u1(dt.u_left),
    rho2(dt.rho_right), p2(dt.p_right), u2(dt.u_right),
    gamma(dt.gamma), t_end(dt.t_end), size(dt.size){
    x_left = 0.5*(dt.x_left - dt.x_right);
    x_right = -x_left;
    shift = -0.5*(dt.x_left + dt.x_right);
    v1 = 1.0/rho1;
    v2 = 1.0/rho2;
    a1 = std::sqrt(gamma*p1/rho1);
    a2 = std::sqrt(gamma*p2/rho2);
    delta_x = (x_right-x_left)/size;
    p_max = std::max(p1, p2);
  }

  void shock_wave_u(double& f, double& f_der,
      double p_i, double u_i, double v_i, double p, bool right_dir);
  void riemann_fan_u(double& f, double& f_der,
      double p_i, double u_i, double a_i, double p, bool right_dir);


  void case_number_choice();
  void f_and_der_calc(double& f, double& f_der, double p);
  double newton_solver();

  void calc_p3();
  void calc_left_sh();
  void calc_right_sh();
  void calc_left_rf();
  void calc_right_rf();
  void calc_left_vac();
  void calc_right_vac();

  void calc_and_print_result(std::ofstream& out);
  double calc_max_velocity_and_waves_velocities();
  void find_values_in_zero(double& rho_disc, double& p_disc, double& u_disc);

};

int Smooth_Davis_Maksimov(data& new_grid, const data& prev_grid);

#endif // DATA_H
