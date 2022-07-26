#ifndef DATA_NODES_HPP
#define DATA_NODES_HPP

#pragma once

#include <vector>
#include <cmath>

class data_node_2d {
public:
  std::vector<double> U;//rho, rho*u, rho*v, e
  std::vector<double> F, G;
  //F: rho*u, p + rho*u^2, rho*u*v, (e + p)*u
  //G: rho*v, rho*u*v, p + rho*v^2, (e + p)*v
  double rho, p, u, v, u_abs, e, a; //u_abs = velocity vector length
  static double gamma;

  virtual void calc_U_when_values_known() = 0;
  virtual void calc_F_when_values_known() = 0;
  virtual void calc_G_when_values_known() = 0;

  virtual void calc_values_from_U() = 0;


  virtual ~data_node_2d(){};

protected:
  data_node_2d() {};

  data_node_2d(double rho_, double p_, double u_, double v_)
    : rho(rho_), p(p_), u(u_), v(v_) {
    u_abs = std::sqrt(u * u + v * v);
    e = p / (gamma - 1) + rho * u_abs * u_abs * 0.5;
    a = std::sqrt(gamma * p / rho);
    U.resize(4);
    F.resize(4);
    G.resize(4);
    /*calc_U_when_values_known();
    calc_F_when_values_known();
    calc_G_when_values_known();*/
  }

};


class data_node_cartesian : public data_node_2d {
public:

  void calc_U_when_values_known() override {
    U[0] = rho;
    U[1] = rho*u;
    U[2] = rho*v;
    U[3] = e;
  }
  void calc_F_when_values_known() override {
    F[0] = rho*u;
    F[1] = p + rho*u*u;
    F[2] = rho*u*v;
    F[3] = (e + p)*u;
  }
  void calc_G_when_values_known() override {
    G[0] = rho*v;
    G[1] = rho*u*v;
    G[2] = p + rho*v*v;
    G[3] = (e + p)*v;
  }

  void calc_values_from_U() override {
    rho = U[0];
    u = U[1]/rho;
    v = U[2]/rho;
    e = U[3];

    u_abs = std::sqrt(u*u + v*v);
    p = (gamma - 1)*(e - rho*u_abs*u_abs*0.5);
  }

  data_node_cartesian(
      double rho_, double p_, double u_, double v_) :
    data_node_2d(rho_, p_, u_, v_) {
    calc_U_when_values_known();
    calc_F_when_values_known();
    calc_G_when_values_known();
  }

  data_node_cartesian(const data_node_2d* node_ptr) :
    data_node_2d(){
    rho = node_ptr->rho;
    p = node_ptr->p;
    u = node_ptr->u;
    v = node_ptr->v;
    u_abs = node_ptr->u_abs;
    e = node_ptr->e;
    a = node_ptr->a;
    U = node_ptr->U;
    F = node_ptr->F;
    G = node_ptr->G;
  }

};

////H ещё нету вообще!!!!
class data_node_axys_symm : public data_node_2d {
public:

  void calc_U_when_values_known() override {
    U[0] = rho * r;
    U[1] = rho * u * r;
    U[2] = rho * v * r;
    U[3] = e * r;
  }
  void calc_F_when_values_known() override {
    F[0] = rho * u * r;
    F[1] = (p + rho * u * u) * r;
    F[2] = rho * u * v * r;
    F[3] = (e + p) * u * r;
  }
  void calc_G_when_values_known() override {
    G[0] = rho * v * r;
    G[1] = rho * u * v * r;
    G[2] = (p + rho * v * v) * r;
    G[3] = (e + p) * v * r;
  }

  void calc_values_from_U() override {
    rho = U[0] / r;
    u = U[1] / U[0];
    v = U[2] / U[0];
    e = U[3] / r;

    u_abs = std::sqrt(u*u + v*v);
    p = (gamma - 1)*(e - rho*u_abs*u_abs*0.5);
  }

  data_node_axys_symm(
      double rho_, double p_, double u_, double v_, double r_) :
    data_node_2d(rho_, p_, u_, v_), r(r_) {
    calc_U_when_values_known();
    calc_F_when_values_known();
    calc_G_when_values_known();
  }

private:
  double r;
};


#endif // DATA_NODES_HPP
