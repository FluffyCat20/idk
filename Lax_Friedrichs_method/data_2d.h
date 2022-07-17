#ifndef DATA_2D_H
#define DATA_2D_H
#pragma once

#include <iomanip>
#include <list>
#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <limits>

#include "data_nodes.h"

void calc_state_after_shock_wave(
  double& rho, double& p, double& u,
    const double rho0, const double p0, const double u0);

struct calculation_params {

  int method_number = -1;
  double x_left, x_right;
  double y_bottom, y_top;
  size_t size_x, size_y;
  double t_end;
  double time_step;
  double delta_x, delta_y, delta_t;
  double D_of_initial_shock = 0;
  size_t x0, y0; //initial discontinuity node number
  double courant_number;
  double gamma;
  double gap_btw_sw_and_bound;
};

template <class node_t>
class data_2d {

public:

  std::vector<std::vector<node_t>> mesh;
  calculation_params par;

  bool stop_now = false;

  std::vector<std::list<std::pair<double, double>>>
    pressure_sensors_on_wall; //<time, pressure> in a wall point


  data_2d(const calculation_params& par_,
          std::unordered_map<std::string, double>& init_data,
          const std::unordered_map<std::string, bool>& init_config) : par(par_) {

    pressure_sensors_on_wall.resize(par.size_y);

    mesh.resize(par.size_y);
    for (size_t i = 0; i < par.size_y; ++i) {
      mesh[i].reserve(par.size_x);
    }

    if (init_config.at("horizontal_flow")) {

      if (init_config.at("shock_wave")) {
        shock_wave_initialization(
          init_data["rho_right"], init_data["p_right"], init_data["u_right"],
          init_data["rho_left"], init_data["p_left"], init_data["u_left"],
          init_data["mach"]);
      }

      size_t i_start = 0;

      if (init_config.at("horizontal_left_contact_disc")) {
        i_start = par.y0;
        for (size_t i = 0; i < i_start; ++i) {
          for (size_t j = 0; j < par.x0; ++j) {
            mesh[i].emplace_back(init_data["rho_left"]*init_data["omega"],
                init_data["p_left"], init_data["u_left"],
                init_data["v_left"]);
          }
          for (size_t j = par.x0; j < par.size_x; ++j) {
            mesh[i].emplace_back(init_data["rho_right"],
                init_data["p_right"], init_data["u_right"],
                init_data["v_right"]);
          }
        }
      }

      for (size_t i = i_start; i < par.size_y; ++i) {
        for (size_t j = 0; j < par.x0; ++j) {
          mesh[i].emplace_back(init_data["rho_left"],
              init_data["p_left"], init_data["u_left"],
              init_data["v_left"]);
        }
        for (size_t j = par.x0; j < par.size_x; ++j) {
          mesh[i].emplace_back(init_data["rho_right"],
              init_data["p_right"], init_data["u_right"],
              init_data["v_right"]);
        }
      }
    }

    if (init_config.at("vertical_flow")) {
      if (init_config.at("shock_wave")) {
        shock_wave_initialization(
          init_data["rho_down"], init_data["p_down"], init_data["v_down"],
          init_data["rho_up"], init_data["p_up"], init_data["v_up"],
          init_data["mach"]);
      }

      for (size_t i = 0; i < par.y0; ++i) {
        for (size_t j = 0; j < par.size_x; ++j) {
          mesh[i].emplace_back(init_data["rho_up"],
              init_data["p_up"], init_data["u_up"],
              init_data["v_up"]);
        }
      }

      for (size_t i = par.y0; i < par.size_y; ++i) {
        for (size_t j = 0; j < par.size_x; ++j) {
          mesh[i].emplace_back(init_data["rho_down"],
              init_data["p_down"], init_data["u_down"],
              init_data["v_down"]);
        }
      }
    }

    if (init_config.at("quadrants")) {

      for (size_t i = 0; i < par.y0; ++i) {
        for (size_t j = 0; j < par.x0; ++j) {
          mesh[i].emplace_back(init_data["rho_down_left"],
              init_data["p_down_left"], init_data["u_down_left"],
              init_data["v_down_left"]);
        }
      }

      for (size_t i = 0; i < par.y0; ++i) {
        for (size_t j = par.x0; j < par.size_x; ++j) {
          mesh[i].emplace_back(init_data["rho_down_right"],
              init_data["p_down_right"], init_data["u_down_right"],
              init_data["v_down_right"]);
        }
      }

      for (size_t i = par.y0; i < par.size_y; ++i) {
        for (size_t j = 0; j < par.x0; ++j) {
          mesh[i].emplace_back(init_data["rho_up_left"],
              init_data["p_up_left"], init_data["u_up_left"],
              init_data["v_up_left"]);
        }
      }

      for (size_t i = par.y0; i < par.size_y; ++i) {
        for (size_t j = par.x0; j < par.size_x; ++j) {
          mesh[i].emplace_back(init_data["rho_up_right"],
              init_data["p_up_right"], init_data["u_up_right"],
              init_data["v_up_right"]);
        }
      }
    }

    if (init_config.at("bubble_near_wall")) {
      double nodes_in_1_x_double = static_cast<double> (par.size_x) / (par.x_right - par.x_left);
      double nodes_in_1_y_double = static_cast<double> (par.size_y) / (par.y_top - par.y_bottom);
      double eps = 1e-06;
      if (std::abs(nodes_in_1_x_double - nodes_in_1_y_double) > eps) {
        std::cout << std::abs(nodes_in_1_x_double - nodes_in_1_y_double) <<
          " < " << eps << std::endl;
        std::cout << "you fucked up: different scales for x, y" << std::endl;
        return;
      }
      size_t nodes_in_1_x = std::round(nodes_in_1_x_double);
      size_t nodes_in_1_y = std::round(nodes_in_1_y_double);
      par.y0 = (init_data["y0"] - par.y_bottom) * nodes_in_1_y;
      par.x0 = (init_data["x0"] - par.x_left) * nodes_in_1_x;
      double rho_after_sw, p_after_sw, u_after_sw;
      double u_bfr_sw_if_sw_stays;
      shock_wave_initialization(rho_after_sw, p_after_sw, u_after_sw,
        init_data["rho_around_bubble"], init_data["p_around_bubble"],
        u_bfr_sw_if_sw_stays, init_data["mach"]);
      par.D_of_initial_shock = u_bfr_sw_if_sw_stays; //we assume there that u before sw = 0 in config!!!
      u_after_sw = - u_after_sw + u_bfr_sw_if_sw_stays;  //we assume there that u before sw = 0 in config!!!
      double v_after_sw = 0.0;
      par.gap_btw_sw_and_bound = init_data["gap_btw_sw_and_bound"];
      size_t shock_wave_initial_x = par.gap_btw_sw_and_bound * nodes_in_1_x;
      double rho_in_bubble = init_data["rho_around_bubble"] * init_data["omega"];
      size_t a = init_data["a"] * nodes_in_1_x;
      double b_div_by_a = init_data["b_div_by_a"];
      size_t b = static_cast<size_t>(a * b_div_by_a);
      for (size_t y = 0; y < par.size_y; ++y) {
        int y_diff = y - par.y0;
        size_t left_bubble_edge = par.x0;
        size_t right_bubble_edge = par.x0;
        if (std::abs(y_diff) <= b) {
          size_t bubble_x = static_cast<size_t>(
            std::sqrt(a*a - y_diff*y_diff/b_div_by_a/b_div_by_a));
          left_bubble_edge -= bubble_x;
          right_bubble_edge += bubble_x;
        }
        for (size_t x = 0; x < shock_wave_initial_x; ++x) {
          mesh[y].emplace_back(rho_after_sw, p_after_sw, u_after_sw,
            v_after_sw);
        }
        for (size_t x = shock_wave_initial_x; x < left_bubble_edge; ++x) {
          mesh[y].emplace_back(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"]);
        }

        if (left_bubble_edge == par.x0 && right_bubble_edge == par.x0) {
          mesh[y].emplace_back(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"]);
        } else {
          for (size_t x = left_bubble_edge; x < right_bubble_edge + 1; ++x) {
            mesh[y].emplace_back(
              rho_in_bubble, init_data["p_around_bubble"],
              init_data["u_around_bubble"], init_data["v_around_bubble"]);
          }
        }

        for (size_t x = right_bubble_edge + 1; x < par.size_x; ++x) {
          mesh[y].emplace_back(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"]);
        }
      }
    }

  }

public:

  void do_one_step(const data_2d<node_t>& prev_grid);
  void shock_wave_initialization(
    double& rho, double& p, double& u,
    double rho1, double p1, double& u1,
    double mach); //1 - before sw, no ind - after sw
  double calc_delta_t();
  void lax_friedrichs(const data_2d& prev_grid);
  void mac_cormack(const data_2d& prev_grid);
  void mac_cormack_predictor_step(data_2d& predictor_grid) const;
  void mac_cormack_corrector_step(
    const data_2d& prev_grid, const data_2d& predictor_grid);
  void mac_cormack_with_davis(const data_2d& prev_grid);
  void calc_davis_artificial_viscosity();

  void calc_values_and_fluxes();
  void boundary_conditions();
  void boundary_conditions_for_bubble_near_wall();
  void boundary_conditions_for_bubble_near_wall_simmetry_on_bottom();
  void boundary_conditions_default();

  void update_pressure_sensors_on_wall(double t);
  void get_pressure_sensors_from_files(const std::string& sensors_folder);


  void calc_pressure_after_reflected_shock_no_bubble(
    const double rho2, const double p2, const double u2,
    double& p3);
  /// 1 - before incident shock
  /// 2 - after incident shock
  /// 3 - between the shock and the wall after the shock's reflection
  /// u1 = u3 = 0 - boundary condition on the wall

  std::vector<double> calc_pressure_impulses_basic(
    double t0, double p0) const;
  /// excess pressure sum comparing to the pressure behind reflected shock
  /// t0 - time of coming of the initial shock on the wall (with no bubble)
  /// p0 - pressure behind reflected shock (with no bubble)


  std::vector<double> calc_pressure_impulses_max_peak_bfr_p0
    (double p0) const;
  /// area of highest peak over the line p = p0
  /// p0 - pressure behind reflected shock (with no bubble)

  std::vector<double>
  calc_pressure_impulses_max_peak_bounded_by_local_min() const;
  ///area of highest peak bounded on the left&right by local minima

  data_2d& operator=(const data_2d& other);
};

#include "data_2d.inl"
#endif // DATA_2D_H
