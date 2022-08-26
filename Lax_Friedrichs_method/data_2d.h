#ifndef DATA_2D_H
#define DATA_2D_H
#pragma once

#include <iomanip>
#include <vector>
#include <list>
#include <unordered_map>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "calculation_info.hpp"

void calc_state_after_shock_wave(
  double& rho, double& p, double& u,
    const double rho0, const double p0, const double u0);

class mesh_and_common_methods {

public:

  mesh_and_common_methods(
      calculation_params par_) : par(par_){
    mesh.resize(par.size_y + 2 * par.method_info.ghost_edges_width);
    for (size_t i = 0; i < mesh.size(); ++i) {
      mesh[i].resize(par.size_x + 2 * par.method_info.ghost_edges_width);
    }
  };

  virtual ~mesh_and_common_methods() {
    for (size_t i = 0; i < mesh.size(); ++i) {
      for (size_t j = 0; j < mesh[i].size(); ++j) {
        delete mesh[i][j];
      }
    }
  }

  const std::vector<std::vector<data_node_2d*>>& get_mesh_const_ref() const {
    return mesh;
  }
  const calculation_params& get_params_const_ref() const {
    return par;
  }
  std::vector<std::vector<data_node_2d*>>::iterator
  y_begin() {
    return mesh.begin() + par.y_begin;
  }
  std::vector<std::vector<data_node_2d*>>::iterator
  y_end() {
    return mesh.begin() + par.y_end;
  }

  std::vector<data_node_2d*>::iterator
  x_begin(std::vector<std::vector<data_node_2d*>>::iterator y_it) {
    return y_it->begin() + par.x_begin;
  }
  std::vector<data_node_2d*>::iterator
  x_end(std::vector<std::vector<data_node_2d*>>::iterator y_it) {
    return y_it->begin() + par.x_end;
  }

  size_t y_begin_ind() {
    return par.y_begin;
  }
  size_t y_end_ind() {
    return par.y_end;
  }
  size_t x_begin_ind() {
    return par.x_begin;
  }
  size_t x_end_ind() {
    return par.x_end;
  }

  void do_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr);


  void shock_wave_initialization(
    double& rho, double& p, double& u,
    double rho1, double p1, double& u1,
    double mach); //1 - before sw, no ind - after sw

  void calc_delta_t(double current_time);
  void set_delta_t(double d_t) {
    par.delta_t = d_t;
  }
  void boundary_conditions();


  virtual void calc_values_and_fluxes() = 0;

  virtual void lax_friedrichs(
    std::shared_ptr<const mesh_and_common_methods> prev_grid) = 0;
  virtual void mac_cormack(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) = 0;
  virtual void mac_cormack_predictor_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) = 0;
  virtual void mac_cormack_corrector_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr,
    std::shared_ptr<const mesh_and_common_methods> predictor_grid_ptr) = 0;
  virtual void mac_cormack_with_davis(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) = 0;
  virtual void calc_davis_artificial_viscosity() = 0;

  virtual void boundary_conditions_for_bubble_near_wall() = 0;
  virtual void boundary_conditions_for_bubble_near_wall_simmetry_on_bottom() = 0;
  virtual void boundary_conditions_default() = 0;

protected:
  std::vector<std::vector<data_node_2d*>> mesh;
  calculation_params par;

};

class mesh_and_methods_cartesian : public mesh_and_common_methods {
public:

  mesh_and_methods_cartesian(const calculation_info& calc_info) :
    mesh_and_common_methods(calc_info.par) {
    auto init_config = calc_info.init_config;
    auto init_data = calc_info.init_data;

    //ghost edges filling
    for (int y = 0; y < par.y_begin; ++y) {
      for (size_t x = 0; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_cartesian(0.0, 0.0, 0.0, 0.0);
      }
    }

    for (int y = par.y_begin; y < par.y_end; ++y) {
      for (size_t x = 0; x < par.x_begin; ++x) {
        mesh[y][x] = new data_node_cartesian(0.0, 0.0, 0.0, 0.0);
      }
      for (size_t x = par.x_end; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_cartesian(0.0, 0.0, 0.0, 0.0);
      }
    }

    for (int y = par.y_end; y < mesh.size(); ++y) {
      for (size_t x = 0; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_cartesian(0.0, 0.0, 0.0, 0.0);
      }
    }

    if (init_config["horizontal_flow"]) {

      if (init_config["shock_wave"]) {
        shock_wave_initialization(
          init_data["rho_right"], init_data["p_right"], init_data["u_right"],
          init_data["rho_left"], init_data["p_left"], init_data["u_left"],
          init_data["mach"]);
      }

      size_t i_start = par.y_begin;

      if (init_config["horizontal_left_contact_disc"]) {
        i_start = par.y0;
        for (size_t i = par.y_begin; i < i_start; ++i) {
          for (size_t j = par.x_begin; j < par.x0; ++j) {
            mesh[i][j] = new data_node_cartesian(
                init_data["rho_left"]*init_data["omega"],
                init_data["p_left"], init_data["u_left"],
                init_data["v_left"]);
          }
          for (size_t j = par.x0; j < par.x_end; ++j) {
            mesh[i][j] = new data_node_cartesian(
                init_data["rho_right"],
                init_data["p_right"], init_data["u_right"],
                init_data["v_right"]);
          }
        }
      }

      for (size_t i = i_start; i < par.y_end; ++i) {
        for (size_t j = par.x_begin; j < par.x0; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_left"],
              init_data["p_left"], init_data["u_left"],
              init_data["v_left"]);
        }
        for (size_t j = par.x0; j < par.x_end; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_right"],
              init_data["p_right"], init_data["u_right"],
              init_data["v_right"]);
        }
      }
    }

    if (init_config["vertical_flow"]) {
      if (init_config["shock_wave"]) {
        shock_wave_initialization(
          init_data["rho_down"], init_data["p_down"], init_data["v_down"],
          init_data["rho_up"], init_data["p_up"], init_data["v_up"],
          init_data["mach"]);
      }

      for (size_t i = par.x_begin; i < par.y0; ++i) {
        for (size_t j = par.x_begin; j < par.x_end; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_up"],
              init_data["p_up"], init_data["u_up"],
              init_data["v_up"]);
        }
      }

      for (size_t i = par.y0; i < par.y_end; ++i) {
        for (size_t j = par.x_begin; j < par.x_end; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_down"],
              init_data["p_down"], init_data["u_down"],
              init_data["v_down"]);
        }
      }
    }

    if (init_config["quadrants"]) {

      for (size_t i = par.y_begin; i < par.y0; ++i) {
        for (size_t j = par.x_begin; j < par.x0; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_down_left"],
              init_data["p_down_left"], init_data["u_down_left"],
              init_data["v_down_left"]);
        }
      }

      for (size_t i = par.y_begin; i < par.y0; ++i) {
        for (size_t j = par.x0; j < par.x_end; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_down_right"],
              init_data["p_down_right"], init_data["u_down_right"],
              init_data["v_down_right"]);
        }
      }

      for (size_t i = par.y0; i < par.y_end; ++i) {
        for (size_t j = par.x_begin; j < par.x0; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_up_left"],
              init_data["p_up_left"], init_data["u_up_left"],
              init_data["v_up_left"]);
        }
      }

      for (size_t i = par.y0; i < par.y_end; ++i) {
        for (size_t j = par.x0; j < par.x_end; ++j) {
          mesh[i][j] = new data_node_cartesian(
              init_data["rho_up_right"],
              init_data["p_up_right"], init_data["u_up_right"],
              init_data["v_up_right"]);
        }
      }
    }

    if (init_config["bubble_near_wall"]) {
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
      par.y0 = par.y_begin + (init_data["y0"] - par.y_bottom) * nodes_in_1_y;//пидорство, теперь в calc_info другие параметры
      par.x0 = par.x_begin + (init_data["x0"] - par.x_left) * nodes_in_1_x;
      double rho_after_sw, p_after_sw, u_after_sw;
      double u_bfr_sw_if_sw_stays;
      shock_wave_initialization(rho_after_sw, p_after_sw, u_after_sw,
        init_data["rho_around_bubble"], init_data["p_around_bubble"],
        u_bfr_sw_if_sw_stays, init_data["mach"]);
      par.D_of_initial_shock = u_bfr_sw_if_sw_stays; //we assume there that u before sw = 0 in config!!!
      u_after_sw = - u_after_sw + u_bfr_sw_if_sw_stays;  //we assume there that u before sw = 0 in config!!!
      double v_after_sw = 0.0;
      size_t shock_wave_initial_x = par.x_begin + par.gap_btw_sw_and_bound * nodes_in_1_x;
      double rho_in_bubble = init_data["rho_around_bubble"] * init_data["omega"];
      size_t a = init_data["a"] * nodes_in_1_x;
      double b_div_by_a = init_data["b_div_by_a"];
      size_t b = static_cast<size_t>(a * b_div_by_a);
      for (size_t y = par.y_begin; y < par.y_end; ++y) {
        int y_diff = y - par.y0;
        size_t left_bubble_edge = par.x0;
        size_t right_bubble_edge = par.x0;
        if (std::abs(y_diff) <= b) {
          size_t bubble_x = static_cast<size_t>(
            std::sqrt(a*a - y_diff*y_diff/b_div_by_a/b_div_by_a));
          left_bubble_edge -= bubble_x;
          right_bubble_edge += bubble_x;
        }
        for (size_t x = par.x_begin; x < shock_wave_initial_x; ++x) {
          mesh[y][x] = new data_node_cartesian(
            rho_after_sw, p_after_sw, u_after_sw, v_after_sw);
        }
        for (size_t x = shock_wave_initial_x; x < left_bubble_edge; ++x) {
          mesh[y][x] = new data_node_cartesian(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"]);
        }

        if (left_bubble_edge == par.x0 && right_bubble_edge == par.x0) {
          mesh[y][par.x0] = new data_node_cartesian(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"]);
        } else {
          for (size_t x = left_bubble_edge; x < right_bubble_edge + 1; ++x) {
            mesh[y][x] = new data_node_cartesian(
              rho_in_bubble, init_data["p_around_bubble"],
              init_data["u_around_bubble"], init_data["v_around_bubble"]);
          }
        }

        for (size_t x = right_bubble_edge + 1; x < par.x_end; ++x) {
          mesh[y][x] = new data_node_cartesian(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"]);
        }
      }
    }

  }

  mesh_and_methods_cartesian(calculation_params par) :
    mesh_and_common_methods(par) {};

  mesh_and_methods_cartesian(
      std::shared_ptr<const mesh_and_common_methods> other_grid) :
    mesh_and_common_methods(other_grid->get_params_const_ref()){

    auto other_mesh = other_grid->get_mesh_const_ref();
    mesh.resize(other_mesh.size());
    for (size_t y =0; y < mesh.size(); ++y) {
      mesh[y].resize(other_mesh[y].size());
      for (size_t x = 0; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_cartesian(other_mesh[y][x]);
      }
    }

  };

  void calc_values_and_fluxes() override;

  virtual void lax_friedrichs(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) override;
  virtual void mac_cormack(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) override;
  virtual void mac_cormack_predictor_step( //must be called from predictor grid
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) override;
  virtual void mac_cormack_corrector_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr,
    std::shared_ptr<const mesh_and_common_methods> predictor_grid_ptr) override;
  virtual void mac_cormack_with_davis(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) override;
  virtual void calc_davis_artificial_viscosity() override;

  virtual void boundary_conditions_for_bubble_near_wall() override;
  virtual void boundary_conditions_for_bubble_near_wall_simmetry_on_bottom() override;
  virtual void boundary_conditions_default() override;
};

class mesh_and_methods_axisymm : public mesh_and_common_methods {
public:
  mesh_and_methods_axisymm(const calculation_info& calc_info) :
    mesh_and_common_methods(calc_info.par) {
    auto init_config = calc_info.init_config;
    auto init_data = calc_info.init_data;

    //ghost edges filling
    for (int y = 0; y < par.y_begin; ++y) {
      double r = (y - par.y_begin) * par.delta_y;
      for (size_t x = 0; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_axisymm(0.0, 0.0, 0.0, 0.0, r);
      }
    }

    for (int y = par.y_begin; y < par.y_end; ++y) {
      double r = (y - par.y_begin) * par.delta_y;
      for (size_t x = 0; x < par.x_begin; ++x) {
        mesh[y][x] = new data_node_axisymm(0.0, 0.0, 0.0, 0.0, r);
      }
      for (size_t x = par.x_end; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_axisymm(0.0, 0.0, 0.0, 0.0, r);
      }
    }

    for (int y = par.y_end; y < mesh.size(); ++y) {
      double r = (y - par.y_begin) * par.delta_y;
      for (size_t x = 0; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_axisymm(0.0, 0.0, 0.0, 0.0, r);
      }
    }


    //initial conditions from config:

    if (init_config["vertical_flow"]) {
      if (init_config["shock_wave"]) {
        shock_wave_initialization(
          init_data["rho_down"], init_data["p_down"], init_data["v_down"],
          init_data["rho_up"], init_data["p_up"], init_data["v_up"],
          init_data["mach"]);
      }

      for (int y = par.y_begin; y < par.y0; ++y) {
        for (size_t x = par.x_begin; x < par.x_end; ++x) {
          mesh[y][x] = new data_node_axisymm(
              init_data["rho_up"],
              init_data["p_up"], init_data["u_up"],
              init_data["v_up"], (y - par.y_begin) * par.delta_y);
        }
      }

      for (int y = par.y0; y < par.y_end; ++y) {
        for (size_t x = par.x_begin; x < par.x_end; ++x) {
          mesh[y][x] = new data_node_axisymm(
              init_data["rho_down"],
              init_data["p_down"], init_data["u_down"],
              init_data["v_down"], (y - par.y_begin) * par.delta_y);
        }
      }
    }

    if (init_config["bubble_near_wall"]) {
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
      par.y0 = par.y_begin + (init_data["y0"] - par.y_bottom) * nodes_in_1_y;//пидорство, теперь в calc_info другие параметры
      par.x0 = par.x_begin + (init_data["x0"] - par.x_left) * nodes_in_1_x;
      double rho_after_sw, p_after_sw, u_after_sw;
      double u_bfr_sw_if_sw_stays;
      shock_wave_initialization(rho_after_sw, p_after_sw, u_after_sw,
        init_data["rho_around_bubble"], init_data["p_around_bubble"],
        u_bfr_sw_if_sw_stays, init_data["mach"]);
      par.D_of_initial_shock = u_bfr_sw_if_sw_stays; //we assume there that u before sw = 0 in config!!!
      u_after_sw = - u_after_sw + u_bfr_sw_if_sw_stays;  //we assume there that u before sw = 0 in config!!!
      double v_after_sw = 0.0;
      size_t shock_wave_initial_x = par.gap_btw_sw_and_bound * nodes_in_1_x;
      double rho_in_bubble = init_data["rho_around_bubble"] * init_data["omega"];
      size_t a = init_data["a"] * nodes_in_1_x;
      double b_div_by_a = init_data["b_div_by_a"];
      size_t b = static_cast<size_t>(a * b_div_by_a);
      for (int y = par.y_begin; y < par.y_end; ++y) {
        int y_diff = y - int(par.y0);
        size_t left_bubble_edge = par.x0;
        size_t right_bubble_edge = par.x0;
        if (std::abs(y_diff) <= b) {
          size_t bubble_x = static_cast<size_t>(
            std::sqrt(a*a - y_diff*y_diff/b_div_by_a/b_div_by_a));
          left_bubble_edge -= bubble_x;
          right_bubble_edge += bubble_x;
        }
        int y_for_r = y - par.y_begin;
        for (size_t x = par.x_begin; x < shock_wave_initial_x; ++x) {
          mesh[y][x] = new data_node_axisymm(
            rho_after_sw, p_after_sw, u_after_sw, v_after_sw, y_for_r * par.delta_y);
        }
        for (size_t x = shock_wave_initial_x; x < left_bubble_edge; ++x) {
          mesh[y][x] = new data_node_axisymm(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"], y_for_r * par.delta_y);
        }

        if (left_bubble_edge == par.x0 && right_bubble_edge == par.x0) {
          mesh[y][par.x0] = new data_node_axisymm(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"], y_for_r * par.delta_y);
        } else {
          for (size_t x = left_bubble_edge; x < right_bubble_edge + 1; ++x) {
            mesh[y][x] = new data_node_axisymm(
              rho_in_bubble, init_data["p_around_bubble"],
              init_data["u_around_bubble"], init_data["v_around_bubble"], y_for_r * par.delta_y);
          }
        }

        for (size_t x = right_bubble_edge + 1; x < par.x_end; ++x) {
          mesh[y][x] = new data_node_axisymm(
            init_data["rho_around_bubble"], init_data["p_around_bubble"],
            init_data["u_around_bubble"], init_data["v_around_bubble"], y_for_r * par.delta_y);
        }
      }
    }

  }

  mesh_and_methods_axisymm(
      std::shared_ptr<const mesh_and_common_methods> other_grid) :
    mesh_and_common_methods(other_grid->get_params_const_ref()){

    auto other_mesh = other_grid->get_mesh_const_ref();
    mesh.resize(other_mesh.size());
    for (size_t y = 0; y < mesh.size(); ++y) {
      mesh[y].resize(other_mesh[y].size());
      for (size_t x = 0; x < mesh[y].size(); ++x) {
        mesh[y][x] = new data_node_axisymm(other_mesh[y][x]);
      }
    }
  }


  void calc_values_and_fluxes() override;

  virtual void lax_friedrichs(
      std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) override {
    throw(std::runtime_error(
      "lax-fridrichs method for axis symmetrical case is not available"));
    return;
  };
  virtual void mac_cormack(
      std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) override;
  virtual void mac_cormack_predictor_step(
      std::shared_ptr<const mesh_and_common_methods> predictor_grid) override;
  virtual void mac_cormack_corrector_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid,
      std::shared_ptr<const mesh_and_common_methods> predictor_grid) override;
  virtual void mac_cormack_with_davis(
      std::shared_ptr<const mesh_and_common_methods> prev_grid) override;
  virtual void calc_davis_artificial_viscosity() override;

  virtual void boundary_conditions_for_bubble_near_wall() override {};
  virtual void boundary_conditions_for_bubble_near_wall_simmetry_on_bottom() override;
  virtual void boundary_conditions_default() override;
};

class statistics_manager {

public:
  statistics_manager(
      std::shared_ptr<mesh_and_common_methods> solver_ptr) :
    mesh_ref(solver_ptr->get_mesh_const_ref()),
    par_ref(solver_ptr->get_params_const_ref()) {

    pressure_sensors_on_wall.resize(par_ref.size_y);

  };

  const std::vector<std::list<std::pair<double, double>>>&
  get_pressure_sensors_on_wall() {
    return pressure_sensors_on_wall;
  }

  void update_pressure_sensors_on_wall(double t);
  void get_pressure_sensors_from_files(
    const std::string& sensors_folder,
    const calculation_params& par);

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

  void impulses_output (double t0, double p0,
    const std::string& output_folder);

private:
  const std::vector<std::vector<data_node_2d* >>& mesh_ref;
  const calculation_params& par_ref;

  std::vector<std::list<std::pair<double, double>>>
    pressure_sensors_on_wall; //<time, pressure> in a wall point

};


#endif // DATA_2D_H
