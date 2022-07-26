#ifndef CALCULATION_INFO_HPP
#define CALCULATION_INFO_HPP
#pragma once

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <functional>
#include <fstream>

#include "data_nodes.hpp"

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

  bool stop_now = false;

};

class calculation_info {

public:
  calculation_params par;
  std::unordered_map<std::string, double> init_data;
  std::unordered_map<std::string, bool> init_config {
    {"horizontal_flow", false},
    {"vertical_flow", false},
    {"shock_wave", false},
    {"horizontal_left_contact_disc", false},
    {"quadrants", false},
    {"bubble_near_wall", false}
  };
  /*using func_t = std::function<std::unique_ptr<base_data_2d>()>;
  std::map<std::string, func_t> coord_types;*/
  std::string method_name = "";
  int coord_type_number = -1;

  calculation_info(
    std::ifstream& config_file_path, std::string& output_path);

private:
  std::map<std::string, int> methods_available = {
    {"lax_friedrichs", 0},
    {"mac_cormack", 1},
    {"mac_cormack+davis", 2}
  };

  std::map<std::string, int> coord_types_available = {
    {"cartesian", 0},
    {"axis_symm", 1}
  };


  /*void prepare_possible_nodes_map() { //must be called after all data getting from config!!!!
    coord_types.emplace("cartesian", [this]{return new data_2d<data_node_cartesian>(par, init_data, init_config);});
    coord_types.emplace("axis_symm", [this]{return new data_2d<data_node_axys_symm>(par, init_data, init_config);});
  }*/
  void get_init_data_map(json &config_data, std::string key);
  void get_data_from_config(
      std::ifstream &input, std::string& output_folder);
};



#endif // CALCULATION_INFO_HPP
