#ifndef CALCULATION_INFO_HPP
#define CALCULATION_INFO_HPP
#pragma once

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <functional>
#include <fstream>

#include "data_nodes.hpp"

struct MethodInfo {
  MethodInfo() :
  method_name(""), method_number(-1), ghost_edges_width(0) {};

  MethodInfo(std::string m_name, int m_num, int gh_e_w) :
    method_name(m_name), method_number(m_num), ghost_edges_width(gh_e_w) {};

  std::string method_name;
  int method_number;
  int ghost_edges_width;
};

struct calculation_params {

  MethodInfo method_info;
  double x_left, x_right; //bounds of calculating area as coordinates
  double y_bottom, y_top;

  int x_begin, x_end, y_begin, y_end; //bounds of calculating area as nodes

  size_t size_x, size_y; //size of calculating area
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
  int coord_type_number = -1;

  calculation_info(
    std::ifstream& config_file_path, std::string& output_path);

private:
  const std::vector<MethodInfo> methods_available = {
    MethodInfo("lax_friedrichs", 0, 1), //name, method number, number of ghost nodes rows
    MethodInfo("mac_cormack", 1, 1),
    MethodInfo("mac_cormack+davis", 2, 2)
  };

  const std::map<std::string, int> coord_types_available = {
    {"cartesian", 0},
    {"axis_symm", 1}
  };


  /*void prepare_possible_nodes_map() { //must be called after all data getting from config!!!!
    coord_types.emplace("cartesian", [this]{return new data_2d<data_node_cartesian>(par, init_data, init_config);});
    coord_types.emplace("axis_symm", [this]{return new data_2d<data_node_axisymm>(par, init_data, init_config);});
  }*/
  void get_init_data_map(json &config_data, std::string key);
  void get_data_from_config(
      std::ifstream &input, std::string& output_folder);
};



#endif // CALCULATION_INFO_HPP
