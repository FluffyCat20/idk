#ifndef CALCULATION_INFO_HPP
#define CALCULATION_INFO_HPP
#pragma once

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <functional>

#include "data_2d.h"

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
  using func_t = std::function<std::unique_ptr<base_data_2d>()>;
  std::map<std::string, func_t> data_creator;
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

  void prepare_possible_nodes_map() { //must be called after all data getting from config!!!!
    data_creator.emplace("cartesian", [this]{return new data_2d<data_node_cartesian>(par, init_data, init_config);});
    data_creator.emplace("axis_symm", [this]{return new data_2d<data_node_axys_symm>(par, init_data, init_config);});
  }
  void get_init_data_map(json &config_data, std::string key);
  void get_data_from_config(
      std::ifstream &input, std::string& output_folder);


};

calculation_info::calculation_info(
    std::ifstream& config_file_path, std::string& output_path) {

  get_data_from_config(config_file_path, output_path);

  par.delta_x = (par.x_right - par.x_left)/par.size_x;
  par.delta_y = (par.y_top - par.y_bottom)/par.size_y;
  par.x0 = par.size_x / 2;
  par.y0 = par.size_y / 2;

  prepare_possible_nodes_map();
}

void calculation_info::get_init_data_map(
    json &config_data, std::string key) {
  if (config_data[key]["on"]) {
    init_config[key] = true;
    for (auto& val : config_data[key].get<json::object_t>()) {
      if (val.second.is_number()) {
        init_data.insert(val);
      }
    }
  }
}

void calculation_info::get_data_from_config(
    std::ifstream &input, std::string& output_folder) {
  json config_data;
  input >> config_data;

  output_folder = config_data["output_folder"];

  method_name = config_data["method_name"];
  auto method_it = methods_available.find(method_name);
  if (method_it != methods_available.end()) {
    par.method_number = method_it->second;
  } else {
    throw (std::invalid_argument("wrong config : unknown method name"));
  }

  auto coord_type_it = coord_types_available.find(config_data["coord_type"]);
  if (coord_type_it != coord_types_available.end()) {
    coord_type_number = coord_type_it->second;
  } else {
    throw (std::invalid_argument("wrong config : unknown coordinates type"));
  }

  par.x_left = config_data["bounds"]["x_left"];
  par.x_right = config_data["bounds"]["x_right"];
  par.y_bottom = config_data["bounds"]["y_bottom"];
  par.y_top = config_data["bounds"]["y_top"];

  for (const auto& pair : init_config) {
    get_init_data_map(config_data, pair.first);
  }

  par.size_x = config_data["size_x"];
  par.size_y = config_data["size_y"];
  par.t_end = config_data["t_end"];
  par.time_step = config_data["time_step"];

  par.courant_number = config_data["courant_number"];

  par.gamma = config_data["gamma"];
  data_node_2d::gamma = par.gamma;

}



#endif // CALCULATION_INFO_HPP
