#include "calculation_info.hpp"

calculation_info::calculation_info(
    std::ifstream& config_file_path, std::string& output_path) {

  get_data_from_config(config_file_path, output_path);

  par.delta_x = (par.x_right - par.x_left)/par.size_x;
  par.delta_y = (par.y_top - par.y_bottom)/par.size_y;
  par.x0 = par.size_x / 2;
  par.y0 = par.size_y / 2;
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

  par.gap_btw_sw_and_bound = init_data["gap_btw_sw_and_bound"];

  par.gamma = config_data["gamma"];
  data_node_2d::gamma = par.gamma;

}
