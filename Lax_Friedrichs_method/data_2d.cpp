#include "data_2d.h"

void data_2d::get_data_from_config(std::ifstream &input) {
  json config_data;
  input >> config_data;

  method_name = config_data["method_name"];

  x_left = config_data["bounds"]["x_left"];
  x_right = config_data["bounds"]["x_right"];
  y_bottom = config_data["bounds"]["y_bottom"];
  y_top = config_data["bounds"]["y_top"];

  rho_left = config_data["rho_left"];
  rho_right = config_data["rho_right"];
  p_left = config_data["p_left"];
  p_right = config_data["p_right"];
  u_left = config_data["u_left"];
  u_right = config_data["u_right"];
  v_left = config_data["v_left"];
  v_right = config_data["v_right"];

  size_x = config_data["size_x"];
  size_y = config_data["size_y"];
  t_end = config_data["t_end"];

  courant_number = config_data["courant_number"];

  gamma = config_data["gamma"];
  mach = config_data["Mach"];
  omega = config_data["omega"];
}

double data_2d::calc_delta_t() {
  double max_u_abs_plus_a = std::numeric_limits<double>::min();
  for (const auto& row : mesh) {
    for (const data_node_2d& node : row) {
      double candidat = node.u_abs + node.a;
      if (candidat > max_u_abs_plus_a) {
        max_u_abs_plus_a = candidat;
      }
    }
  }
  return delta_x/max_u_abs_plus_a*courant_number;
}

void data_2d::boundary_conditions() {
  for (size_t y = 1; y < size_y - 1; ++y) {
    mesh[y][0].U = mesh[y][1].U;
    mesh[y][size_x - 1] = mesh[y][size_x - 2];
  }
  for (size_t x = 1; x < size_x - 1; ++x) {
    mesh[0][x] = mesh[1][x];
    mesh[size_y - 1][x] = mesh[size_y - 2][x];
  }
}

void data_2d::mac_cormack_predictor_step(data_2d& predictor_grid)const{
  for (size_t i = 0; i < size_y - 1; ++i) { //0?? not 1??
    for (size_t j = 0; j < size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        predictor_grid.mesh[i][j].U[k] = mesh[i][j].U[k]
          - delta_t*((mesh[i][j+1].F[k] - mesh[i][j].F[k])/delta_x
                      + (mesh[i+1][j].G[k] - mesh[i][j].G[k])/delta_y);
      }
    }
  }

  for (size_t i = 0; i < size_y - 1; ++i) { //0?? not 1??
    for (size_t j = 0; j < size_x - 1; ++j) {
      data_node_2d& predictor_node = predictor_grid.mesh[i][j];
      predictor_node.calc_values_from_U();
      predictor_node.calc_F_when_values_known();
      predictor_node.calc_G_when_values_known();
    }
  }
}

void data_2d::mac_cormack_corrector_step(
    const data_2d &prev_grid, const data_2d &predictor_grid) {
  for (size_t i = 1; i < size_y - 1; ++i) { //size?? not size - 1??
    for (size_t j = 1; j < size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j].U[k] = 0.5*(prev_grid.mesh[i][j].U[k] + predictor_grid.mesh[i][j].U[k])
          - 0.5*delta_t*((predictor_grid.mesh[i][j].F[k] - predictor_grid.mesh[i][j - 1].F[k])/delta_x
                          + (predictor_grid.mesh[i][j].G[k] - predictor_grid.mesh[i - 1][j].G[k])/delta_y);
      }
    }
  }

  for (size_t i = 1; i < size_y - 1; ++i) { //size?? not size - 1??
    for (size_t j = 1; j < size_x - 1; ++j) {
      mesh[i][j].calc_values_from_U();
      mesh[i][j].calc_F_when_values_known();
      mesh[i][j].calc_G_when_values_known();
    }
  }
}

void data_2d::mac_cormack(const data_2d& prev_grid) {
  data_2d predictor_grid(prev_grid); //TODO: copy common data only
  prev_grid.mac_cormack_predictor_step(predictor_grid);
  mac_cormack_corrector_step(prev_grid, predictor_grid);
  boundary_conditions();
}

void data_2d::write_out_file(std::ofstream &outfile) {
  outfile << "TITLE = \"" << method_name <<
             "; size: " << size_x << "x" << size_y << "\""<< std::endl;
  outfile << "VARIABLES = \"x\", \"y\", \"rho\", \"p\", \"u\", \"v\"" << std::endl;
  outfile << "ZONE T=\"0.0000\", SOLUTIONTIME=0.0000" << std::endl;
  outfile << "I=" << size_y << ", J=" << size_x << ", F=POINT" << std::endl;
  outfile << std::setprecision(4);

  for (size_t i = 0; i < size_x; ++i) {
    for (size_t j = 0; j < size_y; ++j) {
      const data_node_2d& pt = mesh[i][j];
      outfile << j*delta_x << " " << i*delta_y << " " <<
        pt.rho << " " << pt.p << " " << pt.u << " " << pt.v << std::endl;
     }
  }
}

void data_2d::shock_wave_initialization() {

  double a_left = std::sqrt(gamma*p_left/rho_left);
  u_left = mach*a_left;
  p_right = p_left*(2*gamma/(gamma+1)*mach*mach - (gamma-1)/(gamma+1));
  rho_right = 1/(rho_left*((gamma-1)/(gamma+1)) + 2/(gamma+1)/mach/mach);
  u_right = u_left*rho_left/rho_right;

}

data_2d& data_2d::operator=(const data_2d& other) {
  mesh = other.mesh;
  delta_t = other.delta_t;

  return *this;
}
