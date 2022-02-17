#include "data_2d.h"

void data_2d::get_init_data_map(json &config_data, std::string key) {
  if (config_data[key]["on"]) {
    init_config[key] = true;
    for (auto& val : config_data[key].get<json::object_t>()) {
      if (val.second.is_number()) {
        init_data.insert(val);
      }
    }
  }
}

void data_2d::get_data_from_config(std::ifstream &input) {
  json config_data;
  input >> config_data;

  method_name = config_data["method_name"];

  x_left = config_data["bounds"]["x_left"];
  x_right = config_data["bounds"]["x_right"];
  y_bottom = config_data["bounds"]["y_bottom"];
  y_top = config_data["bounds"]["y_top"];

  get_init_data_map(config_data, "horizontal_flow");
  get_init_data_map(config_data, "vertical_flow");
  get_init_data_map(config_data, "shock_wave");
  get_init_data_map(config_data, "horizontal_left_contact_disc");

  size_x = config_data["size_x"];
  size_y = config_data["size_y"];
  t_end = config_data["t_end"];
  time_step = config_data["time_step"];

  courant_number = config_data["courant_number"];

  gamma = config_data["gamma"];
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

void data_2d::calc_values_and_fluxes() {
  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      mesh[i][j].calc_values_from_U();
      mesh[i][j].calc_F_when_values_known();
      mesh[i][j].calc_G_when_values_known();
    }
  }
}

void data_2d::boundary_conditions() {
  for (size_t y = 1; y < size_y - 1; ++y) {
    mesh[y][0] = mesh[y][1];
    mesh[y][size_x - 1] = mesh[y][size_x - 2];
  }
  for (size_t x = 1; x < size_x - 1; ++x) {
    mesh[0][x] = mesh[1][x];
    mesh[size_y - 1][x] = mesh[size_y - 2][x];
  }
}

void data_2d::lax_friedrichs(const data_2d &prev_grid) {
  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      auto& up = prev_grid.mesh[i-1][j];
      auto& down = prev_grid.mesh[i+1][j];
      auto& left = prev_grid.mesh[i][j-1];
      auto& right = prev_grid.mesh[i][j+1];
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j].U[k] = 0.25*(up.U[k] + down.U[k] + left.U[k] + right.U[k]) -
          0.5*delta_t/delta_x * (right.F[k] - left.F[k]) -
          0.5*delta_t/delta_y * (up.G[k] - down.G[k]);
      }
    }
  }
  calc_values_and_fluxes();
  boundary_conditions();
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

  /*for (size_t i = 0; i < size_y - 1; ++i) { //0?? not 1??
    for (size_t j = 0; j < size_x - 1; ++j) {
      data_node_2d& predictor_node = predictor_grid.mesh[i][j];
      predictor_node.calc_values_from_U();
      predictor_node.calc_F_when_values_known();
      predictor_node.calc_G_when_values_known();
    }
  }*/
  predictor_grid.calc_values_and_fluxes();
  predictor_grid.boundary_conditions();
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
}

void data_2d::mac_cormack(const data_2d& prev_grid) {
  data_2d predictor_grid(prev_grid); //TODO: copy common data only
  prev_grid.mac_cormack_predictor_step(predictor_grid);
  mac_cormack_corrector_step(prev_grid, predictor_grid);
  calc_values_and_fluxes();
  boundary_conditions();
}

double inner_product(const std::vector<double>& v1,
                     const std::vector<double>& v2){
  double res = 0;
  for (size_t i = 0; i < v1.size(); ++i){
    res += v1[i]*v2[i];
  }
  return res;
}

inline double phi(double rp, double rm) { //rp = r+, rm = r-
  return std::max ( {
         0.0,
         std::min({2*rp, rm, 1.0}),
         std::min({2*rm, rp, 1.0}) });
}

void data_2d::mac_cormack_with_davis(const data_2d &prev_grid) {
  data_2d predictor_grid(prev_grid); //TODO: copy common data only
  prev_grid.mac_cormack_predictor_step(predictor_grid);
  mac_cormack_corrector_step(prev_grid, predictor_grid);
  boundary_conditions();
  calc_davis_artificial_viscosity();
  calc_values_and_fluxes();
  boundary_conditions();
}

void data_2d::calc_davis_artificial_viscosity() {
  std::vector<std::vector<double>>
    r_x_plus(size_y, std::vector<double>(size_x)),
    r_x_minus(size_y, std::vector<double>(size_x)),
    r_y_plus(size_y, std::vector<double>(size_x)),
    r_y_minus(size_y, std::vector<double>(size_x));
  std::vector<std::vector<std::vector<double>>>
    delta_u_x(size_y, std::vector<std::vector<double>>(size_x, std::vector<double>(4))),
    delta_u_y(size_y, std::vector<std::vector<double>>(size_x, std::vector<double>(4))),
    D_x(size_y, std::vector<std::vector<double>>(size_x, std::vector<double>(4))),
    D_y(size_y, std::vector<std::vector<double>>(size_x, std::vector<double>(4)));

  //delta_u calculation:
  for (size_t i = 0; i < size_y - 1; ++i) {
    for (size_t j = 0; j < size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        delta_u_x[i][j][k] = mesh[i][j+1].U[k] - mesh[i][j].U[k];
        delta_u_y[i][j][k] = mesh[i+1][j].U[k] - mesh[i][j].U[k];
      }
    }
  }

  //r calculation: //TODO: fix: inner products are calculated twice (save squares?)
  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      //x-direction
      double product_mid_x = inner_product(delta_u_x[i][j-1], delta_u_x[i][j]);
      double product_left_x = inner_product(delta_u_x[i][j-1], delta_u_x[i][j-1]);
      double product_right_x = inner_product(delta_u_x[i][j], delta_u_x[i][j]);

      r_x_plus[i][j] = product_mid_x/product_right_x;
      r_x_minus[i][j] = product_mid_x/product_left_x;

      //y_direction
      double product_mid_y = inner_product(delta_u_y[i-1][j], delta_u_y[i][j]);
      double product_up_y = inner_product(delta_u_y[i-1][j], delta_u_y[i-1][j]);
      double product_down_y = inner_product(delta_u_y[i][j], delta_u_y[i][j]);

      r_y_plus[i][j] = product_mid_y/product_down_y;
      r_y_minus[i][j] = product_mid_y/product_up_y;
    }
  }

  //nu, K, D calculation:
  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      double lambda = std::max(std::abs(mesh[i][j].u_abs + mesh[i][j].a),
                               std::abs(mesh[i][j].u_abs - mesh[i][j].a));
      //nu_x, k_x
      double nu_x = delta_t/delta_x*lambda;
      double c_nu_x = nu_x*(1-nu_x);
      if (nu_x > 0.5) {
        c_nu_x = 0.25;
      }
      double k_x = 0.5*c_nu_x*(1 - phi(r_x_plus[i][j], r_x_minus[i][j+1]));

      //nu_y, k_y
      double nu_y = delta_t/delta_y*lambda;
      double c_nu_y = nu_y*(1-nu_y);
      if (nu_y > 0.5) {
        c_nu_y = 0.25;
      }
      double k_y = 0.5*c_nu_y*(1 - phi(r_y_plus[i][j], r_y_minus[i+1][j]));

      //D
      for (size_t k = 0; k < 4; ++k) {
        D_x[i][j][k] = k_x*delta_u_x[i][j][k];
        D_y[i][j][k] = k_y*delta_u_y[i][j][k];
      }

    }
  }

  for (size_t k = 0; k < 4; ++k) {
    for (size_t i = 1; i < size_y - 1; ++i){
      D_x[i][0][k] = 0.0;
      D_x[i][size_x - 1][k] = 0.0;
      D_y[i][0][k] = 0.0;
      D_y[i][size_x - 1][k] = 0.0;
    }
    for (size_t j = 1; j < size_x - 1; ++j){
      D_x[0][j][k] = 0.0;
      D_x[size_y - 1][j][k] = 0.0;
      D_y[0][j][k] = 0.0;
      D_y[size_y - 1][j][k] = 0.0;
    }
  }

  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j].U[k] += D_x[i][j][k] - D_x[i][j-1][k] + D_y[i][j][k] - D_y[i-1][j][k];
      }
    }
  }

}

void data_2d::output_first(std::ofstream &outfile) {
  outfile << "TITLE = \"" << method_name <<
             "; size: " << size_x << "x" << size_y << "\""<< std::endl;
  outfile << "VARIABLES = \"x\", \"y\", \"rho\", \"p\", \"u\", \"v\"" << std::endl;
  outfile << std::setprecision(4);
  output_for_current_time(outfile, 0.0);
}

void data_2d::output_for_current_time(std::ofstream &outfile, double time) {
  outfile << "ZONE T=\"0.0000\", SOLUTIONTIME=" << time << std::endl;
  outfile << "I=" << size_y << ", J=" << size_x << ", F=POINT" << std::endl;

  for (size_t i = 0; i < size_x; ++i) {
    for (size_t j = 0; j < size_y; ++j) {
      const data_node_2d& pt = mesh[j][i];
      outfile << i*delta_x << " " << j*delta_y << " " <<
        pt.rho << " " << pt.p << " " << pt.u << " " << pt.v << std::endl;
     }
  }
}

void data_2d::shock_wave_initialization (
    double& rho, double& p, double& u,
    double rho1, double p1, double& u1) { //1 - before sw, no ind - after sw

  double& mach = init_data["mach"];
  double a1 = std::sqrt(gamma*p1/rho1);
  u1 = mach*a1;
  p = p1*(2*gamma/(gamma+1)*mach*mach - (gamma-1)/(gamma+1));
  rho = 1/(rho1*((gamma-1)/(gamma+1)) + 2/(gamma+1)/mach/mach);
  u = u1*rho1/rho;

}

data_2d& data_2d::operator=(const data_2d& other) {
  mesh = other.mesh;
  delta_t = other.delta_t;

  return *this;
}
