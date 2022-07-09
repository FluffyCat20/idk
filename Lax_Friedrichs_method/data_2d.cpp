#include "utils.hpp"

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

void data_2d::get_data_from_config(std::ifstream &input, std::string& output_folder) {
  json config_data;
  input >> config_data;

  output_folder = config_data["output_folder"];

  method_name = config_data["method_name"];

  x_left = config_data["bounds"]["x_left"];
  x_right = config_data["bounds"]["x_right"];
  y_bottom = config_data["bounds"]["y_bottom"];
  y_top = config_data["bounds"]["y_top"];

  for (const auto& pair : init_config) {
    get_init_data_map(config_data, pair.first);
  }

  size_x = config_data["size_x"];
  size_y = config_data["size_y"];
  t_end = config_data["t_end"];
  time_step = config_data["time_step"];

  courant_number = config_data["courant_number"];

  gamma = config_data["gamma"];
  data_node_2d::gamma = gamma;
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

  //if (init_config["bubble_near_wall"]) {
    //boundary_conditions_for_bubble_near_wall();
    boundary_conditions_for_bubble_near_wall_simmetry_on_bottom();
  //} else {
  //  boundary_conditions_default();
  //}
}

void data_2d::boundary_conditions_default() {
  for (size_t y = 1; y < size_y - 1; ++y) {
    mesh[y][0] = mesh[y][1];
    mesh[y][size_x - 1] = mesh[y][size_x - 2];
  }
  for (size_t x = 1; x < size_x - 1; ++x) {
    mesh[0][x] = mesh[1][x];
    mesh[size_y - 1][x] = mesh[size_y - 2][x];
  }
}

void data_2d::boundary_conditions_for_bubble_near_wall() {
  // solid wall on the right:
  for (size_t y = 1; y < size_y - 1; ++y) {
    data_node_2d& node_to_change = mesh[y][size_x - 1];
    node_to_change = mesh[y][size_x - 2];
    node_to_change.u *= -1;
    node_to_change.U[1] *= -1;
    node_to_change.F[0] *= -1;
    node_to_change.F[2] *= -1;
    node_to_change.F[3] *= -1;
    node_to_change.G[1] *= -1;
  }

  // d/dn = 0 on other bounds:
  for (size_t y = 1; y < size_y - 1; ++y) {
    mesh[y][0] = mesh[y][1];
  }
  for (size_t x = 1; x < size_x - 1; ++x) {
    mesh[0][x] = mesh[1][x];
    mesh[size_y - 1][x] = mesh[size_y - 2][x];
  }
}

void data_2d::boundary_conditions_for_bubble_near_wall_simmetry_on_bottom() {
  // solid wall on the right:
  for (size_t y = 1; y < size_y - 1; ++y) {
    data_node_2d& node_to_change = mesh[y][size_x - 1];
    node_to_change = mesh[y][size_x - 2];
    node_to_change.u *= -1;
    node_to_change.U[1] *= -1;
    node_to_change.F[0] *= -1;
    node_to_change.F[2] *= -1;
    node_to_change.F[3] *= -1;
    node_to_change.G[1] *= -1;
  }

  // d/dn = 0 on left&right bounds and symmetry on bottom:
  for (size_t y = 1; y < size_y - 1; ++y) {
    mesh[y][0] = mesh[y][1];
  }
  for (size_t x = 1; x < size_x - 1; ++x) {
    data_node_2d& node_to_change = mesh[0][x];
    node_to_change = mesh[1][x];
    node_to_change.v *= -1;
    node_to_change.U[2] *= -1;
    node_to_change.F[2] *= -1;
    node_to_change.G[0] *= -1;
    node_to_change.G[1] *= -1;
    node_to_change.G[3] *= -1;

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
      outfile << i*delta_x + x_left << " " << j*delta_y + y_bottom << " " <<
        pt.rho << " " << pt.p << " " << pt.u << " " << pt.v << std::endl;
    }
  }

  std::cout << "Calculations done: time = " << time << std::endl;

}

void data_2d::output_in_wall_point_first(std::ofstream &outfile) {
  outfile << "TITLE = \" in central point \"" << std::endl;
  outfile << "VARIABLES = \"t\", \"rho\", \"p\"" << std::endl;
  outfile << std::setprecision(4);
  output_in_wall_point_for_current_time(outfile, 0.0);
}

void data_2d::output_in_wall_point_for_current_time(std::ofstream &outfile, double time) {
  outfile << time << " " << mesh[1][size_x - 2].rho << " " << mesh[1][size_x - 2].p << std::endl;
}

void data_2d::output_on_symmetry_axis_first(std::ofstream &outfile) {
  outfile << "TITLE = \" on symmetry axis \"" << std::endl;
  outfile << "VARIABLES = \"x\", \"rho\", \"p\"" << std::endl;
  outfile << std::setprecision(4);
  output_on_symmetry_axis_for_current_time(outfile, 0.0);
}

void data_2d::output_on_symmetry_axis_for_current_time(std::ofstream &outfile, double time) {
  outfile << "ZONE I = " << size_x << ", F=POINT, SOLUTIONTIME = "
    << time << " t = " << "\"" << time << "\"" << std::endl;
  for (size_t x = 0; x < size_x; ++x) {
    const data_node_2d& pt = mesh[0][x];
    outfile << x*delta_x + x_left << " " << pt.rho << " " << pt.p << std::endl;
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

void data_2d::update_pressure_sensors_on_wall(double t) {
  size_t y = 0;
  for (auto& sensor : pressure_sensors_on_wall) {
    sensor.emplace_back(t, mesh[y][size_x - 2].p);
    y++;
  }
}

void data_2d::output_pressure_sensors_on_wall(
    const std::string &output_folder) {

  for (size_t i = 0; i < pressure_sensors_on_wall.size(); ++i) {
    std::ofstream fout(output_folder +
      "pressure_on_wall_sensors_norm/pressure_sensor" + std::to_string(i) + ".dat");
    for (const auto& it : pressure_sensors_on_wall[i]) {
      fout << it.first << " " << it.second << std::endl;
    }
    fout.close();
  }
}

void data_2d::get_pressure_sensors_from_files(
    const std::string &sensors_folder) {

  pressure_sensors_on_wall.resize(size_y);
  for (int i = 0; i < size_y; ++i) {
    std::ifstream fin(sensors_folder + "pressure_sensor" +
      std::to_string(i) + ".dat");
    double t, p;
    while (fin >> t >> p) {
      pressure_sensors_on_wall[i].emplace_back(t, p);
    }
    fin.close();
  }

}

inline double F(double p, double rho2, double p2, double u2) {
  //for 1D shock reflection from a wall (for newton solver)
  double gamma = data_node_2d::gamma;

  double res = u2 * u2 - 2.0 * (p - p2) * (p - p2)
    / rho2 / (p2 * (gamma - 1) + p * (gamma + 1));

  return res;
}

inline double F_der(double p, double rho2, double p2) {
  //for 1D shock reflection from a wall (for newton solver)
  double gamma = data_node_2d::gamma;

  double denom = p2 * (gamma - 1) + p * (gamma + 1);

  return 2.0 * (p2 - p) * (p * (gamma + 2) + p2 * (3 * gamma + 2))
      / rho2 / denom / denom;
}

void data_2d::calc_pressure_after_reflected_shock_no_bubble(
    const double rho2, const double p2, const double u2,
    double& p3) {

  p3 = newton_solver(
    std::bind(F, std::placeholders::_1, rho2, p2, u2),
    std::bind(F_der, std::placeholders::_1, rho2, p2));

  //mb add: while p3 == p1, repeat with another newton solver initial state

}

std::vector<double> data_2d::calc_pressure_impulses_basic(
    double t0, double p0) const {

  double eps = 1e-6;

  //t0 - time of coming of initial shock on the wall without bubble
  //t1 - time of coming of initial shock on the wall with bubble
  std::vector<double> impulses(size_y);
  for (size_t i = 0; i < pressure_sensors_on_wall.size(); ++i) {
    const auto& sensor = pressure_sensors_on_wall[i];
    double initial_p = sensor.front().second;
    std::list<std::pair<double, double>>::const_iterator it0 =
      sensor.cbegin();
    std::list<std::pair<double, double>>::const_iterator it1 =
      sensor.cbegin();
    for (; it0 != sensor.cend() && it0->first < t0; ++it0) {}
    for (; it1 != sensor.cend() && std::abs(it1->second - initial_p) < eps; ++it1) {}
    double t1 = t0;
    if (it1 != sensor.cend()) {
      t1 = it1->first;
    }
    for (auto it = t0 < t1 ? it0 : it1; it != sensor.cend(); ++it) {
      //it must not be == sensor.cbegin()
      impulses[i] += (it->second - p0) * (it->first - std::prev(it)->first);
    }
  }

  return impulses;
}

inline std::list<std::pair<double, double>>::const_iterator find_sensor_max(
    const std::list<std::pair<double, double>>& sensor) {

  auto it = sensor.cbegin();
  auto res = it;
  double max = it->second;
  for (; it != sensor.cend(); ++it) {
    if (it->second > max) {
      res = it;
      max = it->second;
    }
  }
  return res;
}

std::vector<double> data_2d::calc_pressure_impulses_max_peak_bfr_p0(
    double p0) const {

  std::vector<double> impulses(size_y);
  for (size_t i = 0; i < pressure_sensors_on_wall.size(); ++i) {
    const auto& sensor = pressure_sensors_on_wall[i];
    auto max_it = find_sensor_max(sensor);
    auto left_peak_bound = max_it;
    left_peak_bound --;
    auto right_peak_bound = max_it;
    while (left_peak_bound != sensor.cbegin() && left_peak_bound->second > p0) {
      double t_next = std::next(left_peak_bound)->first;
      impulses[i] += (left_peak_bound->second - p0) * (
        t_next - left_peak_bound->first);
      left_peak_bound--;
    }
    while (right_peak_bound != sensor.cend() && right_peak_bound->second > p0) {
      double t_prev = std::prev(right_peak_bound)->first;
      impulses[i] += (right_peak_bound->second - p0) * (
        right_peak_bound->first - t_prev);
      right_peak_bound++;
    }
  }
  return impulses;
}

std::vector<double>
data_2d::calc_pressure_impulses_max_peak_bounded_by_local_min() const {

  double diff = 0.5;

  std::vector<double> impulses(size_y);
  for (size_t i = 0; i < pressure_sensors_on_wall.size(); ++i) {
    const auto& sensor = pressure_sensors_on_wall[i];
    auto max_it = find_sensor_max(sensor);
    auto left_peak_bound = max_it;
    left_peak_bound --;
    auto right_peak_bound = max_it;
    while (left_peak_bound != sensor.cbegin() &&
           left_peak_bound->second > std::prev(left_peak_bound)->second - diff) {
      double t_next = std::next(left_peak_bound)->first;
      impulses[i] += left_peak_bound->second * (t_next - left_peak_bound->first);
      left_peak_bound--;
    }
    while (std::next(right_peak_bound) != sensor.cend() &&
           right_peak_bound->second > std::next(right_peak_bound)->second - diff) {
      double t_prev = std::prev(right_peak_bound)->first;
      impulses[i] += right_peak_bound->second * (right_peak_bound->first - t_prev);
      right_peak_bound++;
    }
  }
  return impulses;
}

