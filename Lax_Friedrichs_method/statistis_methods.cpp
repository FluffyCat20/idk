#include "data_2d.h"
#include "utils.hpp"

void statistics_manager::update_pressure_sensors_on_wall(double t) {
  size_t y = par_ref.y_begin;
  for (auto& sensor : pressure_sensors_on_wall) {
    sensor.emplace_back(t, mesh_ref[y][par_ref.x_end - 1]->p);
    y++;
  }
}

void statistics_manager::get_pressure_sensors_from_files(
    const std::string &sensors_folder,
    const calculation_params& par) {

  pressure_sensors_on_wall.resize(par.size_y);
  for (int i = 0; i < par.size_y; ++i) {
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

void statistics_manager::calc_pressure_after_reflected_shock_no_bubble(
    const double rho2, const double p2, const double u2,
    double& p3) {

  p3 = newton_solver(
    std::bind(F, std::placeholders::_1, rho2, p2, u2),
    std::bind(F_der, std::placeholders::_1, rho2, p2));

  //mb add: while p3 == p1, repeat with another newton solver initial state

}

std::vector<double> statistics_manager::calc_pressure_impulses_basic(
  double p0) const {

  double eps = 1e-6;

  //t0 - time of coming of initial shock on the wall without bubble
  //t1 - time of coming of initial shock on the wall with bubble
  std::vector<double> impulses(par_ref.size_y);
  for (size_t i = 0; i < pressure_sensors_on_wall.size(); ++i) {
    const auto& sensor = pressure_sensors_on_wall[i];
    double initial_p = sensor.front().second;
    std::list<std::pair<double, double>>::const_iterator it1 =
      sensor.cbegin();
    for (; it1 != sensor.cend() && std::abs(it1->second - initial_p) < eps; ++it1) {}
    for (auto it = it1; it != sensor.cend(); ++it) {
      impulses[i] += (it->second - p0) * (it->first - std::prev(it)->first);
    }
  }

  return impulses;
}

inline std::list<std::pair<double, double>>::const_iterator
find_sensor_max(const std::list<std::pair<double, double>>& sensor) {

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

std::vector<double> statistics_manager::
calc_pressure_impulses_max_peak_bfr_p0(
    double p0) const {

  std::vector<double> impulses(par_ref.size_y);
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

std::vector<double> statistics_manager::
calc_pressure_impulses_max_peak_bounded_by_local_min() const {

  double diff = 0.5;

  std::vector<double> impulses(par_ref.size_y);
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


void statistics_manager::impulses_output(
    double p0, const std::string& output_folder) {
  std::vector<double> basic_impulses =
    calc_pressure_impulses_basic(p0);
  std::ofstream basic_impulses_file(output_folder
    + "basic_impulses.dat");
  for (size_t i = 0; i < basic_impulses.size(); ++i) {
    basic_impulses_file << i * par_ref.delta_y << " " << basic_impulses[i] << std::endl;
  }
  basic_impulses_file.close();

  std::vector<double> max_peak_bfr_p0_impulses =
    calc_pressure_impulses_max_peak_bfr_p0(p0);
  std::ofstream max_peak_bfr_p0_file(output_folder + "max_peak_bfr_p0_impulses.dat");
  for (size_t i = 0; i < max_peak_bfr_p0_impulses.size(); ++i) {
    max_peak_bfr_p0_file << i * par_ref.delta_y << " " << max_peak_bfr_p0_impulses[i] << std::endl;
  }
  max_peak_bfr_p0_file.close();

  std::vector<double> max_peak_bounded_mins_impulses =
    calc_pressure_impulses_max_peak_bounded_by_local_min();
  std::ofstream max_peak_bounded_mins_file(output_folder +
    "max_peak_bounded_by_local_mins_impulses.dat");
  for (size_t i = 0; i < max_peak_bounded_mins_impulses.size(); ++i) {
    max_peak_bounded_mins_file << i * par_ref.delta_y << " "
      << max_peak_bounded_mins_impulses[i] << std::endl;
  }
  max_peak_bounded_mins_file.close();
}
