#include "data_2d.h"

double data_node_2d::gamma;

void impulses_output (const data_2d& grid,
                      double t0, double p0,
                      const std::string& output_folder) {
  std::vector<double> basic_impulses =
    grid.calc_pressure_impulses_basic(t0, p0);
  std::ofstream basic_impulses_file(output_folder
    + "basic_impulses.dat");
  for (size_t i = 0; i < basic_impulses.size(); ++i) {
    basic_impulses_file << i * grid.delta_y << " " << basic_impulses[i] << std::endl;
  }
  basic_impulses_file.close();

  std::vector<double> max_peak_bfr_p0_impulses =
    grid.calc_pressure_impulses_max_peak_bfr_p0(p0);
  std::ofstream max_peak_bfr_p0_file(output_folder + "max_peak_bfr_p0_impulses.dat");
  for (size_t i = 0; i < max_peak_bfr_p0_impulses.size(); ++i) {
    max_peak_bfr_p0_file << i * grid.delta_y << " " << max_peak_bfr_p0_impulses[i] << std::endl;
  }
  max_peak_bfr_p0_file.close();

  std::vector<double> max_peak_bounded_mins_impulses =
    grid.calc_pressure_impulses_max_peak_bounded_by_local_min();
  std::ofstream max_peak_bounded_mins_file(output_folder +
    "max_peak_bounded_by_local_mins_impulses.dat");
  for (size_t i = 0; i < max_peak_bounded_mins_impulses.size(); ++i) {
    max_peak_bounded_mins_file << i * grid.delta_y << " "
      << max_peak_bounded_mins_impulses[i] << std::endl;
  }
  max_peak_bounded_mins_file.close();
}

int main(int argc, char* argv[]) {

  if (argc < 3 || strcmp(argv[1], "-c")) {
    throw std::runtime_error(
      "Usage: -c config_path");
  }

  std::ifstream input(argv[2]);
  if (input.fail()) {
    throw std::runtime_error("config not found");
  }
  std::string output_folder;
  data_2d grid(input, output_folder);
  input.close();

  /*grid.get_pressure_sensors_from_files(
    output_folder + "pressure_on_wall_sensors/");

  //reflected shock with no bubble
  double p3;
  grid.calc_pressure_after_reflected_shock_no_bubble(
    grid.mesh[0][0].rho, grid.mesh[0][0].p, grid.mesh[0][0].u, p3);
  double t0 = (grid.x_right - grid.x_left -
    grid.init_data["gap_btw_sw_and_bound"]) / grid.D_of_initial_shock;

  //sensors & impulses
  grid.output_pressure_sensors_on_wall(output_folder);
  impulses_output(grid, t0, p3, output_folder);

  return 0;*/





  std::ofstream outfile(output_folder + grid.method_name +
    std::to_string(grid.size_x)
    + "x" + std::to_string(grid.size_y) + ".dat");
  grid.output_first(outfile);

  std::ofstream pressure_diag(output_folder + "pressure_diag" +
    std::to_string(grid.size_x)
    + "x" + std::to_string(grid.size_y) + ".dat");
  grid.output_in_wall_point_first(pressure_diag);

  std::ofstream rho_p_on_symmetry_axis(output_folder + "rho_p_on_symmetry_axis" +
    std::to_string(grid.size_x) + "x" + std::to_string(grid.size_y) + ".dat");
  grid.output_on_symmetry_axis_first(rho_p_on_symmetry_axis);

  std::ofstream peaks(output_folder + "peaks.dat");

  int method_number = -1;

  if (grid.method_name == "lax_friedrichs")
    method_number = 0;

  if (grid.method_name == "mac_cormack")
    method_number = 1;

  if (grid.method_name == "mac_cormack+davis")
    method_number = 2;

  if (method_number == -1){
    std::cout << "Error: unknown method name" << std::endl;
    return -1;
  }

  std::cout << std::scientific;

  data_2d new_grid(grid); //TODO: copy common data only
  int counter = 0;
  double current_t = 0.0;  
  std::vector<double> times;
  for (int i = 0; i*grid.time_step < grid.t_end; ++i) {
    times.push_back(i*grid.time_step);
  }
  double current_time_idx = 0;

  //peaks
  double rho_peak = 0.0, p_peak = 0.0;
  double rho_peak_time = 0.0, p_peak_time = 0.0;
  const auto& corner_pt = grid.mesh[1][grid.size_x - 2];

  while (current_t < new_grid.t_end) {

    new_grid.delta_t = grid.calc_delta_t();
    if (current_t + new_grid.delta_t > new_grid.t_end){
      new_grid.delta_t = new_grid.t_end - current_t;
    }
    grid.delta_t = new_grid.delta_t;

    switch (method_number) {

    case 0: {
      new_grid.lax_friedrichs(grid);
      break;
    }

    case 1: {
      new_grid.mac_cormack(grid);
      break;
    }
    case 2: {
      new_grid.mac_cormack_with_davis(grid);
      break;
    }
    default: {
      std::cout << "????" << std::endl;
      break;
    }
    }

    grid = new_grid;
    current_t += new_grid.delta_t;

    grid.update_pressure_sensors_on_wall(current_t);

    grid.output_in_wall_point_for_current_time(pressure_diag, current_t);

    if (corner_pt.rho > rho_peak) {
      rho_peak = corner_pt.rho;
      rho_peak_time = current_t;
    }

    if (corner_pt.p > p_peak) {
      p_peak = corner_pt.p;
      p_peak_time = current_t;
    }

    grid.output_on_symmetry_axis_for_current_time(
      rho_p_on_symmetry_axis, current_t);

    if (current_time_idx + 1 < times.size() &&
        times[current_time_idx + 1] <= current_t) {
      grid.output_for_current_time(outfile, current_t);
      ++current_time_idx;
    }

    ++counter;
    if (new_grid.stop_now) {
      std::cout << "stop time:" << current_t << " steps: " << counter << std::endl;
      break;
    }

  }

  std::cout << current_t << " steps: " << counter << std::endl;

  grid.output_for_current_time(outfile, current_t);


  //reflected shock with no bubble
  double p3;
  grid.calc_pressure_after_reflected_shock_no_bubble(
    grid.mesh[0][0].rho, grid.mesh[0][0].p, grid.mesh[0][0].u, p3);
  double t0 = (grid.x_right - grid.x_left -
    grid.init_data["gap_btw_sw_and_bound"]) / grid.D_of_initial_shock;

  //sensors & impulses
  grid.output_pressure_sensors_on_wall(output_folder);
  impulses_output(grid, t0, p3, output_folder);

  //peaks:

  peaks << "max rho in corner point:" << std::endl <<
    "time = " << rho_peak_time << " rho = " << rho_peak << std::endl;

  peaks << "max p in corner point:" << std::endl <<
    "time = " << p_peak_time << " p = " << p_peak << std::endl;

  outfile.close();
  pressure_diag.close();
  rho_p_on_symmetry_axis.close();
  peaks.close();

  return 0;
}
