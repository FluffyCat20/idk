#include "data_2d.h"

int main() {

  std::ifstream input("config_2d.json");
  std::string output_folder;
  data_2d grid(input, output_folder);
  input.close();

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

    grid.output_in_wall_point_for_current_time(pressure_diag, current_t);

    if (corner_pt.rho > rho_peak) {
      rho_peak = corner_pt.rho;
      rho_peak_time = current_t;
    }

    if (corner_pt.p > p_peak) {
      p_peak = corner_pt.p;
      p_peak_time = current_t;
    }

    grid.output_on_symmetry_axis_for_current_time(rho_p_on_symmetry_axis, current_t);

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
