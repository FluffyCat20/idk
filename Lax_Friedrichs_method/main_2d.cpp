#include "data_2d_writer.h"

//Someday this shit will be changed...
#define TYPE_NUMBER(type_number) type_number
#if TYPE_NUMBER == 0
  #define NODE_TYPE data_node_cartesian
#elif TYPE_NUMBER == 1
  #define NODE_TYPE data_node_axis_symm
#else
  #define NODE_TYPE error
#endif

double data_node_2d::gamma;

template <typename node_t>
void impulses_output (const data_2d<node_t>& grid,
                      double t0, double p0,
                      const std::string& output_folder) {
  std::vector<double> basic_impulses =
    grid.calc_pressure_impulses_basic(t0, p0);
  std::ofstream basic_impulses_file(output_folder
    + "basic_impulses.dat");
  for (size_t i = 0; i < basic_impulses.size(); ++i) {
    basic_impulses_file << i * grid.par.delta_y << " " << basic_impulses[i] << std::endl;
  }
  basic_impulses_file.close();

  std::vector<double> max_peak_bfr_p0_impulses =
    grid.calc_pressure_impulses_max_peak_bfr_p0(p0);
  std::ofstream max_peak_bfr_p0_file(output_folder + "max_peak_bfr_p0_impulses.dat");
  for (size_t i = 0; i < max_peak_bfr_p0_impulses.size(); ++i) {
    max_peak_bfr_p0_file << i * grid.par.delta_y << " " << max_peak_bfr_p0_impulses[i] << std::endl;
  }
  max_peak_bfr_p0_file.close();

  std::vector<double> max_peak_bounded_mins_impulses =
    grid.calc_pressure_impulses_max_peak_bounded_by_local_min();
  std::ofstream max_peak_bounded_mins_file(output_folder +
    "max_peak_bounded_by_local_mins_impulses.dat");
  for (size_t i = 0; i < max_peak_bounded_mins_impulses.size(); ++i) {
    max_peak_bounded_mins_file << i * grid.par.delta_y << " "
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
  calculation_info calc_info(input, output_folder);
  input.close();

  data_2d<NODE_TYPE> grid(
    calc_info.par, calc_info.init_data, calc_info.init_config);

  data_2d_writer<NODE_TYPE> grid_writer(grid, output_folder);

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

  grid_writer.output_first();

  std::ofstream pressure_diag(output_folder + "pressure_diag" +
    std::to_string(grid.par.size_x)
    + "x" + std::to_string(grid.par.size_y) + ".dat");
  grid_writer.output_in_wall_point_first(pressure_diag);
  grid_writer.output_on_symmetry_axis_first();

  std::ofstream peaks(output_folder + "peaks.dat");

  std::cout << std::scientific;

  data_2d<NODE_TYPE> new_grid(grid); //TODO: copy common data only
  int counter = 0;
  double current_t = 0.0;  
  std::vector<double> times;
  for (int i = 0; i*grid.par.time_step < grid.par.t_end; ++i) {
    times.push_back(i*grid.par.time_step);
  }
  double current_time_idx = 0;

  //peaks
  double rho_peak = 0.0, p_peak = 0.0;
  double rho_peak_time = 0.0, p_peak_time = 0.0;
  const auto& corner_pt = grid.mesh[1][grid.par.size_x - 2];

  while (current_t < new_grid.par.t_end) {

    new_grid.par.delta_t = grid.calc_delta_t();
    if (current_t + new_grid.par.delta_t > new_grid.par.t_end){
      new_grid.par.delta_t = new_grid.par.t_end - current_t;
    }
    grid.par.delta_t = new_grid.par.delta_t;
    new_grid.do_one_step(grid);
    grid = new_grid;
    current_t += new_grid.par.delta_t;

    grid.update_pressure_sensors_on_wall(current_t);
    grid_writer.output_in_wall_point_for_current_time(pressure_diag, current_t);

    if (corner_pt.rho > rho_peak) {
      rho_peak = corner_pt.rho;
      rho_peak_time = current_t;
    }

    if (corner_pt.p > p_peak) {
      p_peak = corner_pt.p;
      p_peak_time = current_t;
    }

    grid_writer.output_on_symmetry_axis_for_current_time(current_t);

    if (current_time_idx + 1 < times.size() &&
        times[current_time_idx + 1] <= current_t) {
      grid_writer.output_for_current_time(current_t);
      ++current_time_idx;
    }

    ++counter;
    if (new_grid.stop_now) {
      std::cout << "stop time:" << current_t << " steps: " << counter << std::endl;
      break;
    }

  }

  std::cout << current_t << " steps: " << counter << std::endl;

  grid_writer.output_for_current_time(current_t);


  //reflected shock with no bubble
  double p3;
  grid.calc_pressure_after_reflected_shock_no_bubble(
    grid.mesh[0][0].rho, grid.mesh[0][0].p, grid.mesh[0][0].u, p3);
  double t0 = (grid.par.x_right - grid.par.x_left -
    grid.par.gap_btw_sw_and_bound) / grid.par.D_of_initial_shock;

  //sensors & impulses
  grid_writer.output_pressure_sensors_on_wall(output_folder);
  impulses_output(grid, t0, p3, output_folder);

  //peaks:

  peaks << "max rho in corner point:" << std::endl <<
    "time = " << rho_peak_time << " rho = " << rho_peak << std::endl;

  peaks << "max p in corner point:" << std::endl <<
    "time = " << p_peak_time << " p = " << p_peak << std::endl;

  pressure_diag.close();
  peaks.close();

  return 0;
}
