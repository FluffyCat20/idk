#include "data_2d_writer.h"

double data_node_2d::gamma;

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

  std::shared_ptr<mesh_and_common_methods> grid_ptr;
  //TODO: refactor with map of functions
  switch (calc_info.coord_type_number) {
  case 0: {
    grid_ptr.reset(new mesh_and_methods_cartesian(calc_info));
    break;
  }
  case 1: {
    grid_ptr.reset(new mesh_and_methods_axisymm(calc_info));
    break;
  }
  default: {
    std::cout << "wrong coordinates type!" << std::endl;
    return -1;
  }
  }

  statistics_manager stat(grid_ptr);

  data_2d_writer grid_writer(grid_ptr->get_mesh_const_ref(),
    calc_info, output_folder);

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

  //TODO: replace into grid_writer or discard
  std::ofstream pressure_diag(output_folder + "pressure_diag" +
    std::to_string(calc_info.par.size_x)
    + "x" + std::to_string(calc_info.par.size_y) + ".dat");
  grid_writer.output_in_wall_point_first(pressure_diag);

  grid_writer.output_on_symmetry_axis_first();

  //TODO: replace into statistic_manager or discard
  std::ofstream peaks(output_folder + "peaks.dat");

  std::cout << std::scientific;

  //TODO: save both previous grid and new grid in one mesh_and_methods
  std::shared_ptr<mesh_and_common_methods> new_grid_ptr;
  switch (calc_info.coord_type_number) {
  case 0: {
    new_grid_ptr.reset(new mesh_and_methods_cartesian(grid_ptr));
    break;
  }
  case 1: {
    new_grid_ptr.reset(new mesh_and_methods_axisymm(grid_ptr));
    break;
  }
  default: {
    std::cout << "wrong coordinates type!" << std::endl;
    return -1;
  }
  }
  int counter = 0;
  double current_t = 0.0;  
  std::vector<double> times;
  for (int i = 0; i*calc_info.par.time_step < calc_info.par.t_end; ++i) {
    times.push_back(i*calc_info.par.time_step);
  }
  double current_time_idx = 0;

  //peaks
  double rho_peak = 0.0, p_peak = 0.0;
  double rho_peak_time = 0.0, p_peak_time = 0.0;
  auto grid_mesh = grid_ptr->get_mesh_const_ref();//HORRIBLE NAMING
  const auto& corner_pt = grid_mesh[1][calc_info.par.size_x - 2];

  const double& delta_t_ref = grid_ptr->get_params_const_ref().delta_t;

  while (current_t < calc_info.par.t_end) {

    grid_ptr->calc_delta_t(current_t);
    new_grid_ptr->set_delta_t(delta_t_ref);
    new_grid_ptr->do_step(grid_ptr);

    std::swap(grid_ptr, new_grid_ptr);
    current_t += delta_t_ref;

    stat.update_pressure_sensors_on_wall(current_t);

    grid_writer.output_in_wall_point_for_current_time(pressure_diag, current_t);

    if (corner_pt->rho > rho_peak) {
      rho_peak = corner_pt->rho;
      rho_peak_time = current_t;
    }

    if (corner_pt->p > p_peak) {
      p_peak = corner_pt->p;
      p_peak_time = current_t;
    }

    grid_writer.output_on_symmetry_axis_for_current_time(current_t);

    if (current_time_idx + 1 < times.size() &&
        times[current_time_idx + 1] <= current_t) {
      grid_writer.output_for_current_time(current_t);
      ++current_time_idx;
    }

    ++counter;
    if (new_grid_ptr->get_params_const_ref().stop_now) {
      std::cout << "stop time:" << current_t << " steps: " << counter << std::endl;
      break;
    }

  }

  std::cout << current_t << " steps: " << counter << std::endl;

  grid_writer.output_for_current_time(current_t);


  //reflected shock with no bubble
  double p3;
  stat.calc_pressure_after_reflected_shock_no_bubble(
    grid_mesh[0][0]->rho, grid_mesh[0][0]->p, grid_mesh[0][0]->u, p3);
  double t0 = (calc_info.par.x_right - calc_info.par.x_left -
    calc_info.par.gap_btw_sw_and_bound) / calc_info.par.D_of_initial_shock;

  //sensors & impulses
  grid_writer.output_pressure_sensors_on_wall(output_folder,
    stat.get_pressure_sensors_on_wall());
  stat.impulses_output(t0, p3, output_folder);

  //peaks:

  peaks << "max rho in corner point:" << std::endl <<
    "time = " << rho_peak_time << " rho = " << rho_peak << std::endl;

  peaks << "max p in corner point:" << std::endl <<
    "time = " << p_peak_time << " p = " << p_peak << std::endl;

  pressure_diag.close();
  peaks.close();

  return 0;
}
