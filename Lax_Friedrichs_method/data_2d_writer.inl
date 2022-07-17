#include "data_2d_writer.h"

template <typename node_t>
void data_2d_writer<node_t>::output_first() {
  grid_ascii_outfile << "TITLE = \"" <<
    "size: " << grid.par.size_x << "x" << grid.par.size_y << "\""<< std::endl;
  grid_ascii_outfile << "VARIABLES = \"x\", \"y\", \"rho\", \"p\", \"u\", \"v\""
    << std::endl;
  grid_ascii_outfile << std::setprecision(4);
  output_for_current_time(0.0);
}

template <typename node_t>
void data_2d_writer<node_t>::output_for_current_time(double time) {
  grid_ascii_outfile << "ZONE T=\"0.0000\", SOLUTIONTIME=" << time << std::endl;
  grid_ascii_outfile << "I=" << grid.par.size_y << ", J=" << grid.par.size_x
          << ", F=POINT" << std::endl;

  for (size_t i = 0; i < grid.par.size_x; ++i) {
    for (size_t j = 0; j < grid.par.size_y; ++j) {
      const data_node_2d& pt = grid.mesh[j][i];
      grid_ascii_outfile << i*grid.par.delta_x + grid.par.x_left << " " <<
        j*grid.par.delta_y + grid.par.y_bottom << " " <<
        pt.rho << " " << pt.p << " " << pt.u << " " << pt.v << std::endl;
    }
  }

  std::cout << "Calculations done: time = " << time << std::endl;

}

template <typename node_t>
void data_2d_writer<node_t>::output_in_wall_point_first(std::ofstream &outfile) {
  outfile << "TITLE = \" in central point \"" << std::endl;
  outfile << "VARIABLES = \"t\", \"rho\", \"p\"" << std::endl;
  outfile << std::setprecision(4);
  output_in_wall_point_for_current_time(outfile, 0.0);
}

template <typename node_t>
void data_2d_writer<node_t>::output_in_wall_point_for_current_time(
    std::ofstream &outfile, double time) {
  outfile << time << " " << grid.mesh[1][grid.par.size_x - 2].rho <<
              " " << grid.mesh[1][grid.par.size_x - 2].p << std::endl;
}

template <typename node_t>
void data_2d_writer<node_t>::output_on_symmetry_axis_first() {
  symmetry_axis_outfile << "TITLE = \" on symmetry axis \"" << std::endl;
  symmetry_axis_outfile << "VARIABLES = \"x\", \"rho\", \"p\"" << std::endl;
  symmetry_axis_outfile << std::setprecision(4);
  output_on_symmetry_axis_for_current_time(0.0);
}

template <typename node_t>
void data_2d_writer<node_t>::output_on_symmetry_axis_for_current_time(
    double time) {
  symmetry_axis_outfile << "ZONE I = " << grid.par.size_x
    << ", F=POINT, SOLUTIONTIME = "
    << time << " t = " << "\"" << time << "\"" << std::endl;
  for (size_t x = 0; x < grid.par.size_x; ++x) {
    const data_node_2d& pt = grid.mesh[0][x];
    symmetry_axis_outfile << x*grid.par.delta_x + grid.par.x_left << " "
            << pt.rho << " " << pt.p << std::endl;
  }
}

template <typename node_t>
void data_2d_writer<node_t>::output_pressure_sensors_on_wall(
    const std::string &output_folder) {

  for (size_t i = 0; i < grid.pressure_sensors_on_wall.size(); ++i) {
    std::ofstream fout(output_folder +
      "pressure_on_wall_sensors_norm/pressure_sensor" + std::to_string(i) + ".dat");
    for (const auto& it : grid.pressure_sensors_on_wall[i]) {
      fout << it.first << " " << it.second << std::endl;
    }
    fout.close();
  }
}
