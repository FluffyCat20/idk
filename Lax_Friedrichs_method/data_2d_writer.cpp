#include "data_2d_writer.h"

void data_2d_writer::output_first() {
  grid_ascii_outfile << "TITLE = \"" << data.method_name <<
    "; size: " << data.size_x << "x" << data.size_y << "\""<< std::endl;
  grid_ascii_outfile << "VARIABLES = \"x\", \"y\", \"rho\", \"p\", \"u\", \"v\""
    << std::endl;
  grid_ascii_outfile << std::setprecision(4);
  output_for_current_time(0.0);
}

void data_2d_writer::output_for_current_time(double time) {
  grid_ascii_outfile << "ZONE T=\"0.0000\", SOLUTIONTIME=" << time << std::endl;
  grid_ascii_outfile << "I=" << data.size_y << ", J=" << data.size_x
          << ", F=POINT" << std::endl;

  for (size_t i = 0; i < data.size_x; ++i) {
    for (size_t j = 0; j < data.size_y; ++j) {
      const data_node_2d& pt = data.mesh[j][i];
      grid_ascii_outfile << i*data.delta_x + data.x_left << " " <<
        j*data.delta_y + data.y_bottom << " " <<
        pt.rho << " " << pt.p << " " << pt.u << " " << pt.v << std::endl;
    }
  }

  std::cout << "Calculations done: time = " << time << std::endl;

}

void data_2d_writer::output_in_wall_point_first(std::ofstream &outfile) {
  outfile << "TITLE = \" in central point \"" << std::endl;
  outfile << "VARIABLES = \"t\", \"rho\", \"p\"" << std::endl;
  outfile << std::setprecision(4);
  output_in_wall_point_for_current_time(outfile, 0.0);
}

void data_2d_writer::output_in_wall_point_for_current_time(
    std::ofstream &outfile, double time) {
  outfile << time << " " << data.mesh[1][data.size_x - 2].rho <<
              " " << data.mesh[1][data.size_x - 2].p << std::endl;
}

void data_2d_writer::output_on_symmetry_axis_first() {
  symmetry_axis_outfile << "TITLE = \" on symmetry axis \"" << std::endl;
  symmetry_axis_outfile << "VARIABLES = \"x\", \"rho\", \"p\"" << std::endl;
  symmetry_axis_outfile << std::setprecision(4);
  output_on_symmetry_axis_for_current_time(0.0);
}

void data_2d_writer::output_on_symmetry_axis_for_current_time(double time) {
  symmetry_axis_outfile << "ZONE I = " << data.size_x
    << ", F=POINT, SOLUTIONTIME = "
    << time << " t = " << "\"" << time << "\"" << std::endl;
  for (size_t x = 0; x < data.size_x; ++x) {
    const data_node_2d& pt = data.mesh[0][x];
    symmetry_axis_outfile << x*data.delta_x + data.x_left << " "
            << pt.rho << " " << pt.p << std::endl;
  }
}

void data_2d_writer::output_pressure_sensors_on_wall(
    const std::string &output_folder) {

  for (size_t i = 0; i < data.pressure_sensors_on_wall.size(); ++i) {
    std::ofstream fout(output_folder +
      "pressure_on_wall_sensors_norm/pressure_sensor" + std::to_string(i) + ".dat");
    for (const auto& it : data.pressure_sensors_on_wall[i]) {
      fout << it.first << " " << it.second << std::endl;
    }
    fout.close();
  }
}
