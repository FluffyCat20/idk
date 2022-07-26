#include "data_2d_writer.h"

void data_2d_writer::output_first() {
  grid_ascii_outfile << "TITLE = \"" << std::to_string(par.method_number) <<
    "; size: " << par.size_x << "x" << par.size_y << "\""<< std::endl;
  grid_ascii_outfile << "VARIABLES = \"x\", \"y\", \"rho\", \"p\", \"u\", \"v\""
    << std::endl;
  grid_ascii_outfile << std::setprecision(4);
  output_for_current_time(0.0);
}

void data_2d_writer::output_for_current_time(double time) {
  grid_ascii_outfile << "ZONE T=\"0.0000\", SOLUTIONTIME=" << time << std::endl;
  grid_ascii_outfile << "I=" << par.size_y << ", J=" << par.size_x
          << ", F=POINT" << std::endl;

  for (size_t i = 0; i < par.size_x; ++i) {
    for (size_t j = 0; j < par.size_y; ++j) {
      const data_node_2d& pt = *mesh[j][i];
      grid_ascii_outfile << i*par.delta_x + par.x_left << " " <<
        j*par.delta_y + par.y_bottom << " " <<
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
  outfile << time << " " << mesh[1][par.size_x - 2]->rho <<
              " " << mesh[1][par.size_x - 2]->p << std::endl;
}

void data_2d_writer::output_on_symmetry_axis_first() {
  symmetry_axis_outfile << "TITLE = \" on symmetry axis \"" << std::endl;
  symmetry_axis_outfile << "VARIABLES = \"x\", \"rho\", \"p\"" << std::endl;
  symmetry_axis_outfile << std::setprecision(4);
  output_on_symmetry_axis_for_current_time(0.0);
}

void data_2d_writer::output_on_symmetry_axis_for_current_time(double time) {
  symmetry_axis_outfile << "ZONE I = " << par.size_x
    << ", F=POINT, SOLUTIONTIME = "
    << time << " t = " << "\"" << time << "\"" << std::endl;
  for (size_t x = 0; x < par.size_x; ++x) {
    const data_node_2d& pt = *mesh[0][x];
    symmetry_axis_outfile << x*par.delta_x + par.x_left << " "
            << pt.rho << " " << pt.p << std::endl;
  }
}

void data_2d_writer::output_pressure_sensors_on_wall(
    const std::string &output_folder,
    const std::vector<std::list<std::pair<double, double>>>&
      pressure_sensors_on_wall) {

  for (size_t i = 0; i < pressure_sensors_on_wall.size(); ++i) {
    std::ofstream fout(output_folder +
      "pressure_on_wall_sensors_norm/pressure_sensor" + std::to_string(i) + ".dat");
    for (const auto& it : pressure_sensors_on_wall[i]) {
      fout << it.first << " " << it.second << std::endl;
    }
    fout.close();
  }
}
