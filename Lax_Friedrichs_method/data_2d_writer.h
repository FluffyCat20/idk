#ifndef DATA2D_WRITER_H
#define DATA2D_WRITER_H
#pragma once

#include "calculation_info.hpp"

class data_2d_writer {

public:

  data_2d_writer(const base_data_2d& grid_,
      const calculation_params& par_,
      const std::string& output_folder)
    : grid(grid_), par(par_){
    std::string grid_ascii_outfile_name =
      output_folder + std::to_string(par.size_x)
      + "x" + std::to_string(par.size_y) + ".dat";
    grid_ascii_outfile.open(grid_ascii_outfile_name);

    std::string symmetry_axis_outfile_name =
      output_folder + "rho_p_on_symmetry_axis" +
      std::to_string(par.size_x) + "x" +
      std::to_string(par.size_y) + ".dat";
    symmetry_axis_outfile.open(symmetry_axis_outfile_name);
  }

  ~data_2d_writer() {
    grid_ascii_outfile.close();
    symmetry_axis_outfile.close();
  }

  void output_first();
  void output_for_current_time(double time);
  void output_in_wall_point_first(std::ofstream& outfile);
  void output_in_wall_point_for_current_time(
    std::ofstream& outfile, double time);
  void output_on_symmetry_axis_first();
  void output_on_symmetry_axis_for_current_time(double time);

  void output_pressure_sensors_on_wall(const std::string& output_folder);

private:
  const base_data_2d& grid;
  const calculation_params& par;
  std::string output_folder;

  std::ofstream grid_ascii_outfile;
  std::ofstream symmetry_axis_outfile;
};

#include "data_2d_writer.inl"
#endif // DATA2D_WRITER_H
