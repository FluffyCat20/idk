#include "data_2d.h"

int main() {

  std::ifstream input("D:/gas_dyn/idk/Lax_Friedrichs_method/config_2d.json");
  data_2d grid(input);
  input.close();

  int method_number = -1;
  if (grid.method_name == "mac_cormack")
    method_number = 0;

  if (method_number == -1){
    std::cout << "Error: unknown method name" << std::endl;
    return -1;
  }

  std::cout << std::scientific;

  data_2d new_grid(grid); //TODO: copy common data only
  double current_t = 0.0;
  int counter = 0;

  while (current_t < new_grid.t_end) {

    switch (method_number) {

    case 0: {
      new_grid.delta_t = grid.calc_delta_t();
      if (current_t + new_grid.delta_t > new_grid.t_end){
        new_grid.delta_t = new_grid.t_end - current_t;
      }
      grid.delta_t = new_grid.delta_t;
      new_grid.mac_cormack(grid);
      break;
    }
    default: {
      std::cout << "????" << std::endl;
      break;
    }
    }

    grid = new_grid;
    current_t += new_grid.delta_t;
    ++counter;
    if (new_grid.stop_now) {
      std::cout << "stop time:" << current_t << " steps: " << counter << std::endl;
      break;
    }

  }

  std::cout << current_t << " steps: " << counter << std::endl;

  std::ofstream outfile(grid.method_name + std::to_string(new_grid.size_x)
    + "x" + std::to_string(new_grid.size_y) + ".dat");
  grid.write_out_file(outfile);
  outfile.close();

  return 0;
}
