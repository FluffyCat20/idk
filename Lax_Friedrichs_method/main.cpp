#include "data.h"
#include "exact_solution.h"
#include "riemann_v2.c"

int main(){

  std::ifstream input("D:/gas_dyn/idk/Lax_Friedrichs_method/config.json");
  data grid(input);
  input.close();

  int method_number = -1;
  if (grid.method_name == "lax_friedrichs")
    method_number = 0;
  if (grid.method_name == "mac_cormack")
    method_number = 1;
  if (grid.method_name == "mac_cormack+davis")
    method_number = 2;
  if (grid.method_name == "godunov")
    method_number = 3;

  if (method_number == -1) {
    std::cout << "Error: unknown method name" << std::endl;
    return -1;
  }

  std::cout << std::scientific;

  data new_grid(grid);
  double current_t = 0.0;
  int counter = 0;
  while (current_t < new_grid.t_end){

    switch (method_number) {

      case 0:
      new_grid.delta_t = grid.calc_delta_t();
      if (current_t + new_grid.delta_t > new_grid.t_end){
        new_grid.delta_t = new_grid.t_end - current_t;
      }
      new_grid.calc_new_U_with_lax_friedrichs(grid);
      new_grid.calc_F_in_nodes();
      break;

      case 1:
      new_grid.delta_t = grid.calc_delta_t();
      if (current_t + new_grid.delta_t > new_grid.t_end){
        new_grid.delta_t = new_grid.t_end - current_t;
      }
      new_grid.calc_new_U_with_mac_cormack(grid);
      new_grid.calc_F_in_nodes();
      break;

      case 2:
      new_grid.delta_t = grid.calc_delta_t();
      if (current_t + new_grid.delta_t > new_grid.t_end) {
        new_grid.delta_t = new_grid.t_end - current_t;
      }
      new_grid.calc_new_U_with_mac_cormack_davis(grid);
      new_grid.calc_F_in_nodes();
      break;

      case 3:
      new_grid.calc_new_U_and_F_with_godunov(grid, current_t); //delta t and fluxes calculated inside
      break;

      default:
      std::cout << "????" << std::endl;
      return -1;

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

  //output
  if (grid.dimension == 1 && grid.method_name != "mac_cormack") {
    std::ofstream fout("out_" + grid.method_name + std::to_string(new_grid.size_x) + ".dat");
    fout << "x     rho     p     u\n";
    for (size_t i = 0; i < new_grid.size_x; ++i) {
      fout << i*new_grid.delta_x << " "
           << new_grid.mesh_1D[i].rho << " "
           << new_grid.mesh_1D[i].p << " "
           << new_grid.mesh_1D[i].u[0] << "\n";
    }
    fout.close();
  }

  if (grid.dimension == 1 && grid.method_name == "mac_cormack") {
    std::ofstream fout("out_" + grid.method_name + std::to_string(new_grid.size_x) + ".dat");
    fout << "x     rho     p     u\n";
    for (size_t i = 0; i < new_grid.size_x; ++i) {
      fout << i*new_grid.delta_x << " "
           << new_grid.mesh[1][i].rho << " "
           << new_grid.mesh[1][i].p << " "
           << new_grid.mesh[1][i].u[0] << "\n";
    }
    fout.close();
  }

  if (grid.dimension == 2 && grid.method_name == "mac_cormack") {
    std::ofstream fout("out_" + grid.method_name + std::to_string(new_grid.size_x)
                       + 'x' + std::to_string(new_grid.size_y)+ ".dat");
    fout << "x     y     rho     p     u     v\n";
    for (size_t i = 0; i < new_grid.size_y; ++i) {
      for (size_t j = 0; j < new_grid.size_x; ++j) {
        fout << j*new_grid.delta_x << " "
             << i*new_grid.delta_y << " "
             << new_grid.mesh[i][j].rho << " "
             << new_grid.mesh[i][j].p << " "
             << new_grid.mesh[i][j].u[0] << " "
             << new_grid.mesh[i][j].u[1] << "\n";
      }
    }
    fout.close();
  }

  std::ofstream fout;
  fout.open("out_exact_solution.dat");
  exact_solution exact_result(grid);
  exact_result.calc_and_print_result(fout);
  fout.close();

  /*std::ofstream output ("error_norms.json");
  output << std::setw(2) << error_norms << std::endl;*/


  return 0;
}
