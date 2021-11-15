#include "data.h"
#include "riemann_v2.c"

int main(){

  std::ifstream input("D:/gas_dyn/Lax_Friedrichs_method/config.json");
  data grid(input);
  input.close();
  grid.delta_t = grid.calc_delta_t();

  int method_number = -1;
  if (grid.method_name == "lax_friedrichs")
    method_number = 0;
  if (grid.method_name == "mac_cormack+davis")
    method_number = 1;

  if (method_number == -1){
    std::cout << "Error: unknown method name" << std::endl;
    return -1;
  }

  std::cout << std::scientific;

  data new_grid (grid);
  double current_t = 0.0;
  int counter = 0;
  while (current_t < new_grid.t_end){
    new_grid.delta_t = grid.calc_delta_t();
    if (current_t + new_grid.delta_t > new_grid.t_end){
      new_grid.delta_t = new_grid.t_end - current_t;
    }

    switch (method_number) {

      case 0:
        new_grid.calc_new_U_with_lax_friedrichs(grid);
        break;
      case 1:
        new_grid.calc_new_U_with_mac_cormack(grid);
        break;
      default:
        std::cout << "????" << std::endl;
        return -1;

    }

    new_grid.calc_F();
    grid = new_grid;
    current_t += new_grid.delta_t;
    ++counter;
    if (new_grid.stop_now){
      std::cout << "stop time:" << current_t << " steps: " << counter << std::endl;
      break;
    }
  }

  //output
  std::ofstream fout("out_" + grid.method_name + ".txt");
  fout << "x     rho     p     u\n";
  for (size_t i = 0; i < new_grid.size; ++i){
    fout << i*new_grid.delta_x << " "
         << new_grid.mesh[i].rho << " "
         << new_grid.mesh[i].p << " "
         << new_grid.mesh[i].u << "\n";
  }
  fout.close();

  /*
  //exact solution calculating
  exact_solution exact(rho1, p1, rho5, p5, gamma, size, shock_x);
  exact.calc_bounds(t_end, delta_x);
  exact.calc_area2(t_end, delta_x);

  //out for exact solution
  fout.open("out_exact_solution.txt");
  exact.output(fout, delta_x);
  fout.close();
  */

  double Zl [4] = {new_grid.rho_left, new_grid.u_left, new_grid.p_left, new_grid.gamma};
  double Zr [4] = {new_grid.rho_right, new_grid.u_right, new_grid.p_right, new_grid.gamma};
  double* R1 = new double [new_grid.size];
  double* U1 = new double [new_grid.size];
  double* P1 = new double [new_grid.size];
  Riemann_Solve(Zl, Zr, R1, U1, P1, current_t/*new_grid.t_end*/, -0.5, 0.5, double(new_grid.size));

  fout.open("out_exact_solution_correct.txt");
  fout << "x     rho     p     u\n";
  double x = 0.0;
  for (size_t i = 0; i < new_grid.size; ++i){
    fout << x << " " << R1[i] << " "
      << P1[i] << " "
      << U1[i] << "\n";
    x += new_grid.delta_x;
  }
  fout.close();

  //error norms calculating
  std::vector<std::vector<double>> norms =
    new_grid.calc_error_norm(R1, P1, U1, new_grid.size);

  std::vector<std::vector<double>> norms_first_half =
    new_grid.calc_error_norm(R1, P1, U1, new_grid.size / 2);

  delete[] R1;
  delete[] U1;
  delete[] P1;

  //writing norms to json
  json error_norms;
  error_norms["whole_grid"]["uniform"] = {
    {"rho", norms[0][0]},
    {"p", norms[0][1]},
    {"u", norms[0][2]}
  };
  error_norms["whole_grid"]["euclid"] = {
    {"rho", norms[1][0]},
    {"p", norms[1][1]},
    {"u", norms[1][2]}
  };
  error_norms["whole_grid"]["deviation"] = {
    {"rho", norms[2][0]},
    {"p", norms[2][1]},
    {"u", norms[2][2]}
  };
  error_norms["whole_grid"]["integral"] = {
    {"rho", norms[3][0]},
    {"p", norms[3][1]},
    {"u", norms[3][2]}
  };

  error_norms["first_half_of_grid"]["uniform"] = {
    {"rho", norms_first_half[0][0]},
    {"p", norms_first_half[0][1]},
    {"u", norms_first_half[0][2]}
  };
  error_norms["first_half_of_grid"]["euclid"] = {
    {"rho", norms_first_half[1][0]},
    {"p", norms_first_half[1][1]},
    {"u", norms_first_half[1][2]}
  };
  error_norms["first_half_of_grid"]["deviation"] = {
    {"rho", norms_first_half[2][0]},
    {"p", norms_first_half[2][1]},
    {"u", norms_first_half[2][2]}
  };
  error_norms["first_half_of_grid"]["integral"] = {
    {"rho", norms_first_half[3][0]},
    {"p", norms_first_half[3][1]},
    {"u", norms_first_half[3][2]}
  };

  std::ofstream output ("error_norms.json");
  output << std::setw(2) << error_norms << std::endl;


  return 0;
}
