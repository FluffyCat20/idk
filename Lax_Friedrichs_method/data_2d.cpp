#include "utils.hpp"

void mesh_and_common_methods::do_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {

  switch (par.method_info.method_number) {

  case 0: {
    lax_friedrichs(prev_grid_ptr);
    break;
  }
  case 1: {
    mac_cormack(prev_grid_ptr);
    break;
  }
  case 2: {
    mac_cormack_with_davis(prev_grid_ptr);
    break;
  }
  case 3: {
    mac_cormack_with_zhmakin_fursenko(prev_grid_ptr);
    break;
  }
  default: {
    std::cout << "????" << std::endl;
    break;
  }
  }
  return;
}


void mesh_and_common_methods::shock_wave_initialization (
    double& rho, double& p, double& u,
    double rho1, double p1, double& u1,
    double mach) { //1 - before sw, no ind - after sw

  double a1 = std::sqrt(par.gamma*p1/rho1);
  u1 = mach*a1;
  p = p1*(2*par.gamma/(par.gamma+1)*mach*mach
          - (par.gamma-1)/(par.gamma+1));
  rho = rho1/(((par.gamma-1)/(par.gamma+1))
           + 2/(par.gamma+1)/mach/mach);
  u = u1*rho1/rho;

}

void mesh_and_common_methods::calc_delta_t(double current_time) {
  double max_u_abs_plus_a = std::numeric_limits<double>::min();
  for (const auto& row : mesh) {
    for (const data_node_2d* node : row) {
      double candidat = node->u_abs + node->a;
      if (candidat > max_u_abs_plus_a) {
        max_u_abs_plus_a = candidat;
      }
    }
  }
  par.delta_t = par.delta_x / max_u_abs_plus_a * par.courant_number;
  if (current_time + par.delta_t > par.t_end){
    par.delta_t = par.t_end - current_time;
  }
  return;
}

void mesh_and_common_methods::boundary_conditions() {

  //if (init_config["bubble_near_wall"]) {
    //boundary_conditions_for_bubble_near_wall();
  //boundary_conditions_for_bubble_near_wall_simmetry_on_bottom();
  boundary_conditions_default();
  //} else {
  //boundary_conditions_default();
  //}
}
