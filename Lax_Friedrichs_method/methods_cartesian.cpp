#include "data_2d.h"



/////////////NUMERICAL METHODS:////////////////////



void mesh_and_methods_cartesian::lax_friedrichs(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {
  auto prev_mesh = prev_grid_ptr->get_mesh_const_ref();
  for (size_t i = 1; i < par.size_y - 1; ++i) {
    for (size_t j = 1; j < par.size_x - 1; ++j) {
      auto up = prev_mesh[i-1][j];
      auto down = prev_mesh[i+1][j];
      auto left = prev_mesh[i][j-1];
      auto right = prev_mesh[i][j+1];
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] = 0.25*(up->U[k] + down->U[k] + left->U[k] + right->U[k]) -
          0.5*par.delta_t/par.delta_x * (right->F[k] - left->F[k]) -
          0.5*par.delta_t/par.delta_y * (up->G[k] - down->G[k]);
      }
    }
  }
  calc_values_and_fluxes();
  boundary_conditions();
}

void mesh_and_methods_cartesian::mac_cormack(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {
  std::shared_ptr<mesh_and_common_methods> predictor_grid_ptr (
    new mesh_and_methods_cartesian(prev_grid_ptr));
  predictor_grid_ptr->mac_cormack_predictor_step(prev_grid_ptr);
  mac_cormack_corrector_step(prev_grid_ptr, predictor_grid_ptr);
  calc_values_and_fluxes();
  boundary_conditions();
}

void mesh_and_methods_cartesian::mac_cormack_predictor_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {
  //must be called from predictor grid

  auto prev_mesh = prev_grid_ptr->get_mesh_const_ref();

  for (size_t i = 0; i < par.size_y - 1; ++i) { //0?? not 1??
    for (size_t j = 0; j < par.size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] = prev_mesh[i][j]->U[k]
          - par.delta_t*((prev_mesh[i][j+1]->F[k]
          - prev_mesh[i][j]->F[k])/par.delta_x
          + (prev_mesh[i+1][j]->G[k] - prev_mesh[i][j]->G[k])/par.delta_y);
      }
    }
  }
  calc_values_and_fluxes();
  boundary_conditions();
}

void mesh_and_methods_cartesian::mac_cormack_corrector_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid,
    std::shared_ptr<const mesh_and_common_methods> predictor_grid) {
  auto prev_mesh = prev_grid->get_mesh_const_ref();
  auto predictor_mesh = predictor_grid->get_mesh_const_ref();
  for (size_t i = 1; i < par.size_y - 1; ++i) { //size?? not size - 1??
    for (size_t j = 1; j < par.size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] =
          0.5*(prev_mesh[i][j]->U[k] + predictor_mesh[i][j]->U[k])
          - 0.5*par.delta_t*((predictor_mesh[i][j]->F[k] - predictor_mesh[i][j - 1]->F[k])/par.delta_x
          + (predictor_mesh[i][j]->G[k] - predictor_mesh[i - 1][j]->G[k])/par.delta_y);
      }
    }
  }
}

void mesh_and_methods_cartesian::mac_cormack_with_davis(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {

  std::shared_ptr<mesh_and_common_methods> predictor_grid_ptr (
    new mesh_and_methods_cartesian(prev_grid_ptr));
  predictor_grid_ptr->mac_cormack_predictor_step(prev_grid_ptr);
  mac_cormack_corrector_step(prev_grid_ptr, predictor_grid_ptr);
  boundary_conditions();
  calc_davis_artificial_viscosity();
  calc_values_and_fluxes();
  boundary_conditions();
}


double inner_product(const std::vector<double>& v1,
                     const std::vector<double>& v2){
  double res = 0;
  for (size_t i = 0; i < v1.size(); ++i){
    res += v1[i]*v2[i];
  }
  return res;
}

inline double phi(double rp, double rm) { //rp = r+, rm = r-
  return std::max ( {
         0.0,
         std::min({2*rp, rm, 1.0}),
         std::min({2*rm, rp, 1.0}) });
}

void mesh_and_methods_cartesian::calc_davis_artificial_viscosity() {
  std::vector<std::vector<double>>
    r_x_plus(par.size_y, std::vector<double>(par.size_x)),
    r_x_minus(par.size_y, std::vector<double>(par.size_x)),
    r_y_plus(par.size_y, std::vector<double>(par.size_x)),
    r_y_minus(par.size_y, std::vector<double>(par.size_x));
  std::vector<std::vector<std::vector<double>>>
    delta_u_x(par.size_y, std::vector<std::vector<double>>(par.size_x, std::vector<double>(4))),
    delta_u_y(par.size_y, std::vector<std::vector<double>>(par.size_x, std::vector<double>(4))),
    D_x(par.size_y, std::vector<std::vector<double>>(par.size_x, std::vector<double>(4))),
    D_y(par.size_y, std::vector<std::vector<double>>(par.size_x, std::vector<double>(4)));

  //delta_u calculation:
  for (size_t i = 0; i < par.size_y - 1; ++i) {
    for (size_t j = 0; j < par.size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        delta_u_x[i][j][k] = mesh[i][j+1]->U[k] - mesh[i][j]->U[k];
        delta_u_y[i][j][k] = mesh[i+1][j]->U[k] - mesh[i][j]->U[k];
      }
    }
  }

  //r calculation: //TODO: fix: inner products are calculated twice (save squares?)
  for (size_t i = 1; i < par.size_y - 1; ++i) {
    for (size_t j = 1; j < par.size_x - 1; ++j) {
      //x-direction
      double product_mid_x = inner_product(delta_u_x[i][j-1], delta_u_x[i][j]);
      double product_left_x = inner_product(delta_u_x[i][j-1], delta_u_x[i][j-1]);
      double product_right_x = inner_product(delta_u_x[i][j], delta_u_x[i][j]);

      r_x_plus[i][j] = product_mid_x/product_right_x;
      r_x_minus[i][j] = product_mid_x/product_left_x;

      //y_direction
      double product_mid_y = inner_product(delta_u_y[i-1][j], delta_u_y[i][j]);
      double product_up_y = inner_product(delta_u_y[i-1][j], delta_u_y[i-1][j]);
      double product_down_y = inner_product(delta_u_y[i][j], delta_u_y[i][j]);

      r_y_plus[i][j] = product_mid_y/product_down_y;
      r_y_minus[i][j] = product_mid_y/product_up_y;
    }
  }

  //nu, K, D calculation:
  for (size_t i = 1; i < par.size_y - 1; ++i) {
    for (size_t j = 1; j < par.size_x - 1; ++j) {
      double lambda = std::max(std::abs(mesh[i][j]->u_abs + mesh[i][j]->a),
                               std::abs(mesh[i][j]->u_abs - mesh[i][j]->a));
      //nu_x, k_x
      double nu_x = par.delta_t/par.delta_x*lambda;
      double c_nu_x = nu_x*(1-nu_x);
      if (nu_x > 0.5) {
        c_nu_x = 0.25;
      }
      double k_x = 0.5*c_nu_x*(1 - phi(r_x_plus[i][j], r_x_minus[i][j+1]));

      //nu_y, k_y
      double nu_y = par.delta_t / par.delta_y*lambda;
      double c_nu_y = nu_y*(1-nu_y);
      if (nu_y > 0.5) {
        c_nu_y = 0.25;
      }
      double k_y = 0.5*c_nu_y*(1 - phi(r_y_plus[i][j], r_y_minus[i+1][j]));

      //D
      for (size_t k = 0; k < 4; ++k) {
        D_x[i][j][k] = k_x*delta_u_x[i][j][k];
        D_y[i][j][k] = k_y*delta_u_y[i][j][k];
      }

    }
  }

  for (size_t k = 0; k < 4; ++k) {
    for (size_t i = 1; i < par.size_y - 1; ++i){
      D_x[i][0][k] = 0.0;
      D_x[i][par.size_x - 1][k] = 0.0;
      D_y[i][0][k] = 0.0;
      D_y[i][par.size_x - 1][k] = 0.0;
    }
    for (size_t j = 1; j < par.size_x - 1; ++j){
      D_x[0][j][k] = 0.0;
      D_x[par.size_y - 1][j][k] = 0.0;
      D_y[0][j][k] = 0.0;
      D_y[par.size_y - 1][j][k] = 0.0;
    }
  }

  for (size_t i = 1; i < par.size_y - 1; ++i) {
    for (size_t j = 1; j < par.size_x - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] += D_x[i][j][k] - D_x[i][j-1][k] + D_y[i][j][k] - D_y[i-1][j][k];
      }
    }
  }

}


///////////////////////BOUNDARY CONDITIONS://///////////////////////////

//in boundary conditions, you need to operate with the whole nodes of mesh
//so don't forget they're pointers! more * motherfucker

void mesh_and_methods_cartesian::boundary_conditions_default() {
  for (size_t y = 1; y < par.size_y - 1; ++y) {
    *mesh[y][0] = *mesh[y][1];
    *mesh[y][par.size_x - 1] = *mesh[y][par.size_x - 2];
  }
  for (size_t x = 1; x < par.size_x - 1; ++x) {
    *mesh[0][x] = *mesh[1][x];
    *mesh[par.size_y - 1][x] = *mesh[par.size_y - 2][x];
  }
}

void mesh_and_methods_cartesian::boundary_conditions_for_bubble_near_wall() {
  // solid wall on the right:
  for (size_t y = 1; y < par.size_y - 1; ++y) {
    data_node_2d& node_to_change = *mesh[y][par.size_x - 1];
    node_to_change = *mesh[y][par.size_x - 2];
    node_to_change.u *= -1;
    node_to_change.U[1] *= -1;
    node_to_change.F[0] *= -1;
    node_to_change.F[2] *= -1;
    node_to_change.F[3] *= -1;
    node_to_change.G[1] *= -1;
  }

  // d/dn = 0 on other bounds:
  for (size_t y = 1; y < par.size_y - 1; ++y) {
    *mesh[y][0] = *mesh[y][1];
  }
  for (size_t x = 1; x < par.size_x - 1; ++x) {
    *mesh[0][x] = *mesh[1][x];
    *mesh[par.size_y - 1][x] = *mesh[par.size_y - 2][x];
  }
}

void mesh_and_methods_cartesian::
boundary_conditions_for_bubble_near_wall_simmetry_on_bottom() {
  // solid wall on the right:
  for (size_t y = 1; y < par.size_y - 1; ++y) {
    data_node_2d& node_to_change = *mesh[y][par.size_x - 1];
    node_to_change = *mesh[y][par.size_x - 2];
    node_to_change.u *= -1;
    node_to_change.U[1] *= -1;
    node_to_change.F[0] *= -1;
    node_to_change.F[2] *= -1;
    node_to_change.F[3] *= -1;
    node_to_change.G[1] *= -1;
  }

  // d/dn = 0 on left&right bounds and symmetry on bottom:
  for (size_t y = 1; y < par.size_y - 1; ++y) {
    *mesh[y][0] = *mesh[y][1];
  }
  for (size_t x = 1; x < par.size_x - 1; ++x) {
    data_node_2d& node_to_change = *mesh[0][x];
    node_to_change = *mesh[1][x];
    node_to_change.v *= -1;
    node_to_change.U[2] *= -1;
    node_to_change.F[2] *= -1;
    node_to_change.G[0] *= -1;
    node_to_change.G[1] *= -1;
    node_to_change.G[3] *= -1;

    *mesh[par.size_y - 1][x] = *mesh[par.size_y - 2][x];

  }
}
