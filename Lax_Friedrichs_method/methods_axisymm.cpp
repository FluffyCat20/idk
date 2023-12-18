#include "data_2d.h"

void mesh_and_methods_axisymm::calc_values_and_fluxes() {
  for (size_t i = y_begin_ind() + 1; i < y_end_ind(); ++i) {
    for (size_t j = x_begin_ind(); j < x_end_ind(); ++j) {
      mesh[i][j]->calc_values_from_U();
      mesh[i][j]->calc_F_when_values_known();
      mesh[i][j]->calc_G_when_values_known();
      mesh[i][j]->calc_H_when_values_known();
    }
  }
}

/////////////NUMERICAL METHODS:////////////////////


void mesh_and_methods_axisymm::mac_cormack(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {
  std::shared_ptr<mesh_and_common_methods> predictor_grid_ptr (
    new mesh_and_methods_axisymm(prev_grid_ptr));
  predictor_grid_ptr->mac_cormack_predictor_step(prev_grid_ptr);
  mac_cormack_corrector_step(prev_grid_ptr, predictor_grid_ptr);
  calc_values_and_fluxes();
  boundary_conditions();
}


void mesh_and_methods_axisymm::mac_cormack_predictor_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {
  //must be called from predictor grid

  auto prev_mesh = prev_grid_ptr->get_mesh_const_ref();

  for (size_t i = y_begin_ind() + 1; i < y_end_ind(); ++i) {
    for (size_t j = x_begin_ind(); j < x_end_ind(); ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] = prev_mesh[i][j]->U[k]
          - par.delta_t*((prev_mesh[i][j+1]->F[k]
          - prev_mesh[i][j]->F[k])/par.delta_x
          + (prev_mesh[i+1][j]->G[k] - prev_mesh[i][j]->G[k])/par.delta_y
          - prev_mesh[i][j]->H[k]);
      }
    }
  }
  calc_values_and_fluxes();
  boundary_conditions();
}


void mesh_and_methods_axisymm::mac_cormack_corrector_step(
    std::shared_ptr<const mesh_and_common_methods> prev_grid,
    std::shared_ptr<const mesh_and_common_methods> predictor_grid) {
  auto prev_mesh = prev_grid->get_mesh_const_ref();
  auto predictor_mesh = predictor_grid->get_mesh_const_ref();
  for (size_t i = y_begin_ind() + 1; i < y_end_ind(); ++i) {
    for (size_t j = x_begin_ind(); j < x_end_ind(); ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] =
          0.5*((prev_mesh[i][j]->U[k] + predictor_mesh[i][j]->U[k])
          - par.delta_t*((predictor_mesh[i][j]->F[k] - predictor_mesh[i][j - 1]->F[k])/par.delta_x
          + (predictor_mesh[i][j]->G[k] - predictor_mesh[i - 1][j]->G[k])/par.delta_y)
          + par.delta_t * predictor_mesh[i][j]->H[k]);
      }
    }
  }
}


void mesh_and_methods_axisymm::mac_cormack_with_davis(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {

  std::shared_ptr<mesh_and_common_methods> predictor_grid_ptr (
    new mesh_and_methods_axisymm(prev_grid_ptr));
  predictor_grid_ptr->mac_cormack_predictor_step(prev_grid_ptr);
  mac_cormack_corrector_step(prev_grid_ptr, predictor_grid_ptr);
  boundary_conditions();
  calc_davis_artificial_viscosity();
  calc_values_and_fluxes();
  fix_small_p_rho();
  boundary_conditions();
}


inline double inner_product(const std::vector<double>& v1,
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

void mesh_and_methods_axisymm::calc_davis_artificial_viscosity() {
  std::vector<std::vector<double>>
    r_x_plus(mesh.size(), std::vector<double>(mesh[0].size())),
    r_x_minus(mesh.size(), std::vector<double>(mesh[0].size())),
    r_y_plus(mesh.size(), std::vector<double>(mesh[0].size())),
    r_y_minus(mesh.size(), std::vector<double>(mesh[0].size()));
  std::vector<std::vector<std::vector<double>>>
    delta_u_x(mesh.size(), std::vector<std::vector<double>>(mesh[0].size(), std::vector<double>(4))),
    delta_u_y(mesh.size(), std::vector<std::vector<double>>(mesh[0].size(), std::vector<double>(4))),
    D_x(mesh.size(), std::vector<std::vector<double>>(mesh[0].size(), std::vector<double>(4))),
    D_y(mesh.size(), std::vector<std::vector<double>>(mesh[0].size(), std::vector<double>(4)));

  //delta_u calculation:
  for (size_t i = 0; i < mesh.size() - 1; ++i) {
    for (size_t j = 0; j < mesh[0].size() - 1; ++j) {
      for (size_t k = 0; k < 4; ++k) {
        delta_u_x[i][j][k] = mesh[i][j+1]->U[k] - mesh[i][j]->U[k];
        delta_u_y[i][j][k] = mesh[i+1][j]->U[k] - mesh[i][j]->U[k];
      }
    }
  }

  //r calculation: //TODO: fix: inner products are calculated twice (save squares?)
  for (size_t i = 1; i < mesh.size() - 1; ++i) {
    for (size_t j = 1; j < mesh[0].size() - 1; ++j) {
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
  for (size_t i = 1; i < mesh.size() - 1; ++i) {
    for (size_t j = 1; j < mesh[0].size() - 1; ++j) {
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

  for (size_t i = y_begin_ind() + 1; i < y_end_ind(); ++i) {
    for (size_t j = x_begin_ind(); j < x_end_ind(); ++j) {
      for (size_t k = 0; k < 4; ++k) {
        mesh[i][j]->U[k] += D_x[i][j][k] - D_x[i][j-1][k] + D_y[i][j][k] - D_y[i-1][j][k];
      }
    }
  }

}

void mesh_and_methods_axisymm::mac_cormack_with_zhmakin_fursenko(
    std::shared_ptr<const mesh_and_common_methods> prev_grid_ptr) {

  std::shared_ptr<mesh_and_common_methods> predictor_grid_ptr (
    new mesh_and_methods_axisymm(prev_grid_ptr));
  predictor_grid_ptr->mac_cormack_predictor_step(prev_grid_ptr);
  mac_cormack_corrector_step(prev_grid_ptr, predictor_grid_ptr);
  //boundary_conditions();
  smooth_2D_zhmakin_fursenko(0.2);
  calc_values_and_fluxes();
  boundary_conditions();
  return;
}

void mesh_and_methods_axisymm::smooth_2D_zhmakin_fursenko(
  const double smooth_intensity) {

  //smooth rows
  for (size_t i = y_begin_ind(); i < y_end_ind(); ++i) {
    Smooth_Array_Zhmakin_Fursenko(mesh[i], smooth_intensity); //by OG
    //smooth_1D_zhmakin_fursenko(mesh[i], smooth_intensity);
  }

  //smooth columns
  for (size_t j = x_begin_ind(); j < x_end_ind(); ++j) {
    std::vector<data_node_2d*> column(mesh.size());
    for (size_t i = 0; i < mesh.size(); ++i) {
      column[i] = mesh[i][j]; //ptr is copied, so nodes from mesh are changed in zhmakin_fursenko
    }
    Smooth_Array_Zhmakin_Fursenko(column, smooth_intensity); //by OG
    //smooth_1D_zhmakin_fursenko(column, smooth_intensity);
  }
}

void mesh_and_methods_axisymm::smooth_1D_zhmakin_fursenko(
    std::vector<data_node_2d*>& vec_to_smooth, const double smooth_intensity) {

  const double equations_amount = vec_to_smooth[0]->U.size();

  for (size_t eq_num = 0; eq_num < equations_amount; ++eq_num) {

    //fluxes:
    std::vector<double> fluxes(vec_to_smooth.size());
    for (size_t i = 0; i < vec_to_smooth.size() - 1; i++) {
      fluxes[i] = smooth_intensity *
        (vec_to_smooth[i + 1]->U[eq_num] - vec_to_smooth[i]->U[eq_num]);
    }

    //diffusion operator (D):
    //diffused = (1+D) applied on vec_to_smooth, D is the difference of fluxes):
    std::vector<double> diffused(vec_to_smooth.size());
    for (size_t i = 1; i < vec_to_smooth.size() - 1; i++) {
      diffused[i] = vec_to_smooth[i]->U[eq_num] +
                    (fluxes[i] - fluxes[i - 1]);
    }

    //difference of diffused:
    std::vector<double> delta_dif(vec_to_smooth.size());
    for (size_t i = 1; i < vec_to_smooth.size() - 2; i++) {
      delta_dif[i] = diffused[i + 1] - diffused[i];
    }

    //flux limiter:
    std::vector<double> fluxes_limited(vec_to_smooth.size());
    for (size_t i = 2; i < vec_to_smooth.size() - 3; i++) {
      int s = fluxes[i] > 0 ? 1 : -1; //sign of the flux
      fluxes_limited[i] = std::min(
        {s * delta_dif[i - 1], std::abs(fluxes[i]), s * delta_dif[i + 1]});
      fluxes_limited[i] = s * std::max(0.0, fluxes_limited[i]);
    }

    //anti-diffusion operator (A):
    //anti_diffused = (1+A) applied on diffused:
    std::vector<double> anti_diffused(vec_to_smooth.size());
    for (size_t i = 3; i < vec_to_smooth.size() - 3; i++) {
      anti_diffused[i] = diffused[i] - fluxes_limited[i] + fluxes_limited[i - 1];
    }

    //writing result:
    for (size_t i = 3; i < vec_to_smooth.size() - 3; i++) {
      vec_to_smooth[i]->U[eq_num] = anti_diffused[i];
    }
  }
}

void mesh_and_methods_axisymm::Smooth_Array_Zhmakin_Fursenko(
  std::vector<data_node_2d*>& vec_to_smooth, const double smooth_intensity)
{
  const int equations_amount = vec_to_smooth[0]->U.size();
  int Nmax = vec_to_smooth.size() - 1;

  for (size_t eq_num = 0; eq_num < equations_amount; ++eq_num) {

    double* Fs = new double[Nmax + 1];
    for (size_t i = 0; i < Nmax + 1; ++i) {
      Fs[i] = vec_to_smooth[i]->U[eq_num];
    }
    double* Fi = new double[Nmax + 1];
    double* Ftd = new double[Nmax + 1];

    int n;
    double sign, f1, f2; // double_smooth_intensity читать из файла
    for (n = 0; n <= Nmax - 1; n++)
      Fi[n] = smooth_intensity*(Fs[n + 1] - Fs[n]);
    for (n = 1; n <= Nmax - 1; n++)
      Ftd[n] = Fs[n] + (Fi[n] - Fi[n - 1]);
    Ftd[0] = Fs[0];
    Ftd[Nmax] = Fs[Nmax];
    for (n = 0; n <= Nmax - 1; n++)
      Fs[n] = Ftd[n + 1] - Ftd[n];
    for (n = 0; n <= Nmax - 1; n++)
    {
      Fi[n] = smooth_intensity*(Ftd[n + 1] - Ftd[n]);
      sign = -1;
      if (Fi[n] >= 0)
        sign = 1;
      Fi[n] = std::abs(Fi[n]);
      if (n == 0)
      {
        f2 = sign*Fs[1];
        if (f2<Fi[0])
          Fi[0] = f2;
      }
      else if (n == Nmax - 1)
      {
        f1 = sign*Fs[Nmax - 2];
        if (f1<Fi[Nmax - 1])
          Fi[Nmax - 1] = f1;
      }
      else
      {
        f1 = sign*Fs[n - 1];
        f2 = sign*Fs[n + 1];
        if (f1<Fi[n])
          Fi[n] = f1;
        if (f2<Fi[n])
          Fi[n] = f2;
      }
      if (Fi[n]<0)
        Fi[n] = 0;
      else
        Fi[n] = sign*Fi[n];
    }
    for (n = 1; n <= Nmax - 1; n++)
      Fs[n] = Ftd[n] - (Fi[n] - Fi[n - 1]);
    Fs[0] = Ftd[0];

    for (size_t i = 0; i < Nmax + 1; ++i) {
      vec_to_smooth[i]->U[eq_num] = Fs[i];
    }

    delete[] Fs;
    delete[] Fi;
    delete[] Ftd;
  }
}

void mesh_and_methods_axisymm::fix_small_p_rho() {

  const double p_min = 0.001, rho_min = 0.001;

  for (size_t i = 0; i < mesh.size(); ++i) {
    for (size_t j = 0; j < mesh[i].size(); ++j) {
      if (mesh[i][j]->p < p_min) {
        std::cout << "x = " << j << " y = " << i << " p = "
         << mesh[i][j]->p << std::endl;
        mesh[i][j]->p = p_min;
        mesh[i][j]->calc_U_when_values_known();
        mesh[i][j]->calc_F_when_values_known();
        mesh[i][j]->calc_G_when_values_known();
        mesh[i][j]->calc_H_when_values_known();
      }
      if (mesh[i][j]->rho < rho_min) {
        std::cout << "x = " << j << " y = " << i << " rho = "
         << mesh[i][j]->rho << std::endl;
        mesh[i][j]->rho = rho_min;
        mesh[i][j]->calc_U_when_values_known();
        mesh[i][j]->calc_F_when_values_known();
        mesh[i][j]->calc_G_when_values_known();
      }
    }
  }

  return;
}

///////////////////////BOUNDARY CONDITIONS://///////////////////////////

//in boundary conditions, you need to operate with the whole nodes of mesh
//so don't forget they're pointers! more * motherfucker

void mesh_and_methods_axisymm::boundary_conditions_default() {

  // d/dn = 0 on top bound (r changes):
  for (size_t y = 1; y < y_begin_ind() + 1; ++y) {
    for (size_t x = 0; x < mesh[0].size(); ++x) {
      data_node_2d& node_to_change = *mesh[y_end_ind() - 1 + y][x];
      double temp_r = node_to_change.r;
      node_to_change = *mesh[y_end_ind() - 1][x];
      node_to_change.r = temp_r;
      //node_to_change.v *= -1;
      node_to_change.calc_U_when_values_known();
      node_to_change.calc_F_when_values_known();
      node_to_change.calc_G_when_values_known();
    }
  }

  // d/dn = 0 on left and right bounds (r doesn't change):
  for (size_t x = 1; x < x_begin_ind() + 1; ++x) {
    for (size_t y = 0; y < mesh.size(); ++y) {
      *mesh[y][x_begin_ind() - x] = *mesh[y][x_begin_ind()];
      *mesh[y][x_end_ind() - 1 + x] = *mesh[y][x_end_ind() - 1];
    }
  }

  // d/dn = 0 & v=0 on symmetry axis:
  for (size_t x = 0; x < mesh[0].size(); ++x) {
    data_node_2d& node_to_change = *mesh[y_begin_ind()][x];
    node_to_change = *mesh[y_begin_ind() + 1][x];
    node_to_change.v = 0.0;
    for (size_t i = 0; i < 4; ++i) {
      node_to_change.U[i] = 0.0;
      node_to_change.F[i] = 0.0;
      node_to_change.G[i] = 0.0;
    }
  }

  // ghost edges on bottom:
  for (size_t row_number = 1; row_number < y_begin_ind() + 1; ++row_number) {
    for (size_t x = 0; x < mesh[row_number].size(); ++x) {
      data_node_2d& node_to_change = *mesh[y_begin_ind() - row_number][x];
      node_to_change = *mesh[y_begin_ind() + row_number][x];
      for (size_t i = 0; i < 4; ++i) {
        node_to_change.U[i] *= -1.0;
        node_to_change.F[i] *= -1.0;
        node_to_change.G[i] *= -1.0;
      }
    }
  }
}


void mesh_and_methods_axisymm::
boundary_conditions_for_bubble_near_wall_simmetry_on_bottom() {
  // solid wall on the right:
  for (size_t x = 1; x < x_begin_ind() + 1; ++x) {
    for (size_t y = 0; y < mesh.size(); ++y) {
      data_node_2d& node_to_change = *mesh[y][x_end_ind() - 1 + x];
      node_to_change = *mesh[y][x_end_ind() - x];
      node_to_change.u *= -1;
      node_to_change.U[1] *= -1;
      node_to_change.F[0] *= -1;
      node_to_change.F[2] *= -1;
      node_to_change.F[3] *= -1;
      node_to_change.G[1] *= -1;
    }
  }


  // d/dn = 0 on top bound (r changes):
  for (size_t y = 1; y < y_begin_ind() + 1; ++y) {
    for (size_t x = 0; x < mesh[0].size(); ++x) {
      data_node_2d& node_to_change = *mesh[y_end_ind() - 1 + y][x];
      double temp_r = node_to_change.r;
      node_to_change = *mesh[y_end_ind() - 1][x];
      node_to_change.r = temp_r;
      node_to_change.v *= -1;
      node_to_change.calc_U_when_values_known();
      node_to_change.calc_F_when_values_known();
      node_to_change.calc_G_when_values_known();
    }
  }

  // d/dn = 0 on left bound (r doesn't change):
  for (size_t x = 1; x < x_begin_ind() + 1; ++x) {
    for (size_t y = 0; y < mesh.size(); ++y) {
      *mesh[y][x_begin_ind() - x] = *mesh[y][x_begin_ind()];
    }
  }

  // d/dn = 0 & v=0 on symmetry axis:
  for (size_t x = 0; x < mesh[0].size(); ++x) {
    data_node_2d& node_to_change = *mesh[y_begin_ind()][x];
    node_to_change = *mesh[y_begin_ind() + 1][x];
    node_to_change.v = 0.0;
    for (size_t i = 0; i < 4; ++i) {
      node_to_change.U[i] = 0.0;
      node_to_change.F[i] = 0.0;
      node_to_change.G[i] = 0.0;
    }
  }

  // ghost edges on bottom:
  for (size_t row_number = 1; row_number < y_begin_ind() + 1; ++row_number) {
    for (size_t x = 0; x < mesh[row_number].size(); ++x) {
      data_node_2d& node_to_change = *mesh[y_begin_ind() - row_number][x];
      node_to_change = *mesh[y_begin_ind() + row_number][x];
      for (size_t i = 0; i < 4; ++i) {
        node_to_change.U[i] *= -1.0;
        node_to_change.F[i] *= -1.0;
        node_to_change.G[i] *= -1.0;
      }
    }
  }
}
