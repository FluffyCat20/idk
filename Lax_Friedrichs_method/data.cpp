#include "data.h"
#include "exact_solution.h"

void data_node::calc_F() {
  F[0] = rho*u[0];
  F[1] = p + rho*u[0]*u[0];
  switch (dimension) {
  case 1: {
    F[2] = (e + p)*u[0];
    break;
  }
  case 2: {
    F[2] = rho*u[0]*u[1];
    F[3] = (e + p)*u[0];
    break;
  }
  default:
    break;
  }
};

void data_node::calc_G() {
  G[0] = rho*u[1];
  G[1] = rho*u[0]*u[1];
  G[2] = p + rho*u[1]*u[1];
  G[3] = (e + p)*u[1];
}

void data_node::init_U_F_1D() {
  U.resize(3);
  F.resize(3);

  U[0] = rho;
  U[1] = rho*u[0];
  U[2] = e;

  calc_F();
}

void data_node::init_U_F_G_2D() {
  U.resize(4);
  F.resize(4);
  G.resize(4);

  U[0] = rho;
  U[1] = rho*u[0];
  U[2] = rho*u[1];
  U[3] = e;

  calc_F();
  calc_G();
}

int data_node::calc_F_in_node(){

  rho = U[0];
  if (rho < eps) {
    std::cout << "density <= 0 :(" << std::endl;
    return -1;
  }

  switch (dimension) {
  case 1: {
    u[0] = U[1]/rho;
    e = U[2];
    p = (e-rho*u[0]*u[0]*0.5)*(gamma-1);

    calc_F();
    break;
  }
  case 2: {
    u[0] = U[1]/rho;
    u[1] = U[2]/rho;
    e = U[3];
    p = (e-rho*u_abs*u_abs*0.5)*(gamma-1);

    calc_F();
    calc_G();
    break;
  }
  default: {
    break;
  }
  }

  return 0;
}

void data::get_data_from_config(std::ifstream& input) {

  json config_data;
  input >> config_data;

  dimension = config_data["dimension"];
  u_left.resize(dimension);
  u_right.resize(dimension);

  method_name = config_data["method_name"];

  x_left = config_data["bounds"]["x_left"];
  x_right = config_data["bounds"]["x_right"];

  if (dimension == 1) {
    y_bottom = 0.0;
    y_top = 1.0;
  } else {
    y_bottom = config_data["bounds"]["y_bottom"];
    y_top = config_data["bounds"]["y_top"];
  }

  rho_left = config_data["initial_values"]["rho_left"];
  rho_right = config_data["initial_values"]["rho_right"];
  p_left = config_data["initial_values"]["p_left"];
  p_right = config_data["initial_values"]["p_right"];

  u_left[0] = config_data["initial_values"]["u_left"];
  u_right[0] = config_data["initial_values"]["u_right"];
  if (dimension > 1) {
    u_left[1] = config_data["initial_values"]["v_left"];
    u_right[1] = config_data["initial_values"]["v_right"];
  }

  size_x = config_data["nodes_number_x"];
  if (dimension > 1) {
    size_y = config_data["nodes_number_y"];
  }
  courant_number = config_data["courant_number"];
  if (config_data["gamma_is_const"]){
    gamma = config_data["gamma"];
  }
  t_end = config_data["t_end"];
}

double data::calc_delta_t(){

  //delete mesh_1D later
  double max_u_abs_plus_a = std::numeric_limits<double>::min();
  switch (dimension) {
  case 1: {
    for (size_t i = 0; i < size_x; ++i) {
      if (std::abs(mesh_1D[i].u[0]) + mesh_1D[i].a > max_u_abs_plus_a)
        max_u_abs_plus_a = std::abs(mesh_1D[i].u[0]) + mesh_1D[i].a;
    }
    break;
  }
  case 2: {
    for (const auto& row : mesh) {
      for (const auto& node : row) {
        if (node.u_abs + node.a > max_u_abs_plus_a)
          max_u_abs_plus_a = node.u_abs + node.a;
      }
    }
    break;
  }
  default: {
    break;
  }
  }

  return delta_x/max_u_abs_plus_a*courant_number;
}

void data::calc_new_U_with_lax_friedrichs(const data& prev_grid){
  for (size_t k = 1; k < size_x - 1; ++k){
    for (size_t i = 0; i < 3; ++i){
      mesh_1D[k].U[i] = 0.5*(prev_grid.mesh_1D[k-1].U[i] + prev_grid.mesh_1D[k+1].U[i]) //delta_t from this* must be calculated before!!
          - delta_t/delta_x*0.5*(prev_grid.mesh_1D[k+1].F[i] - prev_grid.mesh_1D[k-1].F[i]);
    }
  }
  mesh_1D[0].U = mesh_1D[1].U;
  mesh_1D[size_x - 1].U = mesh_1D[size_x - 2].U;
}

double inner_product(const std::vector<double>& v1, const std::vector<double>& v2){
  double res = 0;
  for (size_t i = 0; i < v1.size(); ++i){
    res += v1[i]*v2[i];
  }
  return res;
}

inline double phi(double r){ //flux limiter
  if (r > 0.0)
    return std::min(2*r, 1.0);
  return 0.0;
}

void data::calc_D(){
  std::vector<double> delta_U_prev(3), delta_U_cur(3), delta_U_next(3);
  for (size_t i = 0; i < 3; ++i){
    delta_U_prev[i] = mesh_1D[1].U[i] - mesh_1D[0].U[i];
    delta_U_cur[i] = mesh_1D[2].U[i] - mesh_1D[1].U[i];
  }
  for (size_t k = 1; k < size_x - 2; ++k){
    double nu_k = delta_t/delta_x*std::max(
                  std::abs(mesh_1D[k].u[0] + mesh_1D[k].a),
                  std::abs(mesh_1D[k].u[0] - mesh_1D[k].a)); //  dt/dx*max|lambda|
    double nu_k_plus_1 = delta_t/delta_x*std::max(
                         std::abs(mesh_1D[k+1].u[0] + mesh_1D[k+1].a),
                         std::abs(mesh_1D[k+1].u[0] - mesh_1D[k+1].a));
    double c_k = 0.25, c_k_plus_1 = 0.25;
    if (nu_k <= 0.5) {
      c_k = nu_k*(1 - nu_k);
    }
    if (nu_k_plus_1 <= 0.5){
      c_k_plus_1 = nu_k_plus_1*(1 - nu_k_plus_1);
    }
    for (size_t i = 0; i < 3; ++i){
      delta_U_next[i] = mesh_1D[k+2].U[i] - mesh_1D[k+1].U[i];
    }
    double denom = inner_product(delta_U_cur, delta_U_cur);
    if (std::abs(denom) < eps){
      for (size_t i = 0; i < 3; ++i)
        D[k][i] = 0.0;

      delta_U_prev = delta_U_cur;
      delta_U_cur = delta_U_next;
      continue;
    }
    double r_k_Plus = inner_product(delta_U_prev, delta_U_cur)/denom,
        r_k_plus_1_Minus = inner_product(delta_U_cur, delta_U_next)/denom;

    for (size_t i = 0; i < 3; ++i)
      D[k][i] = (c_k*(1 - phi(r_k_Plus)) +
                 c_k_plus_1*(1 - phi(r_k_plus_1_Minus)))*delta_U_cur[i]; //*0.5 ??

    //for next step:
    delta_U_prev = delta_U_cur;
    delta_U_cur = delta_U_next;
  }
}

void data::mac_cormack_predictor_step(const data& prev_grid){

  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      switch (dimension) {
      case 1: {
        for (size_t k = 0; k < 3; ++k) {
          mesh[i][j].U[k] = prev_grid.mesh[i][j].U[k] -
            delta_t/delta_x*(prev_grid.mesh[i][j+1].F[k] - prev_grid.mesh[i][j].F[k]);
        }
        break;
      }
      case 2: {
        for (size_t k = 0; k < 4; ++k) {
          mesh[i][j].U[k] = prev_grid.mesh[i][j].U[k] -
            delta_t/delta_x*(prev_grid.mesh[i+1][j].F[k] +
            prev_grid.mesh[i][j+1].F[k] - 2.0*prev_grid.mesh[i][j].F[k]);
        }
        break;
      }
      default:
        break;
      }
    }
  }
  calc_F_in_nodes();
}

void data::mac_cormack_corrector_step(const data& prev_grid, const data& predictor_grid) {

  for (size_t i = 1; i < size_y - 1; ++i) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      switch (dimension) {
      case 1: {
        for (size_t k = 0; k < 3; ++k) {
          mesh[i][j].U[k] = 0.5*(prev_grid.mesh[i][j].U[k] + predictor_grid.mesh[i][j].U[k])
            - 0.5*delta_t/delta_x*(predictor_grid.mesh[i][j].F[k] - predictor_grid.mesh[i][j - 1].F[k]);
        }
        break;
      }
      case 2: {
        for (size_t k = 0; k < 4; ++k) {
          mesh[i][j].U[k] = 0.5*(prev_grid.mesh[i][j].U[k] + predictor_grid.mesh[i][j].U[k])
            - 0.5*delta_t/delta_x*(2.0*predictor_grid.mesh[i][j].F[k]
            - predictor_grid.mesh[i][j - 1].F[k] - predictor_grid.mesh[i - 1][j].F[k]);
        }
        break;
      }
      default:
        break;
      }
    }
  }

  calc_F_in_nodes();

}

void data::calc_new_U_with_mac_cormack(data& prev_grid){
  data predictor_grid(*this);
  predictor_grid.mac_cormack_predictor_step(*this);
  mac_cormack_corrector_step(prev_grid, predictor_grid);
  calc_boundary_conditions();
}

void data::calc_new_U_with_mac_cormack_davis(data& prev_grid){
  data predictor_grid(*this);
  predictor_grid.mac_cormack_predictor_step(*this);

  mac_cormack_corrector_step(prev_grid, predictor_grid);

  //Davis artiificial viscosity terms:
  calc_D();
  for (size_t k = 1; k < size_x - 1; ++k) {
    for (size_t i = 0; i < 3; ++i) {
      mesh_1D[k].U[i] += D[k][i] - D[k-1][i];
    }
  }

  //Smooth_Davis_Maksimov(*this, prev_grid);
}

void data::calc_new_U_and_F_with_godunov(data& prev_grid, double current_t) {

  double d_max = 0.0;

  calc_F_btw_nodes(d_max);

  delta_t = courant_number*prev_grid.delta_x/d_max; //delta t but using riemann task exact solution

  if (current_t + delta_t > t_end){
    delta_t = t_end - current_t;
  }

  for (size_t k = 1; k < size_x - 1; ++k) {
    for (size_t i = 0; i < 3; ++i) {
      mesh_1D[k].U[i] = prev_grid.mesh_1D[k].U[i] - delta_t/delta_x*
                     (mesh_1D[k].F[i] - mesh_1D[k-1].F[i]);
    }
  }

  mesh_1D[0].U = mesh_1D[1].U;
  mesh_1D[size_x - 1].U = mesh_1D[size_x - 2].U;

  for (size_t k = 0; k < size_x; ++k) {
    mesh_1D[k].rho = mesh_1D[k].U[0];
    if (mesh_1D[k].rho < eps) {
      stop_now = true;
      std::cout << "rho < 0" << std::endl;
      break;
    }
    mesh_1D[k].u[0] = mesh_1D[k].U[1]/mesh_1D[k].rho;
    mesh_1D[k].e = mesh_1D[k].U[2];
    mesh_1D[k].p = (mesh_1D[k].e -
      mesh_1D[k].rho*mesh_1D[k].u[0]*mesh_1D[k].u[0]*0.5)*(gamma - 1);
  }

}


void data::calc_F_in_nodes() {
  if (dimension == 1) {
    for (size_t k = 0; k < size_x; ++k){
      if (mesh_1D[k].calc_F_in_node() == -1){
        stop_now = true;
        std::cout << k << std::endl;
        break;
      }
    }
  }

  for (size_t i = 0; i < size_y; ++i) {
    for (size_t j = 0; j < size_x; ++j) {
      if (mesh[i][j].calc_F_in_node() == -1) {
        stop_now = true;
        std::cout << i << " " << j << std::endl;
        break;
      }
    }
  }

}

void data::calc_F_btw_nodes(double& max_velocity_abs) {

  max_velocity_abs = 0.0;

  //F[0] is flux btw U[0] and U[1], so last needed F is F[U.size - 2]
  for (size_t k = 0; k < size_x - 1; ++k) {
    const data_node& U_left = mesh_1D[k];
    const data_node& U_right = mesh_1D[k+1];
    exact_solution disc_sln(U_left.p, U_left.rho, U_left.u[0],
                            U_right.p, U_right.rho, U_right.u[0],
                            gamma);
    max_velocity_abs =
        disc_sln.calc_max_velocity_and_waves_velocities();
    std::vector<double>& flux_to_calc = mesh_1D[k].F;
    double rho_disc, p_disc, u_disc;
    disc_sln.find_values_in_zero(rho_disc, p_disc, u_disc);
    flux_to_calc[0] = rho_disc*u_disc;
    flux_to_calc[1] = rho_disc*u_disc*u_disc + p_disc;
    double e_disc = p_disc/(gamma - 1) + 0.5*rho_disc*u_disc*u_disc;
    flux_to_calc[2] = u_disc*(e_disc + p_disc);
  }
}


void data::calc_boundary_conditions() {
  for (size_t i = 0; i < size_y; ++i) {
    mesh[i][0] = mesh[i][1];
    mesh[i][size_x - 1] = mesh[i][size_x - 2];
  }
  if (dimension > 1) {
    for (size_t j = 1; j < size_x - 1; ++j) {
      mesh[0][j] = mesh[1][j];
      mesh[size_y - 1][j] = mesh[size_y - 2][j];
    }
  }
}

data& data::operator=(const data& other) { //only changes mesh and delta_t
  mesh_1D = other.mesh_1D;
  mesh = other.mesh;
  delta_t = other.delta_t;
  return *this;
}

///////////////////////////////////////////////////////////////

double phi2(double x) //взято у Максимова
{
  return std::max(0.0, std::min(2 * std::abs(x), 1.0));
}

int Smooth_Davis_Maksimov(data& new_grid, const data& prev_grid) //алгоритм Андрея Максимова
{
  #define INT_MAX_EQUATION_NUMBER 3
  double numerator, denominator_plus, denominator_minus;

  std::vector<std::vector<double>> DZ (INT_MAX_EQUATION_NUMBER, std::vector<double>(prev_grid.size_x, 0));
  //Z - это U; DZ = delta U - двумерный массив размера [число уравнений]*[размер сетки]
  std::vector<double> RP(prev_grid.size_x, 0); //
  std::vector<double> RM(prev_grid.size_x, 0); //массивы для r+ и r-

  for (int k = 0; k < prev_grid.size_x - 1; k++) //вычисление DZ и сохранение в массив
    for (int int_equation_number = 0; int_equation_number < INT_MAX_EQUATION_NUMBER; int_equation_number++)
      DZ[int_equation_number][k] = new_grid.mesh_1D[k+1].U[int_equation_number] - new_grid.mesh_1D[k].U[int_equation_number];
  for (int k = 1; k < prev_grid.size_x; k++) //вычисление r+ и r- и сохранение в массивы
  {
    numerator = 0;
    denominator_plus = 0;
    denominator_minus = 0;
    for (int int_equation_number = 0; int_equation_number < INT_MAX_EQUATION_NUMBER; int_equation_number++)
    {
      numerator += DZ[int_equation_number][k - 1] * DZ[int_equation_number][k];
      denominator_plus += DZ[int_equation_number][k] * DZ[int_equation_number][k];
      denominator_minus += DZ[int_equation_number][k - 1] * DZ[int_equation_number][k - 1];
    }
    RP[k] = numerator / denominator_plus;
    RM[k] = numerator / denominator_minus;
  }
  for (int k = 1; k < prev_grid.size_x - 1; k++)
  {

    //в версии Максимова double_Courant_number было, видимо, глобальной переменной и никак не определялось,
    //а double_C_nu считалось не в цикле, а снаружи

    double double_Courant_number = prev_grid.delta_t/prev_grid.delta_x*
        std::max(std::abs(prev_grid.mesh_1D[k].u[0] + prev_grid.mesh_1D[k].a),
        std::abs(prev_grid.mesh_1D[k].u[0] - prev_grid.mesh_1D[k].a)); //  dt/dx*max|lambda|

    double double_C_nu;

    if (double_Courant_number > 0.5)
      double_C_nu = 0.25;
    else
      double_C_nu = double_Courant_number*(1 - double_Courant_number);

    //new grid - сетка после МакКормака, прибавляем к ней

    for (int int_equation_number = 0; int_equation_number < INT_MAX_EQUATION_NUMBER; int_equation_number++)
    {
      //множитель 0.5 убран, т.к. с ним остаются некоторые осцилляции
      new_grid.mesh_1D[k].U[int_equation_number] += double_C_nu*((2 - phi2(RP[k]) - phi2(RM[k + 1]))*DZ[int_equation_number][k]
        - (2 - phi2(RP[k - 1]) - phi2(RM[k]))*DZ[int_equation_number][k - 1]);
    }
  }
  #undef INT_MAX_EQUATION_NUMBER
  return 0;
}
