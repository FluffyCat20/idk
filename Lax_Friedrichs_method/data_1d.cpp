#include "data_1d.h"
#include "exact_solution.h"

int data_node::calc_F_in_node(){
  rho = U[0];
  if (rho < eps) {
    std::cout << "density <= 0 :(" << std::endl;
    return -1;
  }
  u = U[1]/rho;
  e = U[2];
  p = (e - rho*u*u*0.5)*(gamma-1);

  F[0] = U[1];
  F[1] = p + rho*u*u;
  F[2] = (e + p)*u;

  return 0;
}

void data::get_data_from_config(std::ifstream& input){

  json config_data;
  input >> config_data;

  method_name = config_data["method_name"];

  x_left = config_data["bounds"]["x_left"];
  x_right = config_data["bounds"]["x_right"];

  rho_left = config_data["initial_values"]["rho_left"];
  rho_right = config_data["initial_values"]["rho_right"];
  p_left = config_data["initial_values"]["p_left"];
  p_right = config_data["initial_values"]["p_right"];
  u_left = config_data["initial_values"]["u_left"];
  u_right = config_data["initial_values"]["u_right"];

  size = config_data["nodes_number"];
  courant_number = config_data["courant_number"];
  if (config_data["gamma_is_const"]){
    gamma = config_data["gamma"];
  }
  t_end = config_data["t_end"];
}

double data::calc_delta_t(){
  //mb abs(u +- a) instead abs(u) +- a ???
  double max_u_abs_plus_a = std::numeric_limits<double>::min();
  for (size_t i = 0; i < size; ++i){
    if (std::abs(mesh[i].u) + mesh[i].a > max_u_abs_plus_a)
      max_u_abs_plus_a = std::abs(mesh[i].u) + mesh[i].a;
  }
  return delta_x/max_u_abs_plus_a*courant_number;
}

void data::calc_new_U_with_lax_friedrichs(const data& prev_grid){
  for (size_t k = 1; k < size - 1; ++k){
    for (size_t i = 0; i < 3; ++i){
      mesh[k].U[i] = 0.5*(prev_grid.mesh[k-1].U[i] + prev_grid.mesh[k+1].U[i]) //delta_t from this* must be calculated before!!
          - delta_t/delta_x*0.5*(prev_grid.mesh[k+1].F[i] - prev_grid.mesh[k-1].F[i]);
    }
  }
  mesh[0].U = mesh[1].U;
  mesh[size - 1].U = mesh[size - 2].U;
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
    delta_U_prev[i] = mesh[1].U[i] - mesh[0].U[i];
    delta_U_cur[i] = mesh[2].U[i] - mesh[1].U[i];
  }
  for (size_t k = 1; k < size - 2; ++k){
    double nu_k = delta_t/delta_x*std::max(std::abs(mesh[k].u + mesh[k].a), std::abs(mesh[k].u - mesh[k].a)); //  dt/dx*max|lambda|
    double nu_k_plus_1 = delta_t/delta_x*std::max(std::abs(mesh[k+1].u + mesh[k+1].a), std::abs(mesh[k+1].u - mesh[k+1].a));
    double c_k = 0.25, c_k_plus_1 = 0.25;
    if (nu_k <= 0.5){
      c_k = nu_k*(1 - nu_k);
    }
    if (nu_k_plus_1 <= 0.5){
      c_k_plus_1 = nu_k_plus_1*(1 - nu_k_plus_1);
    }
    for (size_t i = 0; i < 3; ++i){
      delta_U_next[i] = mesh[k+2].U[i] - mesh[k+1].U[i];
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
      D[k][i] = (c_k*(1 - phi(r_k_Plus)) + c_k_plus_1*(1 - phi(r_k_plus_1_Minus)))*delta_U_cur[i]; //*0.5 ??

    //for next step:
    delta_U_prev = delta_U_cur;
    delta_U_cur = delta_U_next;
  }
}

void data::mac_cormack_predictor_step(const data& prev_grid){
  for (size_t k = 1; k < size - 1; ++k){
    for (size_t i = 0; i < 3; ++i){
      mesh[k].U[i] = prev_grid.mesh[k].U[i] - delta_t/delta_x*(prev_grid.mesh[k+1].F[i] - prev_grid.mesh[k].F[i]);
    }
  }
  calc_F_in_nodes();
}

void data::mac_cormack_corrector_step(const data& prev_grid, const data& predictor_grid){
  for (size_t k = 1; k < size - 1; ++k){
    for (size_t i = 0; i < 3; ++i){
      mesh[k].U[i] = 0.5*(prev_grid.mesh[k].U[i] + predictor_grid.mesh[k].U[i])
        - 0.5*delta_t/delta_x*(predictor_grid.mesh[k].F[i] - predictor_grid.mesh[k-1].F[i]);
    }
  }
  mesh[0].U = mesh[1].U;
  mesh[size - 1].U = mesh[size - 2].U;

}

void data::calc_new_U_with_mac_cormack(data& prev_grid){
  data predictor_grid(*this);
  predictor_grid.mac_cormack_predictor_step(*this);
  mac_cormack_corrector_step(prev_grid, predictor_grid);
}

void data::calc_new_U_with_mac_cormack_davis(data& prev_grid){
  data predictor_grid(*this);
  predictor_grid.mac_cormack_predictor_step(*this);
  mac_cormack_corrector_step(prev_grid, predictor_grid);

  //Davis artiificial viscosity terms:
  calc_D();
  for (size_t k = 1; k < size - 1; ++k) {
    for (size_t i = 0; i < 3; ++i) {
      mesh[k].U[i] += D[k][i] - D[k-1][i];
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

  for (size_t k = 1; k < size - 1; ++k) {
    for (size_t i = 0; i < 3; ++i) {
      mesh[k].U[i] = prev_grid.mesh[k].U[i] - delta_t/delta_x*
                     (mesh[k].F[i] - mesh[k-1].F[i]);
    }
  }

  mesh[0].U = mesh[1].U;
  mesh[size - 1].U = mesh[size - 2].U;

  for (size_t k = 0; k < size; ++k) {
    mesh[k].rho = mesh[k].U[0];
    if (mesh[k].rho < eps) {
      stop_now = true;
      std::cout << "rho < 0" << std::endl;
      break;
    }
    mesh[k].u = mesh[k].U[1]/mesh[k].rho;
    mesh[k].e = mesh[k].U[2];
    mesh[k].p = (mesh[k].e - mesh[k].rho*mesh[k].u*mesh[k].u*0.5)*(gamma - 1);
  }

}


void data::calc_F_in_nodes() {
  for (size_t k = 1; k < size - 1; ++k){
    if (mesh[k].calc_F_in_node() == -1){
      stop_now = true;
      std::cout << k << std::endl;
      break;
    };
  }


  mesh[0].F = mesh[1].F;
  mesh[size - 1].F = mesh[size - 2].F;

}

void data::calc_F_btw_nodes(double& max_velocity_abs) {

  max_velocity_abs = 0.0;

  //F[0] is flux btw U[0] and U[1], so last needed F is F[U.size - 2]
  for (size_t k = 0; k < size - 1; ++k) {
    const data_node& U_left = mesh[k];
    const data_node& U_right = mesh[k+1];
    exact_solution disc_sln(U_left.p, U_left.rho, U_left.u,
                            U_right.p, U_right.rho, U_right.u,
                            gamma);
    max_velocity_abs =
        disc_sln.calc_max_velocity_and_waves_velocities();
    std::vector<double>& flux_to_calc = mesh[k].F;
    double rho_disc, p_disc, u_disc;
    disc_sln.find_values_in_zero(rho_disc, p_disc, u_disc);
    flux_to_calc[0] = rho_disc*u_disc;
    flux_to_calc[1] = rho_disc*u_disc*u_disc + p_disc;
    double e_disc = p_disc/(gamma - 1) + 0.5*rho_disc*u_disc*u_disc;
    flux_to_calc[2] = u_disc*(e_disc + p_disc);
  }
}

data& data::operator=(const data& other){ //only changes field and delta_t
  for (size_t i = 0; i < size; ++i)
    mesh[i] = other.mesh[i];
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

  std::vector<std::vector<double>> DZ (INT_MAX_EQUATION_NUMBER, std::vector<double>(prev_grid.size, 0));
  //Z - это U; DZ = delta U - двумерный массив размера [число уравнений]*[размер сетки]
  std::vector<double> RP(prev_grid.size, 0); //
  std::vector<double> RM(prev_grid.size, 0); //массивы для r+ и r-

  for (int k = 0; k < prev_grid.size - 1; k++) //вычисление DZ и сохранение в массив
    for (int int_equation_number = 0; int_equation_number < INT_MAX_EQUATION_NUMBER; int_equation_number++)
      DZ[int_equation_number][k] = new_grid.mesh[k+1].U[int_equation_number] - new_grid.mesh[k].U[int_equation_number];
  for (int k = 1; k < prev_grid.size; k++) //вычисление r+ и r- и сохранение в массивы
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
  for (int k = 1; k < prev_grid.size - 1; k++)
  {

    //в версии Максимова double_Courant_number было, видимо, глобальной переменной и никак не определялось,
    //а double_C_nu считалось не в цикле, а снаружи

    double double_Courant_number = prev_grid.delta_t/prev_grid.delta_x*
           std::max(std::abs(prev_grid.mesh[k].u + prev_grid.mesh[k].a), std::abs(prev_grid.mesh[k].u - prev_grid.mesh[k].a)); //  dt/dx*max|lambda|

    double double_C_nu;

    if (double_Courant_number > 0.5)
      double_C_nu = 0.25;
    else
      double_C_nu = double_Courant_number*(1 - double_Courant_number);

    //new grid - сетка после МакКормака, прибавляем к ней

    for (int int_equation_number = 0; int_equation_number < INT_MAX_EQUATION_NUMBER; int_equation_number++)
    {
      //множитель 0.5 убран, т.к. с ним остаются некоторые осцилляции - это ваш комментарий
      new_grid.mesh[k].U[int_equation_number] += double_C_nu*((2 - phi2(RP[k]) - phi2(RM[k + 1]))*DZ[int_equation_number][k]
        - (2 - phi2(RP[k - 1]) - phi2(RM[k]))*DZ[int_equation_number][k - 1]);
    }
  }
  #undef INT_MAX_EQUATION_NUMBER
  return 0;
}
