#include "data.h"
int data_node::calc_F(){
  rho = U[0];
  if (rho < eps) {
    std::cout << "density <= 0 :(" << std::endl;
    return -1;
  }
  u = U[1]/rho;
  e = U[2];
  p = (e-rho*u*u*0.5)*(gamma-1);

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
  calc_F();
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
  std::vector<std::vector<double>> fluxes(size - 1, std::vector<double>(3)); //F for prev_grid
  for (size_t k = 0; k < size - 1 ; ++k) {
    double u_disc, p_disc, rho_disc;
    exact_solution disc_sln(prev_grid.mesh[k].p, prev_grid.mesh[k].rho, prev_grid.mesh[k].u,
                       prev_grid.mesh[k + 1].p, prev_grid.mesh[k + 1].rho, prev_grid.mesh[k + 1].u,
                       prev_grid.gamma);
    disc_sln.case_number_choice();
    d_max = std::max(d_max, disc_sln.calc_max_velocity_and_waves_velocities());
    disc_sln.find_values_in_zero(rho_disc, p_disc, u_disc);
    double rho_u_u = rho_disc*u_disc*u_disc;
    double e_disc = p_disc/(gamma - 1) + rho_u_u*0.5;
    fluxes[k][0] = u_disc;
    fluxes[k][1] = p_disc + rho_u_u;
    fluxes[k][2] = (e_disc + p_disc)*u_disc;
  }

  delta_t = courant_number*prev_grid.delta_x/d_max; //delta t but using riemann task exact solution

  if (current_t + delta_t > t_end){
    delta_t = t_end - current_t;
  }

  for (size_t k = 1; k < size - 1; ++k) {  
    for (size_t i = 0; i < 3; ++i) {
      mesh[k].U[i] = prev_grid.mesh[k].U[i] - delta_t/prev_grid.delta_x*
                     (fluxes[k][i] - fluxes[k - 1][i]);
    }
  }

  mesh[0].U = mesh[1].U;
  mesh[size - 1].U = mesh[size - 2].U;

  for (size_t k = 0; k < size; ++k) {
    mesh[k].rho = mesh[k].U[0];
    mesh[k].u = mesh[k].U[1]/mesh[k].rho;
    mesh[k].e = mesh[k].U[2];
    mesh[k].p = (mesh[k].e - mesh[k].rho*mesh[k].u*mesh[k].u*0.5)*(gamma - 1);
  }

}


void data::calc_F(){
  for (size_t k = 1; k < size - 1; ++k){
    if (mesh[k].calc_F() == -1){
      stop_now = true;
      std::cout << k << std::endl;
      break;
    };
  }


  mesh[0].F = mesh[1].F;
  mesh[size - 1].F = mesh[size - 2].F;

}

data& data::operator=(const data& other){ //only changes field and delta_t
  for (size_t i = 0; i < size; ++i)
    mesh[i] = other.mesh[i];
  delta_t = other.delta_t;
  return *this;
}

void exact_solution::shock_wave_u(double& f, double& f_der,
    double p_i, double u_i, double v_i, double p, bool right_dir) {

  double sqrt = std::sqrt((gamma + 1)*p + (gamma - 1)*p_i);
  if (right_dir) {
    f = u_i + (p - p_i)*std::sqrt(2*v_i)/sqrt;
    f_der = std::sqrt(2*v_i)*0.5*((gamma + 1)*p + (3*gamma - 1)*p_i)/std::pow(sqrt, 3);
  }
  else {
    f = u_i - (p - p_i)*std::sqrt(2*v_i)/sqrt;
    f_der = -std::sqrt(2*v_i)*0.5*((gamma + 1)*p + (3*gamma - 1)*p_i)/std::pow(sqrt, 3);
  }
}

void exact_solution::riemann_fan_u(double& f, double& f_der,
    double p_i, double u_i, double a_i, double p, bool right_dir) {
  double pow = std::pow(p/p_i, (gamma - 1)*0.5/gamma);
  double pow_der = std::pow(p/p_i, -(gamma + 1)*0.5/gamma);
  double res = 2*a_i/(gamma - 1)*(1 - pow);
  if (right_dir) {
    f = u_i - res;
    f_der = a_i/gamma*pow_der/p_i;
  }
  else {
    f = u_i + res;
    f_der = -2*a_i/gamma*pow_der/p_i;
  }
}

void exact_solution::f_and_der_calc(double& f, double& f_der, double p){

  double f1, f_der1, //left dir
      f2, f_der2; //right dir
  if (p > p1){
    shock_wave_u(f1, f_der1, p1, u1, v1, p, false);
  }
  else {
    riemann_fan_u(f1, f_der1, p1, u1, a1, p, false);
  }

  if (p > p2){
    shock_wave_u(f2, f_der2, p2, u2, v2, p, true);
  }
  else {
    riemann_fan_u(f2, f_der2, p2, u2, a2, p, true);
  }

  f = f1 - f2;
  f_der = f_der1 - f_der2;

}

void exact_solution::case_number_choice() {
  double f_vac, f_vac_der;

  //vacuum
  f_and_der_calc(f_vac, f_vac_der, 0.0);
  if (f_vac < eps_for_rf_rf) {
    case_number = 4;
    left = 'v';
    right = 'v';
    return;
  }

  bool p1_is_max = true;
  if (p2 > p1)
    p1_is_max = false;

  double f_p1, f_p2, f_der_p1, f_der_p2;

  f_and_der_calc(f_p1, f_der_p1, p1);
  f_and_der_calc(f_p2, f_der_p2, p2);

  const double& f_p_max = (p1_is_max ? f_p1 : f_p2);
  const double& f_p_min = (p1_is_max ? f_p2 : f_p1);

  if (f_p_min <= eps_for_rf_rf){
    case_number = 3;
    left = 'r';
    right = 'r';
    return;
  }

  if (f_p_max > eps_for_rf_rf){
    case_number = 0;
    left = 's';
    right = 's';
    return;
  }

  if (!p1_is_max){
    case_number = 1;
    left = 's';
    right = 'r';
    return;
  }

  case_number = 2;
  left = 'r';
  right = 's';
  return;
}

double exact_solution::newton_solver() {
  const double max_steps_number = 1000;
  double p0 = p_max;
  double f, f_der;
  size_t steps_number = 0;
  f_and_der_calc(f, f_der, p0);

  double delta = p0;

  while (std::abs(delta/p0) > eps && steps_number < max_steps_number){
    delta = f/f_der;
    p0 -= delta;
    f_and_der_calc(f, f_der, p0);
    ++steps_number;
  }
  if (steps_number >= max_steps_number){
    return std::numeric_limits<double>::max();
  }
  //std::cout << "steps: " << steps_number <<std::endl;
  return p0;
}

void exact_solution::calc_p3() {

  switch (case_number) {
  case 3: {
    //p3: analitycal solution
    double power = (1-gamma)*0.5/gamma;
    double nmrtr = (gamma - 1)*0.5*(u1 - u2) + a1 + a2;
    double dnmntr = a1*std::pow(p1, power) + a2*std::pow(p2, power);
    p3 = std::pow(nmrtr/dnmntr, -1/power);
    break;
  }
  case 4: {
    p3 = 0.0;
    break;
  }
  default: {
    p3 = newton_solver();
    break;
  }
  }

}

void exact_solution::calc_left_sh() {

  double tmp = (gamma + 1)*p3+(gamma - 1)*p1;
  rho31 = rho1*tmp/((gamma + 1)*p1 + (gamma - 1)*p3);
  u3 = u1 - (p3 - p1)*std::sqrt(2.0/rho1/tmp);
  d1 = (rho1*u1 - rho31*u3)/(rho1 - rho31);

}

void exact_solution::calc_right_sh() {

  double tmp = (gamma + 1)*p3+(gamma - 1)*p2;
  rho32 = rho2*tmp/((gamma + 1)*p2 + (gamma - 1)*p3);
  u3 = u2 + (p3 - p2)*std::sqrt(2.0/rho2/tmp);
  d2 = (rho2*u2 - rho32*u3)/(rho2 - rho32);

}

void exact_solution::calc_left_rf() {

  rho31 = rho1*std::pow(p3/p1, 1.0/gamma);
  a31 = std::sqrt(gamma*p3/rho31);
  u3 = u1 + 2.0/(gamma - 1)*(a1 - a31);

}

void exact_solution::calc_right_rf() {

  rho32 = rho2*std::pow(p3/p2, 1.0/gamma);
  a32 = std::sqrt(gamma*p3/rho32);
  u3 = u2 - 2.0/(gamma - 1)*(a2 - a32);

}

void exact_solution::calc_left_vac() {

  rho31 = 0.0;
  a31 = 0.0;
  u3 = u1 + 2.0/(gamma - 1)*a1; //a3 = 0

}

void exact_solution::calc_right_vac() {

  rho32 = 0.0;
  a32 = 0.0;
  u3 = u2 + 2.0/(gamma - 1)*a2; //a3 = 0

}

void exact_solution::calc_and_print_result(std::ofstream& out) {

  calc_p3();

  double x = 0;

  switch (left) {

  case 's': {
    double tmp = (gamma+1)*p3+(gamma-1)*p1;
    rho31 = rho1*tmp/((gamma + 1)*p1 + (gamma - 1)*p3);
    u3 = u1 - (p3 - p1)*std::sqrt(2.0/rho1/tmp);
    d1 = (rho1*u1 - rho31*u3)/(rho1-rho31);
    for (x = x_left; x < std::min(x_right, d1*t_end); x+=delta_x){ //steady flow
      out << x - shift << " " << rho1 << " " << p1 << " " << u1 << std::endl;
    }
    for (; x < std::min(x_right, u3*t_end); x+=delta_x){ //btw shock wave and contact discontinuity
      out << x - shift << " " << rho31 << " " << p3 << " " << u3 << std::endl;
    }

    break;
  }

  case 'r': {
    for (x = x_left; x < std::min(x_right, (u1 - a1)*t_end); x+=delta_x){ //steady flow
      out << x - shift << " " << rho1 << " " << p1 << " " << u1 << std::endl;
    }

    rho31 = rho1*std::pow(p3/p1, 1.0/gamma);
    a31 = std::sqrt(gamma*p3/rho31);
    u3 = u1 + 2.0/(gamma - 1)*(a1 - a31);

    for (; x < std::min(x_right, (u3 - a31)*t_end); x+=delta_x){ //riemann fan
      double u = (2.0*(x/t_end + a1) + u1*(gamma - 1))/(gamma + 1);
      double a = u - x/t_end;
      double rho = rho1*std::pow(a/a1, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }

    for (; x < std::min(x_right, u3*t_end); x+=delta_x){ //btw riemann fan and contact discontinuity
      out << x - shift << " " << rho31 << " " << p3 << " " << u3 << std::endl;
    }

    break;
  }

  case 'v': {
    for (x = x_left; x < (u1 - a1)*t_end; x+=delta_x){ //steady flow
      out << x - shift << " " << rho1 << " " << p1 << " " << u1 << std::endl;
    }

    rho31 = 0.0;
    a31 = 0.0;
    u3 = u1 + 2.0/(gamma - 1)*a1; //a3 = 0

    for (; x < u3*t_end; x+=delta_x){ //riemann fan
      double u = (2.0*(x/t_end + a1) + u1*(gamma - 1))/(gamma + 1);
      double a = u - x/t_end;
      double rho = rho1*std::pow(a/a1, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }

    break;
  }
  default: {
    std::cout << "what is going left???" << std::endl;
    break;
  }
  }

  switch (right){

  case 's': {
    double tmp = (gamma+1)*p3+(gamma-1)*p2;
    rho32 = rho2*tmp/((gamma+1)*p2 + (gamma-1)*p3);
    u3 = u2 + (p3 - p2)*std::sqrt(2.0/rho2/tmp);
    d2 = (rho2*u2 - rho32*u3)/(rho2-rho32);

    for (; x < std::min(x_right, d2*t_end); x+=delta_x){ //btw contact discontinuity and shock wave
      out << x - shift << " " << rho32 << " " << p3 << " " << u3 << std::endl;
    }
    for (; x < x_right; x+=delta_x){ //steady flow
      out << x - shift << " " << rho2 << " " << p2 << " " << u2 << std::endl;
    }

    break;
  }

  case 'r': {

    rho32 = rho2*std::pow(p3/p2, 1.0/gamma);
    a32 = std::sqrt(gamma*p3/rho32);
    u3 = u2 - 2.0/(gamma - 1)*(a2 - a32);

    for (; x < std::min((u3 + a32)*t_end, x_right); x+=delta_x){ //btw contact discontinuity and riemann fan
      out << x - shift << " " << rho32 << " " << p3 << " " << u3 << std::endl;
    }
    for (; x < std::min((u2 + a2)*t_end, x_right); x+=delta_x){ //riemann fan
      double u = (2.0*(x/t_end - a2) + u2*(gamma - 1))/(gamma + 1);
      double a = x/t_end - u;
      double rho = rho2*std::pow(a/a2, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }
    for (; x < x_right; x+=delta_x){ //steady flow
      out << x - shift << " " << rho2 << " " << p2 << " " << u2 << std::endl;
    }

    break;
  }

  case 'v': {

    rho32 = 0.0;
    a32 = 0.0;
    u3 = u2 - 2.0/(gamma - 1)*a2; //a3 = 0

    for (; x < std::min(x_right, (u2 + a2)*t_end); x+=delta_x){ //riemann fan
      double u = (2.0*(x/t_end - a2) + u2*(gamma - 1))/(gamma + 1);
      double a = x/t_end - u;
      double rho = rho2*std::pow(a/a2, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }
    for (; x < x_right; x+=delta_x){ //steady flow
      out << x - shift << " " << rho2 << " " << p2 << " " << u2 << std::endl;
    }

    break;
  }

  default: {
    std::cout << "what is going right???" << std::endl;
    break;
  }


  }

}

double exact_solution::calc_max_velocity_and_waves_velocities() {

  calc_p3();

  switch (left) {
    case 's': {
      calc_left_sh();
      waves_velocities.push_back(d1);
      break;
    }

    case 'r': {
      calc_left_rf();
      waves_velocities.push_back(u1 - a1);
      waves_velocities.push_back(u3 - a31);
      break;
    }

    case 'v': {
      calc_left_vac();
      waves_velocities.push_back(u1 - a1);
      break;
    }
  }

  waves_velocities.push_back(u3); //contact discontinuity velocity (or vacuum)

  switch (right) {
    case 's': {
      calc_right_sh();
      waves_velocities.push_back(d2);
      break;
    }

    case 'r': {
      calc_right_rf();
      waves_velocities.push_back(u3 + a32);
      waves_velocities.push_back(u2 + a2);
      break;
    }

    case 'v': {
      calc_right_vac();
      waves_velocities.push_back(u2 + a2);
      break;
    }
  }

  double max_velocity = 0.0;

  for (auto v : waves_velocities) {
    double abs_vel = std::abs(v);
    if (abs_vel > max_velocity) {
      max_velocity = abs_vel;
    }
  }

  return max_velocity;
}

void exact_solution::find_values_in_zero(
    double& rho_disc, double& p_disc, double& u_disc) {
  //waves_velocities must be calculated already

  //if zero in left state:
  if (0.0 < waves_velocities.front()) {
    rho_disc = rho1;
    p_disc = p1;
    u_disc = u1;
    return;
  }

  //if zero in right state:
  if (0.0 > waves_velocities.back()) {
    rho_disc = rho2;
    p_disc = p2;
    u_disc = u2;
    return;
  }

  if (0.0 < u3) { //u3 is velocity of contact discontinuity
    if (left == 's' || left == 'v') {
      p_disc = p3;
      u_disc = u3;
      rho_disc = rho31;
      return;
    }
    //now left == 'r' - riemann fan case
    if (0.0 < waves_velocities[1]) { //inside rf
      u_disc = (2.0*a1 + u1*(gamma - 1))/(gamma + 1);
      rho_disc = rho1*std::pow(u_disc/a1, 2.0/(gamma - 1)); //u = a bc x = 0
      p_disc = u_disc*u_disc*rho_disc/gamma;
      return;
    }
    u_disc = u3;
    p_disc = p3;
    rho_disc = rho31;
    return;
  }

  //now 0 >= u3
  if (right == 's' || right == 'v') {
    p_disc = p3;
    u_disc = u3;
    rho_disc = rho32;
    return;
  }

  //now left == 'r' - riemann fan case
  if (0.0 > *(waves_velocities.end() - 2)) { //inside riemann fan
    u_disc = (-2.0*a2 + u2*(gamma - 1))/(gamma + 1);
    rho_disc = rho2*std::pow(-u_disc/a2, 2.0/(gamma - 1)); //u = -a
    p_disc = u_disc*u_disc*rho_disc/gamma;
    return;
  }

  u_disc = u3;
  p_disc = p3;
  rho_disc = rho32;
  return;
}

void uniform_check(double x1, double x2, double& uniform){
  double diff = std::abs(x1 - x2);
  if (diff > uniform){
    uniform = diff;
  }
}

std::vector<std::vector<double> > data::calc_error_norm(double* R1, double* P1, double* U1, size_t nodes_number){
  std::vector<double> uniform = {0.0, 0.0, 0.0}, //rho, p, u
      euclid = {0.0, 0.0, 0.0},
      deviation = {0.0, 0.0, 0.0},
      integral = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < nodes_number; ++i){
    uniform_check(R1[i], mesh[i].rho, uniform[0]);
    uniform_check(P1[i], mesh[i].p, uniform[1]);
    uniform_check(U1[i], mesh[i].u, uniform[2]);
    euclid[0] += std::pow(R1[i] - mesh[i].rho, 2);
    euclid[1] += std::pow(P1[i] - mesh[i].p, 2);
    euclid[2] += std::pow(U1[i] - mesh[i].u, 2);
    deviation[0] += std::pow(R1[i] - mesh[i].rho, 2);
    deviation[1] += std::pow(P1[i] - mesh[i].p, 2);
    deviation[2] += std::pow(U1[i] - mesh[i].u, 2);
    integral[0] += std::abs(R1[i] - mesh[i].rho);
    integral[1] += std::abs(P1[i] - mesh[i].p);
    integral[2] += std::abs(U1[i] - mesh[i].u);
  }
  for (size_t j = 0; j < 3; ++j){
    euclid[j] = std::sqrt(euclid[j]);
    deviation[j] = std::sqrt(deviation[j]/nodes_number);
    integral[j] *= delta_x;
  }
  std::vector<std::vector<double> > res = {uniform, euclid, deviation, integral};

  return res;
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
