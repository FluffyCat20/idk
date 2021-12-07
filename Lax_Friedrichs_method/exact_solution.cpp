#include "data.h"
#include "exact_solution.h"

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
  double p0 = (p1+p2)*0.5;
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

  case_number_choice();
  calc_p3();

  double x = 0;

  switch (left) {

  case 's': {

    calc_left_sh();

    for (x = x_left; x < std::min(x_right, d1*t_end); x+=delta_x) { //steady flow
      out << x - shift << " " << rho1 << " " << p1 << " " << u1 << std::endl;
    }
    for (; x < std::min(x_right, u3*t_end); x+=delta_x) { //btw shock wave and contact discontinuity
      out << x - shift << " " << rho31 << " " << p3 << " " << u3 << std::endl;
    }

    break;
  }

  case 'r': {
    for (x = x_left; x < std::min(x_right, (u1 - a1)*t_end); x+=delta_x) { //steady flow
      out << x - shift << " " << rho1 << " " << p1 << " " << u1 << std::endl;
    }

    calc_left_rf();

    for (; x < std::min(x_right, (u3 - a31)*t_end); x+=delta_x) { //riemann fan
      double u = (2.0*(x/t_end + a1) + u1*(gamma - 1))/(gamma + 1);
      double a = u - x/t_end;
      double rho = rho1*std::pow(a/a1, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }

    for (; x < std::min(x_right, u3*t_end); x+=delta_x) { //btw riemann fan and contact discontinuity
      out << x - shift << " " << rho31 << " " << p3 << " " << u3 << std::endl;
    }

    break;
  }

  case 'v': {
    for (x = x_left; x < (u1 - a1)*t_end; x+=delta_x) { //steady flow
      out << x - shift << " " << rho1 << " " << p1 << " " << u1 << std::endl;
    }

    rho31 = 0.0;
    a31 = 0.0;
    u3 = u1 + 2.0/(gamma - 1)*a1; //a3 = 0

    for (; x < u3*t_end; x+=delta_x) { //riemann fan
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

    calc_right_sh();

    for (; x < std::min(x_right, d2*t_end); x+=delta_x) { //btw contact discontinuity and shock wave
      out << x - shift << " " << rho32 << " " << p3 << " " << u3 << std::endl;
    }
    for (; x < x_right; x+=delta_x){ //steady flow
      out << x - shift << " " << rho2 << " " << p2 << " " << u2 << std::endl;
    }

    break;
  }

  case 'r': {

    calc_right_rf();

    for (; x < std::min((u3 + a32)*t_end, x_right); x+=delta_x) { //btw contact discontinuity and riemann fan
      out << x - shift << " " << rho32 << " " << p3 << " " << u3 << std::endl;
    }
    for (; x < std::min((u2 + a2)*t_end, x_right); x+=delta_x) { //riemann fan
      double u = (2.0*(x/t_end - a2) + u2*(gamma - 1))/(gamma + 1);
      double a = x/t_end - u;
      double rho = rho2*std::pow(a/a2, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }
    for (; x < x_right; x+=delta_x) { //steady flow
      out << x - shift << " " << rho2 << " " << p2 << " " << u2 << std::endl;
    }

    break;
  }

  case 'v': {

    rho32 = 0.0;
    a32 = 0.0;
    u3 = u2 - 2.0/(gamma - 1)*a2; //a3 = 0

    for (; x < std::min(x_right, (u2 + a2)*t_end); x+=delta_x) { //riemann fan
      double u = (2.0*(x/t_end - a2) + u2*(gamma - 1))/(gamma + 1);
      double a = x/t_end - u;
      double rho = rho2*std::pow(a/a2, 2.0/(gamma - 1));
      double p = a*a*rho/gamma;

      out << x - shift << " " << rho << " " << p << " " << u << std::endl;
    }
    for (; x < x_right; x+=delta_x) { //steady flow
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

  case_number_choice();
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

  //now right == 'r' - riemann fan case
  if (0.0 > u2 + a2) { //inside riemann fan
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
