struct exact_solution{
  double p1, rho1, u1 = 0, a1,
    p3 = 0, rho3 = 0, u3 = 0, a3 = 0, //between rarefaction wave and contact discontinuity
    p4 = 0, rho4 = 0, u4 = 0, a4 = 0, //after shock wave
    p5, rho5, u5 = 0, a5;
  std::vector<double> p2, rho2, u2, a2; //after rarefaction wave
  size_t x1 = 0, x2 = 0, x3 = 0, x4 = 0;
  double shock_velocity = 0;
  const double gamma;
  const size_t size, shock_x;

  exact_solution(double rho1_, double p1_, double rho5_, double p5_, double gamma_, size_t size_, size_t shock_x_)
    : rho1(rho1_), p1(p1_), rho5(rho5_), p5(p5_), gamma(gamma_), size(size_), shock_x(shock_x_){
    a1 = std::sqrt(gamma*p1/rho1);
    a5 = std::sqrt(gamma*p5/rho5);
    shock_velocity = newton_solver();
    p4 = calc_p4();
    rho4 = calc_rho4();
    u4 = calc_u4();
    p3 = p4;
    u3 = u4;
    rho3 = calc_rho3();
    a3 = std::sqrt(gamma*p3/rho3);
    a4 = std::sqrt(gamma*p4/rho4);
  }

  double calc_func_for_newton(double d);
  double calc_func_deriv_for_newton(double d);
  double newton_solver(); //for area 4
  double calc_u4(){return 2/(gamma + 1)*(shock_velocity - 1/shock_velocity);}
  double calc_p4(){return p5*(2*gamma*std::pow(shock_velocity, 2) - gamma + 1)/(gamma + 1);}
  double calc_rho4();
  double calc_rho3(){return rho1*std::pow(p3/p1, 1/gamma);};
  void calc_bounds(double t, double delta_x);
  void calc_area2(double t, double delta_x);
  void output (std::ofstream& fout, double delta_x);

};

double exact_solution::calc_func_for_newton(double d){
  double to_pow = p5/p1/(gamma+1)*(2*gamma*d*d - gamma + 1);
  return d - 1/d - a1*(gamma+1)/(gamma-1)*(1 - std::pow(to_pow, (gamma - 1)*0.5/gamma));
}

double exact_solution::calc_func_deriv_for_newton(double d){
  double to_pow = p5/p1/(gamma+1)*(2*gamma*d*d - gamma + 1);
  return 1 + 1/d/d + a1*p5/p1*2*d*std::pow(to_pow, -(gamma+1)*0.5/gamma);
}

double exact_solution::newton_solver(){
  double x0 = 2.0;
  double f = calc_func_for_newton(x0);
  while (std::abs(f) > eps){
    double f_der = calc_func_deriv_for_newton(x0);
    x0 -= f/f_der;
    f = calc_func_for_newton(x0);
  }
  return x0;
}

double exact_solution::calc_rho4(){
  double rho5_over_rho4 = (2/std::pow(shock_velocity, 2) + gamma - 1)/(gamma + 1);
  return 1/rho5_over_rho4*rho5;
}

void exact_solution::calc_bounds(double t, double delta_x){
  x4 = shock_x + size_t(shock_velocity*t/delta_x);
  x3 = shock_x + size_t(u3*t/delta_x);
  x2 = shock_x + size_t((u3 - a3)*t/delta_x);
  x1 = shock_x - size_t(a1*t/delta_x);
  size_t area2_size = x2-x1;
  rho2.resize(area2_size);
  p2.resize(area2_size);
  u2.resize(area2_size);
  a2.resize(area2_size);
}

void exact_solution::calc_area2(double t, double delta_x){
  for (long x = x1; x < x2; ++x){
    long delta = x - (long)shock_x;
    u2[x - x1] = 2.0/(gamma+1)*(a1 + delta*delta_x/t);
    a2[x - x1] = (1 - gamma)*0.5*u2[x - x1] + a1;
    p2[x - x1] = p1*std::pow(a2[x - x1]/a1, 2*gamma/(gamma - 1));
    rho2[x - x1] = rho1*std::pow(p2[x - x1]/p1, 1/gamma);
  }
}

void exact_solution::output(std::ofstream& fout, double delta_x){
  fout << "x     rho     p     u\n";
  double x = 0.0;
  for (size_t i = 0; i < x1; ++i){
    fout << x << " " << rho1 << " "
      << p1 << " "
      << u1 << "\n";
    x += delta_x;
  }
  for (size_t i = 0; i < rho2.size(); ++i){
    fout << x << " " << rho2[i] << " "
      << p2[i] << " "
      << u2[i] << "\n";
    x += delta_x;
  }
  for (size_t i = 0; i < x3 - x2; ++i){
    fout << x << " " << rho3 << " "
      << p3 << " "
      << u3 << "\n";
    x += delta_x;
  }
  for (size_t i = 0; i < x4 - x3; ++i){
    fout << x << " " << rho4 << " "
      << p4 << " "
      << u4 << "\n";
    x += delta_x;
  }
  for (size_t i = 0; i < size - x4; ++i){
    fout << x << " " << rho5 << " "
      << p5 << " "
      << u5 << "\n";
    x += delta_x;
  }
}
