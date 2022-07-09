#include "data_2d.h"
#include <functional>

void calc_state_after_shock_wave(
    double& rho, double& p, double& u,
    const double rho0, const double p0, const double mach) {
  /// no index - after shock wave
  /// index 0 - before shock wave

  const double gamma = data_node_2d::gamma;

  const double a0 = std::sqrt(gamma*p0/rho0);
  double u0 = mach*a0;
  p = p0*(2*gamma/(gamma+1)*mach*mach - (gamma-1)/(gamma+1));
  rho = 1/(rho0*((gamma-1)/(gamma+1)) + 2/(gamma+1)/mach/mach);
  u = u0*rho0/rho;
}

double newton_solver(
    std::function<double(double)> func,
    std::function<double(double)> func_der) {
  const double max_steps_number = 1000;
  double p0 = 100.0;
  double eps = 1e-6;
  double f, f_der;
  size_t steps_number = 0;
  //f_and_der_calc(f, f_der, p0);
  f = func(p0);
  f_der = func_der(p0);

  double delta = p0;

  while (std::abs(delta/p0) > eps && steps_number < max_steps_number){
    delta = f/f_der;
    p0 -= delta;
    //f_and_der_calc(f, f_der, p0);
    f = func(p0);
    f_der = func_der(p0);
    ++steps_number;
  }
  if (steps_number >= max_steps_number){
    return std::numeric_limits<double>::max();
  }
  //std::cout << "steps: " << steps_number <<std::endl;
  return p0;
}
