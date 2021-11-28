#ifndef EXACT_SOLUTION_H
#define EXACT_SOLUTION_H

class exact_solution {
public:

  //initial values
  double x_left, x_right;

  double p1, rho1, u1,
    p2, rho2, u2;

  double a1, a2, v1, v2;       // v = 1/rho

  double gamma;

  size_t size;
  double t_end;

  //to calculate:
  double u3 = 0, p3,
    rho31, rho32,
    a31 = 0, a32 = 0,
    d1 = 0, d2 = 0;

  double delta_x;

  //left and right from contact discontinuity
  char left, right; //s - shock wave, r - riemann fan, v - vacuum
  int case_number;
  /* 0: Sh-Sh
   * 1: Sh-Rf
   * 2: Rf-Sh
   * 3: Rf-Rf
   * 4: Vacuum */

  double p_max = 1.0;

  double shift = 0.0; //for output: inside solver, interval must be symmetric

  std::vector<double> waves_velocities;

  exact_solution(double p1_, double rho1_, double u1_, double p2_, double rho2_, double u2_,
                 double gamma_, double t_end_, double x_left_, double x_right_, size_t size_)
    : p1(p1_), rho1(rho1_), u1(u1_), p2(p2_), rho2(rho2_), u2(u2_), gamma(gamma_),
      t_end(t_end_), x_left(x_left_), x_right(x_right_), size(size_) {
    //in this constructor we assume that interval is symmetric!
    v1 = 1.0/rho1;
    v2 = 1.0/rho2;
    a1 = std::sqrt(gamma*p1/rho1);
    a2 = std::sqrt(gamma*p2/rho2);
    delta_x = (x_right-x_left)/size;
    p_max = std::max(p1, p2);
  };

  exact_solution(double p1_, double rho1_, double u1_,
                 double p2_, double rho2_, double u2_,
                 double gamma_) :
    p1(p1_), rho1(rho1_), u1(u1_), p2(p2_), rho2(rho2_), u2(u2_), gamma(gamma_) {

    //in this constructor we assume that interval is symmetric!
    v1 = 1.0/rho1;
    v2 = 1.0/rho2;
    a1 = std::sqrt(gamma*p1/rho1);
    a2 = std::sqrt(gamma*p2/rho2);
    p_max = std::max(p1, p2);
    waves_velocities.reserve(5);
  }

  exact_solution(const data& dt) :
    rho1(dt.rho_left), p1(dt.p_left), u1(dt.u_left),
    rho2(dt.rho_right), p2(dt.p_right), u2(dt.u_right),
    gamma(dt.gamma), t_end(dt.t_end), size(dt.size){
    x_left = 0.5*(dt.x_left - dt.x_right);
    x_right = -x_left;
    shift = -0.5*(dt.x_left + dt.x_right);
    v1 = 1.0/rho1;
    v2 = 1.0/rho2;
    a1 = std::sqrt(gamma*p1/rho1);
    a2 = std::sqrt(gamma*p2/rho2);
    delta_x = (x_right-x_left)/size;
    p_max = std::max(p1, p2);
  }

  void shock_wave_u(double& f, double& f_der,
      double p_i, double u_i, double v_i, double p, bool right_dir);
  void riemann_fan_u(double& f, double& f_der,
      double p_i, double u_i, double a_i, double p, bool right_dir);


  void case_number_choice();
  void f_and_der_calc(double& f, double& f_der, double p);
  double newton_solver();

  void calc_p3();
  void calc_left_sh();
  void calc_right_sh();
  void calc_left_rf();
  void calc_right_rf();
  void calc_left_vac();
  void calc_right_vac();

  void calc_and_print_result(std::ofstream& out);
  double calc_max_velocity_and_waves_velocities();
  void find_values_in_zero(double& rho_disc, double& p_disc, double& u_disc);

};


#endif // EXACT_SOLUTION_H
