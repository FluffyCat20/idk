#include <stdio.h>
#include <math.h>

//вспомогательные функции дл€ точного решени€ задачи –имана
double I1(double p, double ul, double al, double pl, double gl)
{
	double alpha = (gl - 1.0) / 2 / gl;
	return ul + 2.0 * al / (gl - 1.0)*(1.0 - pow(p / pl, alpha));
}
double I1prime(double p, double al, double pl, double gl)
{
	double alpha = (gl - 1.0) / 2 / gl;
	return -2.0 * al / (gl - 1.0)*alpha / pl * pow(p / pl, alpha - 1.0);
}
double I3(double p, double ur, double ar, double pr, double gr)
{
	double alpha = (gr - 1.0) / 2 / gr;
	return ur - 2.0 * ar / (gr - 1.0)*(1.0 - pow(p / pr, alpha));
}
double I3prime(double p, double ar, double pr, double gr)
{
	double alpha = (gr - 1.0) / 2 / gr;
	return 2.0 * ar / (gr - 1.0)*alpha / pr * pow(p / pr, alpha - 1.0);
}
double G1(double p, double ul, double al, double pl, double gl)
{
	double beta = (gl + 1.0) / (gl - 1.0);
	return ul + 2.0 * al / sqrt(2.0*gl*(gl - 1.0))*(1.0 - p / pl) / sqrt(1.0 + beta*p / pl);
}
double G1prime(double p, double al, double pl, double gl)
{
	double beta = (gl + 1.0) / (gl - 1.0);
	return 2.0 * al / sqrt(2.0*gl*(gl - 1.0))*((-1.0 / pl*sqrt(1.0 + beta*p / pl) - (1 - p / pl)*beta / pl / 2.0 / sqrt(1.0 + beta*p / pl)) / (1.0 + beta*p / pl));
}
double G3(double p, double ur, double ar, double pr, double gr)
{
	double beta = (gr + 1.0) / (gr - 1.0);
	return ur - 2.0 * ar / sqrt(2.0*gr*(gr - 1.0))*(1.0 - p / pr) / sqrt(1.0 + beta*p / pr);
}
double G3prime(double p, double ar, double pr, double gr)
{
	double beta = (gr + 1.0) / (gr - 1.0);
	return -2.0 * ar / sqrt(2.0*gr*(gr - 1.0))*((-1.0 / pr*sqrt(1.0 + beta*p / pr) - (1 - p / pr)*beta / pr / 2.0 / sqrt(1.0 + beta*p / pr)) / (1.0 + beta*p / pr));
}
double Phi(double p, double ul, double al, double pl, double gl, double ur, double ar, double pr, double gr)
{
	double ret = 0;
	if (p >= pl)
		ret += G1(p, ul, al, pl, gl);
	else
		ret += I1(p, ul, al, pl, gl);
	if (p >= pr)
		ret -= G3(p, ur, ar, pr, gr);
	else
		ret -= I3(p, ur, ar, pr, gr);
	return ret;
}
double Phiprime(double p, double al, double pl, double gl, double ar, double pr, double gr)
{
	double ret = 0;
	if (p >= pl)
		ret += G1prime(p, al, pl, gl);
	else
		ret += I1prime(p, al, pl, gl);
	if (p >= pr)
		ret -= G3prime(p, ar, pr, gr);
	else
		ret -= I3prime(p, ar, pr, gr);
	return ret;
}
double Riemann_Newton(double ul, double al, double pl, double gl, double ur, double ar, double pr, double gr, double p_start)
{
	int n, n_max = 100;
	double p_prev, p_next, eps = 0.000001, phi, phiprime, dp;
	p_prev = p_start;
	for (n = 1; n <= n_max;)
	{
		phi = Phi(p_prev, ul, al, pl, gl, ur, ar, pr, gr);
		if (fabs(phi) < eps)
			return p_prev;
		phiprime = Phiprime(p_prev, al, pl, gl, ar, pr, gr);
		dp = -phi / phiprime;
		p_next = p_prev + dp;
		p_prev = p_next;
		n++;
	}
	//если не сошлось
	printf("Riemann_Newton: n_max = %d reached.\n", n_max);
	return p_next;
}

//Zl - состо€ние слева от разрыва, (rho, u, p, gamma)
//Zr - состо€ние справа от разрыва, (rho, u, p, gamma)
//gamma = 1.4
//предполагаетс€, что разрыв в точке x=0
//R1, U1, P1 - массивы, куда запишетс€ решение
//t1 - момент времени, дл€ которого считать
//x_min, x_max - кра€ расчетной области
//N - число расчетных узлов
int Riemann_Solve(double Zl[4], double Zr[4], double *R1, double *U1, double *P1, double t1, double x_min, double x_max, double N)
{
	double a0, ap1;
	double tmp1, ps, us, as, rsl, rsr, w[5];
	double x, ksi;
	double beta0, betap1;

	#define r0 Zl[0] //rho слева
	#define u0 Zl[1] //u слева
	#define p0 Zl[2] //p слева
	#define g0 Zl[3] //gamma слева
	#define rp1 Zr[0] //rho справа
	#define up1 Zr[1] //u справа
	#define pp1 Zr[2] //p справа
	#define gp1 Zr[3] //gamma справа

	//exact via Newton
	beta0 = (g0 + 1.0) / (g0 - 1.0);
	betap1 = (gp1 + 1.0) / (gp1 - 1.0);
	a0 = sqrt(g0*p0 / r0);
	ap1 = sqrt(gp1*pp1 / rp1);
	if (up1 - u0 + 2 * a0 / (g0 - 1.0) + 2 * ap1 / (gp1 - 1.0)  < 0) //слишком сильный разлет, образуетс€ вакуум
		return 1;
	tmp1 = 0.5*(p0 + pp1);
	//нахождение давлени€ на контактном разрыве методом Ќьютона
	ps = Riemann_Newton(u0, a0, p0, g0, up1, ap1, pp1, gp1, tmp1);
	us = ps >= p0 ? G1(ps, u0, a0, p0, g0) : I1(ps, u0, a0, p0, g0);
	if (ps <= p0)
		rsl = r0 * pow(ps / p0, 1.0 / g0);
	else
		rsl = r0 * ((1.0 + beta0 * ps / p0) / (ps / p0 + beta0));
	if (ps <= pp1)
		rsr = rp1 * pow(ps / pp1, 1.0 / gp1);
	else
		rsr = rp1 * ((1.0 + betap1 * ps / pp1) / (ps / pp1 + betap1));
	//скорости ”¬ и границ ¬–
	if (ps > p0)
	{
		if (fabs(r0 - rsl) > 0.000001)
		{
			w[0] = (r0*u0 - rsl*us) / (r0 - rsl);
			w[1] = w[0];
		}
		else
		{
			as = sqrt(g0*ps / rsl);
			w[0] = us - as*sqrt(1.0 + (g0 + 1.0) / 2.0 / g0*(p0 / ps - 1));
			w[1] = w[0];
		}
	}
	else
	{
		as = sqrt(g0*ps / rsl);
		w[0] = u0 - a0;
		w[1] = us - as;
	}
	w[2] = us;
	if (ps > pp1)
	{
		if (fabs(rp1 - rsr) > 0.000001)
		{
		w[3] = (rp1*up1 - rsr*us) / (rp1 - rsr);
		w[4] = w[3];
		}
		else
		{
			as = sqrt(gp1*ps / rsr);
			w[3] = us + as*sqrt(1.0 + (gp1 + 1.0) / 2.0 / gp1*(pp1 / ps - 1));
			w[4] = w[3];
		}
	}
	else
	{
		as = sqrt(gp1*ps / rsr);
		w[3] = us + as;
		w[4] = up1 + ap1;
	}

	//запись точного решени€ в массивы
	for (int k = 0; k <= N - 1; k++)
	{
		x = x_min + k*(x_max - x_min) / N;
		ksi = x / t1;
		if (w[0] * t1 >= x)
		{
			R1[k] = r0;
			U1[k] = u0;
			P1[k] = p0;
		}
		else if (w[1] * t1 >= x)
		{
			U1[k] = ((g0 - 1.0)*u0 + 2.0*(a0 + ksi)) / (g0 + 1.0);
			R1[k] = pow(pow(r0, g0)*(U1[k] - ksi)*(U1[k] - ksi) / g0 / p0, 1 / (g0 - 1.0));
			P1[k] = p0*pow(R1[k] / r0, g0);
		}
		else if (w[2] * t1 >= x)
		{
			R1[k] = rsl;
			U1[k] = us;
			P1[k] = ps;
		}
		else if (w[3] * t1 >= x)
		{
			R1[k] = rsr;
			U1[k] = us;
			P1[k] = ps;
		}
		else if (w[4] * t1 >= x)
		{
			U1[k] = ((gp1 - 1.0)*up1 - 2.0*(ap1 - ksi)) / (gp1 + 1.0);
			R1[k] = pow(pow(rp1, gp1)*(ksi - U1[k])*(ksi - U1[k]) / gp1 / pp1, 1 / (gp1 - 1.0));
			P1[k] = pp1*pow(R1[k] / rp1, gp1);
		}
		else
		{
			R1[k] = rp1;
			U1[k] = up1;
			P1[k] = pp1;
		}
	}

	#undef r0
	#undef u0
	#undef p0
	#undef g0
	#undef rp1
	#undef up1
	#undef pp1
	#undef gp1
	return 0;
}

