#include "base_poisson_3d.hpp"
// #include <omp.h>
#include <math.h>
#include "parameters.hpp"

void base_poisson_3d::biotsavart_direct_3d(Particle &pi, Particle &pj)
{
	// internal variables
	double dx, dy, dz;
	double rij3, rij2, rij, sij, sij2;
	double gsig;

	// ProgressBar progressBar(pi.num, 35, '|', ' ');

	// #pragma omp parallel for private(j)
	for (int i = 1; i <= pi.num; i++)
	{
		// ++progressBar;		   // record the tick
		// progressBar.display(); // display the bar

		for (int j = 1; j <= pj.num; j++)
		{
			dx = pi.x[i - 1] - pj.x[j - 1];
			dy = pi.y[i - 1] - pj.y[j - 1];
			dz = pi.z[i - 1] - pj.z[j - 1];
			rij2 = std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2);
			rij = std::sqrt(rij2);
			rij3 = rij2 * rij;

			if (rij2 > 0.0e0)
			{
				sij2 = (std::pow(pi.s[i - 1], 2) + std::pow(pj.s[j - 1], 2)) * 0.25e0;
				sij = std::sqrt(sij2);
				regul_func_3d(rij, sij, gsig);

				pi.u[i - 1] = pi.u[i - 1] - (dy * pj.gz[j - 1] - dz * pj.gy[j - 1]) * gsig / rij3;
				pi.v[i - 1] = pi.v[i - 1] - (dz * pj.gx[j - 1] - dx * pj.gz[j - 1]) * gsig / rij3;
				pi.w[i - 1] = pi.w[i - 1] - (dx * pj.gy[j - 1] - dy * pj.gx[j - 1]) * gsig / rij3;
			}
		}
	}
	// progressBar.done();
}

void base_poisson_3d::regul_func_3d(double rij, double sij, double &q)
{
  double qt, qte, rho, rho2, rij2, sij2, sij3;

  rij2 = pow(rij, 2);
  sij2 = pow(sij, 2);
  sij3 = std::sqrt(pow(sij2, 3));
  rho2 = rij2 / sij2;
  rho = std::sqrt(rho2);

  // icutoff =  0.singular ; 1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian
  if (Parameters::icutoff == 0)
  {
    // q = 1.0e0 / (2.0e0 * Parameters::pi);
  }
  else if (Parameters::icutoff == 1)
  {
    // q = ((rho2 * (rho2 + 2.0e0)) / pow((rho2 + 1.0e0), 2)) / (2.0e0 * Parameters::pi);
  }
  else if (Parameters::icutoff == 2)
  {
    qt = -rho2 / (2.0e0);
	qte = std::erf(rho/sqrt(2));
    q = (qte - (rho * sqrt(2/Parameters::pi) * std::exp(qt))) / (4.0e0 * Parameters::pi);
  }
  else if (Parameters::icutoff >= 3)
  {
    // qt = -rho2 / (2.0e0);
    // qte = (1 - qt) * std::exp(qt);
    // q = (1.0e0 - qte) / (2.0e0 * Parameters::pi);
  }
}
