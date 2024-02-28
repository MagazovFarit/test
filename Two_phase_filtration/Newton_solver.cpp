#include "stdafx.h"
#include <iostream>
#include <math.h>
#include "Newton_solver.h"
#include "Gauss_solver.h"

int nonlinear_solver::Newton_solver(std::vector<Cell> &C, std::vector<Cell> &Cn, Parameters &parameters)
{
	double eps = 1e-4;

	// First approximation
	std::copy(C.begin(), C.end(), Cn.begin());

	// J*dX = -F
	// Jacoby matrix J = [dF[i]/dX[j]]
	if (J.empty())
		J.resize(2 * parameters.N, std::vector<double>(2 * parameters.N));
	
	// Right vector F
	if (F.empty())
		F.resize(2 * parameters.N);

	// Increment [dS[0], dP[0], dS[1], dP[1], ...]
	if (dX.empty())
		dX.resize(2 * parameters.N);

	// Iterations of the newton method
	do 
	{
		Jacoby(Cn, C, parameters, J);
		Func(Cn, C, parameters, F);
		gauss::Gauss_solver(J, F, dX);
		for (size_t i = 0; i < parameters.N; i++)
		{
			Cn[i].P -= dX[2 * i + 1];
			Cn[i].S -= dX[2 * i];
		}
	} while (max(dX) > eps);
	return 0;
}

double nonlinear_solver::max(const std::vector<double> &dX) {
	double max = abs(dX[0]);
	for (size_t i = 0; i < dX.size(); i++) {
		if (max < abs(dX[i])) {
			max = abs(dX[i]);
		}
	}
	return max;
}

void nonlinear_solver::Jacoby(const std::vector<Cell> &Cn, const std::vector<Cell> &C, Parameters &parameters,
										//output
										std::vector<std::vector<double>> &J)
{
	// Derivative: (Fr - Fl) / (2*h)
	std::vector<double> Fr(2 * parameters.N), Fl(2 * parameters.N);
	double h = 1e-4;
	std::vector<Cell> Ch(Cn);

	for (size_t j = 0; j < parameters.N; j++) {
		// Derivative dF / dP
		Ch[j].P += h;
		Func(Ch, C, parameters, Fr);
		Ch[j].P -= 2 * h;
		Func(Ch, C, parameters, Fl);
		Ch[j].P += h;
		for (size_t i = 0; i < 2 * parameters.N; i++) {
			J[i][2 * j + 1] = (Fr[i] - Fl[i]) / (2 * h);
		}

		// Derivative dF / dS
		Ch[j].S += h;
		Func(Ch, C, parameters, Fr);
		Ch[j].S -= 2 * h;
		Func(Ch, C, parameters, Fl);
		Ch[j].S += h;
		for (size_t i = 0; i < 2 * parameters.N; i++)
		{
			J[i][2 * j] = (Fr[i] - Fl[i]) / (2 * h);
		}
	}
}

void nonlinear_solver::Func(const std::vector<Cell> &Cn, const std::vector<Cell> &C, Parameters &p,
						   //output
						   std::vector<double> &F) {
	for (size_t i = 1; i < p.N - 1; i++)
	{
		F[2 * i] = Cell::m*Cn[i].S / Cell::Bw(Cn[i].P) - Cell::m*C[i].S / Cell::Bw(C[i].P) - p.dt / pow(p.h, 2) *Cell::k0 / Cell::mw *(
			Cell::kw(Cn[i].S) / Cell::Bw(0.5*(Cn[i + 1].P + Cn[i].P))*(Cn[i + 1].P - Cn[i].P) -
			Cell::kw(Cn[i - 1].S) / Cell::Bw(0.5*(Cn[i - 1].P + Cn[i].P))*(Cn[i].P - Cn[i - 1].P));

		F[2 * i + 1] = Cell::m*(1.0 - Cn[i].S) / Cell::Bo(Cn[i].P) - Cell::m*(1.0 - C[i].S) / Cell::Bo(C[i].P) - p.dt / pow(p.h, 2) *Cell::k0 / Cell::mo *(
			Cell::ko(Cn[i].S) / Cell::Bo(0.5*(Cn[i + 1].P + Cn[i].P))*(Cn[i + 1].P - Cn[i].P) -
			Cell::ko(Cn[i - 1].S) / Cell::Bo(0.5*(Cn[i - 1].P + Cn[i].P))*(Cn[i].P - Cn[i - 1].P));
	}

	// Left bound
	size_t i = 0;
	F[2 * i] = Cell::m*Cn[i].S / Cell::Bw(Cn[i].P) - Cell::m*C[i].S / Cell::Bw(C[i].P) - p.dt / pow(p.h, 2) *Cell::k0 / Cell::mw *(
		Cell::kw(Cn[i].S) / Cell::Bw(0.5*(Cn[i + 1].P + Cn[i].P))*(Cn[i + 1].P - Cn[i].P) -
		Cell::kw(p.S_left(Cn)) / Cell::Bw(0.5*(p.P_left(Cn) + Cn[i].P))*(Cn[i].P - p.P_left(Cn)));

	F[2 * i + 1] = Cell::m*(1.0 - Cn[i].S) / Cell::Bo(Cn[i].P) - Cell::m*(1.0 - C[i].S) / Cell::Bo(C[i].P) - p.dt / pow(p.h, 2) *Cell::k0 / Cell::mo *(
		Cell::ko(Cn[i].S) / Cell::Bo(0.5*(Cn[i + 1].P + Cn[i].P))*(Cn[i + 1].P - Cn[i].P) -
		Cell::ko(p.S_left(Cn)) / Cell::Bo(0.5*(p.P_left(Cn) + Cn[i].P))*(Cn[i].P - p.P_left(Cn)));

	// Right bound
	i = p.N - 1;
	F[2 * i] = Cell::m*Cn[i].S / Cell::Bw(Cn[i].P) - Cell::m*C[i].S / Cell::Bw(C[i].P) - p.dt / pow(p.h, 2) *Cell::k0 / Cell::mw *(
		Cell::kw(Cn[i].S) / Cell::Bw(0.5*(p.P_right(Cn) + Cn[i].P))*(p.P_right(Cn) - Cn[i].P) -
		Cell::kw(Cn[i - 1].S) / Cell::Bw(0.5*(Cn[i - 1].P + Cn[i].P))*(Cn[i].P - Cn[i - 1].P));

	F[2 * i + 1] = Cell::m*(1.0 - Cn[i].S) / Cell::Bo(Cn[i].P) - Cell::m*(1.0 - C[i].S) / Cell::Bo(C[i].P) - p.dt / pow(p.h, 2) *Cell::k0 / Cell::mo *(
		Cell::ko(Cn[i].S) / Cell::Bo(0.5*(p.P_right(Cn) + Cn[i].P))*(p.P_right(Cn) - Cn[i].P) -
		Cell::ko(Cn[i - 1].S) / Cell::Bo(0.5*(Cn[i - 1].P + Cn[i].P))*(Cn[i].P - Cn[i - 1].P));
}