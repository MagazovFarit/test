#pragma once
#include "stdafx.h"
#include "Gauss_solver.h"
#include <iostream>

int gauss::Gauss_solver(std::vector<std::vector<double>> &J, std::vector<double> &F,
	                                    //output
	                                    std::vector<double> &X) {
	if (X.empty())
		X.assign(F.size(), 0.);
	int max_ind;
	double max, tmp;
	for (int i = 0; i < F.size(); i++) {

		max = abs(J[i][i]);
		max_ind = i;
		for (int j = i; j < F.size(); j++)
		{
			if (max <= abs(J[j][i]))
			{
				max = abs(J[j][i]);
				max_ind = j;
			}
		}
		if (i != max_ind)
		{
			for (int j = i; j < F.size(); j++)
			{
				tmp = J[i][j];
				J[i][j] = J[max_ind][j];
				J[max_ind][j] = tmp;
			}
			tmp = F[i];
			F[i] = F[max_ind];
			F[max_ind] = tmp;
		}

		for (int j = i + 1; j < F.size(); j++) {
			tmp = J[j][i] / J[i][i];
			for (int k = 0; k < F.size(); k++) {
				J[j][k] -= tmp*J[i][k];
			}
			F[j] -= tmp*F[i];
		}
	}
	X[F.size() - 1] = F[F.size() - 1] / J[F.size() - 1][F.size() - 1];
	for (int i = F.size() - 2; i >= 0; i--)
	{
		X[i] = F[i];
		for (int j = i + 1; j < F.size(); j++) {
			X[i] -= J[i][j] * X[j];
		}
		X[i] /= J[i][i];
	}
	return 0;
}