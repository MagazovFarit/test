#pragma once
#include "stdafx.h"
#include "Parameters.h"
#include "Cell.h"

class nonlinear_solver {
	//namespace nwt {
public:
		// Newton solver main function
		int Newton_solver(std::vector<Cell> &C, std::vector<Cell> &Cn, Parameters &parameters);
private:
		// Search max(|dX|)
		double max(const std::vector<double> &dX);
		// Jacoby matrix J: J*dX = -F
		void Jacoby(const std::vector<Cell> &Cn, const std::vector<Cell> &C, Parameters &parameters, std::vector<std::vector<double>> &J);
		// Vector function F: J*dX = -F
		void Func(const std::vector<Cell> &Cn, const std::vector<Cell> &C, Parameters &parameters, std::vector<double> &F);

		std::vector<std::vector<double>> J;
		std::vector<double> F;
		std::vector<double> dX;
	//}
};