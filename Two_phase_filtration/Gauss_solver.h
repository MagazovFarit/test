#pragma once
#include "stdafx.h"
#include <vector>

namespace gauss {
	int Gauss_solver(std::vector<std::vector<double>> &J, std::vector<double> &F, std::vector<double> &res);
}