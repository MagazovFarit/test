#include "stdafx.h"
#include "Parameters.h"

double Parameters::P_left(const std::vector<Cell> C) {
	return 22.0e6;
}

double Parameters::S_left(const std::vector<Cell> C) {
	return 1.0;
}

double Parameters::P_right(const std::vector<Cell> C) {
	return 12.5e6;
}

double Parameters::S_right(const std::vector<Cell> C) {
	return 0.0;
}
