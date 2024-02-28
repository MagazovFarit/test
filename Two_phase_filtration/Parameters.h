#pragma once
#include "stdafx.h"
#include <vector>
#include "Cell.h"

// Problem parameters
struct Parameters
{
	int N = 50;		// Cells count
	double L = 850.0;	// Area lenght
	double h = L / N;	// Space step

	double t = 0.0;		// Current time
	double dt = 1000.0;	// Time step
	double tstop = 3.0 * 86400;	// Stop time

	double P = 12.5e6;	// Initial pressure
	double S = 0.0;		// Initial saturation

	// Returns pressure value on the left bound
	double P_left(const std::vector<Cell> C);
	// Returns saturation value on the left bound
	double S_left(const std::vector<Cell> C);
	// Returns pressure value on the right bound
	double P_right(const std::vector<Cell> C);
	// Returns saturation value on the right bound
	double S_right(const std::vector<Cell> C);
};