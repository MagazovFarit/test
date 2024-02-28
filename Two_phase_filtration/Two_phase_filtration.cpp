// Two_phase_filtration.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include "Parameters.h"
#include "Cell.h"
#include "Newton_solver.h"

int write_saturations_to_file(std::ofstream &outfile, std::vector<Cell>& vector)
{
	for (size_t i = 0; i < vector.size(); ++i) {
		outfile << vector[i].S << "\t";
	}
	outfile << "\n";

	return 0;
}

int main()
{
	// Problem parameters
	Parameters parameters;
	nonlinear_solver solver;

	// Initial condition
	std::vector<Cell> C(parameters.N);
	std::vector<Cell> Cn(parameters.N);
	for (size_t i = 0; i < parameters.N; i++)
	{
		C[i].x = i*parameters.h;
		C[i].P = parameters.P;
		C[i].S = parameters.S;
	}

	// Output to file
	std::ofstream fout;
	fout.open("output.txt");
	// Grid output
	for (size_t i = 0; i < parameters.N; i++) {
		fout << C[i].x << "\t";
	}
	fout << "\n";

	while (parameters.t < parameters.tstop) {

		// Newton method
		solver.Newton_solver(C, Cn, parameters);

		// Current time saturation output
		write_saturations_to_file(fout, Cn);

		// Next time step
		parameters.t += parameters.dt;
		Cn.swap(C);
	}

	fout.close();
	system("pause");
    return 0;
}

