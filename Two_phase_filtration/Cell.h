#pragma once
#include "stdafx.h"

// Contains the parameters of the substance in the cell
class Cell
{
public:

	double x = 0.;	// Cell coordinate
	double P = 0.;	// Pressure
	double S = 0.;	// Saturation

	static double m;	// Porousness
	static double mw;	// Water viscosity
	static double mo;	// Oil viscosity
	static double k0;	// Absolute permeability

	// Dependence of relative water permeability on saturation
	static double kw(const double S);
	// Dependence of relative oil permeability on saturation
	static double ko(const double S);
	// Dependence of the volumetric coefficient of water on pressure
	static double Bw(const double P);
	// Dependence of the volumetric coefficient of oil on pressure
	static double Bo(const double P);

	Cell();
	~Cell();

private:

};