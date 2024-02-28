#include "stdafx.h"
#include <math.h>
#include "Cell.h"

double Cell::m = 0.2;
double Cell::mw = 0.003;
double Cell::mo = 0.0003;
double Cell::k0 = 3.0e-11;

double Cell::kw(const double S) {
	return pow(S, 2);
}

double Cell::ko(const double S) {
	return pow(1.0 - S, 2);
}

double Cell::Bw(const double P) {
	return 1.06 - 1e-9*P;
}

double Cell::Bo(const double P) {
	return 1.0 - 6.5e-10*P;
}

Cell::Cell() {
}

Cell::~Cell() {
}

