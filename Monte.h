#pragma once
#include "constants.h"
#include "f.h"
#include <cstdlib>

double Monte(int n) {
	double h = (b - a) / n;
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		double x_i = a + ((double)std::rand() / RAND_MAX) * (b - a);
		sum += f(x_i);
	}
	sum *= h;
	return sum;
}

double MonteErr(int n){
	double sum = 0;
	for (int i = 0; i < 100; ++i) {
		sum += abs(exact - Monte(n));
	}
	sum /= 100.;
	return sum;
}