#pragma once
#include "f.h"
#include "constants.h"

double TrapezoidIntegral(int n) {
	double sum = 0;
	double h = (b - a) / n;
	for (int i = 0; i < n; ++i) {
		double x_i = a + i * h;
		sum += (f(x_i) + f(x_i + h));
	}
	sum *= 0.5 * h;
	return sum;
}