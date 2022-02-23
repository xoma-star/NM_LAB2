#pragma once
#include "f.h"
#include "constants.h"

double SimpsonIntegral(int n) {
	double h = (b - a) / n;
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		double x_i = a + i * h;
		sum += (f(x_i) + 4 * f((x_i + x_i + h) / 2) + f(x_i + h));
	}
	sum *= h / 6;
	return sum;
}