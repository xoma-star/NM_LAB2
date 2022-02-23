#pragma once
#include "f.h"
#include "constants.h"

double RectangleIntegral(int n){
	double h = (b - a) / n;
	double sum = 0;
	//Central rectangles
	for (int i = 0; i < n; ++i) {
		double x_i = a + i * h;
		sum += f((x_i + x_i + h) / 2);
	}
	sum *= h;
	return sum;
}