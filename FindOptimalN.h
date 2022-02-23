#pragma once
#include <iostream>
#include <fstream>
#include "constants.h"

int FindOptimalN(const std::string& name, double (*integral)(int)) {
	int n = 1;
	double value = integral(1);
	double min_error = 10000;
	std::ofstream error_vs_n(name + ".txt");
	while (min_error >= delta && n <= 25) {
		++n;
		value = integral(n);
		double error = abs(exact - value);
		error_vs_n << n << "\t" << error << std::endl;
		if (error < min_error) {
			min_error = error;
		}
	}
	error_vs_n.close();
	return n;
}