#pragma once
#include <math.h>
#include<vector>
#include "constants.h"
#include "f.h"

class Gauss
{
	int n;
	std::vector<double> roots;
	std::vector<double> weights;
	double Legendre(int n, const double& x) {
		if (n == 0)			return 1;
		else if (n == 1)	return x;
		return (2. * n - 1) / n * x * Legendre(n - 1, x) - (n - 1.) / n * Legendre(n - 2, x);
	}
	double LegendreDiff(int n, const double& x) {return n / (1. - pow(x, 2)) * (Legendre(n - 1, x) - x * Legendre(n, x));}
	double ReplacementForRoots(const double& x) {return (a + b) / 2 + (b - a) * x / 2;}

	void CalculateRootsAndWeights(int n) {
		roots.clear();
		weights.clear();
		for (int i = 1; i <= n; ++i) {
			double x = 0;
			double x_next = cos(3.14 * (4. * i - 1) / (4. * n + 2));
			do {
				x = x_next;
				x_next = x - Legendre(n, x) / LegendreDiff(n, x);
			} while (abs(x_next - x) > 1. / 1000000);
			roots.push_back(x);
			double diff = LegendreDiff(n, x);
			weights.push_back(2. / ((1 - x * x) * diff * diff));
		}
	}

public:
	Gauss(int n) {
		this->n = n;
	}
	void setN(int n) {
		this->n = n;
	}
	double getVal(int n) {
		CalculateRootsAndWeights(n);
		double sum = 0;
		for (int i = 0; i < n; ++i) {
			sum += weights[i] * f(ReplacementForRoots(roots[i]));
		}
		sum *= (b - a) / 2;
		return sum;
	}
};

