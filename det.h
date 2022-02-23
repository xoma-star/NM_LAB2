#pragma once
#include <vector>

double Det(std::vector<std::vector<double>> matrix) {
	for (int i = 0; i <= matrix.size() - 1; ++i) {
		for (int j = i + 1; j <= matrix.size() - 1; ++j) {
			double coeff = -matrix[j][i] / matrix[i][i];
			for (int k = 0; k < matrix.size(); ++k) {
				matrix[j][k] += coeff * matrix[i][k];
			}
		}
	}
	double det = 1;
	for (int i = 0; i <= matrix.size() - 1; ++i) {
		det *= matrix[i][i];
	}
	return det;
}