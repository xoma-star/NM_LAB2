#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include "constants.h"
#include "rectangle.h"
#include "FindOptimalN.h"
#include "trapezoid.h"
#include "simpson.h"
#include "det.h"
#include "romberg.h"
#include "Gauss.h"
#include "Monte.h"

std::string file_path = __FILE__;
std::string dir_path = file_path.substr(0, file_path.rfind("\\")) + "\\";

void DrawGraphic(const std::string& data) {
	std::ofstream file("file");
	file << "plot " << data << "; pause mouse keypress" << "\n";
	file.close();
	std::system("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"  -persist file");
}

double GaussWrapper(int n) {
	Gauss x = Gauss(n);
	return x.getVal(n);
}

int main() {
	dir_path.replace(dir_path.begin(), dir_path.end(), "\\", "\\\\");
	setlocale(LC_ALL, "Ru");
	std::cout << "Задание 1\n";
	double rectInt = RectangleIntegral(1);
	int rectangle_n = FindOptimalN("Rectangle", &RectangleIntegral);
	rectInt = RectangleIntegral(rectangle_n);
	std::cout << "Прямоугольный: " << std::endl;
	std::cout << "оптимальный n: " << rectangle_n << std::endl;
	std::cout << "значение: " << std::setprecision(10) << rectInt << std::endl;

	double trapInt = TrapezoidIntegral(1);
	int trapezoid_n = FindOptimalN("Trapezoid", &TrapezoidIntegral);
	trapInt = TrapezoidIntegral(trapezoid_n);
	std::cout << "Трапезоидный: " << std::endl;
	std::cout << "оптимальный n: " << trapezoid_n << std::endl;
	std::cout << "значение: " << std::setprecision(10) << trapInt << std::endl;

	DrawGraphic("\"Rectangle.txt\" using 1:2 title \"n vs Rectangle error\" with lines, \"Trapezoid.txt\" using 1:2 title \"n vs Trapezoid error\" with lines");
	

	std::cout << "Задание 2\n";
	double simpInt = SimpsonIntegral(1);
	int simpson_n = FindOptimalN("Simpson", &SimpsonIntegral);
	simpInt = SimpsonIntegral(simpson_n);
	std::cout << "Прямоугольный: " << std::endl;
	std::cout << "оптимальный n: " << simpson_n << std::endl;
	std::cout << "значение: " << std::setprecision(10) << simpInt << std::endl;

	DrawGraphic("\"Simpson.txt\" using 1:2 title \"n vs Simpson error\" with lines");


	std::cout << "Задание 3\n";
	double h1 = (b - a) / trapezoid_n;
	double h2 = h1 / 2.;
	double h3 = h2 / 2.;
	double trapInt2 = TrapezoidIntegral(trapezoid_n / 2);
	double trapInt3 = TrapezoidIntegral(trapezoid_n / 4);
	std::cout << "Шаг сетки: " << h1 << ", значение: " << trapInt;
	std::cout << "\nШаг сетки: " << h2 << ", значение: " << trapInt2;
	std::cout << "\nШаг сетки: " << h3 << ", значение: " << trapInt3;

	double F = trapInt + pow((trapInt - trapInt2), 2) / (2. * trapInt2 - trapInt - trapInt3);
	std::cout << "\nЭйткин: " << std::setprecision(10) << F << std::endl;
	std::cout << "Погрешность: " << abs(F - exact) << std::endl;

	double p = pow(log(0.5), -1) * log((trapInt3 - trapInt2) / (trapInt2 - trapInt));
	std::cout << "  Порядок точности: " << p << std::endl;
	double sigma = pow(h2, p) / (pow(h2, p) - pow(h1, p));
	double J = sigma * trapInt + (1. - sigma) * trapInt2;
	//double q = (F2 - F1) / (pow(h1, p) - pow(h2, p)) * pow(h1, 2 * p);
	std::cout << "Рунге: " << std::setprecision(10) << J << std::endl;
	std::cout << "Погрешность: " << abs(J - exact) << std::endl;

	//Romberg rule
	std::cout << "Ромберг: " << std::endl;
	double err = 1;
	double new_delta = pow(10, -6);
	double count = 2;
	std::vector<std::vector<double>> up, down;
	do {
		up.resize(count);
		down.resize(count);
		for (int i = 0; i < count; ++i) {
			up[i].resize(count);
			down[i].resize(count);
			up[i][0] = RectangleIntegral(pow(2, i));
			down[i][0] = 1;
		}
		for (int i = 0; i < count; ++i) {
			for (int j = 1; j < count; ++j) {
				double h = (b - a) / (pow(2, i));
				up[i][j] = pow(h, p + j - 1);
				down[i][j] = pow(h, p + j - 1);
			}
		}
		double det = Det(up) / Det(down);
		err = abs(exact - det);
		std::cout << "Сетка: " << count << " Погрешность: " << err << std::endl;
		++count;
	} while (err > new_delta);;

	//Задание 4
	Gauss gauss_integral(1);
	int gauss_n = FindOptimalN("Gauss", &GaussWrapper);
	gauss_integral.setN(gauss_n);
	std::cout << "Гаусс: " << std::endl;
	std::cout << "n: " << gauss_n << std::endl;
	std::cout << "значение: " << std::setprecision(10) << gauss_integral.getVal(gauss_n) << std::endl;

	DrawGraphic("\"Gauss.txt\" using 1:2 title \"n vs Gauss error\" with lines");

	//Задание 6
	std::srand(std::time(nullptr));
	double monte;
	int monte_n = FindOptimalN("MonteCarlo", &Monte);
	monte = Monte(monte_n);
	std::cout << "Монте Карло: " << std::endl;
	std::ofstream error_vs_n("MonteCarlo.txt");
	int N = 1000;
	for (int n = 1; n <= N; ++n) {
		error_vs_n << n << " " << MonteErr(n) << std::endl;
	}
	error_vs_n.close();
	std::cout << "n: " << N << std::endl;
	std::cout << "знач.: " << monte << std::endl;
	DrawGraphic("\"MonteCarlo.txt\" using 1:2 title \"n vs MonteCarlo expected abs error\" with lines");
	return 0;
}