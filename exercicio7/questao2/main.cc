#include <boost/math/quadrature/gauss.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace std;

// Função original para o integrando da função gama
double gamma_integrand(double x, double a) { return pow(x, a - 1) * exp(-x); }

// Alteração de variável para melhorar a precisão
double transformed_integrand(double z, double a, double c) {
  double x = c * z / (1.0 - z);       // Alteração de variável
  double dx_dz = c / pow(1.0 - z, 2); // Derivada de x em relação a z
  return pow(x, a - 1) * exp(-x) * dx_dz;
}

// Função para calcular a função gama usando quadratura gaussiana
double gamma_function(double a, int points = 100) {
  double c = a - 1.0; // Valor de c escolhido para o pico do integrando
  auto integrand = [a, c](double z) { return transformed_integrand(z, a, c); };
  return boost::math::quadrature::gauss<double, 100>::integrate(integrand, 0.0,
                                                                1.0);
}



std::complex<double> intensity(function<double(double)> q, double u, double x,
                               double lambda, double f) {
  double phase = 2.0 * M_PI * x * u / (lambda * f);
  std::complex<double> expTerm = std::polar(1.0, phase); // e^(i*phase)
  return std::sqrt(q(u)) * expTerm;
}

int main() {
  // b) Usando C++ ao inves de Python
  auto alpha = M_PI / (20.0e-6);
  auto q = [alpha](double u) { return sin(alpha*u)*sin(alpha*u); };

  auto lambda = 500.0e-9;
  auto f = 1.0;
  auto w = 0.1;

  return 0;
}
