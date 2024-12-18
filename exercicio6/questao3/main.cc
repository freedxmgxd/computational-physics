#include <boost/math/quadrature/gauss.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>
#include <Eigen/Dense>

using namespace std;

// Regra trapezoidal
double trapezoidal(function<double(double)> f, double a, double b, int N) {
  double h = (b - a) / N;
  auto f_a = (f(a));
  auto f_b = (f(b));
  double sum = 0.5 * (f_a + f_b);
  for (int i = 1; i < N; ++i) {
    sum += f(a + i * h);
  }

  //   double e = 1 / 12 * powf(h, 2) * (derivative(f, a, h) - derivative(f, b,
  //   h));

  return sum * h;
}

// Regra de Simpson
double simpson(function<double(double)> f, double a, double b, int N) {
  if (N % 2 != 0)
    N++; // Garantir que N seja par
  double h = (b - a) / N;
  double sum = f(a) + f(b);
  for (int i = 1; i < N / 2; i++) {
    sum += 4 * f(a + (2 * i - 1) * h);
  }
  for (int i = 2; i < N / 2; i++) {
    sum += 2 * f(a + 2 * i * h);
  }

  //   double e =
  //       1 / 180 * powf(h, 4) * (derivative_3(f, a, h) - derivative_3(f, b,
  //       h));

  return sum * h / 3.0;
}

// Método adaptativo de Simpson
double simpson_adaptativo(function<double(double)> f, double a, double b, int N,
                          double epsilon) {
  if (N % 2 != 0)
    N++; // Garantir que N seja par
  double h = (b - a) / N;
  double sum = f(a) + f(b);
  for (int i = 1; i < N / 2; i++) {
    sum += 4 * f(a + (2 * i - 1) * h);
  }
  for (int i = 2; i < N / 2; i++) {
    sum += 2 * f(a + 2 * i * h);
  }
  auto Ii = sum * h / 3.0;
  auto error = 0.0;

  do {
    Ii = sum * h / 3.0;
    N *= 2;
    h = (b - a) / N;
    sum = f(a) + f(b);
    for (int i = 1; i < N / 2; i++) {
      sum += 4 * f(a + (2 * i - 1) * h);
    }
    for (int i = 2; i < N / 2; i++) {
      sum += 2 * f(a + 2 * i * h);
    }
    error = abs(Ii - sum * h / 3.0) / 15.0;
    cout << "Ii: " << Ii << " N: " << N << " e:  " << error << endl;

  } while (error > epsilon);

  //   double e =
  //       1 / 180 * powf(h, 4) * (derivative_3(f, a, h) - derivative_3(f, b,
  //       h));

  return Ii;
}

// Método de Romberg
double romberg(function<double(double)> f, double a, double b, double epsilon) {
  vector<vector<double>> R(1, vector<double>(1));
  int k = 0;
  R[k][0] = trapezoidal(f, a, b, 1);
  while (true) {
    k++;
    R.emplace_back(k + 1);
    R[k][0] = trapezoidal(f, a, b, 1 << k); // 2^k subintervalos
    // cout << "R[" << k << "][0]: " << R[k][0];
    for (int j = 1; j <= k; j++) {
      R[k][j] = (pow(4, j) * R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) - 1);
      // cout << " R[" << k << "][" << j << "]: " << R[k][j];
      // cout << "e = " << fabs((R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) -
      // 1));
      if (k > 1 &&
          fabs((R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) - 1)) < epsilon) {
        // cout << endl;
        return R[k][j];
      }
    }
    // cout << endl;
  }
}
template <int N>
double gaussian_quadrature(function<double(double)> f, double a, double b) {
  // Convert std::array to Eigen::VectorXd
    auto abscissas = boost::math::quadrature::gauss<double, N>::abscissa();
    auto weights = boost::math::quadrature::gauss<double, N>::weights();

  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    double x = ((b - a) / 2.0) * abscissas[i] + ((b + a) / 2.0);
    double w = ((b - a) / 2.0) * weights[i];
    sum += w * f(x);
  }
  return sum;
}

template <int K>
void print_gaussian_quadrature(std::function<double(double)> func, double a,
                               double b) {
  std::cout << "Quadratura de Gauss: " << K << ": "
            << gaussian_quadrature<K>(func, a, b) << std::endl;
  std::cout << "Quadratura de Gauss boost: " << K << ": "
            << boost::math::quadrature::gauss<double, K>::integrate(func, a, b)
            << std::endl;
}

template <int K>
void unroll_gaussian_quadrature(std::function<double(double)> func, double a,
                                double b) {
  if constexpr (K > 0) {
    unroll_gaussian_quadrature<K - 1>(func, a, b);
    print_gaussian_quadrature<K>(func, a, b);
  }
}

int main() {
  // Definição da função como um lambda
  // integrate ((t^3/(1-t)^5) * exp(-a * (t^2/(1-t)^2)), t, 0, 1)

  auto func = [](double x) {
    return (powf(x, 3.0) / powf(1 - x, 5.0)) *
           exp(-2 * powf((x) / (1 - x), 2.0));
  };

  // auto integral_by_wolfram =
  //     0.5028512714128629086973240134980234735337341059602528410035843691643110437548531357209388898179912121;

  double a = 0;
  double b = 0.99999; // 1.0 não é um ponto válido para a função então
                      // escolhemos um valor muito próximo
  const auto N = 100;
  auto epsilon = 1e-6;

  cout << "Trapezoidal: " << trapezoidal(func, a, b, N) << endl;
  cout << "Simpson: " << simpson(func, a, b, N) << endl;
  cout << "Simpson Adaptativo: " << simpson_adaptativo(func, a, b, N, epsilon)
       << endl;
  cout << "Romberg: " << romberg(func, a, b, epsilon) << endl;

  unroll_gaussian_quadrature<100>(func, a, b);

  cout << "Quadratura de Gauss: " << gaussian_quadrature<N>(func, a, b) << endl;
  cout << "Quadratura de Gauss boost: "
       << boost::math::quadrature::gauss<double, N>::integrate(func, a, b)
       << endl;

  return 0;
}
