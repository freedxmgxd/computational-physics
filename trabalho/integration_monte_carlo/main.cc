#include <cmath>
#include <functional>
#include <iostream>
#include <stdlib.h>

using namespace std;

double integrationMonteCarlo(double (*f)(double), double a, double b, double N,
                             double h, double t) {
  double count = 0;
  double A = abs((a - b) * (h - t));

  for (int i = 0; i < N; ++i) {
    double x = a + (b - a) * rand() / RAND_MAX;
    double y = h + (t - h) * rand() / RAND_MAX;
    if (y <= f(x)) {
      count++;
    }
  }

  return A * count / N;
}

double integrationMonteCarloWeight(double (*f)(double), double a, double b,
                                   double N) {
  double count = 0;
  auto w = [](double x) { return pow(x, -0.5); };

  for (auto i = 0; i < N; ++i) {
    double z = a + (b - a) * rand() / RAND_MAX;
    auto x = z*z;
    count += f(x) / w(x);
  }

  return 2 * count / N; // 2 is the integral of w(x) from 0 to 1
}

pair<double, double> generateGaussianRandomNumber() {
  double sigma = 1;

  auto r = sqrt(-2 * sigma * sigma * log(1 - rand() / RAND_MAX));
  auto theta = 2 * M_PI * rand() / RAND_MAX;

  return make_pair(r * cos(theta), r * sin(theta));
}

int main() {

  auto quarterCircle = [](double x) { return sqrt(1 - x * x); };

  double a = 0;
  double b = 1;
  double h = 0;
  double t = 4;
  double N = 1000000;

  double result = integrationMonteCarlo(quarterCircle, a, b, N, h, t);
  cout << "Result Basic Monte Carlo: " << result << endl;

  cout << "Expected: " << M_PI / 4 << endl;
  cout << "Error: " << abs(result - M_PI / 4) / (M_PI / 4) << endl;

  auto f = [](double x) { return (pow(x, -0.5)) / (exp(x) + 1); };

  double resultWeight = integrationMonteCarloWeight(f, a, b, N);

  cout << "Result Weight Monte Carlo: " << resultWeight << endl;

  cout << "Expected: " << 0.84 << endl;
  cout << "Error: " << abs(resultWeight - 0.84) / (M_PI / 4) << endl;

  return 0;
}