#include <boost/math/quadrature/gauss.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

template <int N>
double gaussian_quadrature(function<double(double)> f, double a, double b) {
  std::vector<double> abscissas =
      boost::math::quadrature::gauss<double, N>::abscissa();
  std::vector<double> weights =
      boost::math::quadrature::gauss<double, N>::weights();

  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    double x = abscissas[i]* ((b - a) / 2.0) + ((b + a) / 2.0);
    sum += ((b - a) / 2.0) * weights[i] * f(x);
  }
  return sum;
}

double hermite(int n, double x) {
  if (n == 0)
    return 1;
  if (n == 1)
    return 2 * x;
  return 2 * x * hermite(n - 1, x) - 2 * (n - 1) * hermite(n - 2, x);
}

double wavefunction(int n, double x) {
  return (1 / sqrt(powf(2.0, n) * tgamma(n + 1) * sqrt(M_PI))) * hermite(n, x) *
         exp(-(x * x) / 2);
}

int main() {

  std::ofstream dataFile("results.dat");
  auto x = -4.0;
  while (x <= 4.0) {
    dataFile << x << " ";
    for (int j = 0; j < 4; j++) {
      dataFile << j << " " << wavefunction(j, x);
    }
    dataFile << "\n";
    x += 0.1;
  }
  dataFile.close();

  // Use your function to make a plot that shows the harmonic oscil-
  // lator wavefunctions for n = 0, 1, 2, and 3, all on the same graph, in the
  // range x = −4 to x = 4.
  // Create GNU Plot script

  std::ofstream gnuplotScript("plot_results.gp");
  gnuplotScript
      << "set terminal qt\n"; // Use qt terminal for interactive window
  gnuplotScript << "set title 'Harmonic Oscillator Wavefunctions'\n";
  gnuplotScript << "set xlabel 'x'\n";
  gnuplotScript << "set ylabel 'Wavefunction'\n";
  gnuplotScript << "set grid\n";
  gnuplotScript << "plot 'results.dat' using 1:2 with lines title 'n=0', \\\n";
  gnuplotScript << "     'results.dat' using 1:3 with lines title 'n=1', \\\n";
  gnuplotScript << "     'results.dat' using 1:4 with lines title 'n=2', \\\n";
  gnuplotScript << "     'results.dat' using 1:5 with lines title 'n=3'\n";
  gnuplotScript.close();

  // Execute GNU Plot script
  system("gnuplot -p plot_results.gp");

  std::ofstream dataFile_b("results_b.dat");
  x = -10.0;
  while (x <= 10.0) {
    dataFile_b << x << " " << wavefunction(30, x) << "\n";
    x += 0.1;
  }
  dataFile_b.close();

  std::ofstream gnuplotScript_b("plot_results_b.gp");
  gnuplotScript_b
      << "set terminal qt\n"; // Use qt terminal for interactive window
  gnuplotScript_b << "set title 'Harmonic Oscillator Wavefunctions'\n";
  gnuplotScript_b << "set xlabel 'x'\n";
  gnuplotScript_b << "set ylabel 'Wavefunction'\n";
  gnuplotScript_b << "set grid\n";
  gnuplotScript_b << "plot 'results_b.dat' using 1:2 with lines title 'n=30'\n";
  gnuplotScript_b.close();

  system("gnuplot -p plot_results_b.gp");

  auto func = [](double x) {
    double psi = wavefunction(5, x);
    return x * x * psi * psi;
  };

  auto a = -10.0;
  auto b = 10.0;
  const auto N = 100;

  cout << "Gaussian Quadrature: " << gaussian_quadrature<N>(func, a, b)
       << endl; // A forma que implementei não está funcionando para o caso
                // então irei usar tambem a biblioteca para o exercicio

  cout << "Gaussian Quadrature lib: "
       << boost::math::quadrature::gauss<double, N>::integrate(func, a, b)
       << endl;

  return 0;
}
