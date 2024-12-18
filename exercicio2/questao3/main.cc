#include <cstdint>
#include <iostream>
#include <matplot/matplot.h>
#include <vector>

auto least_squares_fitting(const std::vector<double> &xs,
                           const std::vector<double> &ys) {
  auto N = xs.size();

  auto sum_x = 0.0;
  for (auto x : xs) {
    sum_x += x;
  }
  auto E_x = sum_x / N;

  auto sum_y = 0.0;
  for (auto y : ys) {
    sum_y += y;
  }
  auto E_y = sum_y / N;

  auto sum_x2 = 0.0;
  for (auto x : xs) {
    sum_x2 += x * x;
  }
  auto E_xx = sum_x2 / N;

  auto sum_xy = 0.0;
  for (std::size_t i = 0; i < N; i++) {
    sum_xy += xs[i] * ys[i];
  }
  auto E_xy = sum_xy / N;

  auto m = (E_xy - E_x * E_y) / (E_xx - E_x * E_x);
  auto c = (E_xx * E_y - E_x * E_xy) / (E_xx - E_x * E_x);

  return std::pair<double, double>{m, c};
}

int main() {

  // // millikan.txt
  // 5.4874e+14 0.5309
  // 6.931e+14 1.0842
  // 7.4307e+14 1.2734
  // 8.2193e+14 1.6598
  // 9.6074e+14 2.19856
  // 1.184e+15 3.10891

  auto xs = std::vector<double>{5.4874e+14, 6.931e+14,  7.4307e+14,
                                8.2193e+14, 9.6074e+14, 1.184e+15};
  auto ys =
      std::vector<double>{0.5309, 1.0842, 1.2734, 1.6598, 2.19856, 3.10891};

  auto [m, c] = least_squares_fitting(xs, ys);


  auto x_ss = matplot::linspace(5.0e+14, 1.2e+15, 100);
  auto y_ss = matplot::transform(x_ss, [m, c](auto x) { return m * x + c; });

  matplot::plot(xs, ys, "-:o");
  matplot::hold(matplot::on);
  matplot::plot(x_ss, y_ss);

  matplot::show();

  const double e = 1.602e-19;

  double h = m * e;
  
  std::cout << "m = " << m << std::endl;
  std::cout << "c = " << c << std::endl;
  std::cout << "h = " << h << std::endl;
  return 0;
}