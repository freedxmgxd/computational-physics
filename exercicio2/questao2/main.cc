#include <iostream>
#include <matplot/matplot.h>
#include <vector>

int main() {


  std::vector<double> xs;
  std::vector<double> ys;

  for (auto r = 1.0; r < 4.0; r += 0.01) {
    auto x = 1.0 / 2.0;
    for (int i = 0; i < 100; i++) {
      x = r * x * (1 - x);
    }
    for (int i = 0; i < 100; i++) {
      x = r * x * (1 - x);
      xs.push_back(r);
      ys.push_back(x);
    }
  }

  matplot::scatter(xs, ys, 2);
  matplot::show();

  return 0;
}
