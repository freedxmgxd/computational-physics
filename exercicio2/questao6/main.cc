#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <matplot/matplot.h>

int main() {
  auto wavelength = 5.0;
  auto k = 2.0 * M_PI / wavelength;
  auto xi0 = 1.0;
  auto separation = 20.0;
  auto side = 100.0;
  auto points = 500;
  auto spacing = side / points;

  auto x1 = side / 2 + separation / 2;
  auto y1 = side / 2;
  auto x2 = side / 2 - separation / 2;
  auto y2 = side / 2;

  std::vector<std::vector<double>> xs(points);

  auto x = 0.0;
  auto y = 0.0;
  for (int i = 0; i < points; i++) {
    y = i * spacing;
    for (int j = 0; j < points; j++) {
      x = j * spacing;

      auto r1 = sqrt(pow(x - x1, 2) + pow(y - y1, 2));
      auto r2 = sqrt(pow(x - x2, 2) + pow(y - y2, 2));

      xs[i].push_back(xi0 * (sin(k * r1) + sin(k * r2)));
    }
  }

  matplot::image(xs);
  matplot::colormap(matplot::palette::gray());
  matplot::show();


  return 0;
}