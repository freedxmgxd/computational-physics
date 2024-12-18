#include <cmath>
#include <iostream>

auto factorial(int x) {
  auto aux = 1;

  for (auto i = 1; i <= x; i++) {
    aux *= i;
  }
  return aux;
}

auto fn_sin(double x, bool debug = false) {

  auto aux = 0.0;

  auto inc = 0.0;
  auto i = 1;

  do {
    inc = (pow(-1, i - 1) * pow(x, (2 * i) - 1)) / factorial((2 * i) - 1);

    aux += inc;

    if (debug) {
      std::cout << "N: " << i << " inc: " << inc << " aux: " << aux
                << std::endl;
    }

    i++;

  } while (fabs(inc) > 0.0000001 * fabs(aux)); // Utilizar o fabs pois o abs é
                                               // valido apenas para inteiros

  return aux;
}

int main() {

  // b) calcule sua série para x≤1 e compare-a com alguma biblioteca que tenha o
  // sen(x) implementado. Pare sua soma no valor de N para o qual o próximo
  // termo da série não seja maior que 10e−7 do somatório até aquele ponto

  std::cout << std::endl << "a) e b)" << std::endl;
  auto x = M_PI / 4;
  auto sen = fn_sin(x);
  auto sen_lib = sin(x);

  std::cout << "sen(" << x << ") = " << sen << std::endl;
  std::cout << "sen(" << x << ") = " << sen_lib << std::endl;

  // c) Examine os termos da série para x≈3π e observe os cancelamentos
  // subtrativos significativos que ocorrem quando termos grandes se somam para
  // dar respostas pequenas.

  std::cout << std::endl << "c)" << std::endl;

  x = 3 * M_PI;

  sen = fn_sin(x, true);
  sen_lib = sin(x);

  return 0;
}