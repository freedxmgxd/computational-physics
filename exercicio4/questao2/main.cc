#include <cmath>
#include <iostream>

unsigned long long int factorial_int(unsigned long long int x) {
  auto aux = 1;

  for (unsigned long long int i = 1; i <= x; i++) {
    aux *= i;
  }
  return aux;
}

long double factorial_double(long double x) {
  long double aux = 1.0;

  for (int i = 1; i <= x; i++) {
    aux *= i;
  }
  return aux;
}

int main() {
  double x = 1.0;
  double y = 1.0 + (1e-14) * std::sqrt(2);

  // Calcula 10^14 * (y - x)
  double result = (1e14) * (y - x);

  // Imprime o resultado e o valor esperado de sqrt(2)
  std::cout << "10^14 * (y - x) = " << result << std::endl;
  std::cout << "sqrt(2) = " << std::sqrt(2) << std::endl;

  // O erro é devido a perda de precisão ao subtrair dois números muito próximos

  // CHAPTER 4 | ACCURACY AND SPEED
  // many digits. Try, for example, doing print(2**1000000) in Python. The cal-
  // culation can be done—it yields a number with 301 030 digits—but it’s so
  // slow that you might as well forget about using your computer for anything
  // else for the next few minutes.3 Exercise 4.1: Write a program to calculate
  // and print the factorial of a number entered by the user. If you wish you
  // can base your program on the user-defined function for facto- rial given in
  // Section 2.6, but write your program so that it calculates the factorial
  // using integer variables, not floating-point ones. Use your program to
  // calculate the factorial of 200. Now modify your program to use
  // floating-point variables instead and again calcu- late the factorial of
  // 200. What do you find? Explain

  int x_int = 10;
  double x_double = 10.0;

  std::cout << "factorial_int(" << x_int << ") = " << factorial_int(x_int)
            << std::endl;
  std::cout << "factorial_double(" << x_double
            << ") = " << factorial_double(x_double) << std::endl;

  x_int = 20;
  x_double = 20.0;

  std::cout << "factorial_int(" << x_int << ") = " << factorial_int(x_int)
            << std::endl;
  std::cout << "factorial_double(" << x_double
            << ") = " << factorial_double(x_double) << std::endl;

  x_int = 200;
  x_double = 200.0;

  std::cout << "factorial_int(" << x_int << ") = " << factorial_int(x_int)
            << std::endl;
  std::cout << "factorial_double(" << x_double
            << ") = " << factorial_double(x_double) << std::endl;

  // Vemos que conforme o valor de x aumenta, o valor de x! calculado começa a
  // divergir devido a perda de precisão ao utilizar variáveis do tipo inteiro
  // para somas muito grandes, e pro caso de valores muito grandes começamos a
  // ter numeros maiores do que é possivel representar com os tipos padrões da
  // linguagem de programação como é o caso para x = 200 tambem pro tipo inteiro

  return 0;
}
