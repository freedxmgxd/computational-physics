#include <cstdint>
#include <iostream>

int main()
{

    // Exercicio 1.4
    //Leia do Teclado N e calcule o valor de S:
    auto N = 0.0;
    auto S = 0.0;

    std::cin >> N;

    for (uint8_t i = 1; i <= N; i++)
    {
        S += i/(N-i+1);
    }

    std::cout << S << std::endl;

        return 0;
}