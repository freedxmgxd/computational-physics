#include <cstdint>
#include <iostream>
#include <list>

bool is_prime(uint64_t n)
{
    if (n == 1)
    {
        return false;
    }

    for (auto i = 2; i < n; i++)
    {
        if (n % i == 0)
        {
            return false;
        }
    }

    return true;
}

int main()
{

    // Exercicio 1.3
    // Escreva um programa que calcule e escreva os 300 primeiros números primos.
    std::list <uint64_t> prime_numbers;
    uint64_t i = 1;

    while (prime_numbers.size() < 300){
        if (is_prime(i)){
            prime_numbers.push_back(i);
        }
        i++;

    }

    for (auto v : prime_numbers)
        std::cout << v << "\n";

        return 0;
}