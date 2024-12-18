#include <cstdint>
#include <iostream>

int main()
{

    // Exercicio 1.2
    uint64_t j = 1;
    uint64_t j_old = 1;

    for (uint8_t i = 0; i < 200; i++)
    {
        std::cout << j << std::endl;

        uint64_t temp = j;
        j_old = j;
        j = temp + j_old;
    }

    

    return 0;
}