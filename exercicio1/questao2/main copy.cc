#include <cstdint>
#include <iostream>
#include <string>
#include <algorithm>

int main()
{

    // Exercicio 1.2
    u_int64_t j = 1;
    u_int64_t j_old = 0;
    u_int64_t temp = 0;

    for (auto i = 0; i < 200; i++)
    {
            std::cout << i + 1 << " " << j << std::endl;

            temp = j;
            j = temp + j_old;
            j_old = temp;
       
    }

    return 0;
}