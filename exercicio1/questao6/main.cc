#include <cstdint>
#include <iostream>
#include <cmath>

int main()
{

    // Exercicio 1.6
    // Um satélite deve ser lançado em uma orbita circular ao redor da Terra para que orbite o planeta uma vez a cada T segundos.

    // a. Apresente que a altitude h acima da superfície da terra que o satélite deve ter é:

    // $$ h = \left( \frac{GMT^2}{4\pi^2} \right)^{1/3} - R $$

    // onde G=6.67 x 10-11 m3kg-1s-2 é a constante gravitacional de Newton, M=5.97 x 10^24 kg é a massa da Terra e R=6371 km é seu raio;

    // b. Escreva um programa que pergunte ao usuário para entrar um valor determinado de T e ele calcula e imprime na tela (ou arquivo) a correta altitude em metros;

    // c. Use este programa para calcular as altitudes dos satélites que orbitam a Terra uma vez ao dia (chamados de geossíncronos), uma vez a cada 90 minutos e uma vez a cada 45 minutos. O que você conclui desses cálculos?

    // d. Tecnicamente, um satélite geossíncrono é aquele que orbita a Terra uma vez no dia persideal, que é de 23,93 horas, não 24 horas. Por que é isso? E quanta diferença fará com a altitude do satélite?
    
    auto G = 6.67e-11; // m3kg-1s-2
    auto M = 5.97e24; // kg
    auto R = 6371e3; // m
    auto T = 0.0;

    std::cin >> T;

    auto h = pow(G*M*pow(T,2)/(4*pow(3.14159265359,2)),1.0/3.0) - R;

    std::cout << h << std::endl;

    return 0;
}