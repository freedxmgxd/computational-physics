#include <cstdint>
#include <iostream>
#include <string>
#include <algorithm>

// Code from https://www.geeksforgeeks.org/how-to-handle-large-numbers-in-cpp/
// Class to handle large numbers
class LargeNumber
{
private:
    // The large number represented as a string
    std::string number;

public:
    LargeNumber() : number("0") {}
    LargeNumber(const std::string &num) : number(num) {}

    // Overloaded operator+ to add two LargeNumber objects
    LargeNumber operator+(const LargeNumber &other) const
    {
        std::string result;
        int carry = 0;
        int maxLength = std::max(number.length(), other.number.length());

        for (int i = 0; i < maxLength || carry; ++i)
        {
            int digit1 = i < number.length()
                             ? number[number.length() - 1 - i] - '0'
                             : 0;
            int digit2 = i < other.number.length()
                             ? other.number[other.number.length() - 1 - i] - '0'
                             : 0;

            int sum = digit1 + digit2 + carry;
            result.push_back(sum % 10 + '0');
            carry = sum / 10;
        }

        // Since the result is reversed, reverse it back to
        // get the correct number
        std::reverse(result.begin(), result.end());
        return LargeNumber(result);
    }

    // Overloaded operator<< to print a LargeNumber object
    friend std::ostream &operator<<(std::ostream &out,
                                    const LargeNumber &num)
    {
        out << num.number;
        return out;
    }
};

int main()
{

    // Exercicio 1.2
    u_int64_t j = 1;
    u_int64_t j_old = 0;
    u_int64_t temp = 0;

    LargeNumber large_j("1");
    LargeNumber large_j_old("0");
    LargeNumber large_temp("0");

    for (auto i = 0; i < 200; i++)
    {

        if (i < 92)
        {
            std::cout << i + 1 << " " << j << std::endl;

            temp = j;
            j = temp + j_old;
            j_old = temp;
        }
        else if (i == 92)
        {
            large_j = std::to_string(j);
            large_j_old = std::to_string(j_old);
            std::cout << i + 1 << " " << large_j << std::endl;

            large_temp = large_j;
            large_j = large_temp + large_j_old;
            large_j_old = large_temp;
        }
        else
        {

            std::cout << i + 1 << " " << large_j << std::endl;

            large_temp = large_j;
            large_j = large_temp + large_j_old;
            large_j_old = large_temp;
        }
    }

    return 0;
}