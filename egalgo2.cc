/*******************************************************************************

    Copyright (C) 2018 Nicholas C. Ham

    This work is licensed under a Creative Commons Attribution-ShareAlike 4.0
    International License. See http://creativecommons.org/licenses/by-sa/4.0/

    For a discussion about this code and the latest version see
    https://gitlab.com/n-ham-paper-files/lattice-path-algorithms

    Compile with:

        g++ -O3 -std=c++11 -Wall -Wextra -pedantic -o egalgo2 egalgo2.cc

*******************************************************************************/

#include "algorithms.h"

int main()
{
    StepSet ss;

    std::ifstream ifs("./input.txt");
    ifs >> ss;
    ifs.close();

    //std::cout << ss << std::endl;

    std::string property;
    Point u;

    algorithm2(ss, 1, property, u);

    //std::cout << property << std::endl;
    /*if(property == "FPP")
        std::cout << "u = <" << u.first << ", " << u.second << ">" << std::endl;*/

    return 0;
}
