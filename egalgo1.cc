/*******************************************************************************

    Copyright (C) 2018 Nicholas C. Ham

    This work is licensed under a Creative Commons Attribution-ShareAlike 4.0
    International License. See http://creativecommons.org/licenses/by-sa/4.0/

    For a discussion about this code and the latest version see
    https://gitlab.com/n-ham-paper-files/lattice-path-algorithms

    Compile with:

        g++ -O3 -std=c++11 -Wall -Wextra -pedantic -o egalgo1 egalgo1.cc

*******************************************************************************/

#include "algorithms.h"

int main(int argc, char* argv[])
{
    int q=10;
    StepSet ss;
    ConstraintVector cv;

    if(argc > 1)
        q = atoi(argv[1]);

    std::ifstream ifs("./input.txt");
    ifs >> ss;
    ifs >> cv;
    ifs.close();

    //std::cout << ss << std::endl;
    //std::cout << cv << std::endl;

    std::vector<EndPointSet> endPointSets;
    std::vector<Graph> graphs;
    algorithm1(ss, cv, q, 1, endPointSets, graphs);

    return 0;
}
