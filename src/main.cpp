/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>
#include <vector>

#include "LUDoolittle.hpp"

int main() {

  // Some example code
/*
  anpi::Matrix<float> A = { {-1,-2,1,2},
                            { 2, 0,1,2},
                            {-1,-1,0,1},
                            { 1, 1,1,1} };
  anpi::Matrix<float> LU;
  
  std::vector<size_t> p;
  anpi::luDoolittle(A,LU,p);
*/
  anpi::Matrix<double> a = { {1, 0, 2},
                             {2, 1, 0},
                             {0, 2, 1} };

  anpi::Matrix<double> b = { {1},
                             {0.5},
                             {2} };

  anpi::Matrix<double> c = a * b;

  std::cout << c[0][0] << std::endl;
  std::cout << c[1][0] << std::endl;
  std::cout << c[2][0] << std::endl;

  return EXIT_FAILURE;
}
  
