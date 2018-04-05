/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: JP
 * @Date  : 27.03.2018
 */

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <functional>
#include <cmath>

#include "Inversion.hpp"

namespace anpi {
  namespace test {

    template<typename T>
    void inversionTest() {

      anpi::Matrix<T> A = { { 3, 4 },
                            { 1, 2 }};

      anpi::Matrix<T> ai;

      anpi::Matrix<T> AI = { {    1,  -2 },
                             { -0.5, 1.5 }};

      anpi::invert<T>(A,ai);

      const T eps = std::numeric_limits<T>::epsilon();

      for(int i = 0; i < A.rows(); i++) {
        for(int j = 0; j < A.rows(); j++) {
          BOOST_CHECK(std::abs(ai[i][j] - AI[i][j]) < eps * 10);
        }
        std::cout << std::endl;
      }

    }

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( Solver )

  BOOST_AUTO_TEST_CASE(inversion) {
    anpi::test::inversionTest<float>();
    anpi::test::inversionTest<double>();
  }

BOOST_AUTO_TEST_SUITE_END()