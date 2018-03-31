/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Gabriel
 * @Date  : 27.03.2018
 */

#include <boost/test/unit_test.hpp>

#include "qr.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
  namespace test {

    /// Test the given closed root finder
    template<typename T>
    void qrTest(const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         Matrix<T>&)>& qr) {


      // The result
      Matrix<T> Q;
      Matrix<T> R;

      // Test que Q sea una matriz ortogonal probando que Q^T * Q = I.
      {
        anpi::Matrix<float> A = { {12,-51, 4},
                                  {6,167,-63},
                                  {-4,24,-41} };

        anpi::Matrix<float> I = { {1,0, 0},
                                  {0,1,0},
                                  {0,0,1} };
        qr(A,Q,R);

        //BOOST_CHECK(gp==p);
      }

      // Test que R sea una matriz triangular superior
      {
        anpi::Matrix<float> A = { {12,-51, 4},
                                  {6,167,-63},
                                  {-4,24,-41} };
        qr(A,Q,R);

        //std::vector<size_t> gp= {1,0,3,2};
        //BOOST_CHECK(gp==p);
      }

      // Test que QR = A
      {
        anpi::Matrix<float> A = { {12,-51, 4},
                                  {6,167,-63},
                                  {-4,24,-41} };
        qr(A,Q,R);


        const T eps = std::numeric_limits<T>::epsilon();

        BOOST_CHECK(Ar.rows()==A.rows());
        BOOST_CHECK(Ar.cols()==A.cols());

        for (size_t i=0;i<Ar.rows();++i) {
          for (size_t j=0;j<Ar.cols();++j) {
            BOOST_CHECK(std::abs(Ar(i,j)-A(i,j)) < eps);
          }
        }
      }
    }

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( QR )

BOOST_AUTO_TEST_CASE(QR)
{
  anpi::test::qrTest<float>(anpi::qr<float>);
  anpi::test::qrTest<double>(anpi::qr<double>);
}


BOOST_AUTO_TEST_SUITE_END()
