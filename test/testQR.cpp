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
        Matrix<T> A = { {12,-51, 4},
                                  {6,167,-63},
                                  {-4,24,-41} };

        Matrix<T> I = { {1,0, 0},
                                  {0,1,0},
                                  {0,0,1} };

        qr(A,Q,R);
          ///METODO transpuesta
          Matrix<T> transpuestaQ(Q.cols(),Q.rows());
          for(int i =0; i < Q.cols(); i++){
              for(int j =0; j < Q.rows(); j++){
                  transpuestaQ[i][j] = Q[j][i];
              }
          }

          Matrix<T> supuestaI = transpuestaQ * Q;
          //const T eps = std::numeric_limits<T>::epsilon();
          const T eps = 1e-6;

          for(int i =0; i < I.cols(); i++){
              for(int j =0; j < I.rows(); j++){
                  if(i==j){
                      BOOST_CHECK(std::abs(1-supuestaI(i,j)) < eps);
                  }else{
                      BOOST_CHECK(std::abs(0-supuestaI(i,j)) < eps);
                  }
              }
          }

      }

      // Test que R sea una matriz triangular superior
      {
        Matrix<T> A = { {12,-51, 4},
                                  {6,167,-63},
                                  {-4,24,-41} };

          Matrix<T> verdaderoR  = { {14, 21, -14},
                                          {0, 175, -70},
                                          {0, 0, -35} };

        qr(A,Q,R);
          //const T eps = std::numeric_limits<T>::epsilon();
          const T eps = 1e-4;

          for(int i =0; i < R.cols(); i++){
              for(int j =0; j < R.rows(); j++){
                  if(i<=j){
                      BOOST_CHECK(std::abs(R(i,j)) > eps);
                  }else{ ///CUIDADO CON ESTE ELSE, REVISAR QUE EL EPS SEA GRANDE
                      BOOST_CHECK(std::abs(R(i,j)) < eps);
                  }
              }
          }

      }

      // Test que QR = A
      {
        Matrix<T> A = { {12,-51, 4},
                                  {6,167,-63},
                                  {-4,24,-41} };
        qr(A,Q,R);
          Matrix<T> supuestaA = Q*R;
          //const T eps = std::numeric_limits<T>::epsilon();
          const T eps = 1e-4;

        for (size_t i=0;i<A.rows();++i) {
          for (size_t j=0;j<A.cols();++j) {
            BOOST_CHECK(std::abs(A(i,j)-supuestaA(i,j)) < eps);
          }
        }

      }



        // Test solveQR
        {
            anpi::Matrix<T> A = { {1,2,-3},
                                  {-3,4,5},
                                  {7,1,2} };
            std::vector<T> b = {1,5,6};
            std::vector<T> verdaderoX = {0.561798, 1.01124, 0.52809};
            std::vector<T> x;

            anpi::solveQR(A,x,b);

            const T eps = 1e-5;

            for (size_t i=0;i<verdaderoX.size();++i) {

                BOOST_CHECK(std::abs(verdaderoX[i] - x[i]) < eps);
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
