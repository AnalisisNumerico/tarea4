/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>

#include "LUCrout.hpp"
#include "LUDoolittle.hpp"

#include <iostream>
#include <complex>

namespace anpi {
  namespace test {

    template<typename T>
    void unpackDoolittleTest() {

      anpi::Matrix<T> A = { {  1, 2, 3, 4 },
                            {  5, 6, 7, 8 },
                            {  9,10,11,12 },
                            { 13,14,15,16 } };

      anpi::Matrix<T> L = { {  1, 0, 0, 0 },
                            {  5, 1, 0, 0 },
                            {  9,10, 1, 0 },
                            { 13,14,15, 1 } };

      anpi::Matrix<T> U = { {  1, 2, 3, 4 },
                            {  0, 6, 7, 8 },
                            {  0, 0,11,12 },
                            {  0, 0, 0,16 } };

      anpi::Matrix<T> l,u;

      anpi::unpackDoolittle<T>(A,l,u);

      BOOST_CHECK(L==l);
      BOOST_CHECK(U==u);

    }

    template<typename T>
    void unpackCroutTest() {

      anpi::Matrix<T> A = { {  1, 2, 3, 4 },
                            {  5, 6, 7, 8 },
                            {  9,10,11,12 },
                            { 13,14,15,16 } };

      anpi::Matrix<T> L = { {  1, 0, 0, 0 },
                            {  5, 6, 0, 0 },
                            {  9,10,11, 0 },
                            { 13,14,15,16 } };

      anpi::Matrix<T> U = { {  1, 2, 3, 4 },
                            {  0, 1, 7, 8 },
                            {  0, 0, 1,12 },
                            {  0, 0, 0, 1 } };

      anpi::Matrix<T> l,u;

      anpi::unpackCrout<T>(A,l,u);

      BOOST_CHECK(L==l);
      BOOST_CHECK(U==u);

    }

    /// Test the given closed root finder
    template<typename T>
    void luTest(const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         std::vector<size_t>&)>& decomp,
                const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         Matrix<T>&)>& unpack) {

      // The result
      Matrix<T> LU;

      // Test if a non-square matrix is successfully detected
      {
        Matrix<T> A = {{ 1, 7, 6, 4},
                       { 2,17,27,17}};
        std::vector<size_t> p;
        try {
          decomp(A,LU,p);
          BOOST_CHECK_MESSAGE(false,"Rectangular matrix not properly catched");
        }
        catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Rectangular matrix properly detected");
        }
      }

      // Test pivoting
      {
        anpi::Matrix<T> A = { {-1,-2, 1, 2},
                              { 2, 0, 1, 2},
                              {-1,-1, 0, 1},
                              { 1, 1, 1, 1} };
        std::vector<size_t> p;
        decomp(A,LU,p);

        std::vector<size_t> gp= {1,0,3,2};
        BOOST_CHECK(gp==p);
      }

      // Test decomposition
      {
        // same matrix as before, but already permuted to force a clean decomposition
        anpi::Matrix<T> A = { { 2, 0, 1, 2},
                              {-1,-2, 1, 2},
                              { 1, 1, 1, 1},
                              {-1,-1, 0, 1} };

        std::vector<size_t> p;
        decomp(A,LU,p);

        Matrix<T> L,U;
        unpack(LU,L,U);
        Matrix<T> Ar=L*U;

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

BOOST_AUTO_TEST_SUITE( LU )

BOOST_AUTO_TEST_CASE(UnpackDoolittle) {

    anpi::test::unpackDoolittleTest<float>();
    anpi::test::unpackDoolittleTest<double>();

}

BOOST_AUTO_TEST_CASE(Doolittle) {

  anpi::test::luTest<float>(anpi::luDoolittle<float>,
                            anpi::unpackDoolittle<float>);
  anpi::test::luTest<double>(anpi::luDoolittle<double>,
                             anpi::unpackDoolittle<double>);
}

BOOST_AUTO_TEST_CASE(UnpackCrout) {

    anpi::test::unpackCroutTest<float>();
    anpi::test::unpackCroutTest<double>();

}

BOOST_AUTO_TEST_CASE(Crout) {

  anpi::test::luTest<float>(anpi::luCrout<float>,anpi::unpackCrout<float>);
  anpi::test::luTest<double>(anpi::luCrout<double>,anpi::unpackCrout<double>);

}

BOOST_AUTO_TEST_SUITE_END()

