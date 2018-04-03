/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Gabriel
 * @Date  : 27.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "qr.hpp"

#ifndef ANPI_QR_HPP
#define ANPI_QR_HPP

namespace anpi {

    template<typename T>
    bool solveQR (const anpi::Matrix <T>&A, std::vector<T>&x, const vector<T>&b){

        Matrix<T> Q;
        Matrix<T> R;
        anpi::qr(A,Q,R);

        ///METODO transpuesta
          anpi::Matrix<T> transpuestaQ(Q.cols(),Q.rows());
          for(int i =0; i < Q.cols(); i++){
              for(int j =0; j < Q.rows(); j++){
                  transpuestaQ[i][j] = Q[j][i];
              }
          }

        ///METODO inversa
          anpi::Matrix<T> inversaR(R.cols(),R.rows());
          for(int i =0; i < R.cols(); i++){
              for(int j =0; j < R.rows(); j++){
                  inversaR[i][j] = R[j][i]; ///Arreglar
              }
          }

        Matrix<T> Temp = inversaR * transpuestaQ;

        Matrix<T> x = Temp * b;

    }
}

#endif