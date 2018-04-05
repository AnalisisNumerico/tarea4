#ifndef ANPI_INVERSION_HPP
#define ANPI_INVERSION_HPP

#include "Solver.hpp"

namespace anpi {

  template<typename T>
  void invert(const anpi::Matrix<T>& A, anpi::Matrix<T>& Ai ){

    anpi::Matrix<T> LU;
    std::vector<size_t> p;
    anpi::lu(A,LU,p);

    anpi::Matrix<T> P;
    anpi::permutationMatrix(p,P);

    int n = A.rows();

    std::vector<T> b(n);
    std::vector<T> x;

    anpi::Matrix<T> ai(n,n);

    for(int i = 0; i < n ; i++) {
      for(int k = 0; k < n; k++) {
        b[k] = P[i][k];
      }
      anpi::solveLU(A,x,b);
      for(int k = 0; k < n; k++) {
        ai[k][i] = x[k];
      }
    }

    Ai = ai;

  }

}//anpi

#endif