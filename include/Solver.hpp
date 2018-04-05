#ifndef ANPI_SOLVER_HPP
#define ANPI_SOLVER_HPP

#include "LUDoolittle.hpp"

namespace anpi {

  /** faster method used for LU decomposition
   */
  template<typename T>
  inline void lu(const anpi::Matrix<T>& A,
                 anpi::Matrix<T> &      LU,
                 std::vector<size_t>&   p) {
    anpi::luDoolittle(A,LU,p);
  }

  /** method used to create  the permutation matrix given a
   * permutation vector
  **/
  template<typename T>
  void permutationMatrix(const std::vector<size_t>& p,
                         anpi::Matrix<T>&     pMatrix) {

    int n = p.size();
    anpi::Matrix<T> matrix(p.size(),p.size());

    for(int i = 0; i < n; i++) {
      matrix[i][p[i]] = T(1);
    }

    pMatrix = matrix;

  }

  /// method used to solve lower triangular matrices
  template<typename T>
  void forwardSubstitution(const anpi::Matrix<T>& L,
                           const std::vector<T>&  b,
                           std::vector<T>&        y) {

    int n = L.rows();

    y = b;

    for(int i = 0; i < n; i++) {
      y[i] =  T(1);
    }

    T sum;

    for(int m = 0; m < n ; m++) {
      sum = T(0);
      for(int i = 0; i < m; i++) {
        sum += L[m][i] * y[i];
      }
      y[m] =  (b[m] - sum)/L[m][m];
    }

  }

  /// method used to solve upper triangular matrices
  template<typename T>
  void backwardSubstitution(const anpi::Matrix<T>& U,
                            const std::vector<T>&  y,
                            std::vector<T>&        x) {
    int n = U.rows();

    x = y;

    for(int i = 0; i < n; i++) {
      x[i] =  T(1);
    }

    T sum;

    x[n-1] = y[n-1] / U[n-1][n-1];

    for(int i = (n-2); i >= 0; i--) {
      sum = T(0);
      for(int j = (n-1); j >= (i+1); j--) {
        sum += U[i][j] * x[j];
      }
      x[i] = (y[i] - sum) / U[i][i];
    }
  }

  template<typename T>
  T vectorNorm(std::vector<T>& b){
    T sum = T(0);
    for(int i = 0; i < b.size(); i++) {
      sum = b[i] * b[i];
    }
    return std::sqrt(sum);
  }

  template<typename T>
  bool solveLU(const anpi::Matrix<T>& A,
               std::vector<T>&        x,
               const std::vector<T>&  b) {

    const T eps = std::numeric_limits<T>::epsilon();

    anpi::Matrix<T> LU;
    std::vector<size_t> p;
    anpi::lu(A,LU,p);

    anpi::Matrix<T> L;
    anpi::Matrix<T> U;
    anpi::unpackDoolittle(LU,L,U);

    anpi::Matrix<T> P;
    anpi::permutationMatrix(p,P);

    std::vector<T>  bprima = b;

    std::vector <T> xPast(A.rows());

    do {

      anpi::Matrix<T>PB = P * bprima;

      std::vector<T> Pb(PB.rows());

      for(int i = 0; i < PB.rows(); i++) {
        Pb[i] = PB[i][0];
      }

      std::vector<T> y;
      anpi::forwardSubstitution(L,Pb,y);

      anpi::backwardSubstitution(U,y,x);

      for(int i = 0; i < A.rows(); i++) {
        x[i] -= xPast[i];
      }

      xPast = x;

      anpi::Matrix<T> BPrima = A * x;

      for(int i = 0; i < PB.rows(); i++) {
        bprima[i] = BPrima[i][0];
      }

      for(int i = 0; i < PB.rows(); i++) {
        bprima[i] -= b[i];
      }

    } while (anpi::vectorNorm<T>(bprima) > eps);

    return 1;

  }

}//anpi

#endif