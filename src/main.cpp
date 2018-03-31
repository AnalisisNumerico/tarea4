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
#include "qr.hpp"

int main() {

  // Some example code
  anpi::Matrix<float> a = { {1,2,3},{ 6, 5, 4} };

  ///METODO imprime matriz
  std::cout << "[";
  for(int i =0; i < a.rows(); i++){
    for(int j =0; j < a.cols(); j++){
      std::cout << a[i][j] << ", " ;
    }
    std::cout << "" << std::endl;
  }
  std::cout << "]"<< std::endl;
  ///FIN METODO imprime matriz

  ///METODO transpuesta
  anpi::Matrix<float> transpuesta(a.cols(),a.rows());
  for(int i =0; i < a.cols(); i++){
    for(int j =0; j < a.rows(); j++){
        transpuesta[i][j] = a[j][i];
    }
  }
  ///FIN METODO transpuesta

    /*anpi::Matrix<float> A = { {2,1},
                              { 2,1},
                              {1,5} };
*/
    anpi::Matrix<float> A = { {12,-51, 4},
                              {6,167,-63},
                              {-4,24,-41} };

    anpi::Matrix<float> Q(A.rows(),A.cols());
    anpi::Matrix<float> R;

    anpi::qr(A,Q,R);

    ///METODO imprime matriz
    std::cout << "Q = [";
    for(int i =0; i < Q.rows(); i++){
        for(int j =0; j < Q.cols(); j++){
            std::cout << Q[i][j] << ", " ;
        }
        std::cout << "" << std::endl;
    }
    std::cout << "]"<< std::endl;

    ///METODO imprime matriz
    std::cout << "R = [";
    for(int i =0; i < R.rows(); i++){
        for(int j =0; j < R.cols(); j++){
            std::cout << R[i][j] << ", " ;
        }
        std::cout << "" << std::endl;
    }
    std::cout << "]"<< std::endl;






    ///METODO identidad
    anpi::Matrix<float> identidad(3,1);
    for(int i =0; i < identidad.cols(); i++){
        for(int j =0; j < identidad.rows(); j++){
            if(i==j){
                identidad[i][j] = 1;
            }else{
                identidad[i][j] = 0;
            }
        }
    }
    ///FIN METODO identidad


    ///METODO vector canonico
    std::vector<float> canonico;
    int limite=3+1;
    for(int i =1; i < limite; i++){
        canonico.push_back(1);
        canonico[i] = 0;
    }
    ///FIN METODO vector canonico


    ///METODO imprime vector
    for(int i =1; i < limite; i++){
        //std::cout << canonico[i] << ", " ;
    }
    ///FIN METODO imprime vector



    std::vector<float> v1T = { -1, 3, -2 };

    ///Metodo producto externo
        std::vector<float> vector= { 0.96,0.28 };
        std::vector<float> vectorTranspuesto= { 0.96,0.28 };
        anpi::Matrix<float> productoExterno1(vector.size(),vectorTranspuesto.size());
        for (int i = 0; i < productoExterno1.rows(); i++) {
            for (int j = 0; j < productoExterno1.cols(); j++) {
                productoExterno1[i][j]=vector[i]*vectorTranspuesto[j];
            }
    }




    //anpi::Matrix<float> quepasa = { 175, 49, 168 };


  return EXIT_FAILURE;
}
  
