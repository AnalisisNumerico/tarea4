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

#ifndef ANPI_QR_HPP
#define ANPI_QR_HPP

namespace anpi {

    template<typename T>
    void qr(const anpi::Matrix<T>& A,anpi::Matrix<T>& Q, anpi::Matrix<T>& R) {
        int k = 0; ///hay que cambiarlo en cada ciclo
        anpi::Matrix<T> Atemp = A;
        while(k<2){
            std::cout << "inicia " << k << std::endl;
            std::vector<T> fila;
             ///PICHA PICHA PICHA
            for (int i = 0+k; i < A.rows(); i++) {
                fila.push_back(Atemp[i][k]);
                std::cout << "fila, deberia ser(-49,168) " << fila[i-k] << std::endl ;
            }



            T norma = T(0);
        for (int i = 0; i < fila.size(); i++) {
            norma = norma + fila[i] * fila[i];
            std::cout << "veo " << fila[i] << std::endl ;

        }
        norma = std::sqrt(norma);
            std::cout << "norma, deberia ser(175) " << norma << std::endl ;


            std::vector<T> uVector;
        uVector.push_back(fila[0] - norma);
            std::cout << "fila, deberia ser(-49) " << fila[0] << std::endl ;
            std::cout << "uVector, deberia ser(-224) " << uVector[0] << std::endl ;
            T uNorma = (fila[0] - norma) * (fila[0] - norma);
        for (int i = 1; i < fila.size(); i++) {
            uVector[i] = fila[i];
            uNorma = uNorma + fila[i] * fila[i];
            std::cout << "uVector, deberia ser(-224,168) " << uVector[i] << std::endl ;

        }
            uNorma = std::sqrt(uNorma);

            std::cout << "uNorma, deberia ser(280) " << uNorma << ", " ;



            std::vector<T> vVector;
        for (int i = 0; i < fila.size(); i++) {
            vVector.push_back(7);
            vVector[i] = (1 / uNorma) * uVector[i];
            std::cout << "vVector, deberia ser(-0.8,0.6) " << vVector[i] << std::endl ;

        }
            //bien

        ///METODO identidad
        anpi::Matrix<T> identidad(A.rows(), A.cols());      ///VER BIEN DONDE SUMAR EL K
        for (int i = 0+k; i < identidad.cols(); i++) {
            for (int j = 0; j < identidad.rows(); j++) {
                if (i == j) {
                    identidad[i][j] = 1;
                } else {
                    identidad[i][j] = 0;
                }
            }
        }


        ///METODO Q = I - 2VV^T
        anpi::Matrix<T> Qtemp(A.rows(), A.cols());
            ///Metodo producto externo
            anpi::Matrix<float> multiVVT(vVector.size(),vVector.size());
            for (int i = 0; i < multiVVT.rows(); i++) {
                for (int j = 0; j < multiVVT.cols(); j++) {
                    multiVVT[i][j]=vVector[i]*vVector[j];
                }
            }

        ///METODO MATRIZ POR ESCALAR
        T escalar = T(2);
        for (int j = 0; j < multiVVT.rows(); j++) {
            for (int i = 0; i < multiVVT.cols(); i++) {
                multiVVT[i][j] = escalar * multiVVT[i][j];
            }
        }
        if(k>0){
            ///agranda matriz
            anpi::Matrix<T> agrandada(multiVVT.rows()+1, multiVVT.cols()+1);      ///VER BIEN DONDE SUMAR EL K
            for (int i = 0; i < agrandada.cols(); i++) {
                for (int j = 0; j < agrandada.rows(); j++) {
                    if(i<k or j<k){
                        if (i == j) {
                            agrandada[i][j] = 1;
                        } else {
                            agrandada[i][j] = 0;
                        }
                    }else{
                        agrandada[i][j] = multiVVT[i-k][j-k];
                    }

                }
            }
            multiVVT = agrandada;
        }



        ///METODO imprime matriz
        std::cout << "multiVVT = [";
        for (int i = 0; i < multiVVT.rows(); i++) {
            for (int j = 0; j < multiVVT.cols(); j++) {
                std::cout << multiVVT[i][j] << ", ";
            }
            std::cout << "" << std::endl;
        }
        std::cout << "]" << std::endl;

            //bien
            ///METODO imprime matriz
            std::cout << k<< "identidad = [";
            for(int i =0; i < identidad.rows(); i++){
                for(int j =0; j < identidad.cols(); j++){
                    std::cout << identidad[i][j] << ", " ;
                }
                std::cout << "" << std::endl;
            }
            std::cout << "]"<< std::endl;

        Qtemp = identidad - multiVVT;

            ///METODO imprime matriz
            std::cout << "Qtemp = [";
            for (int i = 0; i < Qtemp.rows(); i++) {
                for (int j = 0; j < Qtemp.cols(); j++) {
                    std::cout << Qtemp[i][j] << ", ";
                }
                std::cout << "" << std::endl;
            }
            std::cout << "]" << std::endl;
            if(k==0){
                Q = Qtemp;
            }else{
                Q = Q*Qtemp;
            }

            ///METODO imprime matriz
            std::cout << k<< "Q = [";
            for(int i =0; i < Q.rows(); i++){
                for(int j =0; j < Q.cols(); j++){
                    std::cout << Q[i][j] << ", " ;
                }
                std::cout << "" << std::endl;
            }
            std::cout << "]"<< std::endl;

            Atemp = Qtemp * A;
            ///Achica matriz
            for(int i =0; i < Atemp.rows(); i++){
                if(i!=k){
                    Atemp[i][k] = 0;
                }else{
                    Atemp[i][k] = 1;
                }
            }
            for(int i =1; i < Atemp.cols(); i++){
                Atemp[k][i] = 0;
            }///fin achica matriz

            ///METODO imprime matriz
            std::cout << k<< " holi Atemp = [";
            for(int i =0; i < Atemp.rows(); i++){
                for(int j =0; j < Atemp.cols(); j++){
                    std::cout << Atemp[i][j] << ", " ;
                }
                std::cout << "" << std::endl;
            }
            std::cout << "]"<< std::endl;



            k=k+1;

        }

        ///METODO transpuesta
        anpi::Matrix<float> transpuestaQ(Q.cols(), Q.rows());
        for (int i = 0; i < Q.cols(); i++) {
            for (int j = 0; j < Q.rows(); j++) {
                transpuestaQ[i][j] = Q[j][i];
            }
        }




    R = transpuestaQ * A;

    }









/*
    template<typename T>
    void qr(const anpi::Matrix<T>& A,anpi::Matrix<T>& Q, anpi::Matrix<T>& R) {


        std::vector<T> uAntes;
        std::vector<T> uActual;
        for (int j = 0; j < A.cols(); j++) {
            std::vector<T> fila;
            for (int i = 0; i < A.rows(); i++) {
                fila.push_back(7);
                fila[i] = A[i][j];
            }

            if (j == 0) {
                T norma = T(0);
                for (int i = 0; i < fila.size(); i++) {
                    norma = norma + fila[i]*fila[i];
                }
                norma = std::sqrt(norma);
                for (int i = 0; i < fila.size(); i++) {
                    uAntes.push_back(7);
                    uAntes[i] = (1/norma) * fila[i];
                    Q[i][j] = uAntes[i];
                }

            } else {
                T sumaAntes = T(0);
                for (int i = 0; i < Q.rows(); i++) {
                    sumaAntes = sumaAntes + Q[i][j-1];

                }

                for (int i = 0; i < Q.rows(); i++) {
                    uActual.push_back(7);
                    uActual[i] = (fila[i] - sumaAntes * Q[i][j-1]);
                }
                T norma = T(0);

                for (int i = 0; i < uActual.size(); i++) {
                    norma = norma + uActual[i]*uActual[i];
                }
                norma = std::sqrt(norma);
                std::cout <<"norma " <<norma << std::endl;
                for (int i = 0; i < Q.rows(); i++) {
                    uActual[i] = (1/norma) * uActual[i];
                    Q[i][j] = uActual[i];
                }
            }
            uAntes = uActual;
        }

        ///METODO transpuesta
        anpi::Matrix<float> transpuestaQ(Q.cols(),Q.rows());
        for(int i =0; i < Q.cols(); i++){
            for(int j =0; j < Q.rows(); j++){
                transpuestaQ[i][j] = Q[j][i];
            }
        }

        R = transpuestaQ * A;


    }*/



















    /*

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }


    template<typename T>
    void qr(const anpi::Matrix<T>& A,anpi::Matrix<T>& Q, anpi::Matrix<T>& R){
        int i,j,k;
        T scale,sigma,sum,tau;
        for (k=1;k<n;k++) {
            scale= T(0);
            for (i=k;i<=n;i++) scale=fmax(scale,fabs(A[i][k]));
            if (scale == T(0)) {//Singular  case.
                c[k]=d[k]=T(0);
            } else {    //Form Qk and Qk·A.
                for (i=k;i<=n;i++) A[i][k] /= scale;
                for (sum=T(0),i=k;i<=n;i++) sum += (A[i][k])*(A[i][k]);
                sigma= std::sqrt(sum)*anpi::sng(A[k][k]);
                A[k][k] += sigma;
                c[k]=sigma*A[k][k];
                d[k] = -scale*sigma;
                for (j=k+1;j<=n;j++) {
                    for (sum=T(0),i=k;i<=n;i++) sum += A[i][k]*A[i][j];
                    tau=sum/c[k];
                    for (i=k;i<=n;i++) a[i][j] -= tau*A[i][k];
                }
            }
        }
        d[n]=A[n][n];
    }*/



}

#endif

