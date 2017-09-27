/*
lucrout.cpp:Contiene la implementacion
            para la solucion de ecuaciones
            usando descomposion LU Crout. 
Fecha:18/09/2017
Curso: Analisis numerico para ingenieria
Asignacion: Tarea 4
Alumnos:
  Irene
  Gabriel Alfaro Herrera
  Jason Salazar Gonzalez
*/
#ifndef QRHOUSEHOLDER_H
#define QRHOUSEHOLDER_H
#include "Matrix.hpp"
#include <math.h>
#include <vector>
#include <cmath>
template<typename T>
class qrhouseholder{
public:
    //Permite realizar la decomposicion utlizando el metodo de Householder.
    // A es la matriz de coefienciente del sistema
    // En el QR se crea una matrix ortogonal Q y una matriz triangulas superior R
    void qr( anpi::Matrix<T>& A, anpi::Matrix<T>& Q , anpi::Matrix<T>& R){
        anpi::Matrix<T> Aaux(A.rows(),A.rows(),0.0);    //Aaux es una matrix con el mismo tamaño de A, 
                                                              //pero en cada iteración va a disminuir en 1 cada dimensión
        Aaux = A;

        anpi::Matrix<T> Qtr(A.rows(),A.cols(),0.0);
        anpi::Matrix<T> Rn(A.rows(),A.rows(),0.0);      //Rn empleada hasta formar la verdadera Q

        for (int i = 1; i <= A.rows()- 1 ; ++i){        
          
          std::vector <T> X(Aaux.rows());                       //Primera columna de la la matrix A
          anpi::Matrix<T> Qn(Aaux.rows(),Aaux.rows(),0.0);      //Q en cada iteración
          anpi::Matrix<T> QExtendida(A.rows(),A.rows(),0.0);    //Es Qn pero rellno de 1's en la diagonal 
                                                                //y 0's en las demás columnas y filas
          std::vector <T> U(Aaux.rows());                       //Variables necesarias para el calculo de Q       
          std::vector <T> V(Aaux.rows());                       //              ||
          
          T alpha = obtenerPrimerColumnayMagnitud(Aaux, X);     //X y ||X||
          T magnitudU = obtenerUMagnitud(X, alpha, U);          //U y ||U||
          obtenerV(V, U, magnitudU);                            //V = U / ||U||
          obtenerQn(Qn, V);                                     //I -2VVt
          extenderQ(Qn, QExtendida);
          
          if(i == 1){                                    //Se valida se es la primera iteración para almacenar Qext en Q
            Q = QExtendida;
          }
          else{
            Q = Q * QExtendida;
          }

          obtenerRn(QExtendida, A, Rn, Qtr);                   //R = Qt*A
          guardarAaux(Aaux, Rn);
        }
        anpi::Matrix<T> Qt(A.rows(),A.rows(),0.0);
        transponer(Q, Qt);                                      //A partir de Q se obtiene Qt

        R = Qt*A;                                               //Se obtiene R con el Q y la A                                           
    }
    /*Ajusta la matriz Q (de tamaño inferior) al tamaño que Qe (de amaño superior)*/
    void extenderQ(anpi::Matrix<T>& Q,anpi::Matrix<T>& Qe){
      for (int i = Qe.rows()-1; i >= 0; --i){
        for (int j = Qe.rows()-1; j >= 0; --j){
          if(i==j){
            Qe[i][j] = 1;            
          }
        }
      }
      for (int i = Qe.rows()-1,  k = Q.rows()-1; k >= 0; --i,--k){
        for (int j = Qe.rows()-1,  z = Q.rows()-1; z >= 0; --j,--z){
          Qe[i][j] = Q[k][z];
        }
      }
    }
    /*Forma aux a partir de la parte inferior de rn*/
    void guardarAaux(anpi::Matrix<T>& aux,anpi::Matrix<T>& rn){
      anpi::Matrix<T> temp(aux.rows()-1,aux.rows()-1,0.0);
      aux = temp;
      for (int i = rn.rows()-1, k = aux.rows()-1; k>=0 ; --i, --k){
        for (int j = rn.cols()-1, z = aux.cols()-1; z>=0; --j,--z){
          aux[k][z]= rn[i][j];
        }
      }
    }
    /*Obtiene la primer columna de la matriz A, y ademas calcula la magnitud de dicho columna*/
    T obtenerPrimerColumnayMagnitud(anpi::Matrix<T>& A, std::vector<T>& X){
      T alpha = 0.0;
      for (int i = 0; i < A.rows(); ++i){
        X[i] = A[i][0];
        alpha += A[i][0] * A[i][0];
      }
      return std::sqrt(alpha);
    }
    /*Dado el vector X y el valor de alpha se obtiene el vector U*/
    T obtenerUMagnitud(std::vector<T>& X, T & alpha, std::vector<T>& U){
      T mag = 0.0;
      for (int i = 0; i < X.size(); ++i){
        if(i == 0){
          U[i] = X[i] -alpha;
        }
        else{
          U[i] = X[i];
        }
        mag += U[i] * U[i];        
      }
      return std::sqrt(mag);
    }
    /*Obtiene el Vector V a partir del vector U y su magnitud*/
    void obtenerV(std::vector<T>& V, std::vector<T>& U, T & magU){
      for (int i = 0; i < U.size(); ++i)
      {
        V[i] = U[i] / magU;
      }
    }
    /*Dado el vector V, se obtiene Qn*/
    void obtenerQn(anpi::Matrix<T>& Qn, std::vector<T>& V){
      for (int i = 0; i < V.size(); ++i){ 
        for (int j = 0; j < V.size(); ++j){
          if (i == j){
            Qn[i][j] = 1 - 2*(V[i]*V[j]);
          }
          else{
            Qn[i][j] = -2*V[i]*V[j];            
          }
        }
      }
    }
    /*Obtiene Rn a partir de qn y la matriz A*/
    void obtenerRn(anpi::Matrix<T>& qn, anpi::Matrix<T>& A, anpi::Matrix<T>& Rn, anpi::Matrix<T>& Qtr){
        transponer(qn, Qtr);
        Rn = Qtr * A;        
    }
    /*metodo para tranponer una matriz*/
    void transponer(anpi::Matrix<T>& Qn, anpi::Matrix<T>& Qt){
      for (int i = 0; i < Qn.rows(); ++i){
        for (int j = 0; j < Qn.rows(); ++j){
          Qt[j][i] = Qn[i][j];
        }
      }
    }

   double norma(anpi::Matrix<T>& A,anpi::Matrix<T>& QR){
     int n = A.rows();
     anpi::Matrix<T> Dif(n,n,0.0);// Create matrix nxm
     //Se calcula la diferencia entre las matrix y su reconstruccion
      Dif = A-QR;

     //Se calcula la norma matricial Frobenius de la direfecia entre A y su reconstruccion
     T sum = 0;
     for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
              sum = sum + Dif[i][j]*Dif[i][j];
         }
     }
     return std::sqrt(sum);
    }

    double testQR(anpi::Matrix<T>& A, anpi::Matrix<T>& Ar){
      anpi::Matrix<T> Q(A.rows(),A.rows(),0.0);
      anpi::Matrix<T> R(A.rows(),A.rows(),0.0);
      qr(A,Q,R);

      Ar = Q * R;
      /*std::cout<<"Q matrix"<<std::endl;
      Q.printmatrix();
      std::cout<<std::endl;
      std::cout<<"R matrix"<<std::endl;
      R.printmatrix();
      std::cout<<std::endl;*/
      return norma(A, Ar);
    }
    bool solveQR(anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b){
      anpi::Matrix<T> Q(A.rows(),A.rows(),0.0);
      anpi::Matrix<T> Qt(A.rows(),A.rows(),0.0);
      anpi::Matrix<T> R(A.rows(),A.rows(),0.0);

      qr(A,Q,R);

      transponer(Q, Qt);

      std::vector<T> bp = Qt * b;

      backwardSustitution(R, x, bp);
    }

    void backwardSustitution(anpi::Matrix<T>& R, std::vector<T>& x, const std::vector<T>& bp){
      int n = R.rows();
      double sum;
      for (int i=n-1; i >=0; --i){
        sum = bp[i];
        for (int j=i+1; j <n;++j){
          sum -= R[i][j]*x[j];
        }
        x[i] = sum/R[i][i];
      }
    }
};
#endif