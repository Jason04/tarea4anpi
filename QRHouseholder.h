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

        /*int n = A.rows();
        std::vector <T> X(n);        
        std::vector <T> U(n);        
        std::vector <T> V(n);
        anpi::Matrix<T> Qn(n,n,0.0);
        anpi::Matrix<T> Rn(n,n,0.0);

        T alpha = obtenerPrimerColumnayMagnitud(A, X);
        T magnitudU = obtenerUMagnitud(X, alpha, U);
        obtenerV(V, U, magnitudU); //V = U / ||U||
        obtenerQn(Qn, V); //I -2VVt
        obtenerRn(Qn,A, Rn);
        obtenerAPrima(Rn);

        //std::cout<<"Alpha: "<< alpha<<std::endl;
        //std::cout<<"||U||: "<< magnitudU<<std::endl;
        //Qn.printmatrix();
        Rn.printmatrix();*/

        anpi::Matrix<T> Aaux = A;

        for (int i = A.rows(); i >= 2; --i){
        
        std::vector <T> X(Aaux.rows()); 
        anpi::Matrix<T> Qn(Aaux.rows(),Aaux.rows(),0.0);
        anpi::Matrix<T> QExtendida(A.rows(),A.rows(),0.0);
        anpi::Matrix<T> Rn(Aaux.rows(),Aaux.rows(),0.0); 

        std::vector <T> U(Aaux.rows());        
        std::vector <T> V(Aaux.rows());
        
        T alpha = obtenerPrimerColumnayMagnitud(Aaux, X);      //X y ||X||
        T magnitudU = obtenerUMagnitud(X, alpha, U);        //U y ||U||
        obtenerV(V, U, magnitudU);                          //V = U / ||U||
        obtenerQn(Qn, V);                                   //I -2VVt
        extenderQ(Qn, QExtendida);
        if(i == A.rows()){
          Q = QExtendida;
        }
        else{
          Q = Q * QExtendida;
        }

        /*QExtendida.printmatrix();
        std::cout<<"Qextendida"<<std::endl;
        Q.printmatrix();
        std::cout<<"Q"<<std::endl;*/

        obtenerRn(QExtendida,A, Rn);                                //R = Qt*A
        guardarAaux(Aaux, Rn);
        //Aaux.printmatrix();
        //std::cout<<std::endl;
        }
    }
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
    void guardarAaux(anpi::Matrix<T>& aux,anpi::Matrix<T>& rn){
      anpi::Matrix<T> temp(aux.rows()-1,aux.rows()-1,0.0);
      aux = temp;
      for (int i = 1; i < rn.rows(); ++i){
        for (int j = 1; j < rn.rows(); ++j){
          aux[i-1][j-1] = rn[i][j];
        }
      }
    }

    T obtenerPrimerColumnayMagnitud(anpi::Matrix<T>& A, std::vector<T>& X){
      T alpha = 0.0;
      for (int i = 0; i < A.rows(); ++i){
        X[i] = A[i][0];
        alpha += A[i][0] * A[i][0];
      }
      return std::sqrt(alpha);
    }
    
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

    void obtenerV(std::vector<T>& V, std::vector<T>& U, T & magU){
      for (int i = 0; i < U.size(); ++i)
      {
        V[i] = U[i] / magU;
      }
    }

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

    void obtenerRn(anpi::Matrix<T>& qn, anpi::Matrix<T>& A, anpi::Matrix<T>& Rn){
        anpi::Matrix<T> Qt(A.rows(),A.cols(),0.0);
        transponer(qn, Qt);
        Rn = Qt * A;
    }
    void transponer(anpi::Matrix<T>& Qn, anpi::Matrix<T>& Qt){
      for (int i = 0; i < Qn.rows(); ++i){
        for (int j = 0; j < Qn.rows(); ++j){
          Qt[j][i] = Qn[i][j];
        }
      }
    }


    //Permite encontrar la solucion a un sistema de ecuaciones dado.
    bool solveQR(anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b){

        anpi::Matrix<T> LU(A.rows(),A.rows(),0.0);// se crea la metrix LU

        lu(A,LU);// Se descompone la matrix

        solveCrout(A.rows(),LU,b,x);// se realiza la sustitucion de variables

        return 1;
    }


    //Permite hacer sustitucion hacia a delante y hacia atras para encontrar
    //la solucion a un sistema de ecuaciones
   void solveHouseholder(int n, anpi::Matrix<T>& LU,const std::vector<T>& b,std::vector<T>& x){
       T y[n];
       T sum=0;
       for(int i=0;i<n;++i){//Sustitucion hacia adelante
            sum=0;
          for(int k=0;k<i;++k){
              sum+=LU[i][k] * y[k];
          }
          y[i]=(b[i]-sum)/LU[i][i];
       }
       for(int i=n-1;i>=0;--i){//Sustitucion hacia atras
          sum=0;
          for(int k=i+1;k<n;++k){
              sum+=LU[i][k]*x[k];
          }
          x[i]=(y[i]-sum);
       }
   }
   //Permite reeconstruir la matrix A, a partir de su descomposicion LU.
   //Ademas se obtiene la norma de la diferencia entre la matrix A y su reconstruccion
   T testQR(anpi::Matrix<T>& A, anpi::Matrix<T>& LU){
       int n = A.rows();
       //se crea la matrix para reecontruir
       anpi::Matrix<T> AReconstruida(n,n,0.0);// Create matrix nxm
       //Se realiza la recontruccion, utilizando a LU
       for (int i=0;i<n;i++){
           for (int j=0;j<n;j++){
               AReconstruida[i][j]=0;
               for (int k=0;k<=i;k++){
                   if(k<=j){
                       if(k==j){
                           AReconstruida[i][j]=AReconstruida[i][j]+LU[i][k];
                       }else{
                           AReconstruida[i][j]=AReconstruida[i][j]+ LU[i][k]*LU[k][j];
                       }
                   }
                }
             }

          }
       //Se muestra la matrix reeconstruida
      std::cout<<"Reconstruccion matriz A:"<<std::endl;
      AReconstruida.printmatrix();

      return norma(A,AReconstruida);
   }

  //Permite obtener la norma de la diferencia entre la matrix A y su reconstruccion
   T norma(anpi::Matrix<T>& A,anpi::Matrix<T>& R){

   }
};
#endif