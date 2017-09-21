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
#ifndef LUCROUT_H
#define LUCROUT_H
#include "Matrix.hpp"
#include <math.h>
#include <vector>
#include <cmath>
template<typename T>
class lucrout
{
public:
    //Permite realizar la decomposicion utlizando el metodo de Crout.
    // A es la matriz de coefienciente del sistema
    // En el LU se almace la matriz L y la U
    void lu( anpi::Matrix<T>& A, anpi::Matrix<T>& LU ){

        int n = A.rows();
        T sum = 0;    
        for(int k=0;k<n;++k){ //Se forma la parte L
           for(int i=k;i<n;++i){
               sum=0.;
              for(int p=0;p<k;++p){

                  sum+= LU[i][p] * LU[p][k];
              }
              LU[i][k]= A[i][k]-sum;
           }
           for(int j=k+1;j<n;++j){//Se forma la parte U
               sum=0.;
              for(int p=0;p<k;++p)sum+=LU[k][p]*LU[p][j];
              LU[k][j]=(A[k][j]-sum)/LU[k][k];
           }
        }
    }

    //Permite encontrar la solucion a un sistema de ecuaciones dado.
    bool solveLU(anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b){

        anpi::Matrix<T> LU(A.rows(),A.rows(),0.0);// se crea la metrix LU

        lu(A,LU);// Se descompone la matrix

        solveCrout(A.rows(),LU,b,x);// se realiza la sustitucion de variables

        return 1;
    }

    //Permite hacer sustitucion hacia a delante y hacia atras para encontrar
    //la solucion a un sistema de ecuaciones
   void solveCrout(int n, anpi::Matrix<T>& LU,const std::vector<T>& b,std::vector<T>& x){
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
   T testLU(anpi::Matrix<T>& A, anpi::Matrix<T>& LU){
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
       int n = A.rows();
       anpi::Matrix<T> Dif(n,n,0.0);// Create matrix nxm
       //Se calcula la diferencia entre las matrix y su reconstruccion
        Dif = A-R;

       //Se calcula la norma matricial Frobenius de la direfecia entre A y su reconstruccion
       T sum = 0;
       for (int i = 0; i < n; ++i) {
           for (int j = 0; j < n; ++j) {
                sum = sum + Dif[i][j]*Dif[i][j];
           }
       }
       return std::sqrt(sum);
   }
};

#endif
