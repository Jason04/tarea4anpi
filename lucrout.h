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

    std::vector<T> indx;
    
    //Permite realizar la decomposicion utlizando el metodo de Crout.
    // A es la matriz de coefienciente del sistema
    // En el LU se almace la matriz L y la U
    void lu( anpi::Matrix<T>& A, anpi::Matrix<T>& LU ){
       std::vector <float> x(A.rows()); //vector x
       indx = x;                        //Se inicializa el vector indx

      ////////////////////////////************************************/////////////////////
        const double TINY=1.0e-40;
        LU = A;
        int n = A.rows();
        int i,imax,j,k;
        double big,temp;
        std::vector <double> vv(n);
        double d = 1.0; 
        
       // std::vector <double> indx(n); //vector 
        //indx.resize(n);

        for (i=0;i<n;i++) {
          big=0.0;
          for (j=0;j<n;j++){

            if ((temp=std::abs(LU[i][j])) > big){
              big = temp;
            }
          }
          if (big == 0.0){ 
            throw("Singular matrix in LUdcmp");
          }
          vv[i]=1.0/big;
        }

        for (k=0; k<n; k++) { 

          big=0.0;
         // imax=k;          
          for (i = k; i<n; i++) {

            temp = vv[i]*std::abs(LU[i][k]);
            if (temp > big) {
              big=temp;
              imax=i;
            }
          }
          if (k != imax) {
            for (j=0;j<n;j++) {

              temp=LU[imax][j];
              LU[imax][j]=LU[k][j];
              LU[k][j]=temp;
            }
            d = -d;
            vv[imax]=vv[k];
          }
          indx[k]=imax;

          if (LU[k][k] == 0.0){ 
            LU[k][k]=TINY;
          }
          for ( i=k+1;i<n;i++) {

            temp=LU[i][k] /= LU[k][k];
            for (int j=k+1;j<n;j++){
              LU[i][j] -= temp*LU[k][j];
            }
          }
        }

        //std::cout<<"deter: "<<det(LU)<<std::endl;
    

      ////////////////////*************************************///////////////////////////




     

    }

    //Permite encontrar la solucion a un sistema de ecuaciones dado.
   /* bool solveLU(anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b){

        anpi::Matrix<T> LU(A.rows(),A.rows(),0.0);// se crea la metrix LU

        lu(A,LU);// Se descompone la matrix
       //std::cout<<"LU"<<std::endl;
        //LU.printmatrix();


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
   }*/
   //Permite reeconstruir la matrix A, a partir de su descomposicion LU.
   //Ademas se obtiene la norma de la diferencia entre la matrix A y su reconstruccion
   T testLU(anpi::Matrix<T>& A, anpi::Matrix<T>& LU, anpi::Matrix<T>& Ar){
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
      Ar = AReconstruida;

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

  /*double det(anpi::Matrix<T>& LU){
    double dd = d;
    for (int i=0;i<LU.rows();i++){
      dd *= LU[i][i];
    }
    return dd;
  }*/

  //Solucion de sistema de ecuaciones 
  void solveLU(anpi::Matrix<T>& A, std::vector<T>& x, std::vector<T>& b){
    int i,ii=0,ip,j;
    int n = b.size();
    double sum;

    anpi::Matrix<T> LU(A.rows(),A.rows(),0.0);// se crea la metrix LU

    lu(A,LU);// Se descompone la matrix


    if (b.size() != n || x.size() != n){
      throw("LUdcmp::solve bad sizes");
    }
    for (i=0;i<n;i++){
     x[i] = b[i];
    }
    for (i=0;i<n;i++) {
      ip=indx[i];
      sum=x[ip];
      x[ip]=x[i];
      if (ii != 0){
        for (j=ii-1;j<i;j++){
         sum -= LU[i][j]*x[j];
        }
      }
      else if (sum != 0.0){
        ii=i+1;
      }
      x[i]=sum;
    }
    for (i=n-1;i>=0;i--) {
      sum=x[i];
      for (j=i+1;j<n;j++){
        sum -= LU[i][j]*x[j];
      }
      x[i]=sum/LU[i][i];
    }
  }
//solve con matrices
void solve(anpi::Matrix<T>& A,anpi::Matrix<T>& b, anpi::Matrix<T>& x){

  //anpi::Matrix<T> LU(A.rows(),A.rows(),0.0);// se crea la metrix LU

  //lu(A,LU);// Se descompone la matrix

  int n = b.rows();
  int i,j,m=b.cols();
  if (b.rows() != n || x.rows() != n || b.cols() != x.cols()){
    throw("LUdcmp::solve bad sizes");
  }
  std::vector <T> xx(n);
  indx = xx; 
  for (j=0;j<m;j++) {
    for (i=0;i<n;i++){
     xx[i] = b[i][j];
    }
    solveLU(A,xx,xx);
    for (i=0;i<n;i++){ 
      x[i][j] = xx[i];
    }
  }
}

void invert(anpi::Matrix<T>& A, anpi::Matrix<T>& Ai){
  int i,j;
 // ainv.resize(n,n);
  for (i=0;i<A.rows();i++) {
    Ai[i][i] = 1.;
  }

  solve(A,Ai,Ai);
}


void lu_reconstruccion( anpi::Matrix<T>& A, anpi::Matrix<T>& LU ){
     const double tiny = 1.0e-30;
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
              for(int p=0;p<k;++p){
                sum+=LU[k][p]*LU[p][j];
              }
              if(LU[k][k]== 0.0){
                LU[k][k] = tiny;
              }
              LU[k][j]=(A[k][j]-sum)/LU[k][k];
           }
        } 

} 

};

#endif
