#ifndef LUCROUT_H
#define LUCROUT_H
#include "Matrix.hpp"
#include <math.h>
#include <vector>
template<typename T>
class lucrout
{
public:
    //Permite realizar la decomposicion utlizando el metodo de Crout.
    // A es la matriz de coefienciente del sistema
    // En el LU se almace la matriz L y la U
    void lu( anpi::Matrix<T>& A, anpi::Matrix<T>& LU ){

        int n = A.rows();
        float sum = 0;    
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
/*
    //Permite encontrar la solucion a un sistema de ecuaciones dado.
    bool solveLU(anpi::Matrix<T>& A, std::vector<float>& x, const std::vector<float>& b){

        anpi::Matrix<float> LU(A.rows(),A.rows(),0.0);// se crea la metrix LU

        lu(A,LU);// Se descompone la matrix

        solveCrout(A.rows(),LU,b,x);// se realiza la sustitucion de variables

        A.printmatrix(A,"A");//Se muestra la matrix A

        std::cout<<"Vector b"<<std::endl;//Se muestra el vector b

         std::cout<<"b = {";
        for (int i = 0; i < A.rows(); ++i) {
            std::cout<<b[i];
            if(i<A.rows()-1){
                std::cout<<",";
            }
        }
        std::cout<<"}"<<std::endl;
        std::cout<<""<<std::endl;


        LU.printmatrix(LU,"LU");// Se muestra la matriz descompuesta

        std::cout<<"Solutions"<<std::endl;// Se muestra la solucion del sistema
        for (int i = 0; i < A.rows(); ++i) {
            std::cout<<" x"<<i<<"= "<<x[i]<<std::endl;
        }


        return 1;

    }
    //Permite hacer sustitucion hacia a delante y hacia atras para encontrar
    //la solucion a un sistema de ecuaciones
   void solveCrout(int n, anpi::Matrix<T>& LU,const std::vector<float>& b,std::vector<float>& x){
       float y[n];
       float sum=0;
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
   void testLU(anpi::Matrix<T>& A, anpi::Matrix<T>& LU){
       int n = A.rows();
       //se crea la matrix para reecontruir
       anpi::Matrix<float> AReconstruida(n,n,0.0);// Create matrix nxm
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
       AReconstruida.printmatrix(AReconstruida,"A reconstruida");
       std::cout<<"Norma de la diferencia entre A y su reconstruccion:"<<std::endl;
       //Se muestra la norma de la diferencia entre la matrix A y su reconstruccion
       std::cout<<"       Norma = "<<norma(A,AReconstruida)<<std::endl;


   }

   //Permite obtener la norma de la diferencia entre la matrix A y su reconstruccion
   float norma(anpi::Matrix<T>& A,anpi::Matrix<T>& R){
       int n = A.rows();
       anpi::Matrix<float> Dif(n,n,-1.1);// Create matrix nxm
       //Se calcula la diferencia entre las matrix y su reconstruccion
       for (int i = 0; i < n; ++i) {
           for (int j = 0; j < n; ++j) {
               Dif[i][j] = A[i][j]-R[i][j];
           }

       }
       //Se calcula la norma de la direfecia entre A y su reconstruccion
       int sum = 0;
       for (int i = 0; i < n; ++i) {
           for (int j = 0; j < n; ++j) {
                sum = sum + pow(Dif[i][j],2);
           }
       }
       return sqrt(sum);
   }
*/



};

#endif // LUCROUT_H
