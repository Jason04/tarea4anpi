/*
testlu.h:Pruebas para encontrar la solucion
         de ecuaciones usando el metodo 
         descomposion Lu por crout. 
Fecha:20/09/2017
Curso: Analisis numerico para ingenieria
Asignacion: Tarea 4
Alumnos:
  Irene
  Gabriel Alfaro Herrera
  Jason Salazar Gonzalez
*/
#ifndef TESTLU_H
#define TESTLU_H
#include "Matrix.hpp"
#include "lucrout.h"
class testLU
{
public:

    //Menu de pruebas para la descomposion LU: metodo de Croat
    void test(int type){
        std::cout<<"------------Test"<<type<<": Descomposición LU: Método Crout: ";
        switch(type) {
              case 1 :
                 std::cout<<"Matrix 10x10----------------------"<<std::endl;
                 test1();
                 break;
              case 2 :
                 std::cout<<"Matrix 5x5------------------------"<<std::endl;
                 //test2();
                 break;
              case 3 :
                  std::cout<<"Matrix 4x4------------------------"<<std::endl;
                  test3();
                   break;
              default :
                 std::cout<<"Matrix Mal Condicionada------------------------"<<std::endl;
                 test4();
                 break;
        }
        std::cout<<"---------------------Fin Test "<<type<<"------------------------"<<std::endl;
        std::cout<<" "<<std::endl;
    }

private:
    //Prueba del metodo de Croat con una matrix de 10x10
     void test1(){
          //Se crea matrix de 4x4
           anpi::Matrix<float> A = {{-8,5,-6,4,1,-8,2,-7,3,-3},
                            {-7,-5,-5,3,-3,-2,3,-1,2,9},
                            {-3,5,-4,-2,2,-1,4,7,-6,-7},
                            {-2,-4,5,8,-7,-1,2,-6,-5,-4},
                            {-8,-8,2,-2,-7,-5,1,-5,5,-1},
                             {5, 6,7,-5,-6, 4,2,4,-9, 5},
                             {6,-4,9, 8,6, 2,-5,6,-5,-3},
                            {-3,5,-6,-2,8,-4,-3,4,-8,-1},
                            {-1,7,-2,-4,-8,4,-9,-5,-6,-3},
                            {-6,-1,-4,-5,-8,8,8,-5,-7, 4}
                            };
           int n = A.rows();//Para tamano de matriz
           int m = A.cols();//Para tamano de matriz

          std::vector <float> b(n); //vector b

          int aux [] = {84,-109,201,-62,-37,-27,38,258,77,-101};

          b.insert (b.begin(), aux, aux+10);

          std::vector <float> x(n); //vector x

          lucrout<float>CLU;
          //Se resuelve el sistema de ecuaciones
          CLU.solveLU(A,x,b);


           A.printmatrix();//Se muestra la matrix A

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
             std::cout<<"Solutions"<<std::endl;// Se muestra la solucion del sistema
             for (int i = 0; i < A.rows(); ++i) {
              std::cout<<" x"<<i<<"= "<<x[i]<<std::endl;
              }
        //Prueba de resonstruccion       
        anpi::Matrix<float> LU(A.rows(),A.cols(),0.0);// Create LU matrix
         CLU.lu(A,LU);
        //Verifiacion de funciones
       float norma = CLU.testLU( A, LU);
       std::cout<<"Norma de la diferencia entre A y su reconstruccion:"<<std::endl;
       //Se muestra la norma de la diferencia entre la matrix A y su reconstruccion
       std::cout<<"       Norma = "<<norma<<std::endl;

     }
     //Prueba del metodo de Croat con una matrix de 5x5
   /*  void test2(){

         int n = 5; //# of rows
         int m = 5; //# of colums

          //Data to fill matrix A
         float coe_a[n*m]={2,1,1,3,2,
                       1.,2.,2.,1.,1.,
                       1.,2.,9.,1.,5.,
                       3.,1.,1.,7.,1.,
                       2.,1.,5.,1.,8.};


          std::vector <float> b(n); //vector b
          int aux2[] = {-2,4,3,-5,1};
          b.insert (b.begin(), aux2, aux2+5);


          std::vector <float> x(n); //vector x

          anpi::Matrix<float> A(n,m,0.0);// Create A matrix nxm

         //Fill A matrix
         for (int i = 0; i < n; ++i) {
             for (int j = 0; j < m; ++j) {
                 A[i][j] = coe_a[n*i+j] ;// fill matrix with coe_a
             }
         }

        lucrout<float>CLU;
        //Se resuelve el sistema de ecuaciones
        CLU.solveLU(A,x,b);


        anpi::Matrix<float> LU(A.rows(),A.rows(),0.0);// Create LU matrix

        CLU.lu(A,LU);
        //Verifiacion de funciones
        CLU.testLU( A, LU);
     }
     */

     //Prueba del metodo de Croat con una matrix de 4x4
    void test3(){
           //Se crea matrix de 4x4
           anpi::Matrix<float> A = { {1.f, -2.f, 2.f ,-3.f},
                            {3.f, 4.f, -1.f, 1.f},
                            {2.f, -3.f,  2.f,-1.f},
                            {1.f, 1.f,-3.f, -2.f}
                          };
           int n = A.rows();//Para tamano de matriz
           int m = A.cols();//Para tamano de matriz


          std::vector <float> b(n); //vector b
         // float aux2[] = {15,-6,17,-7};
          float aux2[] = {-2.f,4.f,3.f,-5.f};
          b.insert (b.begin(), aux2, aux2+4);

          std::vector <float> x(n); //vector x


          lucrout<float>CLU;
          //Se resuelve el sistema de ecuaciones
          CLU.solveLU(A,x,b);


           A.printmatrix();//Se muestra la matrix A

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
             std::cout<<"Solutions"<<std::endl;// Se muestra la solucion del sistema
             for (int i = 0; i < A.rows(); ++i) {
              std::cout<<" x"<<i<<"= "<<x[i]<<std::endl;
              }
        //Prueba de resonstruccion       
        anpi::Matrix<float> LU(A.rows(),A.cols(),0.0);// Create LU matrix
         CLU.lu(A,LU);
        //Verifiacion de funciones
       float norma = CLU.testLU( A, LU);
       std::cout<<"Norma de la diferencia entre A y su reconstruccion:"<<std::endl;
       //Se muestra la norma de la diferencia entre la matrix A y su reconstruccion
       std::cout<<"       Norma = "<<norma<<std::endl;


        
     }
     //Prueba del metodo de Croat con una matrix mal condicionada
      void test4(){
        //Se crea matrix de 4x4
           anpi::Matrix<float> A = { { 1.f, 0.5f,1.f/3.f},
                                     {0.5f,1.f/3.f,1.f/4.f},
                                     {1.f/3.f,1.f/4.f,1.f/5.f}};
           int n = A.rows();//Para tamano de matriz
           int m = A.cols();//Para tamano de matriz
           std::cout<<"Matirz A:"<<std::endl;
           A.printmatrix();//Se muestra la matrix A

           std::vector <float> b(n); //vector b
           float aux4[] = {3.0,2.1,4.0};
           b.insert (b.begin(), aux4, aux4+3);


          std::vector <float> x(n); //vector x

         lucrout<float>CLU;
         //Se resuelve el sistema de ecuaciones
         CLU.solveLU(A,x,b);

          std::cout<<""<<std::endl;
          std::cout<<"Solutions"<<std::endl;// Se muestra la solucion del sistema
          for (int i = 0; i < A.rows(); ++i) {
            std::cout<<" x"<<i<<"= "<<x[i]<<std::endl;
          }


        anpi::Matrix<float> LU(A.rows(),A.rows(),0.0);// Create LU matrix

        CLU.lu(A,LU);
        
        //Verifiacion de funciones
        float norma = CLU.testLU( A, LU);
        std::cout<<"Norma de la diferencia entre A y su reconstruccion:"<<std::endl;
        //Se muestra la norma de la diferencia entre la matrix A y su reconstruccion
        std::cout<<"       Norma = "<<norma<<std::endl;
     }


};

#endif // TESTLU_H
