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
                 //test1();
                 break;
              case 2 :
                 std::cout<<"Matrix 5x5------------------------"<<std::endl;
                 test2();
                 break;
              case 3 :
                  std::cout<<"Matrix 4x4------------------------"<<std::endl;
                //  test3();
                   break;
              default :
                 std::cout<<"Matrix Mal Condicionada------------------------"<<std::endl;
                // test4();
                 break;
        }
        std::cout<<"---------------------Fin Test "<<type<<"------------------------"<<std::endl;
        std::cout<<" "<<std::endl;
    }

private:
    //Prueba del metodo de Croat con una matrix de 10x10
   /*  void test1(){

         int n = 10; //# of rows
         int m = 10; //# of colums

          //Data to fill matrix A
          const float coe_a[n*m] = {-8,5,-6,4,1,-8,2,-7,3,-3,
                            -7,-5,-5,3,-3,-2,3,-1,2,9,
                            -3,5,-4,-2,2,-1,4,7,-6,-7,
                            -2,-4,5,8,-7,-1,2,-6,-5,-4,
                            -8,-8,2,-2,-7,-5,1,-5,5,-1,
                             5, 6,7,-5,-6, 4,2,4,-9, 5,
                             6,-4,9, 8,6, 2,-5,6,-5,-3,
                            -3,5,-6,-2,8,-4,-3,4,-8,-1,
                            -1,7,-2,-4,-8,4,-9,-5,-6,-3,
                            -6,-1,-4,-5,-8,8,8,-5,-7, 4
                            };
          std::vector <float> b(n); //vector b

          int aux [] = {84,-109,201,-62,-37,-27,38,258,77,-101};

          b.insert (b.begin(), aux, aux+10);

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
        //Verificacion de funciones
        CLU.testLU( A, LU);

     }*/
     //Prueba del metodo de Croat con una matrix de 5x5
     void test2(){
      //Se crea matrix de 4x4
           anpi::Matrix<float> A = { {1, -2, 2 ,-3},
                            {3, 4, -1, 1},
                            {2, -3,  2,-1},
                            {1, 1,-3, -2}
                          };
           int n = A.rows();//Para tamano de matriz
           int m = A.cols();//Para tamano de matriz



          std::vector <float> b(n); //vector b
          int aux2[] = {-2,4,3,-5};
          b.insert (b.begin(), aux2, aux2+5);

          std::vector <float> x(n); //vector x


          lucrout<float>CLU;
          //Se resuelve el sistema de ecuaciones
          CLU.solveLU(A,x,b);


       // anpi::Matrix<float> LU(A.rows(),A.rows(),0.0);// Create LU matrix

       // CLU.lu(A,LU);
        //Verifiacion de funciones
       // CLU.testLU( A, LU);
     }
     //Prueba del metodo de Croat con una matrix de 4x4
    /* void test3(){
         int n = 4; //# of rows
         int m = 4; //# of colums

          //Data to fill matrix A
         float coe_a[n*m]={ 1, -2, 2 ,-3,
                            3, 4, -1., 1,
                            2, -3,  2,-1,
                            1, 1,-3, -2};

          std::vector <float> b(n); //vector b
          float aux3 [] ={15,-6,17,-7};
          b.insert (b.begin(), aux3, aux3+4);


          std::vector <float> x(n); //vector x

          anpi::Matrix<float> A(n,m,0.0);// Create A matrix nxm

         //Fill A matrix
         for (int i = 0; i < n; ++i) {
             for (int j = 0; j < m; ++j) {
                 A[i][j] = coe_a[n*i+j] ;// fill matrix with coe_a
             }
         }

        lucrout<float>CLU;
        //Se resuelte el sistema de ecuaciones
        CLU.solveLU(A,x,b);


        anpi::Matrix<float> LU(A.rows(),A.rows(),0.0);// Create LU matrix

        CLU.lu(A,LU);
        //Verificar funciones
        CLU.testLU( A, LU);
     }
     //Prueba del metodo de Croat con una matrix mal condicionada
     void test4(){
         int n = 2; //# of rows
         int m = 2; //# of colums

          //Data to fill matrix A
         float coe_a[n*m]={ 3, 4,
                            3,4.00001
                          };

          std::vector <float> b(n); //vector b
          float aux4[] = {7,7.00001};
          b.insert (b.begin(), aux4, aux4+2);


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
        //Verificacion de funciones
        CLU.testLU( A, LU);
     }*/


};

#endif // TESTLU_H
