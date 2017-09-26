/*
main.cpp
Fecha:16/09/2017
Curso: Analisis numerico para ingenieria
Asignacion: Tarea 4
Alumnos:
	Irene
	Gabriel Alfaro Herrera
	Jason Salazar Gonzalez
*/

#include <iostream>
#include <exception>
#include <cstdlib>
#include "TEST.h"
int main(){

	/* Pruebas: Descomposicióm LU, Método de Crout*/

    //testLU tl ;//Creación del test Descomposición LU, Método de Crout
    //tl.test(1);//Prueba matrix 10x10
    //tl.test(3);//Prueba matrix 5x5
    //tl.test(3);//Prueba matrix 4x4
   // tl.test(4);//Prueba matrix mal condicionada

    /*Fin Pruebas: Descomposicióm LU, Método de Crout*/

    /*PRUEBAS DE QR POR HOUSEHOLDER*/
    
    test t1;
    //t1.testQR1();
   //t1.testQR2();
    t1.testQR3();
    t1.testSolveQR3();

    /*//prueba mult M x V
    anpi::Matrix<float> A =  {{1, 0, 1},
                            {1, 2, 0},
                            {1, 0, 3}};
    
    std::vector <float> b(3); //vector b
    b[0]= 1;
    b[1]= 3;
    b[2]= 5;

    std::vector<float> mult(3);

    mult = A * b;

    std::cout<<"Resultado de mult: "<<std::endl;
    for (int i = 0; i < mult.size(); ++i){
        std::cout<<mult[i]<<",";
    }
    std::cout<<std::endl;*/

};
