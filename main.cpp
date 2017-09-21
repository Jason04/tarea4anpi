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
#include "Matrix.hpp"
#include "testlu.h"
#include "QRHouseholder.h"
int main(){

	/* Pruebas: Descomposicióm LU, Método de Crout*/

    //testLU tl ;//Creación del test Descomposición LU, Método de Crout
    //tl.test(1);//Prueba matrix 10x10
    //tl.test(3);//Prueba matrix 5x5
    //tl.test(3);//Prueba matrix 4x4
   // tl.test(4);//Prueba matrix mal condicionada

    /*Fin Pruebas: Descomposicióm LU, Método de Crout*/

    /*PRUEBAS DE QR POR HOUSEHOLDER*/

     qrhouseholder<float>QR;
     anpi::Matrix<float> A = {{1, 0, 1},
                            {1, 2, 0},
                            {1, 0, 3}};
    int n = A.rows();//Para tamano de matriz
    int m = A.cols();//Para tamano de matriz
    anpi::Matrix<float> Q(n,m,0.0);
    anpi::Matrix<float> R(n,m,0.0);
    QR.qr(A,Q,R);
};
