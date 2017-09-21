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
int main(){

	/* Pruebas: Descomposicióm LU, Método de Crout*/

    testLU tl ;//Creación del test Descomposición LU, Método de Crout
    //tl.test(1);//Prueba matrix 10x10
    //tl.test(3);//Prueba matrix 5x5
    //tl.test(3);//Prueba matrix 4x4
    tl.test(4);//Prueba matrix mal condicionada

    /*Fin Pruebas: Descomposicióm LU, Método de Crout*/








	/*  Prueba de crecion y acceso de matriz */
	
	// anpi::Matrix<float> A = { {1.f,2.f,3.f,4.f},
	// 						   {5.f,6.f,7.f,8.f},
	// 						   {9.f,10.f,11.f,12.f} };

	//A[0][0] = 1; //Asignar valor a una posicion de la matriz

	//A.printmatrix(); //Imprimir la matrix
	
	//int M00 = A[0][0]; //Acceder al valor de una posicion de la matriz 

	/* Fin prueba creacion y acceso de matriz*/

   //***************************************************************************//

	/* Prueba descomposion LU*/

	/* anpi::Matrix<float> A = {{-8,5,-6,4,1,-8,2,-7,3,-3},
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
    

    anpi::Matrix<float> A = { {1, -2, 2 ,-3},
                            {3, 4, -1, 1},
                            {2, -3,  2,-1},
                            {1, 1,-3, -2}
                        	};


	int n = A.rows();//Para tamano de matriz
	int m = A.cols();//Para tamano de matriz

	anpi::Matrix<float> LU(n,m,0.f);// Creacion de la matrix LU nxm, llena de ceros



    lucrout<float>CLU;//Para usar metodo de crout
	
	CLU.lu(A, LU);//Llamada a la descomposion LU

	LU.printmatrix();//Muestra resultado de la descomposion 

	*/

	/* Fin prueba descomposion LU */



}
