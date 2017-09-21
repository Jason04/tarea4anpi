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
#include "lucrout.h"
int main(){









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

	/* anpi::Matrix<float> A = { {1.f,2.f,3.f,4.f},
	 						   {5.f,6.f,7.f,8.f},
	 						   {9.f,10.f,11.f,12.f} };

	int n = 3;//Para tamano de matriz
	int m = 4;//Para tamano de matriz

	anpi::Matrix<float> LU(n,m,0.f);// Creacion de la matrix LU nxm, llena de ceros



    lucrout<float>CLU;//Para usar metodo de crout
	
	CLU.lu(A, LU);//Llamada a la descomposion LU

	LU.printmatrix();//Muestra resultado de la descomposion */

	/* Fin prueba descomposion LU */



}
