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
int main(){










	//*****************************Prueba de crecion y acceso de matriz*******************
	int n = 4;//Para tamano de matriz
	int m = 4;//Para tamano de matriz

	anpi::Matrix<float> A(n,m,0.f);// Creacion de la matrix A nxm

	A[0][0] = 1; //Asignar valor a una posicion de la matriz

	int M00 = A[0][0]; //Acceder al valor de una posicion de la matriz 

	std::cout<<"Matrix[0][0]= "<<M00<<std::endl;//muestra la posicon 0,0 de la matriz A

	//******************************Fin prueba creacion y acceso de matriz*****************
}
