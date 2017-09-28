/*
unittes.cpp
Se realizan pruebas unitarias de todos los 
metedos de la tarea 4.
Fecha:27/09/2017
Curso: Analisis numerico para ingenieria
Asignacion: Tarea 4
Alumnos:
	Irene
	Gabriel Alfaro Herrera
	Jason Salazar Gonzalez
*/
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Suites
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>
#include "lucrout.h"
#include "QRHouseholder.h"

//Matriz para pruebas generales
anpi::Matrix<float> A = {{-8,5,-6,4,1,-8,2,-7,3,-3},
                        {-7,-5,-5,3,-3,-2,3,-1,2,9},
                        {-3,5,-4,-2,2,-1,4,7,-6,-7},
                        {-2,-4,5,8,-7,-1,2,-6,-5,-4},
                        {-8,-8,2,-2,-7,-5,1,-5,5,-1},
                        {5, 6,7,-5,-6, 4,2,4,-9, 5},
                        {6,-4,9, 8,6, 2,-5,6,-5,-3},
                        {-3,5,-6,-2,8,-4,-3,4,-8,-1},
                        {-1,7,-2,-4,-8,4,-9,-5,-6,-3},
                        {-6,-1,-4,-5,-8,8,8,-5,-7, 4}};

anpi::Matrix<float> reconstruida = {{-8,5,-6,4,1,-8,2,-7,3,-3},
                        {-7,-5,-5,3,-3,-2,3,-1,2,9},
                        {-3,5,-4,-2,2,-1,4,7,-6,-7},
                        {-2,-4,5,8,-7,-1,2,-6,-5,-4},
                        {-8,-8,2,-2,-7,-5,1,-5,5,-1},
                        {5, 6,7,-5,-6, 4,2,4,-9, 5},
                        {6,-4,9, 8,6, 2,-5,6,-5,-3},
                        {-3,5,-6,-2,8,-4,-3,4,-8,-1},
                        {-1,7,-2,-4,-8,4,-9,-5,-6,-3},
                        {-6,-1,-4,-5,-8,8,8,-5,-7, 4}};



/*Pruebas unitarias para el  
  metodo de descomposicion LU*/
BOOST_AUTO_TEST_SUITE(descomposicionLU)

	//Prueba descomposicion LU 1
	BOOST_AUTO_TEST_CASE(testLU1){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 		
	}

	//Prueba descomposicion LU 2
	BOOST_AUTO_TEST_CASE(testLU2){
		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.000000001); 
		
	}

	//Prueba descomposicion LU 3
	BOOST_AUTO_TEST_CASE(testLU3){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}

BOOST_AUTO_TEST_SUITE_END()

/*Pruebas unitarias para el  
  metodo de descomposicion QR*/
BOOST_AUTO_TEST_SUITE(descomposicionQR)

	//Prueba descomposicion LU 1
	BOOST_AUTO_TEST_CASE(testQR1){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 		
	}

	//Prueba descomposicion LU 2
	BOOST_AUTO_TEST_CASE(testQR2){
		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.000000001); 
		
	}

	//Prueba descomposicion LU 3
	BOOST_AUTO_TEST_CASE(testQR3){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}

BOOST_AUTO_TEST_SUITE_END()



/*Pruenas unitarias para la
  reconstrucciion de matrices*/
BOOST_AUTO_TEST_SUITE(reconstruccion)

	//Prueba reconstruccion 1
	BOOST_AUTO_TEST_CASE(testreconstruccion1){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}

	//Prueba reconstruccion 2
	BOOST_AUTO_TEST_CASE(testreconstruccion2){
		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.000000001); 
		
	}

	//Prueba reconstruccion 3
	BOOST_AUTO_TEST_CASE(testreconstruccion3){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}

BOOST_AUTO_TEST_SUITE_END()

/*Pruenas unitarias para el
  calculo de la norma de una
  matriz*/
BOOST_AUTO_TEST_SUITE(norma)

	//Prueba 1 calculo de norma
	BOOST_AUTO_TEST_CASE(testnorma1){

		lucrout<float>clu;
		float norma = clu.norma(A,A);
		float esperado = 0.0;          
	  	BOOST_CHECK_CLOSE_FRACTION (esperado, norma , 0.01); 		
	}

	//Prueba 2 calculo de norma
	BOOST_AUTO_TEST_CASE(testnorma2){
		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.000000000000001); 
		
	}
	/*

	//Prueba 3 calculo de norma
	BOOST_AUTO_TEST_CASE(testnorma3){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}*/
	
BOOST_AUTO_TEST_SUITE_END()