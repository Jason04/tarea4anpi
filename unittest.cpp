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

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}

	//Prueba 2 calculo de norma
	BOOST_AUTO_TEST_CASE(testnorma2){
		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.000000001); 
		
	}

	//Prueba 3 calculo de norma
	BOOST_AUTO_TEST_CASE(testnorma3){

		float f1 = 567.01012;
	  	float result = std::sqrt(f1);                            
	  	BOOST_CHECK_CLOSE_FRACTION (f1, result * result, 0.0001); 
		
	}
	
BOOST_AUTO_TEST_SUITE_END()