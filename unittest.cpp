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
//#define BOOST_TEST_DYN_LINK  
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

anpi::Matrix<float> B = {{1.f/2.f, 1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f},
						{1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f},
						{1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f},
						{1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f, 1.f/9.f},
						{1.f/6.f, 1.f/7.f, 1.f/8.f, 1.f/9.f, 1.f/10.f}};

anpi::Matrix<float> C = {{1,0,1},{1,2,0},{1,0,3}};

/*Pruebas unitarias para el metodo de reconstrucción de matrices por LU y QR*/
BOOST_AUTO_TEST_SUITE(Reconstruccion_norma_Matriz_A)
	anpi::Matrix<float> Ar(A.rows(),A.rows(),0.0);
	float norma;

	//Prueba 1 descomposicion QR
	BOOST_AUTO_TEST_CASE(testQR1){
		qrhouseholder<float>QR;
		norma=QR.testQR(A, Ar);
		for (int i = 0; i < A.rows(); ++i){
			for (int j = 0; j < A.rows(); ++j){
				BOOST_CHECK_CLOSE_FRACTION (A[i][j],Ar[i][j], 0.01);
			}
		}
	}
	//Prueba 1 de Norma con QR
	BOOST_AUTO_TEST_CASE(testNormaQR1){
		BOOST_CHECK(std::abs(norma-0) < 1);
	}

	//Prueba 1 descomposicion LU
	BOOST_AUTO_TEST_CASE(testLU1){
	anpi::Matrix<float> lu(A.rows(),A.rows(),0.0);
		lucrout<float>LU;
		LU.lu_reconstruccion(A, lu);
		norma=LU.testLU(A, lu, Ar);
		for (int i = 0; i < A.rows(); ++i){
			for (int j = 0; j < A.rows(); ++j){
				BOOST_CHECK_CLOSE_FRACTION (A[i][j],Ar[i][j], 0.01);
			}
		}
	}
	//Prueba 1 de Norma con LU
	BOOST_AUTO_TEST_CASE(testNormaLU1){
		BOOST_CHECK(std::abs(norma-0) < 1);
	}
BOOST_AUTO_TEST_SUITE_END()

/*Pruebas unitarias para el metodo de reconstrucción de matrices MAL CONDICIONADA por LU y QR*/
BOOST_AUTO_TEST_SUITE(Reconstruccion_norma_Matriz_B)
	anpi::Matrix<float> Ar(B.rows(),B.rows(),0.0);
	float norma;

	//Prueba 1 descomposicion QR
	BOOST_AUTO_TEST_CASE(testQR1){
		qrhouseholder<float>QR;
		norma=QR.testQR(B, Ar);
		for (int i = 0; i < B.rows(); ++i){
			for (int j = 0; j < B.rows(); ++j){
				BOOST_CHECK_CLOSE_FRACTION (B[i][j],Ar[i][j], 0.01);
			}
		}
	}
	//Prueba 1 de Norma con QR
	BOOST_AUTO_TEST_CASE(testNormaQR1){
		BOOST_CHECK(std::abs(norma-0) < 1);
	}

	//Prueba 1 descomposicion LU
	BOOST_AUTO_TEST_CASE(testLU1){
	anpi::Matrix<float> lu(B.rows(),B.rows(),0.0);
		lucrout<float>LU;
		LU.lu_reconstruccion(B, lu);
		norma=LU.testLU(B, lu, Ar);
		for (int i = 0; i < B.rows(); ++i){
			for (int j = 0; j < B.rows(); ++j){
				BOOST_CHECK_CLOSE_FRACTION (B[i][j],Ar[i][j], 0.01);
			}
		}
	}
	//Prueba 1 de Norma con LU
	BOOST_AUTO_TEST_CASE(testNormaLU1){
		BOOST_CHECK(std::abs(norma-0) < 1);
	}
BOOST_AUTO_TEST_SUITE_END()

//Prubas sobre solucion del sistema de ecuaciones
BOOST_AUTO_TEST_SUITE(Vector_Solucion)

	//Prueba 1 descomposicion QR
	BOOST_AUTO_TEST_CASE(testsolveQR){
		std::vector <float> b ={2,4,6}; //vector b	
	    std::vector <float> x(3); //vector x
	    std::vector <float> xesperado ={0,2,2};
		qrhouseholder<float>QR;
		QR.solveQR(C, x, b);
		    for (int i = 0; i < x.size(); ++i){
		    	BOOST_CHECK(std::abs(xesperado[i]-x[i]) < 0.01);
		    }
	}
	BOOST_AUTO_TEST_CASE(testsolveLU){
		std::vector <float> b ={2,4,6,8,10,12,14,16,18,20}; //vector b	
	    std::vector <float> x(10); //vector x
	    std::vector <float> xesperado ={-4.3816, 0.7757,2.6859,-0.2148,1.2450,1.7639,-1.1757,
	    								0.0109,-1.0824,0.4719};
			lucrout<float>lu;
		lu.solveLU(A, x, b);
		    for (int i = 0; i < x.size(); ++i){
		    	BOOST_CHECK_CLOSE_FRACTION (xesperado[i],x[i], 0.01);
		    }
	}	
BOOST_AUTO_TEST_SUITE_END()

