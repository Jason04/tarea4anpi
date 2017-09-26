/*
main.cpp
Fecha:25/09/2017
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
#include "QRHouseholder.h"
class test{
	public:
		void testQR1(){
			qrhouseholder<float>QR;
		    anpi::Matrix<float> A = {{ 1.f, 0.5f,1.f/3.f},
                                     {0.5f,1.f/3.f,1.f/4.f},
                                     {1.f/3.f,1.f/4.f,1.f/5.f}};
                                     /*{{1, 0, 1},
		                            {1, 2, 0},
		                            {1, 0, 3}};*/
			anpi::Matrix<float> Ar(A.rows(),A.rows(),0.0);
		    double norma = 0.0;
		    norma = QR.testQR(A, Ar);

		    std::cout<<"**********TEST QR1 A4x4 START**********"<<std::endl;
		    std::cout<<"A matrix"<<std::endl;
		    A.printmatrix();
		    std::cout<<std::endl;
		    std::cout<<"A rebuild matrix"<<std::endl;
		    Ar.printmatrix();
		    std::cout<<std::endl;
		    std::cout<<"Norma between A and A rebuild :	"<<norma<<std::endl;
		    std::cout<<"**********TEST QR1 END**********"<<std::endl<<std::endl;
		}

		void testQR2(){
			qrhouseholder<float>QR;
		    anpi::Matrix<float> A = {{1.f/2.f, 1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f},
		    						{1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f},
		    						{1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f},
		    						{1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f, 1.f/9.f},
		    						{1.f/6.f, 1.f/7.f, 1.f/8.f, 1.f/9.f, 1.f/10.f}};
			anpi::Matrix<float> Ar(A.rows(),A.rows(),0.0);
		    double norma = 0.0;
		    norma = QR.testQR(A, Ar);

		    std::cout<<"**********TEST QR2 A5x5 START**********"<<std::endl;
		    std::cout<<"A matrix"<<std::endl;
		    A.printmatrix();
		    std::cout<<std::endl;
		    std::cout<<"A rebuild matrix"<<std::endl;
		    Ar.printmatrix();
		    std::cout<<std::endl;
		    std::cout<<"Norma between A and A rebuild :	"<<norma<<std::endl;
		    std::cout<<"**********TEST QR2 END**********"<<std::endl<<std::endl;		    
		}

		void testQR3(){
			qrhouseholder<float>QR;
		    anpi::Matrix<float> A = {	{-8,5,-6,4,1,-8,2,-7,3,-3},
			                            {-7,-5,-5,3,-3,-2,3,-1,2,9},
			                            {-3,5,-4,-2,2,-1,4,7,-6,-7},
			                            {-2,-4,5,8,-7,-1,2,-6,-5,-4},
			                            {-8,-8,2,-2,-7,-5,1,-5,5,-1},
			                            {5, 6,7,-5,-6, 4,2,4,-9, 5},
			                            {6,-4,9, 8,6, 2,-5,6,-5,-3},
			                            {-3,5,-6,-2,8,-4,-3,4,-8,-1},
			                            {-1,7,-2,-4,-8,4,-9,-5,-6,-3},
			                            {-6,-1,-4,-5,-8,8,8,-5,-7, 4}	};
			anpi::Matrix<float> Ar(A.rows(),A.rows(),0.0);
		    double norma = 0.0;
		    norma = QR.testQR(A, Ar);

		    std::cout<<"**********TEST QR3 A10x10 START**********"<<std::endl;
		    std::cout<<"A matrix"<<std::endl;
		    A.printmatrix();
		    std::cout<<std::endl;
		    std::cout<<"A rebuild matrix"<<std::endl;
		    Ar.printmatrix();
		    std::cout<<std::endl;
		    std::cout<<"Norma between A and A rebuild :	"<<norma<<std::endl;
		    std::cout<<"**********TEST QR3 END**********"<<std::endl<<std::endl;		    
		}

		void testSolveQR1(){
			qrhouseholder<float>QR;
			anpi::Matrix<float> A =  	{{1, 0, 1},
			                            {1, 2, 0},
			                            {1, 0, 3}};    
		    std::vector <float> b(3); //vector b
		    b[0]= 2;
		    b[1]= 4;
		    b[2]= 6;
		    std::vector <float> x(3); //vector x
		    QR.solveQR(A, x, b);
		    std::cout<<"**********TEST SOLVEQR1 A3x3 START**********"<<std::endl;
		    std::cout<<"Vector X: ";
		    for (int i = 0; i < x.size(); ++i){
		    	std::cout<<x[i]<<",";
		    }
		    std::cout<<std::endl;
		    
		    b = A * x;
		    std::cout<<"Vector B recalculado: ";
		    for (int i = 0; i < b.size(); ++i){
		    	std::cout<<b[i]<<",";
		    }
		    std::cout<<std::endl;
		    std::cout<<"**********TEST SOLVEQR1 END**********"<<std::endl<<std::endl;
		}

		void testSolveQR2(){
			qrhouseholder<float>QR;
			anpi::Matrix<float> A = {{1.f/2.f, 1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f},
		    						{1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f},
		    						{1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f},
		    						{1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f, 1.f/9.f},
		    						{1.f/6.f, 1.f/7.f, 1.f/8.f, 1.f/9.f, 1.f/10.f}};   
		    std::vector <float> b(5); //vector b
		    b[0]= 2;
		    b[1]= 4;
		    b[2]= 6;
		    b[3]= 8;
		    b[4]= 10;
		    std::vector <float> x(5); //vector x
		    QR.solveQR(A, x, b);

		    std::cout<<"**********TEST SOLVEQR2 A5x5 START**********"<<std::endl;
		    std::cout<<"Vector X: ";
		    for (int i = 0; i < x.size(); ++i){
		    	std::cout<<x[i]<<",";
		    }
		    std::cout<<std::endl;
		    
		    b = A * x;
		    std::cout<<"Vector B recalculado: ";
		    for (int i = 0; i < b.size(); ++i){
		    	std::cout<<b[i]<<",";
		    }
		    std::cout<<std::endl;
		    std::cout<<"**********TEST SOLVEQR2 END**********"<<std::endl<<std::endl;
		}

		void testSolveQR3(){
			qrhouseholder<float>QR;
			anpi::Matrix<float> A = {	{-8,5,-6,4,1,-8,2,-7,3,-3},
			                            {-7,-5,-5,3,-3,-2,3,-1,2,9},
			                            {-3,5,-4,-2,2,-1,4,7,-6,-7},
			                            {-2,-4,5,8,-7,-1,2,-6,-5,-4},
			                            {-8,-8,2,-2,-7,-5,1,-5,5,-1},
			                            {5, 6,7,-5,-6, 4,2,4,-9, 5},
			                            {6,-4,9, 8,6, 2,-5,6,-5,-3},
			                            {-3,5,-6,-2,8,-4,-3,4,-8,-1},
			                            {-1,7,-2,-4,-8,4,-9,-5,-6,-3},
			                            {-6,-1,-4,-5,-8,8,8,-5,-7, 4}	};    
		    std::vector <float> b(10); //vector b
		    b[0]= 2;
		    b[1]= 4;
		    b[2]= 6;
		    b[3]= 8;
		    b[4]= 10;
		    b[5]= 12;
		    b[6]= 14;
		    b[7]= 16;
		    b[8]= 18;
		    b[9]= 20;
		    std::vector <float> x(10); //vector x
		    QR.solveQR(A, x, b);

		    std::cout<<"**********TEST SOLVEQR2 A10x10 START**********"<<std::endl;
		    std::cout<<"Vector X: ";
		    for (int i = 0; i < x.size(); ++i){
		    	std::cout<<x[i]<<",";
		    }
		    std::cout<<std::endl;
		    
		    b = A * x;
		    std::cout<<"Vector B recalculado: ";
		    for (int i = 0; i < b.size(); ++i){
		    	std::cout<<b[i]<<",";
		    }
		    std::cout<<std::endl;
		    std::cout<<"**********TEST SOLVEQR2 END**********"<<std::endl<<std::endl;	
		}
	private:
};