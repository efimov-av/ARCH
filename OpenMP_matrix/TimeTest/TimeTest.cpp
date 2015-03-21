// TimeTest.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <omp.h>
#include <amp.h>
#include <conio.h>
#include <ctime>
#include <intrin.h>

#include <stdlib.h>

#include <math.h>
#include <iostream>

void Initialization(double* &A, double* &B, double* &C, int &Size);
void Multi(double* A, double* B, double* C, int Size);
void RandInit(double* Matrix, int Size);
void PMulti(double* A, double* B, double* C, int Size);

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	double cache_string=64;
	double t = 2.7*pow(10,-9);//ns
	double b = 12.6*pow(10,9);//mem speed GB/s
	double a = 12.5*pow(10, -9);//latence ns
	double n;
	int numThreads = 4;
	//omp_set_nested(true);
	omp_set_num_threads(numThreads);
		
		double* AM; // Матрица А
		double* BM; // Матрица B
		double* CM; // Матрица С
		
		long time1;
		
		int Size;  // Размер матриц
		// Инициализация данных
		Initialization(AM, BM,CM, Size);
		// Умножение матриц
		
		time1 = clock();
		Multi(AM, BM, CM, Size);
		//PMulti(AM, BM, CM, Size);
		
		time1 = clock() - time1;		
		printf_s("Time - %f s\n", time1 / 1000.0);

		
			double Tcalc;
			double Tmem;
			n = Size;
			/*Tcalc = pow(n, 2)*(2 * n - 1)*t;
			
			Tmem = (2 * pow(n, 3) + pow(n, 2))*(a + cache_string / b);*/

			Tcalc = (pow(n,2)/numThreads)*(2*n-1)*t;

			Tmem = (pow(n,3)+pow(n,3)/numThreads+pow(n,2))*(a+cache_string/b);
			

			//miss rate fromCodeXL
			double missratio = 0.1;//if n 500
			if (n == 1000.0)
			{
				//missratio = 0.3679;
				missratio = 0.19;
			}
			if (n == 1500.0)
			{
				//missratio = 0.3679;
				missratio = 0.21;
			}
			if (n == 2000.0)
			{
				//missratio = 0.3679;
				missratio = 0.15;
			}
			if (n == 2500.0)
			{
				//missratio = 0.3679;
				missratio = 0.2;
			}
			if (n == 3000.0)
			{
				//missratio = 0.3679;
				missratio = 0.19;
			}
			printf_s("Time calc- %f s\n", Tcalc);
			printf_s("Time mem- %f s\n", Tmem*missratio);
			printf_s("Time theory- %f s\n", Tcalc + Tmem*missratio);

			system("pause");
		
		delete[] AM;
		delete[] BM;
		delete[] CM;
		return 0;

}

// Функция выделения памяти и инициализации данных
void Initialization(double* &A,double* &B, double* &C, int &Size)
{
	int i, j;
	do {
		printf("\nВведите размер матриц: ");
		scanf_s("%d", &Size);		
		
	} while (Size <= 0);
	A = new double[Size*Size];
	B = new double[Size*Size];
	C = new double[Size*Size];
	for (i = 0; i<Size; i++)
		for (j = 0; j<Size; j++)
			C[i*Size + j] = 0;
	RandInit(A, Size);
	RandInit(B, Size);
}

void Multi(double* A,double* B, double* C, int Size) {
	int i, j, k;
	for (i = 0; i < Size; i++){
		for (j = 0; j < Size; j++){
			
			for (k = 0; k < Size; k++){
				
				 C[i*Size + j]+= A[i*Size + k] * B[k*Size + j];
			}
			 
		}
	}
}

void RandInit(double* Matrix, int Size) {
	srand(time(0));
	for (int i = 0; i<Size; i++) {
		for (int j = 0; j < Size; j++) Matrix[i*Size + j] = rand() / double(1000);
	}
}

void PMulti(double* A, double* B, double* C, int Size) {
	int row, col, inner;
	
#pragma omp parallel for private (col, inner)
	for (row = 0; row < Size; row++){
		for (col = 0; col < Size; col++){

			for (inner = 0; inner < Size; inner++){

				C[row*Size + col] += A[row*Size + inner] * B[inner*Size + col];
			}

		}
	}
}